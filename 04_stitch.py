#!/usr/bin/env python3
"""
04_stitch.py

Stitch ETFlow-generated modified side chains into an AF2/ESMFold backbone.

For each source residue, verify the molecular graph from script 3 and plan the
atom removals required by its neighbors. Align the cleaned residue to the
canonical backbone and write the stitched PDB.
"""

import argparse
import io
import os
import warnings
from typing import List

from residue_metadata import load_metadata, record_pdb_stage, sha256_file

ROUTING_CONTRACT_VERSION = 1

from Bio import PDB
from Bio.PDB.PDBExceptions import PDBConstructionWarning


# Suppress noisy parser warnings (matches style in older script).
warnings.filterwarnings("ignore", category=PDBConstructionWarning)


def read_pdb_connectivity(pdb_path: str):
    """
    Read undirected atom connectivity from PDB CONECT records.

    Repeated serial numbers used to represent double bonds are collapsed here;
    terminal-atom identification needs adjacency, not bond order.
    """
    adjacency = {}

    with open(pdb_path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            if not raw_line.startswith("CONECT"):
                continue

            serials = []

            for start in range(6, len(raw_line), 5):
                field = raw_line[start:start + 5].strip()
                if field:
                    try:
                        serials.append(int(field))
                    except ValueError:
                        pass

            if len(serials) < 2:
                continue

            source_serial = serials[0]
            adjacency.setdefault(source_serial, set())

            for target_serial in serials[1:]:
                adjacency.setdefault(source_serial, set()).add(target_serial)
                adjacency.setdefault(target_serial, set()).add(source_serial)

    return adjacency


def clean_pdb_text_in_memory(pdb_path: str) -> str:
    """
    Read a PDB file and remove lines starting with CONECT or END in-memory.

    This is equivalent to your clean_et_flow_mod_pdb.py filtering logic,
    but avoids writing temporary files.
    """
    if not os.path.exists(pdb_path):
        raise FileNotFoundError(f"Modification PDB not found: {pdb_path}")

    cleaned_lines: List[str] = []

    with open(pdb_path, "r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("CONECT") or line.startswith("END"):
                continue
            cleaned_lines.append(line)

    return "".join(cleaned_lines)


def load_source_structure_from_cleaned_text(parser: PDB.PDBParser, structure_id: str, cleaned_pdb_text: str):
    """
    Load a cleaned PDB string into BioPython using StringIO (in-memory file handle).
    """
    return parser.get_structure(structure_id, io.StringIO(cleaned_pdb_text))


def clean_source_termini(source_res, graph, pdb_connectivity, *, incoming, outgoing):
    """Plan and apply the removals for the ordinary peptide joins we support."""
    atoms = {int(serial): data for serial, data in graph["atoms"].items()}
    pdb_atoms = {atom.get_serial_number(): atom for atom in source_res}
    if atoms.keys() != pdb_atoms.keys() or any(
        (data["name"], data["element"]) != (pdb_atoms[serial].name.strip(), pdb_atoms[serial].element)
        for serial, data in atoms.items()
    ):
        raise ValueError("Standalone PDB atom names or serials disagree with its molecular graph.")

    bonds = {serial: {} for serial in atoms}
    for first, second, order in graph["bonds"]:
        bonds[first][second] = bonds[second][first] = order
    if {serial: set(neighbors) for serial, neighbors in bonds.items() if neighbors} != pdb_connectivity:
        raise ValueError("Standalone PDB CONECT bonds disagree with its molecular graph.")

    ids = {data["name"]: serial for serial, data in atoms.items()}
    for name, element in (("N", "N"), ("CA", "C"), ("C", "C")):
        if name not in ids or atoms[ids[name]]["element"] != element:
            raise ValueError(f"Cannot identify the backbone atom {name} in the molecular graph.")
    n, ca, c = (ids[name] for name in ("N", "CA", "C"))
    cb = ids.get("CB")
    # Accept the existing named alpha anchors and the N-CB-CA-C placement.
    n_to_ca = bonds[n].get(ca) == 1 or (
        bonds[n].get(cb) == 1 and bonds[ca].get(cb) == 1 and atoms[cb]["element"] == "C"
    )
    if not n_to_ca or bonds[ca].get(c) != 1:
        raise ValueError("Named N, CA, and C do not form a supported backbone in the molecular graph.")

    remove = set()
    if incoming:
        if atoms[n]["charge"] != 0 or sorted(bonds[n].values()) != [1, 1, 1]:
            raise ValueError("Incoming peptide bond requires a neutral backbone N with three single bonds.")
        n_hydrogens = [
            serial for serial in bonds[n]
            if atoms[serial]["element"] == "H" and len(bonds[serial]) == 1
        ]
        if not n_hydrogens:
            raise ValueError("Incoming peptide bond requires a removable hydrogen bonded to backbone N.")
        # Preserve the previous hybrid tie-break for multiple N-bound hydrogens.
        remove.add(max(n_hydrogens))

    if outgoing:
        oxygen = ids.get("O")
        leaving_oxygens = [
            serial for serial, order in bonds[c].items()
            if order == 1 and atoms[serial]["element"] == "O"
        ]
        if (
            bonds[c].get(oxygen) != 2 or atoms[oxygen]["element"] != "O"
            or len(leaving_oxygens) != 1 or len(bonds[c]) != 3
            or atoms[c]["charge"] != 0 or atoms[oxygen]["charge"] != 0
            or len(bonds[oxygen]) != 1
        ):
            raise ValueError("Outgoing peptide bond requires a supported free carboxyl group on backbone C.")
        leaving = leaving_oxygens[0]
        attached = set(bonds[leaving]) - {c}
        is_acid = atoms[leaving]["charge"] == 0 and len(attached) == 1 and all(
            atoms[serial]["element"] == "H" and bonds[serial] == {leaving: 1}
            for serial in attached
        )
        is_carboxylate = atoms[leaving]["charge"] == -1 and not attached
        if not (is_acid or is_carboxylate):
            raise ValueError("Outgoing peptide bond requires a free carboxyl group; ester or blocked oxygen detected.")
        remove.update(attached | {leaving})

    def reachable(excluded):
        visited = set()
        todo = [ca]
        while todo:
            serial = todo.pop()
            if serial not in visited and serial not in excluded:
                visited.add(serial)
                todo.extend(bonds[serial])
        return visited

    if remove and reachable(set()) - remove != reachable(remove):
        raise ValueError("Terminal cleanup would leave a disconnected fragment.")

    removed_names = [atoms[serial]["name"] for serial in sorted(remove)]
    for serial in remove:
        source_res.detach_child(pdb_atoms[serial].id)
    print(f"    Bond-based cleanup: removed {removed_names if removed_names else 'nothing'}")


def perform_swap(target_struct, target_res_id: int, source_struct, chain_length: int,
                 source_graph, source_connectivity) -> bool:
    """Clean the source termini, align by the backbone anchors, and replace the residue."""
    try:
        target_model = target_struct[0]

        # --- BLIND CHAIN DETECTION START (from older logic) ---
        chains = list(target_model.get_chains())

        if len(chains) == 0:
            print("    [!] Error: Target structure has no chains.")
            return False

        # Assume peptide is in the first chain for this single-sequence pipeline.
        target_chain = chains[0]

        if len(chains) > 1:
            print(
                f"    [!] Warning: Target has multiple chains ({len(chains)}). "
                f"Operating on first one: {target_chain.id}"
            )
        # --- BLIND CHAIN DETECTION END ---

        # Chain membership/indexing accepts int residue IDs in BioPython.
        if target_res_id not in target_chain:
            print(f"    [!] Error: Residue {target_res_id} not found in chain {target_chain.id}.")
            return False

        target_res = target_chain[target_res_id]

        source_residues = list(source_struct.get_residues())
        if not source_residues:
            print("    [!] Error: Source structure has no residues.")
            return False
        source_res = source_residues[0]

    except (KeyError, IndexError) as exc:
        print(f"    [!] Error accessing residues: {exc}")
        return False

    clean_source_termini(
        source_res, source_graph, source_connectivity,
        incoming=target_res_id > 1, outgoing=target_res_id < chain_length,
    )

    # --- Alignment logic (same anchors as old script) ---
    potential_anchors = ["N", "CA", "C", "CB"]
    fixed_atoms = []   # Atoms in target residue (reference frame)
    moving_atoms = []  # Corresponding atoms in source residue (to move)

    for atom_name in potential_anchors:
        if atom_name in target_res and atom_name in source_res:
            if not target_res[atom_name].is_disordered() and not source_res[atom_name].is_disordered():
                fixed_atoms.append(target_res[atom_name])
                moving_atoms.append(source_res[atom_name])

    # Need at least 3 atoms to define a stable 3D superimposition.
    if len(fixed_atoms) < 3:
        print(f"    [!] Critical: Fewer than 3 common backbone atoms at pos {target_res_id}. Skipping.")
        return False

    # Apply superimposition to move the source residue into target frame.
    superimposer = PDB.Superimposer()
    superimposer.set_atoms(fixed_atoms, moving_atoms)
    superimposer.apply(source_res.get_atoms())

    # --- Swapping logic ---
    old_id = target_res.id
    target_chain.detach_child(old_id)

    # Reuse target residue ID so sequence numbering is preserved.
    source_res.id = old_id
    target_chain.add(source_res)

    # Keep residues ordered by sequence number.
    target_chain.child_list.sort(key=lambda x: x.id[1])

    print(
        f"    -> Swapped pos {target_res_id} on Chain {target_chain.id} "
        f"with {source_res.get_resname()} (RMSD: {superimposer.rms:.4f} A)"
    )
    return True


def main() -> None:
    """CLI entry point for single-structure stitching."""
    parser = argparse.ArgumentParser(
        description="Stitch ETFlow modification residues into backbone PDB using BioPython superimposition."
    )
    parser.add_argument(
        "--mods",
        required=True,
        help="Path to modifications mapping text file",
    )
    parser.add_argument(
        "--backbone",
        required=True,
        help="Path to backbone PDB to modify",
    )
    parser.add_argument(
        "--mod_dir",
        required=True,
        help="Directory containing side-chain PDB files mod_{position}_{code}.pdb",
    )
    parser.add_argument(
        "--out",
        required=True,
        help="Output stitched PDB path",
    )

    args = parser.parse_args()

    # Read modification definitions (position + code).
    modifications = load_metadata(args.mods)["records"]
    if not modifications:
        print("No modifications found in mapping file; saving backbone unchanged.")

    # Load target backbone structure.
    if not os.path.exists(args.backbone):
        raise FileNotFoundError(f"Backbone PDB not found: {args.backbone}")

    pdb_parser = PDB.PDBParser(QUIET=True)
    target_struct = pdb_parser.get_structure("target", args.backbone)

    # Count residues in first chain so perform_swap knows about terminal positions.
    target_chain = list(target_struct[0].get_chains())[0]
    chain_length = len(list(target_chain.get_residues()))
    print(f"[Stitch] Chain has {chain_length} residues.")

    # Apply each swap in order from the modifications file.
    swaps_applied = 0
    failed_swaps: List[str] = []

    for record in modifications:
        position, mod_code = record["position"], record["code"]
        mod_filename = record["modification_file"]
        mod_path = os.path.join(args.mod_dir, mod_filename)

        print(f"\n[Swap] position={position}, code={mod_code}")

        if not os.path.exists(mod_path):
            print(f"    [!] Missing modification file: {mod_path}")
            failed_swaps.append(
                f"position {position} ({mod_code}): missing side-chain PDB"
            )
            continue

        try:
            source_graph = record.get("stitch_graph")
            if source_graph is None or source_graph.get("pdb_sha256") != sha256_file(mod_path):
                raise ValueError(
                    "Missing or stale standalone residue graph. "
                    "Run the updated 03_run_sidechains.py to regenerate it."
                )
            # Read chemical adjacency before removing CONECT records.
            source_connectivity = read_pdb_connectivity(mod_path)

            # Clean source PDB text in-memory for Bio.PDB parsing.
            cleaned_text = clean_pdb_text_in_memory(mod_path)

            # 2) Parse cleaned string directly via StringIO.
            source_struct = load_source_structure_from_cleaned_text(
                pdb_parser,
                structure_id=f"source_{position}_{mod_code}",
                cleaned_pdb_text=cleaned_text,
            )

            # 3) Perform alignment + replacement.
            if perform_swap(
                target_struct,
                position,
                source_struct,
                chain_length,
                source_graph,
                source_connectivity,
            ):
                swaps_applied += 1
            else:
                failed_swaps.append(
                    f"position {position} ({mod_code}): swap was rejected"
                )

        except Exception as exc:
            print(f"    [!] Failed to process {mod_filename}: {exc}")
            failed_swaps.append(f"position {position} ({mod_code}): {exc}")

    if failed_swaps:
        details = "\n  - ".join(failed_swaps)
        raise RuntimeError(
            "Stitching aborted; no stitched PDB was written because one or more "
            f"modifications failed:\n  - {details}"
        )

    # Ensure output directory exists.
    out_dir = os.path.dirname(args.out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    # Save final stitched structure (or unchanged backbone if no swaps succeeded).
    io_writer = PDB.PDBIO()
    io_writer.set_structure(target_struct)
    io_writer.save(args.out)

    try:
        record_pdb_stage(args.mods, args.out, "stitched")
    except Exception:
        if os.path.exists(args.out):
            os.remove(args.out)
        raise

    print("\n" + "=" * 60)
    print(f"Stitching complete. Applied {swaps_applied}/{len(modifications)} swaps.")
    print(f"Saved stitched structure: {args.out}")


if __name__ == "__main__":
    main()
