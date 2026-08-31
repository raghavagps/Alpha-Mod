#!/usr/bin/env python3
"""
04_stitch.py

Stitch ETFlow-generated modified side chains into an AF2/ESMFold backbone.

Pipeline behavior:
1) Read modification positions/codes from a text file:
      position : modification_code : SMILES
2) Load target backbone structure (output/backbone.pdb).
3) For each modification, load output/mod_{position}_{code}.pdb,
   but first clean it in-memory by removing lines starting with CONECT or END.
4) Apply the same swap logic used in your older dist_batch_2.py:
   - sanitize source residue (remove OXT + all H atoms)
   - align using common anchors N, CA, C, CB
   - replace target residue at requested position
5) Save stitched structure as output/stitched.pdb.
"""

import argparse
import io
import os
import warnings
from typing import List, Tuple

from Bio import PDB
from Bio.PDB.PDBExceptions import PDBConstructionWarning


# Suppress noisy parser warnings (matches style in older script).
warnings.filterwarnings("ignore", category=PDBConstructionWarning)


def parse_modifications_file(mods_path: str) -> List[Tuple[int, str, bool]]:
    """
    Return (position, modification_code, is_novel) records.

    Generated formats:
        Legacy/RCSB: position : code : SMILES
        Novel:       position : code : SMILES : parent_letter

    The parser writes separators as " : ", so colons inside a SMILES string
    are not confused with field separators.
    """
    if not os.path.exists(mods_path):
        raise FileNotFoundError(f"Modifications file not found: {mods_path}")

    parsed: List[Tuple[int, str, bool]] = []

    with open(mods_path, "r", encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            parts = [part.strip() for part in line.split(" : ", 3)]

            if len(parts) not in (3, 4):
                raise ValueError(
                    f"Invalid line {line_no}: {line!r}. Expected either "
                    "'position : code : SMILES' or "
                    "'position : code : SMILES : parent'."
                )

            pos_text, mod_code = parts[0], parts[1]

            if not pos_text.isdigit():
                raise ValueError(
                    f"Invalid position on line {line_no}: {pos_text!r}"
                )

            if not mod_code:
                raise ValueError(
                    f"Missing modification code on line {line_no}."
                )

            is_novel = len(parts) == 4 and bool(parts[3])
            parsed.append((int(pos_text), mod_code, is_novel))

    return parsed

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


def perform_swap(target_struct, target_res_id: int, source_struct, chain_length: int, is_novel: bool, source_connectivity) -> bool:
    """
    Performs alignment and residue swap using logic adapted from dist_batch_2.py.

    POSITION-AWARE SANITIZATION:
    - Internal NCAA: remove OXT, HXT, H2, H3 (terminal atoms only)
    - N-terminal NCAA: remove OXT, HXT (keep H2, H3 — they belong at chain start)
    - C-terminal NCAA: remove H2, H3 (keep OXT, HXT — they belong at chain end)
    - ALL side-chain hydrogens are always kept.

    Alignment anchors:
    - N, CA, C, CB (requires at least 3 common atoms).
    """
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

    # =========================================================================
    # NOVEL_CONNECTIVITY_CLEANUP_V2
    #
    # Legacy/RCSB residues retain their established OXT/HXT/H2/H3 cleanup.
    #
    # Novel SMILES residues must not use those names to infer chemistry.
    # Names such as H3, O1 and O2 are systematic and do not necessarily mean
    # amino-terminal H3 or carboxyl-terminal OXT.
    #
    # For novel residues, identify terminal atoms using CONECT bonds:
    #   - backbone carbonyl oxygen = atom named O bonded to atom C
    #   - terminal carboxyl oxygen = the other oxygen bonded to C
    #   - amino hydrogen = a hydrogen directly bonded to N
    # =========================================================================
    is_n_terminal = target_res_id == 1
    is_c_terminal = target_res_id == chain_length

    pos_type = (
        "N-terminal" if is_n_terminal
        else ("C-terminal" if is_c_terminal else "internal")
    )

    def atom_element(atom):
        element = getattr(atom, "element", "")
        element = element.strip().upper() if element else ""

        if element:
            return element

        name = atom.name.strip().upper()
        return "H" if name.startswith("H") else name[:1]

    if is_novel:
        required_names = {"N", "CA", "C", "O"}
        source_names = {atom.name.strip() for atom in source_res}
        missing = required_names - source_names

        if missing:
            raise ValueError(
                f"Novel residue is missing required atom(s): {sorted(missing)}"
            )

        if not source_connectivity:
            raise ValueError(
                "Novel residue PDB contains no usable CONECT records; "
                "connectivity-based terminal cleanup cannot proceed safely."
            )

        atoms_by_serial = {
            int(atom.get_serial_number()): atom
            for atom in source_res
        }

        n_atom = source_res["N"]
        c_atom = source_res["C"]
        backbone_o_atom = source_res["O"]

        n_serial = int(n_atom.get_serial_number())
        c_serial = int(c_atom.get_serial_number())
        backbone_o_serial = int(backbone_o_atom.get_serial_number())

        detach_atom_ids = set()
        cleanup_messages = []

        # -------------------------------------------------------------
        # Remove the free carboxyl -OH group unless this is C-terminal.
        # -------------------------------------------------------------
        if not is_c_terminal:
            terminal_oxygen_serials = [
                serial
                for serial in source_connectivity.get(c_serial, set())
                if serial in atoms_by_serial
                and atom_element(atoms_by_serial[serial]) == "O"
                and serial != backbone_o_serial
            ]

            if len(terminal_oxygen_serials) != 1:
                found_names = [
                    atoms_by_serial[serial].name.strip()
                    for serial in terminal_oxygen_serials
                ]

                raise ValueError(
                    "Could not identify exactly one non-backbone carboxyl "
                    f"oxygen bonded to C; found {found_names!r}."
                )

            terminal_o_serial = terminal_oxygen_serials[0]
            terminal_o_atom = atoms_by_serial[terminal_o_serial]
            detach_atom_ids.add(terminal_o_atom.id)

            attached_hydrogens = [
                serial
                for serial in source_connectivity.get(
                    terminal_o_serial, set()
                )
                if serial in atoms_by_serial
                and atom_element(atoms_by_serial[serial]) == "H"
            ]

            for hydrogen_serial in attached_hydrogens:
                detach_atom_ids.add(atoms_by_serial[hydrogen_serial].id)

            removed_group = [terminal_o_atom.name.strip()] + [
                atoms_by_serial[serial].name.strip()
                for serial in attached_hydrogens
            ]

            cleanup_messages.append(
                f"terminal carboxyl group {removed_group}"
            )

        # -------------------------------------------------------------
        # Remove one hydrogen actually BONDED to N unless N-terminal.
        # Never remove a hydrogen merely because it is named H2 or H3.
        # -------------------------------------------------------------
        if not is_n_terminal:
            n_hydrogen_serials = sorted(
                serial
                for serial in source_connectivity.get(n_serial, set())
                if serial in atoms_by_serial
                and atom_element(atoms_by_serial[serial]) == "H"
            )

            if not n_hydrogen_serials:
                raise ValueError(
                    "Internal novel residue has no hydrogen bonded to N; "
                    "cannot form the expected peptide bond safely."
                )

            # Deterministic selection. For a normal free primary amino group,
            # the two N-bound hydrogens are chemically equivalent here.
            removed_n_h_serial = n_hydrogen_serials[-1]
            removed_n_h_atom = atoms_by_serial[removed_n_h_serial]
            detach_atom_ids.add(removed_n_h_atom.id)

            cleanup_messages.append(
                f"N-bound hydrogen {removed_n_h_atom.name.strip()}"
            )

        detached_names = []

        for atom_id in sorted(detach_atom_ids):
            if atom_id in source_res:
                detached_names.append(source_res[atom_id].name.strip())
                source_res.detach_child(atom_id)

        remaining_names = {atom.name.strip() for atom in source_res}
        missing_after_cleanup = required_names - remaining_names

        if missing_after_cleanup:
            raise ValueError(
                "Novel cleanup accidentally removed required backbone "
                f"atom(s): {sorted(missing_after_cleanup)}"
            )

        print(
            f"    Novel connectivity cleanup ({pos_type}): "
            f"removed {detached_names if detached_names else 'nothing'}"
        )

        if cleanup_messages:
            print(
                "    Chemical roles removed: "
                + "; ".join(cleanup_messages)
            )

        print(
            "    Side-chain hydrogens were not removed by H-number names."
        )

    else:
        # Preserve the existing RCSB/CCD naming behavior.
        carboxyl_terminal_atoms = {"OXT", "HXT"}
        amino_terminal_atoms = {"H2", "H3"}
        atoms_to_detach: List[str] = []

        for atom in source_res:
            atom_name = atom.name.strip()

            if (
                atom_name in carboxyl_terminal_atoms
                and not is_c_terminal
            ):
                atoms_to_detach.append(atom.name)
                continue

            if (
                atom_name in amino_terminal_atoms
                and not is_n_terminal
            ):
                atoms_to_detach.append(atom.name)
                continue

        for atom_name in atoms_to_detach:
            if atom_name in source_res:
                source_res.detach_child(atom_name)

        print(
            f"    Legacy RCSB cleanup ({pos_type}): removed "
            f"{atoms_to_detach if atoms_to_detach else 'nothing'}"
        )
    # =========================================================================

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
    modifications = parse_modifications_file(args.mods)
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

    for position, mod_code, is_novel in modifications:
        mod_filename = f"mod_{position}_{mod_code}.pdb"
        mod_path = os.path.join(args.mod_dir, mod_filename)

        print(f"\n[Swap] position={position}, code={mod_code}")

        if not os.path.exists(mod_path):
            print(f"    [!] Missing modification file: {mod_path}")
            failed_swaps.append(
                f"position {position} ({mod_code}): missing side-chain PDB"
            )
            continue

        try:
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
                is_novel,
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

    print("\n" + "=" * 60)
    print(f"Stitching complete. Applied {swaps_applied}/{len(modifications)} swaps.")
    print(f"Saved stitched structure: {args.out}")


if __name__ == "__main__":
    main()
