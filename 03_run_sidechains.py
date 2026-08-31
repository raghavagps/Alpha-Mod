#!/usr/bin/env python3
"""
03_run_sidechains.py

Generate one 3D side-chain structure per modified residue using ETFlow,
then standardize atom names to RCSB/PDB style.

Input format (required):
    position : modification_code : SMILES : parent_code (optional)
"""

import argparse
import os
from typing import List, Tuple

import requests
from rdkit import Chem
from rdkit.Chem import rdFMCS
from etflow import BaseFlow
from etflow.commons.covmat import set_rdmol_positions
from etflow.commons.featurization import get_mol_from_smiles

AA_1_TO_3 = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}

class PDBStandardizer:
    def __init__(self, res_name: str, parent_letter: str = ""):
        self.res_name = res_name.upper()
        self.parent_letter = parent_letter.upper() if parent_letter else ""
        
        if self.parent_letter and self.parent_letter in AA_1_TO_3:
            # New pathway: Novel NCAA anchored by parent structure
            self.fetch_code = AA_1_TO_3[self.parent_letter]
            self.is_novel = True
        else:
            # Legacy pathway: CCD lookup
            self.fetch_code = self.res_name
            self.is_novel = False
            
        self.ref_mol = self._fetch_structure_sdf(self.fetch_code)
        self.atom_names = self._fetch_names_cif(self.fetch_code)
        self._apply_ref_names()

    def _fetch_structure_sdf(self, res_name: str):
        url = f"https://files.rcsb.org/ligands/view/{res_name}_ideal.sdf"
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        mol = Chem.MolFromMolBlock(response.text, removeHs=False)
        if mol is None:
            raise ValueError(f"Could not parse SDF for residue '{res_name}'.")
        return mol

    def _fetch_names_cif(self, res_name: str) -> List[str]:
        url = f"https://files.rcsb.org/ligands/view/{res_name}.cif"
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        lines = response.text.splitlines()
        in_atom_loop = False
        names: List[str] = []
        for line in lines:
            line = line.strip()
            if line.startswith("loop_"):
                in_atom_loop = False
            if line.startswith("_chem_comp_atom."):
                in_atom_loop = True
            if in_atom_loop and line and not line.startswith("_") and not line.startswith("loop_"):
                parts = line.split()
                if len(parts) > 1:
                    names.append(parts[1].replace('"', ""))
        return names

    def _apply_ref_names(self) -> None:
        for idx, atom in enumerate(self.ref_mol.GetAtoms()):
            if idx < len(self.atom_names):
                atom.SetProp("pdb_name", self.atom_names[idx])

    def _trim_parent(self, mol):
        """Structurally pre-trim the free termini from the parent reference."""
        atoms_to_remove = []
        for atom in mol.GetAtoms():
            name = atom.GetProp("pdb_name").strip() if atom.HasProp("pdb_name") else ""
            if name in ["H", "H2", "H3", "OXT", "HXT"]:
                atoms_to_remove.append(atom.GetIdx())
        
        mw = Chem.RWMol(mol)
        # Sort descending to maintain index integrity during removal
        for idx in sorted(atoms_to_remove, reverse=True):
            mw.RemoveAtom(idx)
        return mw.GetMol()

    def standardize(self, target_mol):
        if self.is_novel:
            return self._standardize_novel(target_mol)
        else:
            return self._standardize_legacy(target_mol)

    def _standardize_legacy(self, target_mol):
        """Original direct CCD lookup standardization."""
        if self.ref_mol.GetNumAtoms() > target_mol.GetNumAtoms():
            target_mol = Chem.AddHs(target_mol)
            
        mcs = rdFMCS.FindMCS(
            [self.ref_mol, target_mol],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareOrder,
            ringMatchesRingOnly=True,
        )
        common = Chem.MolFromSmarts(mcs.smartsString)
        if common is None:
            return target_mol

        ref_match = self.ref_mol.GetSubstructMatch(common)
        target_match = target_mol.GetSubstructMatch(common)

        if not ref_match or not target_match:
            return target_mol

        idx_map = dict(zip(target_match, ref_match))
        for target_idx, ref_idx in idx_map.items():
            target_atom = target_mol.GetAtomWithIdx(target_idx)
            ref_atom = self.ref_mol.GetAtomWithIdx(ref_idx)

            if ref_atom.HasProp("pdb_name"):
                info = Chem.AtomPDBResidueInfo(
                    atomName=f"{ref_atom.GetProp('pdb_name'):<4}",
                    residueName=self.res_name,
                    residueNumber=1,
                    isHeteroAtom=True,
                )
                target_atom.SetMonomerInfo(info)
        return target_mol

    def _standardize_novel(self, target_mol):
        """CONSENSUS_MCS_FALLBACK_V2_ZERO_SAFE: safely name ambiguous novel NCAAs."""
        trimmed_parent = self._trim_parent(self.ref_mol)

        if target_mol.GetNumAtoms() < trimmed_parent.GetNumAtoms():
            target_mol = Chem.AddHs(target_mol)

        Chem.GetSSSR(trimmed_parent)
        Chem.GetSSSR(target_mol)

        mcs = rdFMCS.FindMCS(
            [trimmed_parent, target_mol],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareOrder,
            ringMatchesRingOnly=True,
        )
        common = Chem.MolFromSmarts(mcs.smartsString)
        if common is None:
            raise ValueError("MCS SMARTS construction failed.")

        ref_matches = trimmed_parent.GetSubstructMatches(common, uniquify=True)
        target_matches = target_mol.GetSubstructMatches(common, uniquify=True)

        bb_indices = {}
        for atom in trimmed_parent.GetAtoms():
            if atom.HasProp("pdb_name"):
                name = atom.GetProp("pdb_name").strip()
                if name in ["N", "CA", "C", "CB"]:
                    bb_indices[name] = atom.GetIdx()

        required_anchors = ["N", "CA", "C"]
        for anchor in required_anchors:
            if anchor not in bb_indices:
                raise ValueError(
                    f"Parent {self.fetch_code} is missing critical backbone atom {anchor}."
                )

        valid_mappings = []
        for ref_match in ref_matches:
            if not all(bb_indices[anchor] in set(ref_match) for anchor in required_anchors):
                continue

            for target_match in target_matches:
                valid_mappings.append(
                    {
                        target_match[q_idx]: ref_match[q_idx]
                        for q_idx in range(len(common.GetAtoms()))
                    }
                )

        if not valid_mappings:
            # The global maximum common substructure can prefer a side-chain
            # fragment and omit the complete parent N/CA/C backbone. This does
            # not mean the novel molecule lacks an amino-acid backbone.
            #
            # Do not fail here. The existing `else` branch below will identify
            # N--CA--C(=O)O directly from the novel molecular graph.
            print(
                f"  [WARNING] Zero backbone-containing global MCS mappings for "
                f"{self.res_name}; using deterministic direct-backbone fallback."
            )

        def ref_atom_name(ref_idx):
            atom = trimmed_parent.GetAtomWithIdx(ref_idx)
            return atom.GetProp("pdb_name").strip() if atom.HasProp("pdb_name") else ""

        assigned_names = {}

        if len(valid_mappings) == 1:
            # Keep the previous behavior for an unambiguous novel residue.
            for target_idx, ref_idx in valid_mappings[0].items():
                name = ref_atom_name(ref_idx)
                if name:
                    assigned_names[target_idx] = name

        else:
            print(
                f"  [WARNING] {len(valid_mappings)} MCS mappings for {self.res_name}; "
                "using consensus names plus deterministic backbone naming."
            )

            # Transfer a parent name only when every valid mapping agrees on
            # exactly the same target atom.
            for ref_atom in trimmed_parent.GetAtoms():
                ref_idx = ref_atom.GetIdx()
                name = ref_atom_name(ref_idx)
                if not name:
                    continue

                mapped_target_indices = []

                for mapping in valid_mappings:
                    matches = [
                        target_idx
                        for target_idx, mapped_ref_idx in mapping.items()
                        if mapped_ref_idx == ref_idx
                    ]

                    if len(matches) != 1:
                        mapped_target_indices = []
                        break

                    mapped_target_indices.append(matches[0])

                if (
                    len(mapped_target_indices) == len(valid_mappings)
                    and len(set(mapped_target_indices)) == 1
                ):
                    assigned_names[mapped_target_indices[0]] = name

            def is_carbonyl_carbon(atom):
                return any(
                    neighbor.GetAtomicNum() == 8
                    and bond.GetBondType() == Chem.BondType.DOUBLE
                    for bond in atom.GetBonds()
                    for neighbor in [bond.GetOtherAtom(atom)]
                )

            # Identify one alpha-amino-acid backbone directly from the novel SMILES:
            # N -- CA -- C(=O)O
            alpha_candidates = []

            for atom in target_mol.GetAtoms():
                if atom.GetAtomicNum() != 6:
                    continue

                nitrogen_neighbors = [
                    neighbor for neighbor in atom.GetNeighbors()
                    if neighbor.GetAtomicNum() == 7
                ]

                carbonyl_neighbors = [
                    neighbor for neighbor in atom.GetNeighbors()
                    if neighbor.GetAtomicNum() == 6 and is_carbonyl_carbon(neighbor)
                ]

                if len(nitrogen_neighbors) == 1 and len(carbonyl_neighbors) == 1:
                    alpha_candidates.append(
                        (atom, nitrogen_neighbors[0], carbonyl_neighbors[0])
                    )

            if len(alpha_candidates) != 1:
                raise ValueError(
                    f"Ambiguous MCS fallback could not identify one alpha-amino-acid "
                    f"backbone for {self.res_name} "
                    f"(found {len(alpha_candidates)} candidates)."
                )

            alpha_atom, nitrogen_atom, carbonyl_atom = alpha_candidates[0]

            double_oxygen_neighbors = [
                bond.GetOtherAtom(carbonyl_atom)
                for bond in carbonyl_atom.GetBonds()
                if bond.GetBondType() == Chem.BondType.DOUBLE
                and bond.GetOtherAtom(carbonyl_atom).GetAtomicNum() == 8
            ]

            if len(double_oxygen_neighbors) != 1:
                raise ValueError(
                    f"Ambiguous MCS fallback could not identify one carbonyl oxygen "
                    f"for {self.res_name}."
                )

            def assign_certain_name(target_idx, name):
                # Keep PDB atom names unique within this residue.
                for old_target_idx, old_name in list(assigned_names.items()):
                    if old_name == name and old_target_idx != target_idx:
                        del assigned_names[old_target_idx]

                assigned_names[target_idx] = name

            # These atoms are chemically unambiguous in an alpha amino acid.
            assign_certain_name(nitrogen_atom.GetIdx(), "N")
            assign_certain_name(alpha_atom.GetIdx(), "CA")
            assign_certain_name(carbonyl_atom.GetIdx(), "C")
            assign_certain_name(double_oxygen_neighbors[0].GetIdx(), "O")

            # CB is assigned only when CA has exactly one heavy side-chain neighbor.
            # AIB has two such carbons, so neither is arbitrarily called CB.
            sidechain_neighbors = [
                neighbor
                for neighbor in alpha_atom.GetNeighbors()
                if neighbor.GetIdx() not in {
                    nitrogen_atom.GetIdx(),
                    carbonyl_atom.GetIdx(),
                }
                and neighbor.GetAtomicNum() != 1
            ]

            if len(sidechain_neighbors) == 1:
                assign_certain_name(sidechain_neighbors[0].GetIdx(), "CB")
                print("  [INFO] Fallback anchors: N, CA, C, O, CB.")

            else:
                for old_target_idx, old_name in list(assigned_names.items()):
                    if old_name == "CB":
                        del assigned_names[old_target_idx]

                print(
                    "  [INFO] Fallback anchors: N, CA, C, O; "
                    "CB is absent or not uniquely defined."
                )

        used_names = set()

        for target_idx, pdb_name in assigned_names.items():
            target_atom = target_mol.GetAtomWithIdx(target_idx)
            used_names.add(pdb_name)

            info = Chem.AtomPDBResidueInfo(
                atomName=f"{pdb_name:<4}",
                residueName=self.res_name,
                residueNumber=1,
                isHeteroAtom=True,
            )
            target_atom.SetMonomerInfo(info)

        # Any atom whose parent identity is not proven receives a deterministic
        # PDB-safe name: C1, C2, O1, H1, and so on.
        element_counts = {}

        for atom in target_mol.GetAtoms():
            if not atom.GetMonomerInfo() or not atom.GetMonomerInfo().GetName().strip():
                sym = atom.GetSymbol().upper()
                element_counts[sym] = element_counts.get(sym, 0) + 1

                new_name = f"{sym}{element_counts[sym]}"

                while new_name in used_names:
                    element_counts[sym] += 1
                    new_name = f"{sym}{element_counts[sym]}"

                used_names.add(new_name)

                info = Chem.AtomPDBResidueInfo(
                    atomName=f"{new_name:<4}",
                    residueName=self.res_name,
                    residueNumber=1,
                    isHeteroAtom=True,
                )
                atom.SetMonomerInfo(info)

        return target_mol

def parse_modifications_file(mods_path: str) -> List[Tuple[int, str, str, str]]:
    if not os.path.exists(mods_path):
        raise FileNotFoundError(f"Modifications file not found: {mods_path}")

    records: List[Tuple[int, str, str, str]] = []
    with open(mods_path, "r", encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            parts = [p.strip() for p in line.split(":")]
            if len(parts) < 3:
                raise ValueError(f"Invalid format at line {line_no}: {line!r}")

            position = int(parts[0])
            mod_code = parts[1]
            smiles = parts[2]
            parent_code = parts[3] if len(parts) >= 4 else ""

            records.append((position, mod_code, smiles, parent_code))
    return records


def validate_sidechain_pdb(pdb_path: str, mod_code: str) -> None:
    """Reject raw or improperly named PDB output before it reaches the stitcher."""
    with open(pdb_path, "r", encoding="utf-8") as handle:
        atom_lines = [
            line for line in handle
            if line.startswith(("ATOM  ", "HETATM"))
        ]

    if not atom_lines:
        raise ValueError("PDB contains no ATOM/HETATM records.")

    residue_names = {line[17:20].strip() for line in atom_lines}
    expected_residue_name = mod_code.upper()[:3]

    if residue_names != {expected_residue_name}:
        raise ValueError(
            f"Expected residue name {expected_residue_name!r}, "
            f"found {sorted(residue_names)!r}."
        )

    if "UNL" in residue_names:
        raise ValueError("Raw RDKit residue name 'UNL' is not acceptable.")

    atom_names = {line[12:16].strip() for line in atom_lines}
    missing_anchors = {"N", "CA", "C"} - atom_names

    if missing_anchors:
        raise ValueError(
            f"Missing peptide-backbone anchor atom(s): "
            f"{sorted(missing_anchors)!r}."
        )


def run_etflow_and_save_pdbs(records: List[Tuple[int, str, str, str]], output_dir: str, cache_dir: str) -> None:
    if not records:
        print("No modification records found. Nothing to generate.")
        return

    os.makedirs(output_dir, exist_ok=True)

    print("Loading ETFlow model...")
    model = BaseFlow.from_default(model="drugs-o3", cache=cache_dir)

    all_smiles = [rec[2] for rec in records]
    print("Running ETFlow batch prediction with num_samples=1...")
    results = model.predict(all_smiles, num_samples=1, as_mol=False)

    for position, mod_code, smiles, parent_code in records:
        print(f"Processing position={position}, code={mod_code}...")

        if smiles not in results or len(results[smiles]) == 0:
            raise RuntimeError(f"No ETFlow output for {mod_code} ({smiles}).")

        conformer_pos = results[smiles][0]
        mol_2d = get_mol_from_smiles(smiles)

        if mol_2d is None:
            raise RuntimeError(f"RDKit could not parse SMILES for {mod_code}: {smiles}")

        mol_3d = set_rdmol_positions(mol_2d, conformer_pos)

        try:
            standardizer = PDBStandardizer(mod_code, parent_code)
            final_mol = standardizer.standardize(mol_3d)
        except Exception as exc:
            raise RuntimeError(
                f"PDB standardization failed for {mod_code}; "
                "raw RDKit output will not be used."
            ) from exc

        out_path = os.path.join(output_dir, f"mod_{position}_{mod_code}.pdb")
        tmp_path = f"{out_path}.tmp"

        try:
            Chem.MolToPDBFile(final_mol, tmp_path)
            validate_sidechain_pdb(tmp_path, mod_code)
            os.replace(tmp_path, out_path)
        except Exception:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
            raise

        print(f"  -> Saved: {out_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate standardized 3D side-chain PDBs from modifications.txt using ETFlow.")
    parser.add_argument("--mods", required=True)
    parser.add_argument("--cache", default="./cache")
    parser.add_argument("--out_dir", required=True)

    args = parser.parse_args()
    records = parse_modifications_file(args.mods)
    run_etflow_and_save_pdbs(records, output_dir=args.out_dir, cache_dir=args.cache)


if __name__ == "__main__":
    main()
