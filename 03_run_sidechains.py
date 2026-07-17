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
        """New multi-match guardrail standardizer for novel NCAAs."""
        trimmed_parent = self._trim_parent(self.ref_mol)
        
        if target_mol.GetNumAtoms() < trimmed_parent.GetNumAtoms():
            target_mol = Chem.AddHs(target_mol)

        # --- ADD THESE TWO LINES TO FIX THE RINGINFO ERROR ---
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
            
        # 1. Enumerate ALL candidate matches
        ref_matches = trimmed_parent.GetSubstructMatches(common, uniquify=True)
        target_matches = target_mol.GetSubstructMatches(common, uniquify=True)
        
        # 2. Identify known backbone indices in trimmed parent
        bb_indices = {}
        for atom in trimmed_parent.GetAtoms():
            if atom.HasProp("pdb_name"):
                name = atom.GetProp("pdb_name").strip()
                if name in ["N", "CA", "C", "CB"]:
                    bb_indices[name] = atom.GetIdx()
        
        required_anchors = ["N", "CA", "C"]
        for a in required_anchors:
            if a not in bb_indices:
                raise ValueError(f"Parent {self.fetch_code} is missing critical backbone atom {a}.")
                
        # 3. Backbone-anchor verification and decision logic
        valid_mappings = []
        for r_match in ref_matches:
            r_match_set = set(r_match)
            # Check if this specific r_match shape covers N, CA, and C
            has_bb = all(bb_indices[a] in r_match_set for a in required_anchors)
            if not has_bb:
                continue
                
            for t_match in target_matches:
                mapping = {}
                for q_idx in range(len(common.GetAtoms())):
                    mapping[t_match[q_idx]] = r_match[q_idx]
                valid_mappings.append(mapping)
                
        if len(valid_mappings) == 0:
            raise ValueError(f"Zero candidate matches recovered the backbone (N, CA, C) for {self.res_name}.")
        elif len(valid_mappings) > 1:
            raise ValueError(f"Ambiguous matches! {len(valid_mappings)} candidates recovered the backbone for {self.res_name}.")
            
        final_mapping = valid_mappings[0]
        used_names = set()
        
        # 4. Transfer mapped names
        for target_idx, ref_idx in final_mapping.items():
            target_atom = target_mol.GetAtomWithIdx(target_idx)
            ref_atom = trimmed_parent.GetAtomWithIdx(ref_idx)
            
            if ref_atom.HasProp("pdb_name"):
                p_name = ref_atom.GetProp("pdb_name").strip()
                used_names.add(p_name)
                info = Chem.AtomPDBResidueInfo(
                    atomName=f"{p_name:<4}",
                    residueName=self.res_name,
                    residueNumber=1,
                    isHeteroAtom=True,
                )
                target_atom.SetMonomerInfo(info)
                
        # 5. Assign systematic names to unmatched unique NCAA atoms
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
            print(f"  -> Warning: No ETFlow output for {mod_code}. Skipping.")
            continue

        conformer_pos = results[smiles][0]
        mol_2d = get_mol_from_smiles(smiles)
        if mol_2d is None:
            continue

        mol_3d = set_rdmol_positions(mol_2d, conformer_pos)

        try:
            standardizer = PDBStandardizer(mod_code, parent_code)
            final_mol = standardizer.standardize(mol_3d)
        except Exception as exc:
            print(f"  -> Warning: Standardization failed for {mod_code} ({exc}). Saving raw ETFlow output.")
            final_mol = mol_3d

        out_path = os.path.join(output_dir, f"mod_{position}_{mod_code}.pdb")
        Chem.MolToPDBFile(final_mol, out_path)
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