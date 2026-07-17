#!/usr/bin/env python3
"""
03_run_sidechains.py

Generate one 3D side-chain structure per modified residue using ETFlow,
then standardize atom names to RCSB/PDB style using the notebook's
PDBStandardizer logic.

Input format (required):
    position : modification_code : SMILES

Example:
    4 : 5PG : CN[C@H](C(O)=O)c1ccc(O)cc1

Outputs:
    <out_dir>/mod_{position}_{modification_code}.pdb
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


# -----------------------------------------------------------------------------
# PDB standardization class copied/adapted from notebook logic.
# It:
# 1) Downloads the reference ligand structure (ideal SDF) from RCSB.
# 2) Downloads official atom names from the ligand CIF.
# 3) Uses MCS atom mapping to transfer official PDB atom names onto ETFlow output.
# -----------------------------------------------------------------------------
class PDBStandardizer:
    def __init__(self, res_name: str):
        self.res_name = res_name.upper()

        # Fetch reference geometry and official atom names from RCSB.
        self.ref_mol = self._fetch_structure_sdf(self.res_name)
        self.atom_names = self._fetch_names_cif(self.res_name)

        # Attach names to reference atoms so they can be transferred by mapping.
        self._apply_ref_names()

    def _fetch_structure_sdf(self, res_name: str):
        """Download the RCSB ideal SDF for this residue code."""
        url = f"https://files.rcsb.org/ligands/view/{res_name}_ideal.sdf"
        print(f"Fetching reference structure for {res_name}...")
        response = requests.get(url, timeout=60)
        response.raise_for_status()

        mol = Chem.MolFromMolBlock(response.text, removeHs=False)
        if mol is None:
            raise ValueError(f"Could not parse SDF for residue '{res_name}'.")
        return mol

    def _fetch_names_cif(self, res_name: str) -> List[str]:
        """Download ligand CIF and extract official atom names."""
        url = f"https://files.rcsb.org/ligands/view/{res_name}.cif"
        print(f"Fetching official names for {res_name}...")
        response = requests.get(url, timeout=60)
        response.raise_for_status()

        lines = response.text.splitlines()
        in_atom_loop = False
        names: List[str] = []

        for line in lines:
            line = line.strip()

            # CIF uses loop_ blocks; when a new loop starts, reset state.
            if line.startswith("loop_"):
                in_atom_loop = False

            # Start collecting atom table when chem_comp_atom fields appear.
            if line.startswith("_chem_comp_atom."):
                in_atom_loop = True

            # Collect rows from the atom loop and grab the atom_id column.
            if in_atom_loop and line and not line.startswith("_") and not line.startswith("loop_"):
                parts = line.split()
                if len(parts) > 1:
                    names.append(parts[1].replace('"', ""))

        return names

    def _apply_ref_names(self) -> None:
        """Attach extracted atom names onto reference atoms as 'pdb_name'."""
        for idx, atom in enumerate(self.ref_mol.GetAtoms()):
            if idx < len(self.atom_names):
                atom.SetProp("pdb_name", self.atom_names[idx])

    def standardize(self, target_mol):
        """
        Map ETFlow-generated target molecule onto reference molecule via MCS,
        then transfer PDB atom names + residue metadata.
        """
        # If target lacks hydrogens relative to reference, add Hs for better mapping.
        if self.ref_mol.GetNumAtoms() > target_mol.GetNumAtoms():
            target_mol = Chem.AddHs(target_mol)

        # Maximum common substructure mapping between reference and target.
        mcs = rdFMCS.FindMCS(
            [self.ref_mol, target_mol],
            atomCompare=rdFMCS.AtomCompare.CompareElements,
            bondCompare=rdFMCS.BondCompare.CompareOrder,
            ringMatchesRingOnly=True,
        )

        common = Chem.MolFromSmarts(mcs.smartsString)
        if common is None:
            print("Warning: MCS SMARTS construction failed; returning unstandardized molecule.")
            return target_mol

        ref_match = self.ref_mol.GetSubstructMatch(common)
        target_match = target_mol.GetSubstructMatch(common)

        # If mapping fails, keep raw ETFlow naming rather than crashing the run.
        if not ref_match or not target_match:
            print("Warning: Could not perfectly map atoms. Names might be incomplete.")
            return target_mol

        idx_map = dict(zip(target_match, ref_match))

        # Transfer atom names and residue information onto target atoms.
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


def parse_modifications_file(mods_path: str) -> List[Tuple[int, str, str]]:
    """
    Parse the modifications text file.

    Expected strict format per non-empty line:
        position : modification_code : SMILES

    Returns:
        List of tuples: (position, modification_code, smiles)
    """
    if not os.path.exists(mods_path):
        raise FileNotFoundError(f"Modifications file not found: {mods_path}")

    records: List[Tuple[int, str, str]] = []

    with open(mods_path, "r", encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()

            # Skip empty lines quietly.
            if not line:
                continue

            # Split into exactly three fields (position, code, smiles).
            parts = [p.strip() for p in line.split(":", 2)]
            if len(parts) != 3:
                raise ValueError(
                    f"Invalid format at line {line_no}: {line!r}. "
                    "Expected: position : modification_code : SMILES"
                )

            pos_text, mod_code, smiles = parts

            if not pos_text.isdigit():
                raise ValueError(f"Invalid position at line {line_no}: {pos_text!r}")

            position = int(pos_text)
            if position <= 0:
                raise ValueError(f"Position must be >= 1 at line {line_no}: {position}")

            if not mod_code:
                raise ValueError(f"Missing modification code at line {line_no}.")
            if not smiles:
                raise ValueError(f"Missing SMILES at line {line_no}.")

            records.append((position, mod_code, smiles))

    return records


def run_etflow_and_save_pdbs(records: List[Tuple[int, str, str]], output_dir: str, cache_dir: str) -> None:
    """
    Run ETFlow with num_samples=1, standardize each output molecule, and write PDB files.
    """
    if not records:
        print("No modification records found. Nothing to generate.")
        return

    os.makedirs(output_dir, exist_ok=True)

    # Initialize ETFlow model exactly as notebook style (drugs-o3) with local cache.
    print("Loading ETFlow model...")
    model = BaseFlow.from_default(model="drugs-o3", cache=cache_dir)

    # Batch all SMILES in one ETFlow call for efficiency.
    all_smiles = [rec[2] for rec in records]
    print("Running ETFlow batch prediction with num_samples=1...")
    results = model.predict(all_smiles, num_samples=1, as_mol=False)

    # Process each modification independently so one failure does not kill the whole batch.
    for position, mod_code, smiles in records:
        print(f"Processing position={position}, code={mod_code}...")

        # ETFlow output dictionary should contain an entry keyed by original SMILES.
        # if smiles not in results or not results[smiles]:
        # 
        # 
        # 
        if smiles not in results or len(results[smiles]) == 0:
            print(f"  -> Warning: No ETFlow output for {mod_code} ({smiles}). Skipping.")
            continue

        # Use the first (and only requested) conformer sample.
        conformer_pos = results[smiles][0]

        # Build RDKit molecule from SMILES and attach generated 3D coordinates.
        mol_2d = get_mol_from_smiles(smiles)
        if mol_2d is None:
            print(f"  -> Warning: RDKit could not parse SMILES for {mod_code}. Skipping.")
            continue

        mol_3d = set_rdmol_positions(mol_2d, conformer_pos)

        # Apply RCSB/PDB naming with residue code; fallback to raw molecule if it fails.
        try:
            standardizer = PDBStandardizer(mod_code)
            final_mol = standardizer.standardize(mol_3d)
        except Exception as exc:
            print(
                f"  -> Warning: Standardization failed for {mod_code} ({exc}). "
                "Saving raw ETFlow output."
            )
            final_mol = mol_3d

        # Required output naming convention for downstream stitching.
        out_path = os.path.join(output_dir, f"mod_{position}_{mod_code}.pdb")
        Chem.MolToPDBFile(final_mol, out_path)
        print(f"  -> Saved: {out_path}")


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Generate standardized 3D side-chain PDBs from modifications.txt using ETFlow."
    )
    parser.add_argument(
        "--mods",
        required=True,
        help="Path to modifications text file (e.g., input/modifications.txt)",
    )
    parser.add_argument(
        "--cache",
        default="./cache",
        help="ETFlow cache directory (default: ./cache)",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Directory where side-chain PDBs will be written",
    )

    args = parser.parse_args()

    records = parse_modifications_file(args.mods)
    run_etflow_and_save_pdbs(records, output_dir=args.out_dir, cache_dir=args.cache)


if __name__ == "__main__":
    main()
