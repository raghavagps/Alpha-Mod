#!/usr/bin/env python3
"""
05_minimize.py (MACE-OFF Version)

Takes the stitched.pdb file, adds missing hydrogens using OpenMM (safely 
handling NCAAs), and relaxes the structure using the MACE-OFF machine 
learning force field via ASE.
"""

import argparse
import json
import os
import sys
import time
import tempfile
import numpy as np
import torch

from ase.io import read, write
from ase.optimize import LBFGS
from mace.calculators import mace_off
from pdbfixer import PDBFixer
from openmm.app import PDBFile

STANDARD_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
    'THR', 'TRP', 'TYR', 'VAL',
    'HID', 'HIE', 'HIP', 'CYX', 'ASH', 'GLH',
    'ACE', 'NME', 'NMA',
    'HOH', 'WAT', 'NA', 'CL', 'K', 'MG', 'CA', 'ZN',
}

def read_expected_modifications(mods_path):
    """Return {position: expected_three_letter_PDB_name} from modifications.txt."""
    if not os.path.exists(mods_path):
        raise FileNotFoundError(f"Modifications file not found: {mods_path}")

    expected = {}
    with open(mods_path, "r", encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue

            # Split only twice: the rest may contain SMILES and an optional parent field.
            parts = [part.strip() for part in line.split(":", 2)]
            if len(parts) != 3 or not parts[0].isdigit() or not parts[1]:
                raise ValueError(
                    f"Invalid modifications.txt line {line_no}: {line!r}"
                )

            position = int(parts[0])
            if position < 1:
                raise ValueError(
                    f"Invalid residue position on line {line_no}: {position}"
                )

            if position in expected:
                raise ValueError(
                    f"Multiple modifications at residue {position} are not supported "
                    "by this single-residue stitching pipeline."
                )

            # PDB residue-name field is exactly three characters wide.
            expected[position] = parts[1].upper()[:3]

    return expected


def validate_stitched_residue_names(stitched_pdb, mods_path):
    """Ensure every requested modification is present before MACE starts."""
    expected = read_expected_modifications(mods_path)

    if not expected:
        print("[Validate] No modifications declared; residue-name validation passed.")
        return

    observed = {}
    with open(stitched_pdb, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            if not raw_line.startswith(("ATOM  ", "HETATM")):
                continue

            try:
                position = int(raw_line[22:26])
            except ValueError:
                continue

            if position not in expected:
                continue

            key = (raw_line[21], position)
            entry = observed.setdefault(key, {"names": set(), "atoms": set()})
            entry["names"].add(raw_line[17:20].strip())
            entry["atoms"].add(raw_line[12:16].strip())

    errors = []

    for position, expected_name in expected.items():
        matches = [
            (chain_id, data)
            for (chain_id, residue_number), data in observed.items()
            if residue_number == position
        ]

        if len(matches) == 0:
            errors.append(
                f"position {position}: expected {expected_name}, "
                "but no residue was found"
            )
            continue

        if len(matches) > 1:
            chains = [chain_id or "<blank>" for chain_id, _ in matches]
            errors.append(
                f"position {position}: found residues in multiple chains {chains!r}"
            )
            continue

        chain_id, data = matches[0]
        actual_names = data["names"]

        if actual_names != {expected_name}:
            errors.append(
                f"position {position}, chain {chain_id or '<blank>'}: "
                f"expected {expected_name}, found {sorted(actual_names)!r}"
            )

        missing_anchors = {"N", "CA", "C"} - data["atoms"]
        if missing_anchors:
            errors.append(
                f"position {position}, chain {chain_id or '<blank>'}: "
                f"missing backbone atom(s) {sorted(missing_anchors)!r}"
            )

    if errors:
        details = "\n  - ".join(errors)
        raise ValueError(
            "Pre-minimization validation failed; MACE will not run:\n  - "
            + details
        )

    summary = ", ".join(
        f"{name}{position}" for position, name in sorted(expected.items())
    )
    print(f"[Validate] Confirmed stitched modified residues: {summary}")


# ==============================================================================
# 1. CORE FUNCTIONS (From your MACE Notebook)
# ==============================================================================

def preprocess(pdb_path, output_dir):
    def _is_h(line):
        atom_name = line[12:16].strip()
        element   = line[76:78].strip() if len(line) > 76 else ''
        if element == 'H': return True
        if not element:
            return atom_name.startswith('H') or (len(atom_name) > 1 and atom_name[0].isdigit() and atom_name[1] == 'H')
        return False

    def _res_key(line):
        return (line[17:20].strip(), line[21], line[22:26].strip())

    # --- STEP 1: save NCAA H from init, write H-stripped init to temp file ---
    ncaa_h_saved = {}   # res_key -> list of original ATOM lines for H atoms
    stripped_lines = []

    with open(str(pdb_path)) as fh:
        for raw in fh:
            if raw.startswith(('ATOM', 'HETATM')):
                res_name = raw[17:20].strip()
                key = _res_key(raw)
                if res_name not in STANDARD_RESIDUES and _is_h(raw):
                    ncaa_h_saved.setdefault(key, []).append(raw)
                else:
                    stripped_lines.append(raw)
            else:
                stripped_lines.append(raw)

    tmp_fd, tmp_path = tempfile.mkstemp(suffix='_stripped.pdb')
    os.close(tmp_fd)
    with open(tmp_path, 'w') as fh:
        fh.writelines(stripped_lines)

    # --- STEP 2: run PDBFixer exactly as before, on the stripped file ---
    fixer = PDBFixer(filename=tmp_path)
    os.unlink(tmp_path)

    n_initial = sum(1 for _ in fixer.topology.atoms()) + sum(len(v) for v in ncaa_h_saved.values())

    ncaa_list = []
    for res in fixer.topology.residues():
        if res.name not in STANDARD_RESIDUES:
            key = (res.name, '', str(res.id))
            n_h = len(ncaa_h_saved.get(key, []))
            ncaa_list.append((res.name, res.index + 1, n_h))

    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    tmp_fd2, tmp_proto = tempfile.mkstemp(suffix='_pdbfixer_out.pdb')
    os.close(tmp_fd2)
    with open(tmp_proto, 'w') as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh)

    # --- STEP 3: in PDBFixer's output, remove NCAA H and paste back original H ---
    pdbfixer_lines = open(tmp_proto).readlines()
    os.unlink(tmp_proto)

    last_idx = {}
    for i, raw in enumerate(pdbfixer_lines):
        if raw.startswith(('ATOM', 'HETATM')):
            res_name = raw[17:20].strip()
            if res_name not in STANDARD_RESIDUES:
                last_idx[_res_key(raw)] = i

    out_lines = []
    serial = 1
    inserted = set()
    for i, raw in enumerate(pdbfixer_lines):
        if raw.startswith(('ATOM', 'HETATM')):
            res_name = raw[17:20].strip()
            key = _res_key(raw)
            if res_name not in STANDARD_RESIDUES and _is_h(raw):
                continue
            out_lines.append(f"{raw[:6]}{serial:5d}{raw[11:]}")
            serial += 1
            if res_name not in STANDARD_RESIDUES and i == last_idx.get(key) and key not in inserted:
                for h_raw in ncaa_h_saved.get(key, []):
                    out_lines.append(f"{h_raw[:6]}{serial:5d}{h_raw[11:]}")
                    serial += 1
                inserted.add(key)
        else:
            out_lines.append(raw)

    protonated_path = os.path.join(output_dir, 'protonated.pdb')
    with open(protonated_path, 'w') as fh:
        fh.writelines(out_lines)

    final_atom_lines = [l for l in out_lines if l.startswith(('ATOM', 'HETATM'))]
    n_final   = len(final_atom_lines)
    elements  = sorted(set(l[76:78].strip() for l in final_atom_lines if len(l) > 76 and l[76:78].strip()))

    return protonated_path, n_initial, n_final, ncaa_list, elements

class ExplosionError(Exception):
    pass

class ExplosionDetector:
    def __init__(self, atoms, initial_energy, initial_positions, check_every=50):
        self.atoms = atoms
        self.e0 = initial_energy
        self.pos0 = initial_positions.copy()
        self.check_every = check_every
        self.step_count = 0

    def __call__(self):
        self.step_count += 1
        if self.step_count % self.check_every != 0:
            return

        e_now = self.atoms.get_potential_energy()
        delta_e = e_now - self.e0

        if delta_e < -1000:
            raise ExplosionError(f"Energy dropped {delta_e:.0f} eV. Structure is exploding.")
        if e_now > 0:
            raise ExplosionError(f"Energy is positive ({e_now:.0f} eV). Atoms are overlapping catastrophically.")

        pos_now = self.atoms.get_positions()
        max_disp = np.max(np.linalg.norm(pos_now - self.pos0, axis=1))
        if max_disp > 10.0:
            raise ExplosionError(f"Max atomic displacement is {max_disp:.1f} Å. Structure is disintegrating.")

def minimize(protonated_pdb, output_dir, model_size, fmax, max_steps, device):
    atoms = read(protonated_pdb)
    n_atoms = len(atoms)

    calc = mace_off(model=model_size, device=device, default_dtype="float64")
    atoms.calc = calc

    e0 = atoms.get_potential_energy()
    f0_max = np.max(np.linalg.norm(atoms.get_forces(), axis=1))
    pos0 = atoms.get_positions().copy()

    if e0 > 0:
        print(f"    ⚠ Initial energy is POSITIVE ({e0:.0f} eV) — input structure has severe overlaps")

    logfile = os.path.join(output_dir, 'opt.log')
    trajectory = os.path.join(output_dir, 'opt.traj')
    ase_raw = os.path.join(output_dir, '_ase_raw.pdb')

    opt = LBFGS(atoms, logfile=logfile, trajectory=trajectory)
    # detector = ExplosionDetector(atoms, e0, pos0, check_every=50)
    # opt.attach(detector)

    t0 = time.time()
    exploded = False
    explosion_msg = ""

    try:
        converged = opt.run(fmax=fmax, steps=max_steps)
    except ExplosionError as ex:
        exploded = True
        explosion_msg = str(ex)
        converged = False

    elapsed = time.time() - t0
    e1 = atoms.get_potential_energy()
    f1_max = np.max(np.linalg.norm(atoms.get_forces(), axis=1))
    write(ase_raw, atoms)

    return {
        'n_atoms': n_atoms,
        'e0': e0, 'e1': e1, 'delta_e': e1 - e0,
        'f0_max': f0_max, 'f1_max': f1_max,
        'steps': opt.nsteps, 'converged': converged,
        'exploded': exploded, 'explosion_msg': explosion_msg,
        'time_sec': elapsed,
        'time_per_step': elapsed / max(opt.nsteps, 1),
        'ase_raw_path': ase_raw,
    }

def reformat_pdb(protonated_pdb, ase_raw_pdb, output_pdb):
    prot_lines = []
    with open(protonated_pdb) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                prot_lines.append(line)
    min_coords = []
    with open(ase_raw_pdb) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                min_coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    if len(prot_lines) != len(min_coords):
        import shutil
        shutil.copy(ase_raw_pdb, output_pdb)
        return False
    with open(output_pdb, 'w') as f:
        for i in range(len(prot_lines)):
            line = prot_lines[i]
            x, y, z = min_coords[i]
            f.write(line[:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + line[54:])
        f.write("TER\n")
        f.write("END\n")
    return True

def validate(protonated_pdb, minimized_pdb):
    """
    Compute structural drift between pre-minimization and post-minimization.
    Returns RMSD for all atoms, heavy atoms only, and max single-atom displacement.
    """
    before = read(protonated_pdb)
    after = read(minimized_pdb)
    pos_b = before.get_positions()
    pos_a = after.get_positions()

    # All-atom RMSD
    rmsd_all = float(np.sqrt(np.mean(np.sum((pos_b - pos_a)**2, axis=1))))

    # Heavy-atom RMSD (exclude hydrogen)
    heavy = np.array([s != 'H' for s in after.get_chemical_symbols()])
    rmsd_heavy = float(np.sqrt(np.mean(np.sum((pos_b[heavy] - pos_a[heavy])**2, axis=1))))

    # Max displacement of any single atom
    max_disp = float(np.max(np.linalg.norm(pos_b - pos_a, axis=1)))

    return {'rmsd_all': rmsd_all, 'rmsd_heavy': rmsd_heavy, 'max_disp': max_disp}

# ==============================================================================
# 2. THE PIPELINE ORCHESTRATOR 
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(description="Minimize stitched PDB using MACE-OFF.")
    parser.add_argument("--input", required=True, help="Path to input stitched.pdb")
    parser.add_argument("--mods", required=True, help="Path to modifications.txt")
    parser.add_argument("--out", required=True, help="Path to output final_minimized.pdb")
    args = parser.parse_args()

    input_pdb = os.path.abspath(args.input)
    mods_path = os.path.abspath(args.mods)
    final_pdb = os.path.abspath(args.out)
    work_dir = os.path.dirname(input_pdb)

    if not os.path.exists(input_pdb):
        print(f"[ERROR] Minimizer input not found: {input_pdb}")
        sys.exit(1)

    if not os.path.exists(mods_path):
        print(f"[ERROR] Modifications file not found: {mods_path}")
        sys.exit(1)

    print(f"[MACE] Starting Machine Learning Minimization...")

    # Determine Device
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[MACE] PyTorch Device: {device.upper()}")

    try:
        print("[Validate] Checking stitched residue identities before protonation...")
        validate_stitched_residue_names(input_pdb, mods_path)

        # Step 1: Preprocess
        print("[MACE] Preprocessing and protonating structure...")
        protonated_path, n_initial, n_final, ncaa_list, elements = preprocess(input_pdb, work_dir)

        # Step 2: Minimize
        print(f"[MACE] Running MACE-OFF (medium model, float64)...")
        stats = minimize(
            protonated_pdb=protonated_path, 
            output_dir=work_dir, 
            model_size="medium", 
            fmax=0.1, 
            max_steps=2000, 
            device=device
        )

        # Step 3: Reformat (runs whether it exploded or converged)
        if stats['exploded']:
            print(f"[MACE] EXPLOSION DETECTED: {stats['explosion_msg']}")
            print(f"[MACE] Saving exploded coordinates to {final_pdb} for inspection...")
        else:
            status_word = "Converged" if stats['converged'] else "Max steps reached"
            print(f"[MACE] {status_word} in {stats['steps']} steps. ΔE = {stats['delta_e']:.2f} eV")

        reformat_pdb(protonated_path, stats['ase_raw_path'], final_pdb)
        print(f"[MACE] Final PDB saved to: {final_pdb}")

        # Step 4: Validate (compute structural drift)
        print("[MACE] Computing structural drift...")
        val = validate(protonated_path, final_pdb)
        print(f"[MACE] RMSD(heavy): {val['rmsd_heavy']:.4f} Å, Max disp: {val['max_disp']:.4f} Å")

        # Step 5: Write status JSON so main.py knows exactly what happened
        status = {
            'converged': bool(stats['converged']),
            'exploded': bool(stats['exploded']),
            'explosion_msg': stats.get('explosion_msg', ''),
            'n_atoms': stats['n_atoms'],
            'steps': stats['steps'],
            'e0': float(stats['e0']),
            'e1': float(stats['e1']),
            'delta_e': float(stats['delta_e']),
            'time_sec': float(stats['time_sec']),
            'time_per_step': float(stats['time_per_step']),
            'rmsd_all': val['rmsd_all'],
            'rmsd_heavy': val['rmsd_heavy'],
            'max_disp': val['max_disp'],
        }
        status_path = os.path.join(work_dir, 'minimize_status.json')
        with open(status_path, 'w') as f:
            json.dump(status, f, indent=2)
        print(f"[MACE] Status written to: {status_path}")

        # Exit with code 1 if exploded (so main.py knows via exit code too)
        if stats['exploded']:
            sys.exit(1)

    except Exception as e:
        print(f"[ERROR] MACE minimization failed: {e}")
        sys.exit(1)
        
    finally:
        # Clean up only the intermediate files (protonated.pdb is kept for reference)
        temp_files = ["_ase_raw.pdb", "opt.traj"]
        for f in temp_files:
            path = os.path.join(work_dir, f)
            if os.path.exists(path):
                os.remove(path)

if __name__ == "__main__":
    main()
