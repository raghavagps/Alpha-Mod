#!/usr/bin/env python3
"""
02_run_backbone.py

Command-line wrapper to predict peptide backbone structure from a canonical
amino acid FASTA sequence.

Behavior:
- Reads sequence from FASTA (ignores header lines starting with '>').
- Runs prediction using ESMFold API by default.
- Falls back to local ColabFold/AlphaFold2 if ESMFold fails.
- Saves final backbone to a user-specified output path (--out_pdb).
"""

import argparse
import glob
import os
import shutil
import time
from typing import List

import requests


def read_fasta_sequence(fasta_path: str) -> str:
    """
    Read a FASTA file and return the concatenated sequence string.

    Notes:
    - Header lines ('>...') are ignored.
    - Multiple sequence lines are concatenated.
    - Whitespace is stripped.
    """
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")

    seq_chunks: List[str] = []
    with open(fasta_path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq_chunks.append(line)

    sequence = "".join(seq_chunks).strip()
    if not sequence:
        raise ValueError(f"No sequence found in FASTA file: {fasta_path}")

    return sequence


def setup_colabfold_params() -> None:
    """
    Ensure ColabFold parameter files are available in local cache.

    The function tries two likely source locations for 'colabfold_params':
    1) Next to this script (project root use-case).
    2) One directory above this script (if script is inside a subfolder).
    """
    # Absolute path to the current script file.
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Candidate locations for the model parameter source folder.
    candidate_dirs = [
        os.path.join(script_dir, "colabfold_params"),
        os.path.join(os.path.dirname(script_dir), "colabfold_params"),
    ]

    source_dir = None
    for path in candidate_dirs:
        if os.path.exists(path):
            source_dir = path
            break

    if source_dir is None:
        raise FileNotFoundError(
            "[ERROR] Could not find 'colabfold_params' folder. "
            f"Checked: {candidate_dirs}"
        )

    # Standard ColabFold cache target used by local runs.
    target_dir = os.path.join(os.path.expanduser("~"), ".cache", "colabfold", "params")
    os.makedirs(target_dir, exist_ok=True)

    copied = False
    for fname in os.listdir(source_dir):
        src = os.path.join(source_dir, fname)
        dst = os.path.join(target_dir, fname)

        # Only copy files that are missing in cache.
        if os.path.isfile(src) and not os.path.exists(dst):
            shutil.copy(src, dst)
            copied = True

    if not copied:
        print("[PEPstrMOD2] All model files already installed.")
    print(f"[PEPstrMOD2] Parameters available in: {target_dir}")


def predict_structure(
    sequence: str,
    out_dir: str,
    jobname: str = "test_pred",
    method: str = "alphafold2",
    retries: int = 3,
    delay: int = 5,
) -> str:
    """
    Predict structure for a sequence.

    Parameters:
    - sequence: canonical amino acid sequence string.
    - out_dir: directory where intermediate/final outputs are stored.
    - jobname: base name for output PDB.
    - method: 'alphafold2' (default) or 'esmfold'.
    - retries/delay: ESMFold API retry behavior.

    Returns:
    - Absolute/relative path to final PDB file (out_dir/jobname.pdb).
    """
    os.makedirs(out_dir, exist_ok=True)
    pdb_path = os.path.join(out_dir, f"{jobname}.pdb")

    def run_esmfold() -> str:
        """Run ESMFold API request with retry logic."""
        headers = {"Content-Type": "application/x-www-form-urlencoded"}

        for attempt in range(1, retries + 1):
            print(f"\nSubmitting to ESMFold (attempt {attempt})...")
            try:
                response = requests.post(
                    "https://api.esmatlas.com/foldSequence/v1/pdb/",
                    headers=headers,
                    data=sequence,
                    timeout=300,
                )

                if response.status_code == 200:
                    with open(pdb_path, "w", encoding="utf-8") as handle:
                        handle.write(response.text)
                    print(f"\n[✓] ESMFold structure saved to {pdb_path}")
                    return pdb_path

                print(
                    f"[Attempt {attempt}] API error: {response.status_code} "
                    f"— {response.text.strip()}"
                )
            except requests.RequestException as exc:
                print(f"[Attempt {attempt}] Request failed: {exc}")

            if attempt < retries:
                print(f"Retrying in {delay} seconds...")
                time.sleep(delay)

        raise RuntimeError("ESMFold API failed after multiple attempts")

    def run_alphafold2() -> str:
        """
        Run local ColabFold/AlphaFold2 prediction.

        It first tries MSA mode, then falls back to single-sequence mode
        if the MSA run fails.
        """
        # Import here so ESMFold-only usage does not require ColabFold install.
        try:
            from colabfold.batch import run  # type: ignore
        except ImportError as exc:
            raise ImportError(
                "colabfold is not installed in this environment. "
                "Please run this script in your AlphaFold2/ColabFold conda environment."
            ) from exc

        setup_colabfold_params()
        queries = [(jobname, sequence, None, None)]

        try:
            print("\nRunning alphafold2 with MSA mode...")
            run(
                queries=queries,
                result_dir=out_dir,
                msa_mode="mmseqs2_uniref_env",
                num_models=5,
                is_complex=False,
                user_agent="PEPstrMOD/2.0",
            )
        except Exception as exc:
            print("\n[WARNING] MSA-based alphafold2 run failed! Error:", exc)
            print("\nRetrying in single-sequence mode...")
            run(
                queries=queries,
                result_dir=out_dir,
                msa_mode="single_sequence",
                model_type="alphafold2_ptm",
                num_models=5,
                is_complex=False,
                user_agent="PEPstrMOD/2.0",
            )

        # Prefer the top-ranked unrelaxed model, matching senior script behavior.
        pattern = os.path.join(out_dir, f"{jobname}_unrelaxed_rank_001*.pdb")
        found = glob.glob(pattern)

        if found:
            shutil.copy(found[0], pdb_path)
            print(f"[✓] Top alphafold2 model copied → {pdb_path}")
            return pdb_path

        # Final fallback: use any PDB found in result directory.
        candidates = [f for f in os.listdir(out_dir) if f.endswith(".pdb")]
        if candidates:
            shutil.copy(os.path.join(out_dir, candidates[0]), pdb_path)
            print(f"[✓] Fallback alphafold2 model used → {pdb_path}")
            return pdb_path

        raise RuntimeError("alphafold2 did not produce a PDB file.")

    # Method control flow: ESMFold first with fallback, or direct alphafold2.
    # method = method.lower().strip()
    # if method == "esmfold":
    #     try:
    #         return run_esmfold()
    #     except Exception as exc:
    #         print("\n[WARNING] ESMFold prediction failed! Error:", exc)
    #         print("\n→ Falling back to alphafold2...")
    #         return run_alphafold2()

    # if method == "alphafold2":
    #     return run_alphafold2()

    # raise ValueError("Invalid method. Choose either 'esmfold' or 'alphafold2'.")
    # Method control flow: alphafold2 first with fallback, or direct esmfold.
    method = method.lower().strip()
    if method == "alphafold2":
        try:
            return run_alphafold2()
        except Exception as exc:
            print("\n[WARNING] AlphaFold2 prediction failed! Error:", exc)
            print("\n→ Falling back to ESMFold...")
            return run_esmfold()

    if method == "esmfold":
        return run_esmfold()

    raise ValueError("Invalid method. Choose either 'esmfold' or 'alphafold2'.")


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Predict peptide backbone PDB from a FASTA sequence."
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Path to input FASTA file from step 01 (e.g., input/parsed_sequence.fasta)",
    )
    parser.add_argument(
        "--method",
        default="alphafold2",
        choices=["esmfold", "alphafold2"],
        help="Prediction method (default: alphafold2)",
    )
    parser.add_argument(
        "--out_pdb",
        required=True,
        help="Final output backbone PDB path",
    )

    args = parser.parse_args()

    # Read the canonical amino acid sequence from FASTA.
    sequence = read_fasta_sequence(args.fasta)

    # Derive ColabFold working directory and jobname from requested output path.
    # predict_structure writes to: <out_dir>/<jobname>.pdb
    out_pdb = os.path.abspath(args.out_pdb)
    out_dir = os.path.dirname(out_pdb) or "."
    jobname = os.path.splitext(os.path.basename(out_pdb))[0]

    if not jobname:
        raise ValueError("Invalid --out_pdb path: could not derive job name.")

    final_pdb = predict_structure(
        sequence=sequence,
        out_dir=out_dir,
        jobname=jobname,
        method=args.method,
    )

    print(f"\n[✓] Final backbone written to: {final_pdb}")


if __name__ == "__main__":
    main()
