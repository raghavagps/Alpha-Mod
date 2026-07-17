# Alpha-Mod: Batch Peptide Pipeline

**Alpha-Mod** is a modular, hybrid computational pipeline for predicting the 3D structures of peptides containing **non-canonical amino acids (NCAAs)** — chemically modified residues that fall outside the 20 standard amino acids.

Most structure predictors either ignore modified residues entirely or only support a fixed vocabulary of them. Alpha-Mod instead predicts the peptide's canonical backbone with **AlphaFold2**, generates a 3D conformer for *any* NCAA directly from its **SMILES** string using **ET-Flow**, stitches the two together with a Kabsch/BioPython superimposition, and relaxes the resulting structure with the **MACE-OFF23** machine-learning force field. Because NCAAs are modeled independently as small molecules rather than via residue-specific force-field parameters, Alpha-Mod can in principle handle **any chemically modified residue expressible as a SMILES string**, on standard (non-HPC) hardware.

Alpha-Mod was benchmarked against AlphaFold3, PEPstrMOD, and PEPstrMOD2 across three peptide datasets (ModPep_257, ModPep_16, PEP_SOLO), and is competitive with — and in several metrics outperforms — these existing tools. Full methodology, benchmarking, and results are described in the accompanying manuscript (preprint link to be added).

## How it works

Alpha-Mod runs each peptide through five sequential steps:

| Step | Script | What it does |
|---|---|---|
| 1. Parse | `01_parse_input.py` | Parses the input sequence, resolves each modification against `modifications.json`, and writes a canonical-sequence FASTA plus a `modifications.txt` map (position : code : SMILES : parent residue). |
| 2. Backbone | `02_run_backbone.py` | Predicts the canonical backbone with AlphaFold2 (via ColabFold, MMseqs2 MSA). Falls back to single-sequence mode if MSA generation fails, and to the ESMFold API if AlphaFold2 fails outright. |
| 3. Side chains | `03_run_sidechains.py` | Generates a 3D conformer for each NCAA from its SMILES string using ET-Flow, then standardizes atom naming against the RCSB Chemical Component Dictionary (for known codes) or against the parent canonical residue (for novel/user-defined NCAAs) using maximum-common-substructure matching. |
| 4. Stitch | `04_stitch.py` | Superimposes each NCAA conformer onto the backbone at the corresponding position using shared anchor atoms (N, Cα, C, Cβ) and swaps it in, with position-aware handling of terminal atoms (N- vs C- vs internal residue). |
| 5. Minimize | `05_minimize.py` | Adds missing hydrogens with PDBFixer (NCAA-aware), then relaxes the stitched structure with the MACE-OFF23 ML force field (ASE + L-BFGS) to remove steric clashes introduced by stitching, with an explosion detector to catch and report failed minimizations. |

`main.py` orchestrates all five steps for one or many sequences, isolating each job in its own output folder, running each conda-environment-specific step as a subprocess, and — in the Colab setting — sorting completed jobs into `converged/` or `exploded/` folders on Google Drive and writing a `batch_summary.csv`.

## Repository structure

```
.
├── main.py                  # Orchestrator: runs the 5-step pipeline per sequence
├── 01_parse_input.py         # Step 1 — sequence parsing & modification resolution
├── 02_run_backbone.py         # Step 2 — AlphaFold2 / ESMFold backbone prediction (af2_env)
├── 03_run_sidechains.py       # Step 3 — ET-Flow NCAA conformer generation (etflow_env)
├── 04_stitch.py               # Step 4 — Kabsch-based stitching of side chains into backbone
├── 05_minimize.py             # Step 5 — MACE-OFF23 energy minimization (mace_env)
├── modifications.json        # Lookup table of ~1,650 known NCAA modifications (CCD code, SMILES, parent AA, etc.)
└── alpha_mod_pipeline.ipynb  # Google Colab notebook — one-click setup + single/batch prediction UI
```

## Input formats

A peptide sequence can mix canonical amino acids with modified residues in three ways:

1. **Legacy PDB three-letter code**, e.g. `APG(5PG)APG` — the code in parentheses is looked up in `modifications.json`.
2. **MAP format**, e.g. `APGA{ptm:chloro}APG`, or terminal tags `{nt:...}` / `{ct:...}` for N-/C-terminal modifications. To learn more about the MAP format itself, refer to [this paper](https://arxiv.org/abs/2505.03403).
3. **Novel/explicit format** for residues not in `modifications.json`, using pipe delimiters: `APG|SMILES,ParentLetter|APG` — e.g. `APG|CC(C)C(N)C(O)=O,A|APG` defines a custom NCAA by its own SMILES string and parent canonical residue.

## Getting started

### Option A — Google Colab (recommended)

Open `alpha_mod_pipeline.ipynb` in Google Colab. It contains four cells:

1. **Smart Setup** — installs three isolated conda environments (`etflow_env`, `af2_env`, `mace_env`), mounts Google Drive, clones this repository, and downloads/caches ET-Flow and AlphaFold2 model weights to Drive so subsequent sessions skip the download. Run once per session.
2. **Single Sequence Prediction** — enter a sequence and job name, run the full pipeline, and auto-download the result.
3. **Batch Prediction from FASTA** — upload a multi-sequence FASTA file and run two jobs in parallel; outputs and a summary CSV are auto-downloaded as a zip.
4. **Recovery Cell** — if a batch run crashes partway through, this rescans Google Drive for completed jobs and rebuilds the summary CSV without re-running anything.

No local installation or GPU is required — everything runs in the Colab environment.

### Option B — Local installation

The pipeline expects three separate conda environments (mirroring the Colab setup), since AlphaFold2/ColabFold, ET-Flow, and MACE-OFF have overlapping/conflicting dependencies:

```bash
# af2_env — AlphaFold2 backbone prediction
conda create -n af2_env python=3.10 -y
conda run -n af2_env pip install "colabfold[alphafold] @ git+https://github.com/sokrypton/ColabFold"

# etflow_env — NCAA conformer generation
conda create -n etflow_env python=3.10 -y
conda run -n etflow_env pip install torch etflow rdkit biopython requests

# mace_env — energy minimization
conda create -n mace_env python=3.10 -y
conda run -n mace_env pip install mace-torch ase openmm
conda run -n mace_env pip install git+https://github.com/openmm/pdbfixer.git
```

Then run a single sequence:

```bash
python main.py --sequence "APG(5PG)APG" --json modifications.json --job_name my_peptide
```

or a batch from a multi-FASTA file:

```bash
python main.py --fasta_in sequences.fasta --json modifications.json
```

`main.py` writes results to `output/{job_name}/`, with the final relaxed structure at `output/{job_name}/final_minimized.pdb` and a `minimize_status.json` describing convergence, energy change, and structural drift from minimization.

## Benchmarking

Alpha-Mod was evaluated against AlphaFold3, PEPstrMOD, and PEPstrMOD2 on three curated datasets (ModPep_257, ModPep_16, PEP_SOLO) using Cα-RMSD, backbone RMSD, all-atom RMSD, a local NCAA-specific RMSD, and DSSP-based secondary structure recovery (Q3/Q8). On the primary ModPep_257 benchmark, Alpha-Mod's all-atom RMSD (4.30 Å) outperformed AlphaFold3 (4.89 Å) and PEPstrMOD (4.96 Å), while trailing PEPstrMOD2 (3.94 Å). Full results, including per-dataset breakdowns and the effect of MACE-OFF23 minimization, are reported in the manuscript.

## Citation

If you use Alpha-Mod, please cite the accompanying preprint (citation details to be added once available).

## License

Add your chosen license here (e.g. MIT, Apache-2.0).
