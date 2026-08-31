# Alpha-Mod: Modified-Peptide Structure Pipeline

Alpha-Mod predicts three-dimensional structures of peptides containing canonical and non-canonical amino acids (NCAAs). It combines an AlphaFold2/ColabFold peptide backbone, ET-Flow conformers for modified residues, residue stitching, and MACE-OFF23 geometry minimization.

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/raghavagps/Alpha-Mod/blob/main/alpha_mod_pipeline_latest.ipynb)

## Pipeline overview

Alpha-Mod runs five sequential stages:

| Stage | Script | Purpose |
|---|---|---|
| 1. Parse | `01_parse_input.py` | Converts the input into a canonical parent sequence and writes `modifications.txt`. Known CCD residues and user-defined novel residues are routed separately. |
| 2. Backbone | `02_run_backbone.py` | Predicts the canonical parent backbone with AlphaFold2 through ColabFold. It can fall back to single-sequence AlphaFold2 and then ESMFold if the configured AlphaFold2 route fails. |
| 3. Side chains | `03_run_sidechains.py` | Generates one ET-Flow conformer for every modified residue and assigns PDB-compatible atom and residue names. Known residues use their RCSB CCD reference; novel residues use the declared canonical parent and a consensus/direct-backbone fallback. |
| 4. Stitch | `04_stitch.py` | Aligns each modified residue to the predicted parent residue and replaces it. Novel-residue terminal atoms are handled from molecular connectivity rather than names such as `H2`, `H3`, or `OXT`. |
| 5. Minimize | `05_minimize.py` | Validates the stitched residue identities, adds hydrogens with PDBFixer, and minimizes the structure with MACE-OFF23 for at most 2,000 L-BFGS steps. |

`main.py` orchestrates the five stages and reports whether minimization converged, reached the 2,000-step limit, or failed before a usable minimization result was produced.

## Repository files

```text
.
├── main.py
├── 01_parse_input.py
├── 02_run_backbone.py
├── 03_run_sidechains.py
├── 04_stitch.py
├── 05_minimize.py
├── modifications.json
├── alpha_mod_pipeline_latest.ipynb
└── README.md
```

Files whose names contain `_old_` are retained only as historical versions. The files without `_old_` are the active pipeline.

## Input formats

### 1. Canonical sequence

```text
ACDEFGHIK
```

### 2. Known PDB/CCD modification code

Put a known modification code in parentheses:

```text
KETAAAK(NVA)ERQH(NLE)DS
```

The code must be resolvable through `modifications.json`. Known residues remain on the direct RCSB CCD naming path.

### 3. MAP notation

MAP-style blocks supported by `modifications.json` may be mixed with the canonical sequence, for example:

```text
APGA{ptm:chloro}APG
```

See the [MAP specification paper](https://arxiv.org/abs/2505.03403) for the notation itself.

### 4. Explicit novel-residue format

For a residue that is not available through the current modification dictionary, provide its isomeric SMILES and canonical parent letter between pipe characters:

```text
ACDEK|N[C@@](C)(C)C(=O)O,A|FGHIK
```

The fields inside the pipe block are:

```text
SMILES,ParentLetter
```

Multiple novel residues can occur in one sequence:

```text
A|N[C@@](C)(C)C(=O)O,A|G|N[C@@H](CCO)C(=O)O,S|K|N[C@@H](CCC)C(=O)O,V|L
```

They receive unique three-character PDB residue names in sequence order:

```text
N_1, N_2, ... N_9, N10, N11, ... N99
```

The parent letter tells Alpha-Mod which canonical residue occupies that position during backbone prediction and which canonical CCD structure should be used as a naming/alignment reference. It does not convert an arbitrary molecule into an amino acid.

## Google Colab usage

Open `alpha_mod_pipeline_latest.ipynb` and use a GPU runtime. The notebook contains five cells:

1. **Smart Setup** — mounts Google Drive, clones or updates this repository, creates the three conda environments, restores/downloads model weights, verifies the required pipeline files, and warms the MACE model cache before parallel work.
2. **Single Sequence Prediction** — runs one sequence and downloads either the final PDB or all job files.
3. **Batch Prediction from FASTA** — runs two sequences in parallel, writes `batch_summary.csv`, preserves hard-failure diagnostics, and downloads a ZIP archive.
4. **Visualize Final Minimized Structures** — displays completed structures found on Google Drive.
5. **Recovery** — rebuilds a summary and ZIP from jobs already present on Drive after an interrupted batch session.

The notebook stores results under:

```text
MyDrive/try_pipeline/
├── converged/
├── max_steps/
├── failure_logs/
└── batch_summary.csv
```

Classification means:

- `converged`: MACE met the configured force threshold within 2,000 steps.
- `max_steps`: MACE completed 2,000 steps without reaching that threshold. Metrics and the final coordinates are still saved for diagnosis, but the structure must not be assumed to be chemically valid merely because a PDB was written.
- `failure_logs`: the pipeline stopped before producing a normal minimization status. `pipeline.log` and `failure_status.json` record the error.

The batch cell clears previous `converged`, `max_steps`, and `failure_logs` contents before starting. Download or move results that must be retained before running a new batch.

## FASTA batch input

Use standard multi-FASTA headers followed by one Alpha-Mod sequence per record:

```fasta
>known_nva_nle
KETAAAK(NVA)ERQH(NLE)DS
>novel_aib
ACDEK|N[C@@](C)(C)C(=O)O,A|FGHIK
```

FASTA record identifiers are used as job-folder names, so use short filesystem-safe identifiers.

## Output files

A completed job normally contains:

```text
parsed_sequence.fasta
modifications.txt
backbone.pdb
mod_<position>_<code>.pdb
stitched.pdb
protonated.pdb
final_minimized.pdb
minimize_status.json
```

`minimize_status.json` and `batch_summary.csv` report the convergence flag, step count, initial and final energies, energy change, elapsed time, RMSD values, and maximum atomic displacement when minimization completed normally.

## Important scope and limitations

- The explicit novel-residue pathway currently recognizes an **alpha-amino-acid** backbone with the connectivity `N-CA-C(=O)O`. A beta-amino acid such as `N-CA-CB-C(=O)O` changes the peptide backbone itself and is not currently supported by the side-chain-replacement strategy.
- Novel residues must contain one unambiguous alpha-amino-acid backbone. Atom identities that cannot be transferred safely from the parent receive deterministic element-based names such as `C1`, `C2`, `O1`, and `H1`.
- For an internal novel residue, the stitcher removes only a hydrogen actually bonded to the amino nitrogen. Heavy N-substituents such as a methyl group are preserved. If no removable N-bound hydrogen exists, the stitcher fails instead of deleting a heavy atom.
- Reaching `max_steps` is not equivalent to convergence. Inspect geometry, peptide-bond distances, energies, and displacement metrics before using such a structure.
- AlphaFold2 and ET-Flow predictions, proton placement, and MACE optimization can produce small run-to-run differences. A peptide may converge in a different number of steps between otherwise equivalent runs.

## Advanced local setup

The supplied notebook is the maintained execution path. The scripts use three isolated environments because their dependencies can conflict:

```bash
conda create -n af2_env python=3.10 -y
conda run -n af2_env pip install "colabfold[alphafold] @ git+https://github.com/sokrypton/ColabFold"

conda create -n etflow_env python=3.10 -y
conda run -n etflow_env pip install requests rdkit biopython torch etflow

conda create -n mace_env python=3.10 -y
conda run -n mace_env pip install mace-torch ase openmm
conda run -n mace_env pip install git+https://github.com/openmm/pdbfixer.git
```

The current delivery paths in `main.py` and the notebook are Colab/Google-Drive oriented. A local deployment may need its output-delivery paths adjusted.

## Citation

If you use Alpha-Mod, cite the accompanying Alpha-Mod manuscript when its final citation information is available.

## License

Add the project license and redistribution terms in a repository `LICENSE` file before publishing a formal release.
