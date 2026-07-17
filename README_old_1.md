# Modified Peptide Structure Pipeline

This pipeline predicts the 3D structure of modified peptides.

In simple terms, it does the following:
1. Reads your peptide sequence
2. Converts the sequence into a clean canonical FASTA form
3. Predicts the peptide backbone
4. Builds the modified side chains
5. Stitches the backbone and side chains together
6. Minimizes the final structure

The backbone is predicted using AlphaFold2 first. If that does not work, the pipeline falls back to ESMFold. The modified side chains are generated using the ETFlow model.

For a finished run, the main final structure file is usually:

```text
final_minimized.pdb
```

## Supported Input Styles

This pipeline supports:
- normal amino acid sequences
- legacy format such as `(SEP)` or `(5PG)`
- MAP format such as `{ptm:Code}`, `{nnr:Code}`, `{nt:Code}`, `{ct:Code}`

If you want to understand the MAP format in more detail, you can read the paper here:

[MAP format paper](https://arxiv.org/abs/2505.03403)

Examples:

```text
AKTA
APG(5PG)APG
MKA{ptm:6CV}A
MK{nnr:0FL}TA
{nt:DiMe1}AKTA
MKTA{ct:0NC}
```

## Main Ways To Run The Pipeline

### Run one sequence directly

```bash
python main.py --sequence "MKA{ptm:6CV}A" --json modifications.json --job_name my_test
```

### Run many sequences from a FASTA file

```bash
python main.py --fasta_in batch_input.fasta --json modifications.json
```

## Output

Each run gets its own folder inside:

```text
output/<job_name>/
```

Important files in that folder are usually:
- `parsed_sequence.fasta`
- `modifications.txt`
- `backbone.pdb`
- `stitched.pdb`
- `final_minimized.pdb`

If you want the final result for a sequence, usually look at:

```text
final_minimized.pdb
```

## How To Use `final_batch.ipynb`

The notebook `final_batch.ipynb` is the easiest way to try the pipeline, especially if you are not very technical.

### What the notebook does

It helps you:
- install the required software
- connect Google Drive
- set up the environments
- run either a single sequence or a batch FASTA file

### Simple notebook workflow

1. Open `final_batch.ipynb`
2. Run the setup cells from top to bottom
3. Wait for the installations to finish
4. Choose either the single-sequence run cell or the batch FASTA run cell
5. Run the cell you want
6. Wait for the pipeline to finish
7. Check the results in the output folder linked to your Google Drive

By default, the notebook is set to save results in a folder named `try_pipeline` in your Google Drive. You can change that folder name if you wish.

## Single Sequence In `final_batch.ipynb`

For single-sequence prediction, the notebook has placeholders where you can enter:
- `sequence`
- `job_name`

As shown in the notebook cell, you simply type your sequence in the sequence box and type the folder name you want in the job name box, then run the cell.

You can enter the sequence in:
- MAP format
- or the older legacy format like `MKA(SEP)A`

Examples:

```text
MKA{ptm:6CV}A
APG(5PG)APG
```

The `job_name` you enter will be used as the output folder name inside the `try_pipeline` folder in Google Drive.

## Batch FASTA In `final_batch.ipynb`

For batch prediction, the notebook has a separate cell for FASTA upload.

When you run cell 4, it will ask you to upload a `.fasta` file. After you upload the file, it will do the rest of the work automatically.

So for batch mode, you mainly need to:
- run the batch cell
- upload the FASTA file when asked
- wait for the pipeline to finish

Each FASTA entry will be processed separately.

## Good First-Time Advice

If this is your first time using the pipeline:
1. Start with `final_batch.ipynb`
2. First try one small single sequence
3. After that, try your real sequence
4. Then try batch FASTA mode if needed

This is usually the easiest order.

## Summary

Use `final_batch.ipynb` if you want the easiest guided workflow.

Use `main.py` directly if you want command-line control.

Single sequence:

```bash
python main.py --sequence "MKA{ptm:6CV}A" --json modifications.json --job_name my_test
```

Batch FASTA:

```bash
python main.py --fasta_in batch_input.fasta --json modifications.json
```
