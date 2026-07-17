#!/usr/bin/env python3
"""
01_parse_input.py

Adapter parser for peptide sequences containing:
- Legacy modification blocks in parentheses, e.g. APG(5PG)APG
- MAP-style blocks in curly braces, e.g. MKT{ptm:DiMe1}A

Outputs remain unchanged:
1) FASTA with the translated canonical amino acid sequence.
2) TXT with one modification per line in this exact format:
   position : code : SMILES
"""

import argparse
import json
import os
from typing import Dict, List, Tuple


def load_modification_index(json_path: str) -> Dict[str, dict]:
    """Load JSON and index entries by both 'Three letter code' and 'User Code'."""
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    if not isinstance(data, list):
        raise ValueError("JSON must be a list of modification records.")

    index: Dict[str, dict] = {}
    for entry in data:
        if not isinstance(entry, dict):
            continue
        three_letter_code = str(entry.get("Three letter code", "")).strip()
        user_code = str(entry.get("User Code", "")).strip()

        if three_letter_code:
            index[three_letter_code] = entry
        if user_code:
            index[user_code] = entry
    return index


def get_modification_record(mod_code: str, mod_index: Dict[str, dict]) -> dict:
    """Fetch one modification record by either three-letter code or user code."""
    if mod_code not in mod_index:
        raise ValueError(f"Modification code '{mod_code}' not found in JSON.")
    return mod_index[mod_code]


def get_three_letter_code(record: dict, input_code: str) -> str:
    """Extract the official PDB three-letter code for downstream logging."""
    three_letter_code = str(record.get("Three letter code", "")).strip()
    if not three_letter_code:
        raise ValueError(
            f"Missing 'Three letter code' in JSON record for modification '{input_code}'."
        )
    return three_letter_code


def get_smiles(record: dict, mod_code: str) -> str:
    """Extract SMILES from one JSON record."""
    smiles = str(record.get("SMILES", "")).strip()
    if not smiles:
        raise ValueError(f"Missing SMILES for modification code '{mod_code}'.")
    return smiles


def extract_parent_one_letter(natural_aa_field: str) -> str:
    """
    Extract parent one-letter code from strings like 'Alanine/Ala/A'.
    We use the last slash-delimited token as the canonical 1-letter residue.
    """
    parts = [p.strip() for p in str(natural_aa_field).split("/") if p.strip()]
    if not parts:
        raise ValueError(f"Invalid 'Natural Amino Acid' field: {natural_aa_field!r}")
    one_letter = parts[-1]
    if len(one_letter) != 1 or not one_letter.isalpha():
        raise ValueError(f"Could not extract 1-letter code from: {natural_aa_field!r}")
    return one_letter.upper()


def parse_sequence(
    sequence: str, mod_index: Dict[str, dict]
) -> Tuple[str, List[Tuple[int, str, str]]]:
    """
    Parse input sequence while maintaining 1-based residue position.

    Rules:
    - Alphabetic characters are treated as canonical residues.
    - Parenthesized text '(XXX)' is legacy format and behaves like a standalone
      non-natural residue: map to the parent canonical residue and insert it.
    - Curly-brace MAP tags support:
      * {ptm:Code}: modify the immediately preceding residue
      * {nnr:Code}: standalone residue, insert parent canonical residue
      * {nt:Code}: N-terminal modification, log at position 1
      * {ct:Code}: C-terminal modification, log at final canonical position
      * {Code}: fallback shorthand, treated like standalone non-natural residue

    Returns:
    - canonical sequence string
    - list of tuples: (position, modification_code, smiles)
    """
    caa_chars: List[str] = []
    modification_specs: List[Tuple[object, str, str]] = []

    residue_position = 1  # 1-based position in final translated sequence
    i = 0
    n = len(sequence)

    while i < n:
        ch = sequence[i]

        # Canonical residue case: single alphabetic character.
        if ch.isalpha():
            caa_chars.append(ch.upper())
            residue_position += 1
            i += 1
            continue

        # Legacy modification case: text inside parentheses behaves like a
        # standalone modified residue, so we insert the parent canonical letter.
        if ch == "(":
            end_idx = sequence.find(")", i + 1)
            if end_idx == -1:
                raise ValueError(f"Unclosed parenthesis starting at index {i}.")

            mod_code = sequence[i + 1 : end_idx].strip()
            if not mod_code:
                raise ValueError(f"Empty modification code at index {i}.")

            record = get_modification_record(mod_code, mod_index)
            official_code = get_three_letter_code(record, mod_code)
            parent_one_letter = extract_parent_one_letter(record.get("Natural Amino Acid", ""))
            smiles = get_smiles(record, mod_code)

            # Insert mapped parent residue into canonical sequence.
            caa_chars.append(parent_one_letter)

            # Record position for downstream side-chain generation/stitching.
            modification_specs.append((residue_position, official_code, smiles))

            residue_position += 1
            i = end_idx + 1
            continue

        # MAP-format case: text inside curly braces.
        if ch == "{":
            end_idx = sequence.find("}", i + 1)
            if end_idx == -1:
                raise ValueError(f"Unclosed curly brace starting at index {i}.")

            block_text = sequence[i + 1 : end_idx].strip()
            if not block_text:
                raise ValueError(f"Empty MAP block at index {i}.")

            # Fallback shorthand: {CODE} behaves like a standalone modified residue.
            if ":" not in block_text:
                mod_code = block_text
                record = get_modification_record(mod_code, mod_index)
                official_code = get_three_letter_code(record, mod_code)
                parent_one_letter = extract_parent_one_letter(record.get("Natural Amino Acid", ""))
                smiles = get_smiles(record, mod_code)

                caa_chars.append(parent_one_letter)
                modification_specs.append((residue_position, official_code, smiles))

                residue_position += 1
                i = end_idx + 1
                continue

            prefix, mod_code = [part.strip() for part in block_text.split(":", 1)]
            if not prefix or not mod_code:
                raise ValueError(f"Invalid MAP block at index {i}: {block_text!r}")

            prefix = prefix.lower()
            record = get_modification_record(mod_code, mod_index)
            official_code = get_three_letter_code(record, mod_code)
            smiles = get_smiles(record, mod_code)

            if prefix == "ptm":
                # PTM modifies the immediately preceding canonical residue.
                if residue_position == 1:
                    raise ValueError(f"PTM block '{block_text}' has no preceding residue to modify.")
                modification_specs.append((residue_position - 1, official_code, smiles))
                i = end_idx + 1
                continue

            if prefix == "nnr":
                # Non-natural residue is a standalone residue in the sequence.
                parent_one_letter = extract_parent_one_letter(record.get("Natural Amino Acid", ""))
                caa_chars.append(parent_one_letter)
                modification_specs.append((residue_position, official_code, smiles))
                residue_position += 1
                i = end_idx + 1
                continue

            if prefix == "nt":
                # N-terminal modifications are always logged at position 1.
                modification_specs.append((1, official_code, smiles))
                i = end_idx + 1
                continue

            if prefix == "ct":
                # C-terminal modifications are resolved after the full canonical
                # sequence length is known.
                modification_specs.append(("CT", official_code, smiles))
                i = end_idx + 1
                continue

            raise ValueError(f"Unsupported MAP prefix '{prefix}' in block '{block_text}'.")

        # Any other character is invalid for this minimal parser.
        raise ValueError(
            f"Unexpected character '{ch}' at index {i}. "
            "Only letters, '(CODE)' blocks, and '{PREFIX:CODE}' blocks are supported."
        )

    caa_sequence = "".join(caa_chars)
    final_position = len(caa_sequence)

    modifications: List[Tuple[int, str, str]] = []
    for position_spec, mod_code, smiles in modification_specs:
        if position_spec == "CT":
            if final_position == 0:
                raise ValueError(f"C-terminal modification '{mod_code}' cannot be applied to an empty sequence.")
            position = final_position
        else:
            position = int(position_spec)
        modifications.append((position, mod_code, smiles))

    return caa_sequence, modifications


def write_outputs(caa_sequence: str, modifications: List[Tuple[int, str, str]], fasta_path: str, mods_path: str) -> None:
    """Write FASTA and modification mapping text files."""
    fasta_dir = os.path.dirname(fasta_path)
    mods_dir = os.path.dirname(mods_path)

    if fasta_dir:
        os.makedirs(fasta_dir, exist_ok=True)
    if mods_dir:
        os.makedirs(mods_dir, exist_ok=True)

    # Standard FASTA output for AF2/local colabfold backbone generation.
    with open(fasta_path, "w", encoding="utf-8") as f:
        f.write(">parsed_sequence\n")
        f.write(f"{caa_sequence}\n")

    # Exact line format required by downstream ETFlow + stitching logic.
    with open(mods_path, "w", encoding="utf-8") as f:
        for position, code, smiles in modifications:
            f.write(f"{position} : {code} : {smiles}\n")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse peptide sequence with legacy '(CODE)' or MAP '{prefix:Code}' modifications into CAA FASTA + modification map."
    )
    parser.add_argument(
        "--sequence",
        required=True,
        help="Input peptide sequence, e.g. APG(5PG)APG",
    )
    parser.add_argument(
        "--json",
        required=True,
        help="Path to modifications JSON file",
    )
    parser.add_argument(
        "--fasta_out",
        required=True,
        help="Output FASTA path",
    )
    parser.add_argument(
        "--mods_out",
        required=True,
        help="Output modification map path",
    )

    args = parser.parse_args()

    mod_index = load_modification_index(args.json)
    caa_sequence, modifications = parse_sequence(args.sequence, mod_index)
    write_outputs(caa_sequence, modifications, args.fasta_out, args.mods_out)


if __name__ == "__main__":
    main()
