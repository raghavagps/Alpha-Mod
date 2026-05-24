#!/usr/bin/env python3
"""
01_parse_input.py  (updated approach)

Adapter parser for peptide sequences containing:
- Legacy modification blocks in parentheses, e.g. APG(5PG)APG
- MAP-style blocks in curly braces, e.g. MKT{ptm:meth}A

Two separate lookup strategies are used depending on the tag type:

  Strategy A  —  nt, ct, ptm tags
      Index key: (user_input_string_lowercased, parent_one_letter_AA)
      Example:   ("{ct:amid}", "L")  →  L02 entry (Leucine amide)
      Resolves the overwrite problem and selects the correct residue-specific
      entry automatically. Case-insensitive on both dimensions.

  Strategy B  —  nnr tag and legacy (TLC) format
      Index key: three_letter_code_uppercased
      Example:   "NVA"  →  Norvaline entry
      TLCs are unique identifiers so no overwrite problem exists here.

Outputs (unchanged, fully compatible with downstream scripts):
  1) FASTA  — canonical amino acid sequence for AlphaFold/ESMFold backbone
  2) TXT    — one modification per line:  position : TLC : SMILES
"""

import argparse
import json
import os
from typing import Dict, List, Tuple, Any


# ---------------------------------------------------------------------------
# Type aliases
# ---------------------------------------------------------------------------
IndexA = Dict[Tuple[str, str], dict]   # (user_input_lower, one_letter_AA) → entry
IndexB = Dict[str, dict]               # tlc_upper → entry


# ---------------------------------------------------------------------------
# Index building
# ---------------------------------------------------------------------------

def load_modification_index(json_path: str) -> Tuple[IndexA, IndexB]:
    """
    Build two separate lookup indexes from the modifications JSON file.

    Sub-index A  —  keyed by (user_input.lower(), one_letter_AA.upper())
                    Used for nt / ct / ptm MAP tags.
                    Eliminates overwrite problem and selects the correct
                    residue-specific entry for each (modification, amino acid)
                    combination automatically.

    Sub-index B  —  keyed by three_letter_code.upper()
                    Used for nnr MAP tag and legacy (TLC) parenthesis format.
                    TLCs are unique PDB identifiers so a flat index is correct.

    Both lookups are case-insensitive by normalising keys at build time.
    """
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    if not isinstance(data, list):
        raise ValueError("Modifications JSON must be a top-level list of records.")

    sub_index_a: IndexA = {}
    sub_index_b: IndexB = {}

    for entry in data:
        if not isinstance(entry, dict):
            continue

        # --- Sub-index B: keyed by Three Letter Code (nnr + legacy format) ---
        tlc = str(entry.get("Three letter code") or "").strip()
        if tlc:
            sub_index_b[tlc.upper()] = entry

        # --- Sub-index A: keyed by (User Input lower, parent AA) ---
        user_input = str(entry.get("User Input") or "").strip()
        natural_aa = str(entry.get("Natural Amino Acid") or "").strip()

        if not user_input:
            # Entry has no User Input → only reachable via TLC (sub-index B).
            continue

        one_letter = _extract_one_letter_safe(natural_aa)
        if one_letter is None:
            # Natural AA is Unknown or missing → cannot build a valid A key.
            # Entry is still reachable via sub-index B if it has a TLC.
            continue

        key_a = (user_input.lower(), one_letter.upper())
        sub_index_a[key_a] = entry

    return sub_index_a, sub_index_b


def _extract_one_letter_safe(natural_aa_field: str) -> str | None:
    """
    Silently return None instead of raising when the Natural Amino Acid field
    is missing, 'Unknown', or in an unexpected format.

    Used only during index build.  Parsing raises a proper error if a
    looked-up entry turns out to have an unusable Natural AA field.
    """
    parts = [p.strip() for p in natural_aa_field.split("/") if p.strip()]
    if not parts:
        return None
    candidate = parts[-1]
    if len(candidate) == 1 and candidate.isalpha():
        return candidate.upper()
    return None


# ---------------------------------------------------------------------------
# Record retrieval helpers
# ---------------------------------------------------------------------------

def _lookup_index_a(
    full_block: str,
    amino_acid: str,
    sub_index_a: IndexA,
    context: str,
) -> dict:
    """
    Look up an entry in sub-index A using the full MAP block and amino acid.

    full_block  — the complete block as it appeared in the sequence,
                  e.g. '{ct:Amid}' or '{ptm:meth}'.  Lowercased here.
    amino_acid  — one-letter canonical AA code (uppercase), e.g. 'L'.
    context     — human-readable description for error messages.
    """
    key = (full_block.lower(), amino_acid.upper())
    entry = sub_index_a.get(key)
    if entry is None:
        raise ValueError(
            f"No JSON entry found for {context}.\n"
            f"  Looked up: User Input='{full_block.lower()}' + "
            f"Amino Acid='{amino_acid.upper()}'\n"
            f"  Check that your modifications JSON has an entry with:\n"
            f"    \"User Input\": \"{full_block.lower()}\"  and  "
            f"\"Natural Amino Acid\": \".../{amino_acid.upper()}\""
        )
    return entry


def _lookup_index_b(
    tlc: str,
    sub_index_b: IndexB,
    context: str,
) -> dict:
    """
    Look up an entry in sub-index B using the Three Letter Code.

    tlc     — the TLC as written in the sequence, e.g. 'NVA' or '5PG'.
               Uppercased here.
    context — human-readable description for error messages.
    """
    entry = sub_index_b.get(tlc.upper())
    if entry is None:
        raise ValueError(
            f"No JSON entry found for {context}.\n"
            f"  Looked up Three Letter Code: '{tlc.upper()}'\n"
            f"  Check that your modifications JSON has an entry with "
            f"\"Three letter code\": \"{tlc.upper()}\""
        )
    return entry


def _get_tlc(record: dict, context: str) -> str:
    """Extract and validate the Three Letter Code from a JSON record."""
    tlc = str(record.get("Three letter code") or "").strip()
    if not tlc:
        raise ValueError(f"Missing 'Three letter code' in JSON record for {context}.")
    return tlc


def _get_smiles(record: dict, context: str) -> str:
    """Extract and validate the SMILES string from a JSON record."""
    smiles = str(record.get("SMILES") or "").strip()
    if not smiles:
        raise ValueError(f"Missing 'SMILES' in JSON record for {context}.")
    return smiles


def _extract_parent_one_letter(record: dict, context: str) -> str:
    """
    Extract the parent canonical one-letter amino acid code from a JSON record.

    Reads the 'Natural Amino Acid' field which is formatted as:
        'Leucine/Leu/L'   →  returns 'L'
        'Glycine/Gly/G'   →  returns 'G'

    Raises a clear error if the field is missing, Unknown, or malformed.
    Used by nnr tag and legacy (TLC) format — both of which insert a
    canonical letter into the FASTA sequence.
    """
    natural_aa = str(record.get("Natural Amino Acid") or "").strip()
    parts = [p.strip() for p in natural_aa.split("/") if p.strip()]
    if not parts:
        raise ValueError(
            f"Missing or empty 'Natural Amino Acid' field for {context}.\n"
            f"  This entry cannot be used as a standalone residue."
        )
    one_letter = parts[-1]
    if one_letter.lower() == "unknown":
        raise ValueError(
            f"'Natural Amino Acid' is 'Unknown' for {context}.\n"
            f"  This entry cannot be used as a standalone residue (nnr or legacy format).\n"
            f"  Fix the JSON entry or use a different modification code."
        )
    if len(one_letter) != 1 or not one_letter.isalpha():
        raise ValueError(
            f"Could not extract a valid one-letter code from "
            f"'Natural Amino Acid'={natural_aa!r} for {context}."
        )
    return one_letter.upper()


# ---------------------------------------------------------------------------
# Core sequence parser
# ---------------------------------------------------------------------------

def parse_sequence(
    sequence: str,
    sub_index_a: IndexA,
    sub_index_b: IndexB,
) -> Tuple[str, List[Tuple[int, str, str]]]:
    """
    Parse an input peptide sequence and return the canonical amino acid
    sequence plus a list of modification records.

    Supported input formats
    -----------------------
    Plain letters         : treated as canonical amino acids, added directly.

    (TLC)                 : Legacy format.  TLC is the PDB Three Letter Code.
                            Inserts the parent canonical residue into the FASTA
                            and logs a modification at that position.
                            Lookup via sub-index B.

    {nnr:TLC}             : Non-natural residue in MAP format.
                            Identical behaviour to legacy (TLC).
                            Lookup via sub-index B.

    {ptm:Code}            : Post-translational modification on the immediately
                            preceding canonical residue.  Does NOT insert a new
                            residue.  Amino acid is known immediately (caa_chars[-1]).
                            Lookup via sub-index A using (full_block, preceding_AA).

    {nt:Code}             : N-terminal modification.  Logged at position 1.
                            Amino acid is not known until the full sequence is
                            built → DEFERRED.  Lookup via sub-index A after loop.

    {ct:Code}             : C-terminal modification.  Logged at final position.
                            Amino acid is not known until the full sequence is
                            built → DEFERRED.  Lookup via sub-index A after loop.

    Returns
    -------
    caa_sequence   : str
        Canonical amino acid sequence (letters only, no modification blocks).

    modifications  : List[Tuple[int, str, str]]
        Each tuple is (position, three_letter_code, smiles).
        Position is 1-based.
        This is written directly to modifications.txt by write_outputs().
    """
    caa_chars: List[str] = []

    # Resolved during parsing (legacy, nnr, ptm).
    resolved_mods: List[Tuple[int, str, str]] = []

    # Deferred until full sequence is known (nt, ct).
    # Each item: {"tag": "NT"|"CT", "full_block": str}
    deferred: List[Dict[str, Any]] = []

    residue_position = 1   # 1-based, tracks next position to be assigned
    i = 0
    n = len(sequence)

    while i < n:
        ch = sequence[i]

        # ── Plain canonical residue ──────────────────────────────────────────
        if ch.isalpha():
            caa_chars.append(ch.upper())
            residue_position += 1
            i += 1
            continue

        # ── Legacy format: (TLC) ─────────────────────────────────────────────
        if ch == "(":
            end_idx = sequence.find(")", i + 1)
            if end_idx == -1:
                raise ValueError(
                    f"Unclosed parenthesis starting at index {i} in sequence."
                )

            tlc = sequence[i + 1 : end_idx].strip()
            if not tlc:
                raise ValueError(f"Empty TLC inside parentheses at index {i}.")

            context = f"legacy block '({tlc})'"
            record = _lookup_index_b(tlc, sub_index_b, context)
            official_tlc = _get_tlc(record, context)
            smiles = _get_smiles(record, context)
            parent_letter = _extract_parent_one_letter(record, context)

            caa_chars.append(parent_letter)
            resolved_mods.append((residue_position, official_tlc, smiles))
            residue_position += 1
            i = end_idx + 1
            continue

        # ── MAP format: {prefix:Code} or {TLC} ──────────────────────────────
        if ch == "{":
            end_idx = sequence.find("}", i + 1)
            if end_idx == -1:
                raise ValueError(
                    f"Unclosed curly brace starting at index {i} in sequence."
                )

            block_text = sequence[i + 1 : end_idx].strip()
            if not block_text:
                raise ValueError(f"Empty MAP block at index {i}.")

            # Full block including braces, used as lookup key for sub-index A.
            full_block = sequence[i : end_idx + 1]   # e.g. '{ct:Amid}'

            # ── No colon: fallback shorthand {TLC} ──────────────────────────
            # Treat exactly like legacy (TLC): insert parent letter + log mod.
            if ":" not in block_text:
                tlc = block_text
                context = f"shorthand block '{{{tlc}}}'"
                record = _lookup_index_b(tlc, sub_index_b, context)
                official_tlc = _get_tlc(record, context)
                smiles = _get_smiles(record, context)
                parent_letter = _extract_parent_one_letter(record, context)

                caa_chars.append(parent_letter)
                resolved_mods.append((residue_position, official_tlc, smiles))
                residue_position += 1
                i = end_idx + 1
                continue

            # ── Split prefix and code ────────────────────────────────────────
            prefix_raw, mod_code = [p.strip() for p in block_text.split(":", 1)]
            if not prefix_raw or not mod_code:
                raise ValueError(
                    f"Invalid MAP block at index {i}: {block_text!r}. "
                    "Expected format: {{prefix:Code}}"
                )
            prefix = prefix_raw.lower()

            # ── {nnr:TLC} — non-natural residue ─────────────────────────────
            # Uses sub-index B (TLC lookup). Inserts parent canonical letter.
            if prefix == "nnr":
                context = f"MAP block '{full_block}'"
                record = _lookup_index_b(mod_code, sub_index_b, context)
                official_tlc = _get_tlc(record, context)
                smiles = _get_smiles(record, context)
                parent_letter = _extract_parent_one_letter(record, context)

                caa_chars.append(parent_letter)
                resolved_mods.append((residue_position, official_tlc, smiles))
                residue_position += 1
                i = end_idx + 1
                continue

            # ── {ptm:Code} — post-translational modification ─────────────────
            # Uses sub-index A. Amino acid = immediately preceding residue.
            # Known right now so lookup happens immediately, no deferral.
            if prefix == "ptm":
                if not caa_chars:
                    raise ValueError(
                        f"PTM block '{full_block}' at index {i} has no "
                        "preceding residue to modify."
                    )
                preceding_aa = caa_chars[-1]
                context = (
                    f"MAP block '{full_block}' modifying preceding "
                    f"residue '{preceding_aa}'"
                )
                record = _lookup_index_a(full_block, preceding_aa, sub_index_a, context)
                official_tlc = _get_tlc(record, context)
                smiles = _get_smiles(record, context)

                # PTM does NOT insert a new residue. Modifies the preceding one.
                resolved_mods.append((residue_position - 1, official_tlc, smiles))
                i = end_idx + 1
                continue

            # ── {nt:Code} — N-terminal modification ─────────────────────────
            # Uses sub-index A. Amino acid = first residue of full sequence.
            # Not known yet → DEFERRED until after the parsing loop.
            if prefix == "nt":
                deferred.append({"tag": "NT", "full_block": full_block})
                i = end_idx + 1
                continue

            # ── {ct:Code} — C-terminal modification ─────────────────────────
            # Uses sub-index A. Amino acid = last residue of full sequence.
            # Not known yet → DEFERRED until after the parsing loop.
            if prefix == "ct":
                deferred.append({"tag": "CT", "full_block": full_block})
                i = end_idx + 1
                continue

            raise ValueError(
                f"Unsupported MAP prefix '{prefix}' in block '{full_block}'. "
                "Supported prefixes: nt, ct, ptm, nnr."
            )

        # ── Unexpected character ─────────────────────────────────────────────
        raise ValueError(
            f"Unexpected character '{ch}' at index {i}. "
            "Only letters, '(TLC)' blocks, and '{prefix:Code}' blocks are supported."
        )

    # ── Build canonical sequence ─────────────────────────────────────────────
    caa_sequence = "".join(caa_chars)
    final_position = len(caa_sequence)

    if not caa_sequence:
        raise ValueError("Parsed canonical sequence is empty.")

    # ── Resolve deferred nt and ct lookups ───────────────────────────────────
    # Now that the full sequence is known, we can identify the first and last
    # residues and perform the sub-index A lookup for nt and ct modifications.
    for item in deferred:
        tag = item["tag"]
        full_block = item["full_block"]

        if tag == "NT":
            amino_acid = caa_sequence[0]
            position = 1
            context = (
                f"MAP block '{full_block}' (N-terminal modification on "
                f"first residue '{amino_acid}')"
            )

        elif tag == "CT":
            amino_acid = caa_sequence[-1]
            position = final_position
            context = (
                f"MAP block '{full_block}' (C-terminal modification on "
                f"last residue '{amino_acid}')"
            )

        else:
            raise RuntimeError(f"Unknown deferred tag type: {tag!r}")

        record = _lookup_index_a(full_block, amino_acid, sub_index_a, context)
        official_tlc = _get_tlc(record, context)
        smiles = _get_smiles(record, context)

        resolved_mods.append((position, official_tlc, smiles))

    # ── Combine and sort all modifications by position ───────────────────────
    # Sorting ensures modifications.txt is always in sequence order regardless
    # of where nt/ct tags appeared in the input string.
    modifications: List[Tuple[int, str, str]] = sorted(
        resolved_mods, key=lambda x: x[0]
    )

    return caa_sequence, modifications


# ---------------------------------------------------------------------------
# Output writer  (unchanged — fully compatible with downstream scripts)
# ---------------------------------------------------------------------------

def write_outputs(
    caa_sequence: str,
    modifications: List[Tuple[int, str, str]],
    fasta_path: str,
    mods_path: str,
) -> None:
    """
    Write the two output files consumed by downstream pipeline steps.

    parsed_sequence.fasta
        Standard single-entry FASTA consumed by 02_run_backbone.py
        (AlphaFold2 / ESMFold).

    modifications.txt
        One line per modification in the exact format required by
        03_run_sidechains.py (ETFlow) and 04_stitch.py:
            position : three_letter_code : SMILES
    """
    fasta_dir = os.path.dirname(fasta_path)
    mods_dir = os.path.dirname(mods_path)

    if fasta_dir:
        os.makedirs(fasta_dir, exist_ok=True)
    if mods_dir:
        os.makedirs(mods_dir, exist_ok=True)

    with open(fasta_path, "w", encoding="utf-8") as f:
        f.write(">parsed_sequence\n")
        f.write(f"{caa_sequence}\n")

    with open(mods_path, "w", encoding="utf-8") as f:
        for position, code, smiles in modifications:
            f.write(f"{position} : {code} : {smiles}\n")


# ---------------------------------------------------------------------------
# CLI entry point  (unchanged — fully compatible with main.py)
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Parse a peptide sequence containing legacy '(TLC)' or MAP "
            "'{prefix:Code}' modification blocks into a canonical FASTA file "
            "and a modifications mapping text file."
        )
    )
    parser.add_argument(
        "--sequence",
        required=True,
        help=(
            "Input peptide sequence string. Examples:\n"
            "  Legacy format : APG(5PG)APG\n"
            "  MAP format    : MKT{ptm:meth}A{ct:Amid}"
        ),
    )
    parser.add_argument(
        "--json",
        required=True,
        help="Path to the modifications JSON file (e.g. merged_final_data.json).",
    )
    parser.add_argument(
        "--fasta_out",
        required=True,
        help="Output path for the canonical FASTA file.",
    )
    parser.add_argument(
        "--mods_out",
        required=True,
        help="Output path for the modifications mapping text file.",
    )

    args = parser.parse_args()

    sub_index_a, sub_index_b = load_modification_index(args.json)
    caa_sequence, modifications = parse_sequence(args.sequence, sub_index_a, sub_index_b)
    write_outputs(caa_sequence, modifications, args.fasta_out, args.mods_out)

    print(f"[parse_input] Canonical sequence : {caa_sequence}")
    print(f"[parse_input] Modifications found: {len(modifications)}")
    for pos, code, _ in modifications:
        print(f"  position {pos:3d}  →  {code}")


if __name__ == "__main__":
    main()
