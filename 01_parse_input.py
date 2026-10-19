#!/usr/bin/env python3
"""
01_parse_input.py  (updated approach)

Adapter parser for peptide sequences containing:
- Legacy modification blocks in parentheses, e.g. APG(5PG)APG
- MAP-style blocks in curly braces, e.g. MKT{ptm:meth}A
- D-amino-acid tags after canonical residues, e.g. A{d}CG
- Novel NCAA blocks in pipe delimiters, e.g. APG|SMILES,F|APG
"""

import argparse
import json
import os
from typing import Dict, List, Tuple, Any

from residue_metadata import (
    make_record, normalize_actual_code, text_fields, write_metadata,
)

ROUTING_CONTRACT_VERSION = 1

IndexA = Dict[Tuple[str, str], dict]
IndexB = Dict[str, dict]

def _normalize_map_block(value: str) -> str:
    block = value.strip()

    if not (block.startswith("{") and block.endswith("}")):
        return block.lower()

    inner = block[1:-1].strip()

    if ":" in inner:
        prefix, code = inner.split(":", 1)
        inner = f"{prefix.strip()}:{code.strip()}"

    return "{" + inner.lower() + "}"

def load_modification_index(json_path: str) -> Tuple[IndexA, IndexB]:
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    if not isinstance(data, list):
        raise ValueError("Modifications JSON must be a top-level list of records.")

    sub_index_a: IndexA = {}
    sub_index_b: IndexB = {}

    for entry_number, entry in enumerate(data, 1):
        if not isinstance(entry, dict):
            continue

        entry["_source_entry_number"] = entry_number
        tlc = str(entry.get("Three letter code") or "").strip()
        if tlc:
            sub_index_b[tlc.upper()] = entry

        user_input = str(entry.get("User Input") or "").strip()
        natural_aa = str(entry.get("Natural Amino Acid") or "").strip()

        if not user_input:
            continue

        one_letter = _extract_one_letter_safe(natural_aa)
        if one_letter is None:
            continue

        key_a = (_normalize_map_block(user_input), one_letter.upper())
        sub_index_a[key_a] = entry

    # Register D-residue entries from MAP notation, e.g. A{d}.
    # These entries can lack User Input, so examine all JSON records.
    canonical_letters = set("ACDEFGHIKLMNPQRSTVWY")

    for entry in data:
        if not isinstance(entry, dict):
            continue

        natural_aa = str(entry.get("Natural Amino Acid") or "").strip()
        parent_letter = _extract_one_letter_safe(natural_aa)
        map_notation = str(entry.get("MAP notation") or "").strip()

        if parent_letter not in canonical_letters or not map_notation:
            continue
        if map_notation[0].upper() != parent_letter:
            continue
        if _normalize_map_block(map_notation[1:]) != "{d}":
            continue

        key_d = ("{d}", parent_letter)
        previous = sub_index_a.get(key_d)
        if previous is not None and (
            {k: v for k, v in previous.items() if k != "_source_entry_number"}
            != {k: v for k, v in entry.items() if k != "_source_entry_number"}
        ):
            raise ValueError(
                f"Conflicting D-residue entries for parent '{parent_letter}'."
            )

        sub_index_a[key_d] = entry

    return sub_index_a, sub_index_b


def _extract_one_letter_safe(natural_aa_field: str) -> str | None:
    parts = [p.strip() for p in natural_aa_field.split("/") if p.strip()]
    if not parts:
        return None
    candidate = parts[-1]
    if len(candidate) == 1 and candidate.isalpha():
        return candidate.upper()
    return None

def _lookup_index_a(full_block: str, amino_acid: str, sub_index_a: IndexA, context: str) -> dict:
    key = (_normalize_map_block(full_block), amino_acid.upper())
    entry = sub_index_a.get(key)
    if entry is None:
        raise ValueError(f"No JSON entry found for {context}.")
    return entry

def _lookup_index_b(tlc: str, sub_index_b: IndexB, context: str) -> dict:
    entry = sub_index_b.get(tlc.upper())
    if entry is None:
        raise ValueError(f"No JSON entry found for {context}.")
    return entry

def _get_tlc(record: dict, context: str) -> str:
    tlc = str(record.get("Three letter code") or "").strip()
    if not tlc:
        raise ValueError(f"Missing 'Three letter code' in JSON record for {context}.")
    return tlc

def _get_smiles(record: dict, context: str) -> str:
    smiles = str(record.get("SMILES") or "").strip()
    if not smiles:
        raise ValueError(f"Missing 'SMILES' in JSON record for {context}.")
    return smiles

def _extract_parent_one_letter(record: dict, context: str) -> str:
    natural_aa = str(record.get("Natural Amino Acid") or "").strip()
    parts = [p.strip() for p in natural_aa.split("/") if p.strip()]
    if not parts:
        raise ValueError(f"Missing or empty 'Natural Amino Acid' field for {context}.")
    one_letter = parts[-1]
    if one_letter.lower() == "unknown":
        raise ValueError(f"'Natural Amino Acid' is 'Unknown' for {context}.")
    if len(one_letter) != 1 or not one_letter.isalpha():
        raise ValueError(f"Could not extract a valid one-letter code for {context}.")
    return one_letter.upper()


def parse_sequence(sequence: str, sub_index_a: IndexA, sub_index_b: IndexB,
                   metadata_records=None) -> Tuple[str, List[Tuple[int, str, str, str]]]:
    caa_chars: List[str] = []
    resolved_mods: List[Tuple[int, str, str, str]] = []
    deferred: List[Dict[str, Any]] = []
    identities = []

    def append_json_modification(position, record, source_format, input_block, context):
        identity = make_record(
            position, _get_tlc(record, context), _get_smiles(record, context),
            _extract_parent_one_letter(record, context),
            normalize_actual_code(record, context), source_format, input_block,
            record.get("_source_entry_number"),
        )
        identities.append(identity)
        resolved_mods.append(text_fields(identity))

    residue_position = 1
    new_ncaa_counter = 1
    i = 0
    n = len(sequence)

    # while i < n:
    #     ch = sequence[i]
    
    while i < n:
        ch = sequence[i]

        if ch.isspace():
            i += 1
            continue

        if ch.isalpha():

        # if ch.isalpha():
            caa_chars.append(ch.upper())
            residue_position += 1
            i += 1
            continue

        # ── Legacy format: (TLC) 
        if ch == "(":
            end_idx = sequence.find(")", i + 1)

            if end_idx == -1:
                raise ValueError(
                    f"Unclosed PDB-code block: missing ')' "
                    f"for '(' at character {i + 1}."
                )
        
            tlc = sequence[i + 1 : end_idx].strip()
            context = f"legacy block '({tlc})'"
            record = _lookup_index_b(tlc, sub_index_b, context)
            official_tlc = _get_tlc(record, context)
            smiles = _get_smiles(record, context)
            parent_letter = _extract_parent_one_letter(record, context)

            caa_chars.append(parent_letter)
            append_json_modification(residue_position, record, "pdb", sequence[i:end_idx + 1], context)
            residue_position += 1
            i = end_idx + 1
            continue

        # ── Novel NCAA format: |SMILES,Parent| 
        if ch == "|":
            end_idx = sequence.find("|", i + 1)
            if end_idx == -1:
                raise ValueError(f"Unclosed pipe delimiter starting at index {i}.")
            
            block_text = sequence[i + 1 : end_idx].strip()
            if "," not in block_text:
                raise ValueError(f"Novel NCAA block must contain a comma: |SMILES,Parent|. Got: {block_text}")
            
            # smiles_part, parent_raw = [p.strip() for p in block_text.rsplit(",", 1)]
            # parent_letter = parent_raw.upper()
            smiles_part, parent_raw = [p.strip() for p in block_text.rsplit(",", 1)]

            if not smiles_part:
                raise ValueError(
                    f"Hybrid residue at position {residue_position}: SMILES is empty. "
                    "Provide SMILES before the comma in |SMILES,Parent|."
                )

            parent_letter = parent_raw.upper()

            canonical_letters = set("ACDEFGHIKLMNPQRSTVWY")

            if parent_letter not in canonical_letters:
                raise ValueError(
                    f"Hybrid residue at position {residue_position}: "
                    f"parent '{parent_letter}' is not one of the "
                    "20 canonical amino acid letters."
                )
                
            # PDB residue names are exactly three characters wide.
            # Keep every pipe-defined novel residue uniquely identifiable.
            if 1 <= new_ncaa_counter <= 9:
                mod_code = f"N_{new_ncaa_counter}"      # N_1 ... N_9
            elif 10 <= new_ncaa_counter <= 99:
                mod_code = f"N{new_ncaa_counter:02d}"   # N10 ... N99
            else:
                raise ValueError(
                    "A maximum of 99 pipe-defined novel residues is supported "
                    "in one sequence because PDB residue names have only 3 characters."
                )

            new_ncaa_counter += 1
            
            caa_chars.append(parent_letter)
            identity = make_record(
                residue_position, mod_code, smiles_part, parent_letter,
                None, "hybrid", sequence[i:end_idx + 1],
            )
            identities.append(identity)
            resolved_mods.append(text_fields(identity))
            
            residue_position += 1
            i = end_idx + 1
            continue

        # ── MAP format: {prefix:Code} or {TLC} 
        if ch == "{":
            end_idx = sequence.find("}", i + 1)

            if end_idx == -1:
                raise ValueError(
                    f"Unclosed MAP block: missing '}}' "
                    f"for '{{' at character {i + 1}."
                )

            block_text = sequence[i + 1 : end_idx].strip()
            full_block = sequence[i : end_idx + 1]

            if block_text.lower() == "d":
                # Find the character immediately before the tag, ignoring spaces.
                previous_index = i - 1
                while previous_index >= 0 and sequence[previous_index].isspace():
                    previous_index -= 1

                if previous_index < 0:
                    raise ValueError(
                        "D-amino-acid tag '{d}' requires a preceding "
                        "canonical amino acid letter."
                    )

                parent_letter = sequence[previous_index].upper()
                if parent_letter not in set("ACDEFGHIKLMNPQRSTVWY"):
                    raise ValueError(
                        "D-amino-acid tag '{d}' must directly follow a "
                        "canonical amino acid letter, allowing whitespace. "
                        "It cannot follow another modification or residue block."
                    )

                position = residue_position - 1

                if parent_letter == "G":
                    print(
                        f"[parse_input] Position {position}: glycine has no "
                        "D/L distinction; {d} leaves it unchanged."
                    )
                else:
                    context = f"D-amino-acid tag '{{d}}' on '{parent_letter}' at position {position}"
                    record = _lookup_index_a(
                        "{d}", parent_letter, sub_index_a, context
                    )
                    official_tlc = _get_tlc(record, context)
                    smiles = _get_smiles(record, context)
                    append_json_modification(position, record, "d", full_block, context)

                i = end_idx + 1
                continue

            if ":" not in block_text:
                tlc = block_text
                context = f"shorthand block '{{{tlc}}}'"
                record = _lookup_index_b(tlc, sub_index_b, context)
                official_tlc = _get_tlc(record, context)
                smiles = _get_smiles(record, context)
                parent_letter = _extract_parent_one_letter(record, context)

                caa_chars.append(parent_letter)
                append_json_modification(residue_position, record, "shorthand", full_block, context)
                residue_position += 1
                i = end_idx + 1
                continue

            prefix_raw, mod_code = [p.strip() for p in block_text.split(":", 1)]
            prefix = prefix_raw.lower()

            if prefix == "nnr":
                context = f"MAP block '{full_block}'"
                record = _lookup_index_b(mod_code, sub_index_b, context)
                official_tlc = _get_tlc(record, context)
                smiles = _get_smiles(record, context)
                parent_letter = _extract_parent_one_letter(record, context)

                caa_chars.append(parent_letter)
                append_json_modification(residue_position, record, "nnr", full_block, context)
                residue_position += 1
                i = end_idx + 1
                continue

            # if prefix == "ptm":
                # preceding_aa = caa_chars[-1]
            if prefix == "ptm":
                if not caa_chars:
                    raise ValueError(
                        f"PTM tag '{full_block}' requires a preceding residue."
                    )

                preceding_aa = caa_chars[-1]
                context = f"MAP block '{full_block}' modifying preceding residue '{preceding_aa}'"
                record = _lookup_index_a(full_block, preceding_aa, sub_index_a, context)
                official_tlc = _get_tlc(record, context)
                smiles = _get_smiles(record, context)
                append_json_modification(residue_position - 1, record, "ptm", full_block, context)
                i = end_idx + 1
                continue

            if prefix == "nt":
                deferred.append({"tag": "NT", "full_block": full_block})
                i = end_idx + 1
                continue

            if prefix == "ct":
                deferred.append({"tag": "CT", "full_block": full_block})
                i = end_idx + 1
                continue

            raise ValueError(f"Unsupported MAP prefix '{prefix}'.")

        raise ValueError(f"Unexpected character '{ch}' at index {i}.")

    caa_sequence = "".join(caa_chars)
    final_position = len(caa_sequence)

    # for item in deferred:
    #     tag = item["tag"]
    #     full_block = item["full_block"]

    #     if tag == "NT":
    
    for item in deferred:
        tag = item["tag"]
        full_block = item["full_block"]

        if not caa_sequence:
            terminal_name = "N-terminal" if tag == "NT" else "C-terminal"
            raise ValueError(
                f"{terminal_name} tag '{full_block}' "
                "requires a nonempty peptide sequence."
            )

        if tag == "NT":
            amino_acid = caa_sequence[0]
            position = 1
            context = f"MAP block '{full_block}' (N-terminal on '{amino_acid}')"
        elif tag == "CT":
            amino_acid = caa_sequence[-1]
            position = final_position
            context = f"MAP block '{full_block}' (C-terminal on '{amino_acid}')"
        
        record = _lookup_index_a(full_block, amino_acid, sub_index_a, context)
        official_tlc = _get_tlc(record, context)
        smiles = _get_smiles(record, context)
        append_json_modification(position, record, tag.lower(), full_block, context)
    
    modifications: List[Tuple[int, str, str, str]] = sorted(
        resolved_mods, key=lambda x: x[0]
    )

    seen_positions = {}

    for position, code, smiles, parent in modifications:
        if position in seen_positions:
            previous_code = seen_positions[position]
            raise ValueError(
                f"Multiple modification definitions target residue position "
                f"{position}: '{previous_code}' and '{code}'. "
                "Only one modification definition per position is supported."
            )

        seen_positions[position] = code

    if metadata_records is not None:
        metadata_records.extend(sorted(identities, key=lambda item: item["position"]))
    return caa_sequence, modifications

def write_outputs(caa_sequence: str, modifications: List[Tuple[int, str, str, str]],
                  fasta_path: str, mods_path: str, metadata_records=None,
                  input_sequence="", database_path=None) -> None:
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
        for position, code, smiles, parent in modifications:
            if parent:
                # Direct hybrid inputs and JSON entries with null actual codes use parent-based naming.
                f.write(f"{position} : {code} : {smiles} : {parent}\n")
            else:
                # Actual-code entries use the explicit reference recorded in metadata.
                f.write(f"{position} : {code} : {smiles}\n")

    if metadata_records is not None:
        write_metadata(mods_path, metadata_records, caa_sequence,
                       input_sequence=input_sequence, database_path=database_path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Parse a peptide sequence into a canonical FASTA file and modifications map.")
    parser.add_argument("--sequence", required=True)
    parser.add_argument("--json", required=True)
    parser.add_argument("--fasta_out", required=True)
    parser.add_argument("--mods_out", required=True)
    args = parser.parse_args()

    sub_index_a, sub_index_b = load_modification_index(args.json)
    metadata_records = []
    caa_sequence, modifications = parse_sequence(
        args.sequence, sub_index_a, sub_index_b, metadata_records=metadata_records,
    )
    write_outputs(
        caa_sequence, modifications, args.fasta_out, args.mods_out,
        metadata_records=metadata_records, input_sequence=args.sequence,
        database_path=args.json,
    )

    print(f"[parse_input] Canonical sequence : {caa_sequence}")
    print(f"[parse_input] Modifications found: {len(modifications)}")
    print("[parse_input] Wrote companion metadata and residue_identity.csv")

if __name__ == "__main__":
    main()