"""Shared residue identity and reference-routing contract; standard library only.

The text file retains its three/four-field layout. Its companion JSON carries
the explicit backend and actual CCD reference. Consumers require the companion
file, so an old text file cannot silently trigger a lookup using a made-up TLC.
"""

import ast
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile

ROUTING_CONTRACT_VERSION = 1
AA_1_TO_3 = dict(zip(
    "ACDEFGHIKLMNPQRSTVWY",
    "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR".split(),
))


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def metadata_path(mods_path):
    return Path(mods_path).with_suffix(".metadata.json")


def validate_code(code, label="residue code"):
    if not isinstance(code, str) or not re.fullmatch(r"[A-Za-z0-9_]+", code):
        raise ValueError(f"Invalid {label}: {code!r}. Use letters, digits or underscores.")
    return code


def normalize_actual_code(record, context):
    if "Actual PDB code" not in record:
        raise ValueError(f"Missing 'Actual PDB code' for {context}; explicit null is required for hybrid routing.")
    value = record["Actual PDB code"]
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (str, int)):
        raise ValueError(f"Invalid 'Actual PDB code' for {context}: {value!r}.")
    code = str(value).strip().upper()
    if not re.fullmatch(r"[A-Z0-9]+", code) or code in {"NULL", "NONE"}:
        raise ValueError(f"Empty or invalid 'Actual PDB code' for {context}: {value!r}; this is not explicit JSON null.")
    return code


def make_record(position, code, smiles, parent, actual_code, source_format,
                input_block, json_entry_number=None):
    validate_code(code)
    if parent not in AA_1_TO_3:
        raise ValueError(f"Parent {parent!r} at position {position} is not one of the 20 canonical amino acid letters.")
    if not isinstance(smiles, str) or not smiles.strip():
        raise ValueError(f"SMILES is empty at position {position}.")
    backend = "hybrid" if actual_code is None else "actual_pdb"
    return {
        "position": position, "code": code, "pdb_residue_name": code.upper()[:3],
        "smiles": smiles, "parent_letter": parent, "actual_pdb_code": actual_code,
        "backend": backend,
        "reference_code": AA_1_TO_3[parent] if backend == "hybrid" else actual_code,
        "source_format": source_format, "input_block": input_block,
        "json_entry_number": json_entry_number,
        "modification_file": f"mod_{position}_{code}.pdb", "observations": {},
    }


def text_fields(record):
    parent = record["parent_letter"] if record["backend"] == "hybrid" else ""
    return (record["position"], record["code"], record["smiles"], parent)


def read_text_records(mods_path):
    records = []
    seen = set()
    for number, raw in enumerate(Path(mods_path).read_text(encoding="utf-8").splitlines(), 1):
        if not raw.strip():
            continue
        parts = [part.strip() for part in raw.strip().split(" : ")]
        if len(parts) not in (3, 4) or not parts[0].isdigit() or int(parts[0]) < 1:
            raise ValueError(f"Invalid modifications record on line {number}: {raw!r}.")
        position = int(parts[0])
        if position in seen:
            raise ValueError(f"Multiple modification definitions target residue position {position}.")
        seen.add(position)
        validate_code(parts[1])
        if not parts[2]:
            raise ValueError(f"SMILES is empty on line {number}.")
        records.append((position, parts[1], parts[2], parts[3] if len(parts) == 4 else ""))
    return records


def _atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(dir=path.parent, prefix=f".{path.name}.")
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            handle.write(text)
        os.replace(temp_name, path)
    finally:
        if os.path.exists(temp_name):
            os.unlink(temp_name)


def save_metadata(mods_path, document):
    _atomic_text(metadata_path(mods_path), json.dumps(document, indent=2, ensure_ascii=False) + "\n")
    write_identity_csv(mods_path, document)


def write_metadata(mods_path, records, canonical_sequence, input_sequence="",
                   database_path=None):
    if [text_fields(record) for record in records] != read_text_records(mods_path):
        raise ValueError("Residue metadata and modifications.txt disagree; no metadata was written.")
    document = {
        "schema_version": ROUTING_CONTRACT_VERSION,
        "modifications_sha256": sha256_file(mods_path),
        "canonical_sequence": canonical_sequence, "input_sequence": input_sequence,
        "database_file": str(Path(database_path).resolve()) if database_path else None,
        "database_sha256": sha256_file(database_path) if database_path else None,
        "records": records,
    }
    save_metadata(mods_path, document)
    return document


def load_metadata(mods_path):
    path = metadata_path(mods_path)
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing residue metadata: {path}. Rerun the updated 01_parse_input.py; "
            "old modifications.txt files do not contain the actual reference code."
        )
    document = json.loads(path.read_text(encoding="utf-8"))
    if document.get("schema_version") != ROUTING_CONTRACT_VERSION:
        raise ValueError("Unsupported residue-metadata version. Use all updated scripts together.")
    if document.get("modifications_sha256") != sha256_file(mods_path):
        raise ValueError("Residue metadata is stale: modifications.txt has changed. Rerun parsing.")
    records = document.get("records")
    if not isinstance(records, list):
        raise ValueError("Residue metadata has no valid records list.")
    for record in records:
        rebuilt = make_record(
            record["position"], record["code"], record["smiles"],
            record["parent_letter"], record["actual_pdb_code"],
            record["source_format"], record["input_block"], record.get("json_entry_number"),
        )
        for key in ("backend", "reference_code", "pdb_residue_name", "modification_file"):
            if record.get(key) != rebuilt[key]:
                raise ValueError(f"Inconsistent {key} in metadata at position {record['position']}.")
        actual = record["actual_pdb_code"]
        if actual is not None and normalize_actual_code({"Actual PDB code": actual}, "metadata") != actual:
            raise ValueError("Actual reference code in metadata is not normalized.")
        position = record["position"]
        sequence = document.get("canonical_sequence", "")
        if not isinstance(position, int) or not 1 <= position <= len(sequence):
            raise ValueError(f"Invalid metadata position {position!r}.")
        if sequence[position - 1] != record["parent_letter"]:
            raise ValueError(f"Metadata parent disagrees with canonical sequence at position {position}.")
    if [text_fields(record) for record in records] != read_text_records(mods_path):
        raise ValueError("Residue metadata and modifications.txt disagree. Rerun parsing.")
    return document


def write_identity_csv(mods_path, document):
    import io
    fields = [
        "sequence_position", "chain", "pdb_residue_number", "insertion_code",
        "full_residue_code", "pdb_residue_name", "actual_pdb_code", "parent_letter",
        "smiles", "input_format", "input_block", "json_entry_number", "backend",
        "reference_code", "modification_file", "stitched_verified",
        "protonated_verified", "final_verified",
    ]
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=fields)
    writer.writeheader()
    for record in document["records"]:
        observations = record.get("observations", {})
        observed = observations.get("final") or observations.get("protonated") or observations.get("stitched") or {}
        writer.writerow({
            "sequence_position": record["position"], "chain": observed.get("chain", ""),
            "pdb_residue_number": observed.get("residue_number", ""),
            "insertion_code": observed.get("insertion_code", ""),
            "full_residue_code": record["code"], "pdb_residue_name": record["pdb_residue_name"],
            "actual_pdb_code": "null" if record["actual_pdb_code"] is None else record["actual_pdb_code"],
            "parent_letter": record["parent_letter"], "smiles": record["smiles"],
            "input_format": record["source_format"], "input_block": record["input_block"],
            "json_entry_number": record.get("json_entry_number"), "backend": record["backend"],
            "reference_code": record["reference_code"], "modification_file": record["modification_file"],
            **{f"{stage}_verified": stage in observations for stage in ("stitched", "protonated", "final")},
        })
    _atomic_text(Path(mods_path).parent / "residue_identity.csv", buffer.getvalue())


def record_pdb_stage(mods_path, pdb_path, stage):
    """Verify labels and record observed identifiers; this does not verify chemistry."""
    if stage not in {"stitched", "protonated", "final"}:
        raise ValueError(f"Unsupported structure stage {stage!r}.")
    document = load_metadata(mods_path)
    observed = {}
    model_count = 0
    for line in Path(pdb_path).read_text(encoding="utf-8").splitlines():
        if line.startswith("MODEL "):
            model_count += 1
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        key = (line[21:22], int(line[22:26]), line[26:27].strip())
        data = observed.setdefault(key, {"names": set(), "atoms": set()})
        data["names"].add(line[17:20].strip())
        data["atoms"].add(line[12:16].strip())
    if model_count > 1:
        raise ValueError("Residue identity reporting requires a single-model PDB.")
    for record in document["records"]:
        previous = record.get("observations", {}).get("stitched")
        if previous and stage != "stitched":
            key = (previous["chain"], previous["residue_number"], previous["insertion_code"])
            candidates = [(key, observed[key])] if key in observed else []
        else:
            candidates = [(key, value) for key, value in observed.items() if key[1] == record["position"]]
        if len(candidates) != 1:
            raise ValueError(f"Cannot verify {stage} identity at position {record['position']}: expected one residue, found {len(candidates)}.")
        key, data = candidates[0]
        if data["names"] != {record["pdb_residue_name"]} or not {"N", "CA", "C"} <= data["atoms"]:
            raise ValueError(f"Residue identity/anchor mismatch in {stage} at position {record['position']}: {sorted(data['names'])}.")
        record.setdefault("observations", {})[stage] = {
            "chain": key[0], "residue_number": key[1], "insertion_code": key[2],
            "pdb_residue_name": record["pdb_residue_name"],
        }
    save_metadata(mods_path, document)
    return document


def check_pipeline_versions(project_root):
    """Detect mixed old/new stages without importing prediction dependencies."""
    for name in ("01_parse_input.py", "03_run_sidechains.py", "04_stitch.py", "05_minimize.py", "main.py"):
        version = None
        for node in ast.parse((Path(project_root) / name).read_text(encoding="utf-8")).body:
            if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "ROUTING_CONTRACT_VERSION" for t in node.targets):
                version = ast.literal_eval(node.value)
        if version != ROUTING_CONTRACT_VERSION:
            raise ValueError(f"{name} is not from this reference-routing update. Replace all updated files together.")
