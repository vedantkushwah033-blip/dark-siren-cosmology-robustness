"""Fail-closed inventory of historical random seeds before training/holdout partitioning.

This utility inventories extracted, retained experiment artifacts. It does not fetch
GitHub Actions artifacts and cannot establish completeness unless every relevant
historical artifact has first been supplied under the declared roots.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path
from typing import Any

SUPPORTED = {".json", ".csv"}
SEED_KEY = re.compile(r"seed", re.IGNORECASE)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _seed_int(value: Any, *, source: str, field: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{source}: boolean seed in {field}")
    if isinstance(value, int):
        return value
    if isinstance(value, float) and value.is_integer():
        return int(value)
    if isinstance(value, str) and re.fullmatch(r"[+-]?\\d+", value.strip()):
        return int(value.strip())
    raise ValueError(f"{source}: non-integer seed value in {field}: {value!r}")


def _walk_json(value: Any, source: str, prefix: str = "$") -> list[dict[str, Any]]:
    found: list[dict[str, Any]] = []
    if isinstance(value, dict):
        for key, child in value.items():
            field = f"{prefix}.{key}"
            if SEED_KEY.search(str(key)):
                if isinstance(child, (list, tuple)):
                    for index, item in enumerate(child):
                        if isinstance(item, (dict, list)):
                            found.extend(_walk_json(item, source, f"{field}[{index}]"))
                        else:
                            found.append({"field": f"{field}[{index}]", "seed": _seed_int(item, source=source, field=field)})
                elif isinstance(child, dict):
                    # Preserve nested seed maps by recursively inspecting them.
                    found.extend(_walk_json(child, source, field))
                else:
                    found.append({"field": field, "seed": _seed_int(child, source=source, field=field)})
            else:
                found.extend(_walk_json(child, source, field))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            found.extend(_walk_json(child, source, f"{prefix}[{index}]"))
    return found


def inventory_file(path: Path, root: Path) -> dict[str, Any]:
    suffix = path.suffix.lower()
    relative = path.relative_to(root).as_posix()
    digest = sha256_file(path)
    if suffix == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        seeds = _walk_json(payload, relative)
    elif suffix == ".csv":
        seeds = []
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise ValueError(f"{relative}: CSV has no header")
            for row_number, row in enumerate(reader, start=2):
                for column, value in row.items():
                    if SEED_KEY.search(column or "") and value not in (None, ""):
                        seeds.append({
                            "field": f"{column}@row{row_number}",
                            "seed": _seed_int(value, source=relative, field=column),
                        })
    else:
        raise ValueError(f"{relative}: unsupported artifact format {suffix}")
    return {"path": relative, "bytes": path.stat().st_size, "sha256": digest, "seeds": seeds}


def inventory_roots(roots: list[Path]) -> tuple[list[dict[str, Any]], list[str]]:
    records: list[dict[str, Any]] = []
    errors: list[str] = []
    seen: set[tuple[str, str]] = set()
    for root in roots:
        root = root.resolve()
        if not root.is_dir():
            errors.append(f"artifact root does not exist or is not a directory: {root}")
            continue
        files = sorted(path for path in root.rglob("*") if path.is_file())
        for path in files:
            if path.suffix.lower() in {".yaml", ".yml"}:
                errors.append(f"unsupported YAML artifact requires explicit conversion/review: {path}")
                continue
            if path.suffix.lower() not in SUPPORTED:
                continue
            record = inventory_file(path, root)
            identity = (record["sha256"], record["path"])
            if identity not in seen:
                records.append({"root": str(root), **record})
                seen.add(identity)
    return records, errors


def _schedule_pairs(schedule: dict[str, Any], block: str) -> list[dict[str, int]]:
    rows = schedule.get(block)
    if not isinstance(rows, list) or not rows:
        raise ValueError(f"schedule.{block} must be a non-empty list of seed-pair objects")
    pairs = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict) or not {"population_seed", "validation_seed"} <= row.keys():
            raise ValueError(f"schedule.{block}[{index}] must contain population_seed and validation_seed")
        population = _seed_int(row["population_seed"], source=f"schedule.{block}", field="population_seed")
        validation = _seed_int(row["validation_seed"], source=f"schedule.{block}", field="validation_seed")
        if population == validation:
            raise ValueError(f"schedule.{block}[{index}] reuses population and validation seed")
        pairs.append({"population_seed": population, "validation_seed": validation})
    all_values = [value for pair in pairs for value in pair.values()]
    if len(all_values) != len(set(all_values)):
        raise ValueError(f"schedule.{block} contains a seed collision across pairs or streams")
    return pairs


def audit(roots: list[Path], schedule_path: Path) -> dict[str, Any]:
    schedule = json.loads(schedule_path.read_text(encoding="utf-8"))
    training = _schedule_pairs(schedule, "training")
    holdout = _schedule_pairs(schedule, "holdout")
    train_values = {value for pair in training for value in pair.values()}
    holdout_values = {value for pair in holdout for value in pair.values()}
    if train_values & holdout_values:
        raise ValueError("training and holdout schedules share one or more seed values")

    records, errors = inventory_roots(roots)
    historical: dict[int, list[dict[str, str]]] = {}
    for record in records:
        for item in record["seeds"]:
            historical.setdefault(item["seed"], []).append({
                "root": record["root"], "path": record["path"], "field": item["field"]
            })
    overlap_training = sorted(train_values & historical.keys())
    overlap_holdout = sorted(holdout_values & historical.keys())
    status = "PASS" if records and not errors and not overlap_training and not overlap_holdout else "NOT_EVALUABLE"
    return {
        "schema_version": 1,
        "status": status,
        "scope_warning": "Completeness is not established by scanning alone; the operator must supply every relevant historical artifact and document its provenance.",
        "source_files": records,
        "source_file_count": len(records),
        "historical_seed_count": len(historical),
        "historical_seed_values": sorted(historical),
        "historical_seed_occurrences": historical,
        "training_pairs": training,
        "holdout_pairs": holdout,
        "training_historical_overlap": overlap_training,
        "holdout_historical_overlap": overlap_holdout,
        "inventory_errors": errors,
        "schedule_sha256": sha256_file(schedule_path),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-root", action="append", required=True, type=Path)
    parser.add_argument("--schedule", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        report = audit(args.artifact_root, args.schedule)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        report = {"schema_version": 1, "status": "NOT_EVALUABLE", "fatal_error": str(exc)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\\n", encoding="utf-8")
    print(f"{report['status']}: {args.output}")
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
