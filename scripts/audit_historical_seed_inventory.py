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

import numpy as np
from pathlib import Path
from typing import Any

SUPPORTED = {".json", ".csv", ".npz"}
SEED_KEY = re.compile(r"(?:^seed$|_seed$|^seeds$|_seeds$|^seed_start$|_seed_start$)", re.IGNORECASE)

def is_seed_field(name: Any) -> bool:
    """Match concrete seed fields, not aggregate metrics such as *_across_seeds."""
    normalized = str(name).lower()
    if (
        "_across_seeds" in normalized
        or normalized in {"n_seeds", "replicates_per_seed", "seed_count", "seeds_count"}
        or normalized.startswith(("n_", "num_", "total_")) and normalized.endswith("_seeds")
    ):
        return False
    return SEED_KEY.search(normalized) is not None



def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _seed_int(value: Any, *, source: str, field: str) -> int:
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, bool):
        raise ValueError(f"{source}: boolean seed in {field}")
    if isinstance(value, int):
        return value
    if isinstance(value, float) and np.isfinite(value) and value.is_integer():
        return int(value)
    if isinstance(value, str) and re.fullmatch(r"[+-]?\d+", value.strip()):
        return int(value.strip())
    raise ValueError(f"{source}: non-integer seed value in {field}: {value!r}")


def _walk_json(value: Any, source: str, prefix: str = "$") -> list[dict[str, Any]]:
    found: list[dict[str, Any]] = []
    if isinstance(value, dict):
        for key, child in value.items():
            field = f"{prefix}.{key}"
            # These fields describe deterministic seed-generation rules; they are
            # metadata, not consumed integer RNG seeds. Their concrete values are
            # inventoried from the master/stride fields and replicate ledgers.
            if str(key).lower() in {"seed_schedule", "seed_design", "seed_description"} and isinstance(child, str):
                continue
            if is_seed_field(key):
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


def _ledger_seed_reuse_warnings(path: Path, relative: str) -> list[dict[str, str]]:
    """Flag seed reuse across replicates within design strata.

    Repeated rows for nested validation sizes are collapsed. Scenario columns
    separate independent design strata, while gamma remains a stratum so both
    common-random-number and independently seeded gamma designs are admissible.
    """
    warnings: list[dict[str, str]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "replicate" not in reader.fieldnames:
            return warnings
        seed_columns = [
            column for column in reader.fieldnames
            if column and is_seed_field(column)
            and column.lower() not in {"seed", "seeds", "seed_start", "master_seed", "base_seed"}
        ]
        if len(seed_columns) < 2:
            return warnings
        context_columns = [
            column for column in (
                "experiment", "gamma", "catalogue_floor", "target_retention",
                "floor", "scenario", "scenario_id", "dgp", "dgp_id",
                "cell_type", "condition", "replicate_group",
            )
            if column in reader.fieldnames
        ]
        by_stratum_replicate: dict[tuple[tuple[str, str], ...], dict[str, dict[str, set[int]]]] = {}
        for row_number, row in enumerate(reader, start=2):
            replicate = row.get("replicate")
            if replicate in (None, ""):
                continue
            stratum = tuple((column, str(row.get(column, ""))) for column in context_columns)
            streams = by_stratum_replicate.setdefault(stratum, {}).setdefault(
                str(replicate), {column: set() for column in seed_columns}
            )
            for column in seed_columns:
                raw = row.get(column)
                if raw not in (None, ""):
                    streams[column].add(
                        _seed_int(raw, source=relative, field=f"{column}@row{row_number}")
                    )

    for stratum, by_replicate in by_stratum_replicate.items():
        stratum_label = ",".join(f"{key}={value}" for key, value in stratum) or "all_rows"
        stable_by_replicate: dict[str, dict[str, int]] = {}
        for replicate, streams in by_replicate.items():
            stable: dict[str, int] = {}
            for column, values in streams.items():
                if len(values) > 1:
                    warnings.append({
                        "path": relative,
                        "kind": "seed_stream_varies_within_replicate_stratum",
                        "detail": f"stratum={stratum_label} replicate={replicate} stream={column} has {len(values)} distinct values",
                    })
                elif len(values) == 1:
                    stable[column] = next(iter(values))
            stable_by_replicate[replicate] = stable

        for column in seed_columns:
            by_value: dict[int, list[str]] = {}
            for replicate, streams in stable_by_replicate.items():
                if column in streams:
                    by_value.setdefault(streams[column], []).append(replicate)
            for seed_value, replicate_ids in by_value.items():
                if len(replicate_ids) > 1:
                    warnings.append({
                        "path": relative,
                        "kind": f"{column}_reused_across_replicates",
                        "detail": f"stratum={stratum_label} seed={seed_value} appears in replicate IDs {','.join(sorted(replicate_ids))}",
                    })

        # A seed used by different stream columns in different replicates is
        # still a collision, even if no single column repeats it by itself.
        cross_stream_uses: dict[int, list[tuple[str, str]]] = {}
        for replicate, streams in stable_by_replicate.items():
            for column, seed_value in streams.items():
                cross_stream_uses.setdefault(seed_value, []).append((replicate, column))
        for seed_value, uses in cross_stream_uses.items():
            replicate_ids = sorted({replicate for replicate, _ in uses})
            stream_names = sorted({column for _, column in uses})
            if len(replicate_ids) > 1 and len(stream_names) > 1:
                warnings.append({
                    "path": relative,
                    "kind": "seed_reused_across_replicates",
                    "detail": f"stratum={stratum_label} seed={seed_value} is reused across streams {','.join(stream_names)} and replicate IDs {','.join(replicate_ids)}",
                })
    return warnings


def inventory_file(path: Path, root: Path) -> dict[str, Any]:
    suffix = path.suffix.lower()
    relative = path.relative_to(root).as_posix()
    digest = sha256_file(path)
    parse_warnings: list[dict[str, str]] = []
    if suffix == ".json":
        raw_text = path.read_text(encoding="utf-8")
        try:
            payload = json.loads(raw_text)
        except json.JSONDecodeError:
            # Repair only one known terminal serialization defect from older artifacts.
            stripped = raw_text.rstrip()
            if not stripped.endswith("\\n"):
                raise
            payload = json.loads(stripped[:-2])
            parse_warnings.append({
                "path": relative,
                "kind": "terminal_literal_backslash_n_removed",
                "detail": "Removed one literal backslash-n suffix after a valid JSON document; original bytes remain represented by sha256.",
            })
        seeds = _walk_json(payload, relative)
    elif suffix == ".csv":
        seeds = []
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise ValueError(f"{relative}: CSV has no header")
            for row_number, row in enumerate(reader, start=2):
                for column, value in row.items():
                    if is_seed_field(column or "") and value not in (None, ""):
                        seeds.append({
                            "field": f"{column}@row{row_number}",
                            "seed": _seed_int(value, source=relative, field=column),
                        })
    elif suffix == ".npz":
        seeds = []
        # Never enable pickle while inspecting research artifacts.
        with np.load(path, allow_pickle=False) as archive:
            for key in archive.files:
                if not is_seed_field(key):
                    continue
                values = np.asarray(archive[key]).reshape(-1)
                for index, value in enumerate(values):
                    seeds.append({
                        "field": f"{key}[{index}]",
                        "seed": _seed_int(value, source=relative, field=key),
                    })
    else:
        raise ValueError(f"{relative}: unsupported artifact format {suffix}")
    ledger_warnings = _ledger_seed_reuse_warnings(path, relative) if suffix == ".csv" else []
    return {"path": relative, "bytes": path.stat().st_size, "sha256": digest, "seeds": seeds, "ledger_warnings": ledger_warnings, "parse_warnings": parse_warnings}


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
        raise ValueError(f"schedule.{block} must be a non-empty list of seed-stream objects")
    pairs: list[dict[str, int]] = []
    for index, row in enumerate(rows):
        if not isinstance(row, dict) or not {"population_seed", "validation_seed"} <= row.keys():
            raise ValueError(f"schedule.{block}[{index}] must contain population_seed and validation_seed")
        seed_fields = [key for key in row if is_seed_field(key)]
        if set(seed_fields) != set(row):
            raise ValueError(f"schedule.{block}[{index}] contains non-seed fields; use explicit *_seed stream names only")
        values: dict[str, int] = {}
        for key in seed_fields:
            raw = row[key]
            if isinstance(raw, (list, tuple)):
                if not raw:
                    raise ValueError(f"schedule.{block}[{index}].{key} must not be empty")
                for stream_index, item in enumerate(raw):
                    values[f"{key}[{stream_index}]"] = _seed_int(
                        item, source=f"schedule.{block}", field=f"{key}[{stream_index}]"
                    )
            else:
                values[str(key)] = _seed_int(
                    raw, source=f"schedule.{block}", field=str(key)
                )
        if len(values.values()) != len(set(values.values())):
            raise ValueError(f"schedule.{block}[{index}] reuses population and validation seed values across random streams")
        pairs.append(values)
    all_values = [value for pair in pairs for value in pair.values()]
    if len(all_values) != len(set(all_values)):
        raise ValueError(f"schedule.{block} contains a seed collision across replicates or streams")
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
    ledger_warnings: list[dict[str, str]] = []
    inventory_warnings: list[dict[str, str]] = []
    for record in records:
        ledger_warnings.extend(record.get("ledger_warnings", []))
        inventory_warnings.extend(record.get("parse_warnings", []))
        for item in record["seeds"]:
            historical.setdefault(item["seed"], []).append({
                "root": record["root"], "path": record["path"], "field": item["field"]
            })
    overlap_training = sorted(train_values & historical.keys())
    overlap_holdout = sorted(holdout_values & historical.keys())
    # Historical seed reuse prevents certification of a disjoint candidate partition.
    status = (
        "PASS"
        if records and historical and not errors and not ledger_warnings
        and not overlap_training and not overlap_holdout
        else "NOT_EVALUABLE"
    )
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
        "historical_seed_warnings": ledger_warnings,
        "inventory_warnings": inventory_warnings,
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
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"{report['status']}: {args.output}")
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
