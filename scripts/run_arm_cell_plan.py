"""Run a fully planned arm-cell design sequentially with audit-gated checkpoints.

Plan JSON schema:
{"schema_version": 1, "cells": [
  {"replicate": 1, "gamma": 0.0, "validation_size": 100,
   "arm": "A_catalogue_only", "bundle_dir": "bundles/r1_g0",
   "output_dir": "cells/r1_g0_n100_A"}
]}

All cells are checkpointed as NOT_RUN before the first cell executes. This
coordinator is deliberately serial. Use --recover-interrupted only after the
previous coordinator process has stopped; it finalizes every started but
unfinished cell in the checkpoint before resuming.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path

from dark_siren_robustness.cell_ledger_store import (
    initialize_checkpoint,
    read_checkpoint,
    finalize_interrupted_cells,
    requeue_not_scheduled_cells,
)
from dark_siren_robustness.lifecycle_ledger import validate_ledger
from dark_siren_robustness.arm_cell_execution import ARM_NAMES


def _load_runner():
    path = Path(__file__).with_name("run_audited_arm_cell.py")
    spec = importlib.util.spec_from_file_location("_audit_gated_arm_cell_runner", path)
    if spec is None or spec.loader is None:
        raise ImportError("could not load audit-gated arm-cell runner")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.run_audited_cell


run_audited_cell = _load_runner()


def _paths_overlap(left: str, right: str) -> bool:
    """Return whether two filesystem paths are equal or one contains the other."""
    a = Path(left).resolve()
    b = Path(right).resolve()
    return a == b or a in b.parents or b in a.parents


def _key(row: dict) -> tuple:
    return (int(row["replicate"]), float(row["gamma"]),
            int(row["validation_size"]), str(row["arm"]))


def _normalize_plan(path: Path) -> list[dict]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if raw.get("schema_version") != 1 or not isinstance(raw.get("cells"), list):
        raise ValueError("plan must have schema_version=1 and a cells list")
    rows = []
    seen = set()
    output_paths = set()
    audit_paths = set()
    for source in raw["cells"]:
        row = dict(source)
        required = {"replicate", "gamma", "validation_size", "arm", "bundle_dir", "output_dir"}
        missing = required - set(row)
        if missing:
            raise ValueError(f"planned cell is missing fields: {sorted(missing)}")
        row["replicate"] = int(row["replicate"])
        row["gamma"] = float(row["gamma"])
        row["validation_size"] = int(row["validation_size"])
        row["arm"] = str(row["arm"])
        key = _key(row)
        if key in seen:
            raise ValueError(f"duplicate cell in plan: {key}")
        seen.add(key)
        if row["replicate"] < 1 or row["validation_size"] < 5:
            raise ValueError(f"replicate and validation_size must be positive/plausible: {key}")
        if not math.isfinite(row["gamma"]):
            raise ValueError(f"gamma must be finite: {key}")
        if row["arm"] not in ARM_NAMES:
            raise ValueError(f"unknown arm in plan: {row['arm']}")
        for field in ("bundle_dir", "output_dir", "audit_report_path"):
            if row.get(field) is not None:
                candidate = Path(row[field])
                if not candidate.is_absolute():
                    candidate = path.parent / candidate
                row[field] = str(candidate.resolve())
        output_path = row["output_dir"]
        if output_path in output_paths:
            raise ValueError(f"multiple cells cannot share an output_dir: {output_path}")
        output_paths.add(output_path)
        audit_path = row.get("audit_report_path") or output_path + ".independent_audit.json"
        if audit_path in audit_paths:
            raise ValueError(f"multiple cells cannot share an audit_report_path: {audit_path}")
        audit_paths.add(audit_path)
        if row["bundle_dir"] == output_path:
            raise ValueError("bundle_dir and output_dir cannot be the same path")
        rows.append(row)
    if not rows:
        raise ValueError("plan must contain at least one cell")

    output_paths = [row["output_dir"] for row in rows]
    audit_paths = [
        row.get("audit_report_path") or row["output_dir"] + ".independent_audit.json"
        for row in rows
    ]
    bundle_paths = [row["bundle_dir"] for row in rows]
    for i, left in enumerate(output_paths):
        for right in output_paths[i + 1:]:
            if _paths_overlap(left, right):
                raise ValueError(f"cell output directories overlap: {left} and {right}")
    for i, left in enumerate(audit_paths):
        for right in audit_paths[i + 1:]:
            if _paths_overlap(left, right):
                raise ValueError(f"cell audit-report paths overlap: {left} and {right}")
    for output in output_paths:
        for audit_path in audit_paths:
            if _paths_overlap(output, audit_path):
                raise ValueError(f"cell output and audit-report paths overlap: {output} and {audit_path}")
        for bundle in bundle_paths:
            if _paths_overlap(output, bundle):
                raise ValueError(f"cell output and shared-bundle paths overlap: {output} and {bundle}")
    for audit_path in audit_paths:
        for bundle in bundle_paths:
            if _paths_overlap(audit_path, bundle):
                raise ValueError(f"cell audit-report and shared-bundle paths overlap: {audit_path} and {bundle}")
    return rows


def run_plan(plan_path: Path, checkpoint_path: Path, *, recover_interrupted: bool = False) -> dict:
    plan_path = Path(plan_path).resolve()
    checkpoint_path = Path(checkpoint_path).resolve()
    planned = _normalize_plan(plan_path)

    if checkpoint_path.exists():
        checkpoint = read_checkpoint(checkpoint_path)
        if checkpoint["planned_cells"] != planned:
            raise ValueError("checkpoint plan differs from supplied plan; immutable plan cannot change")
        if recover_interrupted:
            finalize_interrupted_cells(
                checkpoint_path,
                reason_detail="Coordinator restart after prior process stopped.",
            )
            requeue_not_scheduled_cells(
                checkpoint_path,
                reason_detail="Never-started cell requeued after coordinator restart.",
            )
    else:
        if recover_interrupted:
            raise ValueError("--recover-interrupted requires an existing checkpoint")
        initialize_checkpoint(checkpoint_path, planned)

    results = []
    for cell in planned:
        checkpoint = read_checkpoint(checkpoint_path)
        key = _key(cell)
        row = next(row for row in checkpoint["cell_ledger"] if _key(row) == key)
        if row["status"] in {"SUCCESS", "NON_EVALUABLE", "FAILED"}:
            results.append({"cell_key": list(key), "status": row["status"], "skipped": True})
            continue
        result = run_audited_cell(
            checkpoint_path,
            cell["bundle_dir"],
            cell["output_dir"],
            replicate=cell["replicate"],
            gamma=cell["gamma"],
            validation_size=cell["validation_size"],
            arm=cell["arm"],
            audit_report_path=cell.get("audit_report_path"),
        )
        results.append(result)

    final = read_checkpoint(checkpoint_path)
    audit = validate_ledger(final["planned_cells"], final["cell_ledger"])
    counts = audit["status_counts"]
    if counts["NOT_RUN"]:
        run_status = "INCOMPLETE"
    elif counts["FAILED"] or counts["NON_EVALUABLE"]:
        run_status = "COMPLETED_WITH_NON_SUCCESS"
    else:
        run_status = "PASS"
    return {
        "status": run_status,
        "ledger_validation_status": audit["status"],
        "checkpoint": str(checkpoint_path),
        "plan_sha256": final["plan_sha256"],
        "planned_cells": audit["planned_cells"],
        "status_counts": counts,
        "cell_results": results,
        "execution_mode": "SERIAL",
        "scope_limit": (
            "Sequential synthetic arm-cell execution only. SUCCESS requires independent audit PASS. "
            "No confirmatory calibration or multi-worker concurrency claim is implied."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument(
        "--recover-interrupted", action="store_true",
        help="After all previous workers have stopped, mark unfinished started cells interrupted before resuming.",
    )
    args = parser.parse_args()
    report = run_plan(args.plan, args.checkpoint, recover_interrupted=args.recover_interrupted)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
