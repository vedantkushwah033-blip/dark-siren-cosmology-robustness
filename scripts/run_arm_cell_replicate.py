"""Restartable coordinator for independent arm × validation-size cells.

This coordinator is scoped to one frozen replicate/gamma shared-input bundle.
It creates one lifecycle row per requested (validation_size, arm), claims cells
atomically, runs one worker and one independent auditor per cell, and persists
only that cell's terminal outcome.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd

from dark_siren_robustness.cell_ledger_store import (
    finalize_interrupted_cells,
    initialize_checkpoint,
    mark_cell_started,
    read_checkpoint,
    record_cell_outcome,
)
from dark_siren_robustness.shared_input_bundle import load_shared_input_manifest
from dark_siren_robustness.arm_cell_execution import ARM_NAMES


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="microseconds")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _cell_id(row: dict) -> str:
    gamma = f"{float(row['gamma']):+g}".replace("+", "p").replace("-", "m")
    return f"rep{int(row['replicate']):04d}_g{gamma}_n{int(row['validation_size'])}_{row['arm']}"


def _execute_one(task: dict) -> dict:
    root = Path(task["input_dir"])
    output = Path(task["output_dir"])
    checkpoint = Path(task["checkpoint_path"])
    row = task["plan_row"]
    key = (row["replicate"], row["gamma"], row["validation_size"], row["arm"])
    started = _now()
    mark_cell_started(checkpoint, key, started)

    cell_root = output / "cells" / row["cell_id"]
    cell_root.mkdir(parents=True, exist_ok=True)
    prior_attempts = sorted(cell_root.glob("attempt_*"))
    attempt_number = len(prior_attempts) + 1
    attempt_dir = cell_root / f"attempt_{attempt_number:04d}"
    attempt_dir.mkdir()
    attempt_record = {
        "cell_id": row["cell_id"], "attempt": attempt_number,
        "attempt_started_at": started, "first_cell_started_at": started,
        "worker_status": "NOT_RUN", "audit_status": "NOT_RUN",
    }
    (attempt_dir / "attempt_started.json").write_text(
        json.dumps(attempt_record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    scripts = Path(__file__).resolve().parent
    worker = scripts / "run_arm_cell.py"
    auditor = scripts / "audit_arm_cell.py"
    worker_cmd = [
        sys.executable, str(worker), "--input-dir", str(root),
        "--output-dir", str(attempt_dir / "result"),
        "--replicate", str(row["replicate"]),
        "--validation-size", str(row["validation_size"]),
        "--arm", str(row["arm"]),
    ]
    audit_cmd = [
        sys.executable, str(auditor), "--input-dir", str(root),
        "--cell-dir", str(attempt_dir / "result"),
        "--output", str(attempt_dir / "audit_report.json"),
    ]
    audit_status = "NOT_RUN"
    reason_code = "CHILD_PROCESS_ERROR"
    reason_detail = ""
    try:
        worker_result = subprocess.run(worker_cmd, capture_output=True, text=True)
        (attempt_dir / "worker.stdout.txt").write_text(worker_result.stdout, encoding="utf-8")
        (attempt_dir / "worker.stderr.txt").write_text(worker_result.stderr, encoding="utf-8")
        if worker_result.returncode != 0:
            error = worker_result.stderr[-5000:] or worker_result.stdout[-5000:]
            if any(token in error.lower() for token in ("sha-256", "manifest", "provenance", "byte count")):
                reason_code = "DATA_INTEGRITY_FAILURE"
            raise RuntimeError(f"cell worker failed with exit code {worker_result.returncode}: {error}")

        audit_result = subprocess.run(audit_cmd, capture_output=True, text=True)
        (attempt_dir / "audit.stdout.txt").write_text(audit_result.stdout, encoding="utf-8")
        (attempt_dir / "audit.stderr.txt").write_text(audit_result.stderr, encoding="utf-8")
        if audit_result.returncode != 0:
            audit_status = "NOT_PASS"
            reason_code = "CHILD_AUDIT_FAILED"
            error = audit_result.stderr[-5000:] or audit_result.stdout[-5000:]
            raise RuntimeError(f"independent cell audit failed: {error}")
        report = json.loads((attempt_dir / "audit_report.json").read_text(encoding="utf-8"))
        if report.get("status") != "PASS":
            audit_status = "NOT_PASS"
            reason_code = "CHILD_AUDIT_FAILED"
            raise RuntimeError("independent cell auditor did not report PASS")
        audit_status = "PASS"
        completed = _now()
        attempt_record.update({
            "attempt_completed_at": completed,
            "worker_status": "SUCCESS", "audit_status": "PASS",
            "audit_report": "audit_report.json",
        })
        (attempt_dir / "attempt_result.json").write_text(
            json.dumps(attempt_record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        report_path = str((attempt_dir / "audit_report.json").relative_to(output))
        record_cell_outcome(
            checkpoint, key, status="SUCCESS", audit_status="PASS",
            completed_at=completed, audit_report_path=report_path,
        )
        return {"cell_id": row["cell_id"], "status": "SUCCESS", "audit_status": "PASS",
                "attempt": attempt_number, "audit_report_path": report_path}
    except Exception as exc:
        completed = _now()
        reason_detail = f"{type(exc).__name__}: {exc}"[:6000]
        attempt_record.update({
            "attempt_completed_at": completed, "worker_status": "FAILED",
            "audit_status": audit_status, "reason_code": reason_code,
            "reason_detail": reason_detail,
        })
        (attempt_dir / "attempt_result.json").write_text(
            json.dumps(attempt_record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        record_cell_outcome(
            checkpoint, key, status="FAILED", audit_status=audit_status,
            reason_code=reason_code, reason_detail=reason_detail,
            completed_at=completed,
        )
        return {"cell_id": row["cell_id"], "status": "FAILED",
                "audit_status": audit_status, "attempt": attempt_number,
                "reason_code": reason_code, "reason_detail": reason_detail}


def run(args: argparse.Namespace) -> dict:
    root, output = args.input_dir.resolve(), args.output_dir.resolve()
    config = json.loads((root / "config.json").read_text(encoding="utf-8"))
    shared = load_shared_input_manifest(root, expected_validation_sizes=config["validation_sizes"])
    if config.get("status") != "shared_inputs_only_no_inference":
        raise ValueError("coordinator requires an inference-free shared-input bundle")
    if args.replicate < 0 or args.max_workers < 1:
        raise ValueError("replicate must be nonnegative and max_workers positive")
    arms = list(args.arms)
    if not arms or len(set(arms)) != len(arms) or any(arm not in ARM_NAMES for arm in arms):
        raise ValueError("arms must be a nonempty unique subset of the five declared arms")

    output.mkdir(parents=True, exist_ok=True)
    checkpoint_path = output / "cell_lifecycle_checkpoint.json"
    input_manifest_hash = _sha256(root / "shared_input_manifest.json")
    config_hash = _sha256(root / "config.json")
    provenance = shared["provenance"]
    planned = []
    for nval in config["validation_sizes"]:
        for arm in arms:
            row = {
                "cell_id": None,
                "replicate": int(args.replicate),
                "gamma": float(config["curvature_gamma"]),
                "validation_size": int(nval),
                "arm": arm,
                "input_dir": str(root),
                "shared_input_manifest_sha256": input_manifest_hash,
                "config_sha256": config_hash,
                "analysis_seed": int(provenance["seed"]),
                "event_seed": provenance.get("event_seed"),
                "validation_seed": int(provenance["validation_seed"]),
                "validation_sample_seed": provenance.get("validation_sample_seed"),
                "h0_true": float(config["h0_true"]),
                "inclusion_model_degree": int(config["inclusion_model_degree"]),
            }
            row["cell_id"] = _cell_id(row)
            planned.append(row)

    if args.resume:
        if not checkpoint_path.is_file():
            raise ValueError("--resume requires an existing checkpoint")
        existing = read_checkpoint(checkpoint_path)
        if existing["planned_cells"] != planned:
            raise ValueError("resume plan differs from immutable checkpoint plan")
        finalize_interrupted_cells(checkpoint_path, "Coordinator restart recovered unfinished cells.")
    else:
        if any(path.name != ".gitkeep" for path in output.iterdir()):
            raise ValueError("output directory must be empty unless --resume is specified")
        initialize_checkpoint(checkpoint_path, planned)

    checkpoint = read_checkpoint(checkpoint_path)
    terminal = {"SUCCESS", "NON_EVALUABLE", "FAILED"}
    pending_keys = {
        (row["replicate"], row["gamma"], row["validation_size"], row["arm"])
        for row in checkpoint["cell_ledger"] if row["status"] not in terminal
    }
    tasks = []
    for row in planned:
        key = (row["replicate"], row["gamma"], row["validation_size"], row["arm"])
        if key in pending_keys:
            tasks.append({
                "input_dir": str(root), "output_dir": str(output),
                "checkpoint_path": str(checkpoint_path), "plan_row": row,
            })
    outcomes = []
    if args.max_workers == 1:
        outcomes = [_execute_one(task) for task in tasks]
    elif tasks:
        with ProcessPoolExecutor(max_workers=args.max_workers) as pool:
            futures = [pool.submit(_execute_one, task) for task in tasks]
            for future in as_completed(futures):
                outcomes.append(future.result())

    checkpoint = read_checkpoint(checkpoint_path)
    ledger = checkpoint["cell_ledger"]
    ledger_frame = pd.DataFrame(ledger)
    summaries = []
    for row in ledger:
        if row["status"] != "SUCCESS":
            continue
        report_path = output / row["audit_report_path"]
        cell_dir = report_path.parent / "result"
        summary_path = cell_dir / "summary.json"
        if not summary_path.is_file():
            raise RuntimeError(f"audited cell summary missing: {summary_path}")
        summaries.append(json.loads(summary_path.read_text(encoding="utf-8")))
    summary_frame = pd.DataFrame(summaries, columns=[
        "replicate", "gamma", "validation_size", "arm", "status", "h0_true",
        "median_error", "posterior_median", "q025", "q05", "q16", "q84", "q95", "q975",
        "weight_diagnostics", "realized_bandwidths", "input_sha256",
    ])
    ledger_frame.to_csv(output / "cell_ledger.csv", index=False)
    summary_frame.to_csv(output / "cell_summaries.csv", index=False)
    status_counts = {
        status: int((ledger_frame["status"] == status).sum())
        for status in ("FAILED", "NON_EVALUABLE", "NOT_RUN", "SUCCESS")
    }
    run_status = "PASS" if len(ledger) and (ledger_frame["status"] == "SUCCESS").all() and (ledger_frame["audit_status"] == "PASS").all() else "INCOMPLETE"
    run_summary = {
        "status": run_status, "replicate": int(args.replicate),
        "gamma": float(config["curvature_gamma"]),
        "planned_cells": len(planned), "status_counts": status_counts,
        "completed_this_invocation": len(outcomes),
        "successful_this_invocation": sum(item["status"] == "SUCCESS" for item in outcomes),
        "failed_this_invocation": sum(item["status"] == "FAILED" for item in outcomes),
        "max_workers": int(args.max_workers),
        "claim_limit": "Single replicate/gamma cell orchestration; not confirmatory calibration.",
    }
    (output / "run_summary.json").write_text(
        json.dumps(run_summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    files = sorted(
        path for path in output.rglob("*")
        if path.is_file() and path.name not in {"run_manifest.json"}
        and not path.name.endswith(".lock")
    )
    run_manifest = {
        "status": run_status, "file_count": len(files),
        "files": {
            str(path.relative_to(output)): {
                "bytes": path.stat().st_size, "sha256": _sha256(path),
            } for path in files
        },
    }
    (output / "run_manifest.json").write_text(
        json.dumps(run_manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return run_summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--arms", type=lambda value: [item.strip() for item in value.split(",")],
                        default=list(ARM_NAMES))
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    print(json.dumps(run(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
