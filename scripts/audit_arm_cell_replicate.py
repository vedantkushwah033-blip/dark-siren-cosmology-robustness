"""Independent aggregate audit for the restartable arm-cell coordinator."""
from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

STATUSES = {"SUCCESS", "NON_EVALUABLE", "FAILED", "NOT_RUN"}
AUDITS = {"PASS", "NOT_PASS", "NOT_RUN"}
REASONS = {
    "SUCCESS": {""},
    "NON_EVALUABLE": {
        "NONFINITE_LIKELIHOOD", "INVALID_POSTERIOR_NORMALIZATION",
        "SUPPORT_BOUNDARY_MODE", "PREDECLARED_WEIGHT_INSTABILITY",
    },
    "FAILED": {"CHILD_PROCESS_ERROR", "CHILD_AUDIT_FAILED", "DATA_INTEGRITY_FAILURE"},
    "NOT_RUN": {"RUN_NOT_STARTED", "RUN_INTERRUPTED", "UPSTREAM_FAILURE", "NOT_SCHEDULED"},
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical_hash(value) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _timestamp(value: str | None, label: str):
    if value is None:
        return None
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise AssertionError(f"{label} lacks an explicit UTC offset")
    return parsed


def audit(output_dir: Path, input_dir: Path) -> dict:
    output, root = output_dir.resolve(), input_dir.resolve()
    checkpoint_path = output / "cell_lifecycle_checkpoint.json"
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    plan, ledger = checkpoint.get("planned_cells"), checkpoint.get("cell_ledger")
    if not isinstance(plan, list) or not isinstance(ledger, list) or not plan:
        raise AssertionError("checkpoint must contain a nonempty plan and ledger")
    if checkpoint.get("schema_version") != 1 or checkpoint.get("plan_sha256") != _canonical_hash(plan):
        raise AssertionError("checkpoint schema or immutable plan hash mismatch")
    key_fields = ("replicate", "gamma", "validation_size", "arm")
    plan_keys = [tuple(row[field] for field in key_fields) for row in plan]
    ledger_keys = [tuple(row[field] for field in key_fields) for row in ledger]
    if len(set(plan_keys)) != len(plan_keys) or len(set(ledger_keys)) != len(ledger_keys):
        raise AssertionError("plan or ledger contains duplicate cell keys")
    if set(plan_keys) != set(ledger_keys):
        raise AssertionError("ledger keys do not exactly cover the immutable plan")
    planned_by_key = {tuple(row[field] for field in key_fields): row for row in plan}
    counts = {status: 0 for status in sorted(STATUSES)}
    for row in ledger:
        key = tuple(row[field] for field in key_fields)
        planned = planned_by_key[key]
        for field, value in planned.items():
            if row.get(field) != value:
                raise AssertionError(f"ledger provenance mismatch for {field}")
        status, audit_status, reason = row.get("status"), row.get("audit_status"), row.get("reason_code")
        if status not in STATUSES or audit_status not in AUDITS or reason not in REASONS[status]:
            raise AssertionError("invalid lifecycle status, audit status, or reason code")
        if status in {"SUCCESS", "NON_EVALUABLE"} and audit_status != "PASS":
            raise AssertionError("audited terminal status lacks audit PASS")
        if status == "FAILED" and audit_status == "PASS":
            raise AssertionError("FAILED row cannot claim audit PASS")
        if status == "NOT_RUN" and audit_status != "NOT_RUN":
            raise AssertionError("NOT_RUN row must retain audit NOT_RUN")
        if status == "SUCCESS" and reason:
            raise AssertionError("SUCCESS row has a failure reason")
        if status != "SUCCESS" and not reason:
            raise AssertionError("non-success row lacks a reason")
        started = _timestamp(row.get("started_at"), "started_at")
        completed = _timestamp(row.get("completed_at"), "completed_at")
        if status in {"SUCCESS", "NON_EVALUABLE", "FAILED"} and (started is None or completed is None):
            raise AssertionError("terminal row lacks lifecycle timestamps")
        if started is not None and completed is not None and completed < started:
            raise AssertionError("cell completion precedes its start")
        if status in {"SUCCESS", "NON_EVALUABLE"}:
            path = row.get("audit_report_path")
            if not isinstance(path, str) or not path:
                raise AssertionError("audited cell has no audit report path")
            report_path = output / path
            report = json.loads(report_path.read_text(encoding="utf-8"))
            if report.get("status") != "PASS":
                raise AssertionError("cell audit report is not PASS")
            if tuple(report.get(field) for field in ("replicate", "gamma", "validation_size", "arm")) != key:
                raise AssertionError("cell audit report identity differs from ledger key")
        counts[status] += 1

    # The current project has no approved arm-specific NON_EVALUABLE thresholds.
    if counts["NON_EVALUABLE"]:
        raise AssertionError("NON_EVALUABLE rows are prohibited until their policy is preregistered")
    if pd.read_csv(output / "cell_ledger.csv").shape[0] != len(ledger):
        raise AssertionError("CSV ledger export row count differs from checkpoint")
    summaries = pd.read_csv(output / "cell_summaries.csv")
    if len(summaries) != counts["SUCCESS"]:
        raise AssertionError("cell summary row count differs from audited SUCCESS count")
    if set(summaries["status"].astype(str)) != ({"SUCCESS"} if counts["SUCCESS"] else set()):
        raise AssertionError("cell summary contains a non-success status")
    for row in ledger:
        if row["status"] == "SUCCESS":
            match = summaries.loc[
                (summaries["replicate"] == row["replicate"])
                & (summaries["gamma"] == row["gamma"])
                & (summaries["validation_size"] == row["validation_size"])
                & (summaries["arm"] == row["arm"])
            ]
            if len(match) != 1:
                raise AssertionError("SUCCESS ledger row does not have exactly one summary")

    run_summary = json.loads((output / "run_summary.json").read_text(encoding="utf-8"))
    expected_status = "PASS" if counts["SUCCESS"] == len(plan) else "INCOMPLETE"
    if run_summary.get("status") != expected_status or run_summary.get("status_counts") != counts:
        raise AssertionError("run summary does not reconcile to the checkpoint lifecycle counts")

    manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    if manifest.get("status") != expected_status:
        raise AssertionError("run manifest status differs from lifecycle outcomes")
    listed = set(manifest.get("files", {}))
    actual = {
        str(path.relative_to(output)) for path in output.rglob("*")
        if path.is_file() and path.name != "run_manifest.json" and not path.name.endswith(".lock")
    }
    if listed != actual:
        raise AssertionError("run manifest file set differs from actual run outputs")
    for relative, metadata in manifest["files"].items():
        path = output / relative
        if path.stat().st_size != int(metadata["bytes"]) or _sha256(path) != metadata["sha256"]:
            raise AssertionError(f"run manifest hash mismatch: {relative}")

    # Re-run the independent cell auditor for every SUCCESS, not merely trust
    # the cell report created by the coordinator.
    auditor = Path(__file__).resolve().with_name("audit_arm_cell.py")
    audited_again = 0
    for row in ledger:
        if row["status"] != "SUCCESS":
            continue
        report_path = output / row["audit_report_path"]
        cell_dir = report_path.parent / "result"
        command = [
            sys.executable, str(auditor), "--input-dir", str(root),
            "--cell-dir", str(cell_dir),
        ]
        check = subprocess.run(command, capture_output=True, text=True)
        if check.returncode != 0:
            raise AssertionError(
                "repeat independent cell reconstruction failed: "
                + (check.stderr[-4000:] or check.stdout[-4000:])
            )
        repeated = json.loads(check.stdout)
        if repeated.get("status") != "PASS":
            raise AssertionError("repeat independent cell audit did not PASS")
        audited_again += 1

    return {
        "status": "PASS",
        "planned_cells": len(plan),
        "status_counts": counts,
        "successful_cells_reaudited": audited_again,
        "run_manifest_files_verified": len(listed),
        "claim_limit": "Lifecycle and independent per-cell artifact audit; not calibration.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    report = audit(args.output_dir, args.input_dir)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
