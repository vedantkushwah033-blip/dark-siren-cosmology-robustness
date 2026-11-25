"""Independent aggregate audit of a four-state arm-cell checkpoint.

This module intentionally does not import the checkpoint store, lifecycle
ledger, coordinator, or cell evaluator. It reconstructs plan coverage, hashes,
state invariants, and successful-cell audit-report bindings from JSON.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

STATUSES = {"SUCCESS", "NON_EVALUABLE", "FAILED", "NOT_RUN"}
REASONS = {
    "SUCCESS": {""},
    "NON_EVALUABLE": {
        "NONFINITE_LIKELIHOOD", "INVALID_POSTERIOR_NORMALIZATION",
        "SUPPORT_BOUNDARY_MODE", "PREDECLARED_WEIGHT_INSTABILITY",
    },
    "FAILED": {"CHILD_PROCESS_ERROR", "CHILD_AUDIT_FAILED", "DATA_INTEGRITY_FAILURE"},
    "NOT_RUN": {"RUN_NOT_STARTED", "RUN_INTERRUPTED", "UPSTREAM_FAILURE", "NOT_SCHEDULED"},
}
AUDIT = {"PASS", "NOT_PASS", "NOT_RUN"}
KEY_FIELDS = ("replicate", "gamma", "validation_size", "arm")


def _key(row):
    if not isinstance(row, dict) or any(k not in row for k in KEY_FIELDS):
        raise ValueError("cell row has an incomplete key")
    return (int(row["replicate"]), float(row["gamma"]),
            int(row["validation_size"]), str(row["arm"]))


def _canonical_hash(value):
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _timestamp(value, label):
    if not isinstance(value, str) or not value:
        raise ValueError(f"{label} must be a non-empty ISO-8601 timestamp")
    from datetime import datetime
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError(f"{label} is not ISO-8601") from exc
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise ValueError(f"{label} must include a UTC offset")
    return parsed


def audit_checkpoint(checkpoint_path: Path) -> dict:
    path = Path(checkpoint_path).resolve()
    checkpoint = json.loads(path.read_text(encoding="utf-8"))
    if checkpoint.get("schema_version") != 1:
        raise ValueError("unsupported checkpoint schema")
    plan, ledger = checkpoint.get("planned_cells"), checkpoint.get("cell_ledger")
    if not isinstance(plan, list) or not plan or not isinstance(ledger, list):
        raise ValueError("checkpoint must contain non-empty plan and ledger lists")
    if checkpoint.get("plan_sha256") != _canonical_hash(plan):
        raise ValueError("immutable plan SHA-256 mismatch")

    plan_by_key = {}
    for cell in plan:
        key = _key(cell)
        if key in plan_by_key:
            raise ValueError(f"duplicate planned cell: {key}")
        plan_by_key[key] = cell
    ledger_by_key = {}
    for row in ledger:
        key = _key(row)
        if key in ledger_by_key:
            raise ValueError(f"duplicate ledger cell: {key}")
        ledger_by_key[key] = row
    if set(plan_by_key) != set(ledger_by_key):
        raise ValueError("ledger does not exactly cover the immutable plan")

    counts = {name: 0 for name in sorted(STATUSES)}
    verified_reports = []
    for key, row in ledger_by_key.items():
        planned = plan_by_key[key]
        for field, value in planned.items():
            if row.get(field) != value:
                raise ValueError(f"planned provenance mismatch for {field}: {key}")
        status, reason, audit = row.get("status"), row.get("reason_code"), row.get("audit_status")
        if status not in STATUSES or reason not in REASONS[status] or audit not in AUDIT:
            raise ValueError(f"invalid state/reason/audit tuple for {key}")
        started, completed = row.get("started_at"), row.get("completed_at")
        report_path = row.get("audit_report_path")
        if status == "NOT_RUN":
            if audit != "NOT_RUN" or completed is not None or report_path is not None:
                raise ValueError(f"NOT_RUN has terminal metadata: {key}")
            if reason == "RUN_INTERRUPTED" and started is None:
                raise ValueError(f"interrupted cell lacks first-start timestamp: {key}")
            if reason in {"UPSTREAM_FAILURE", "NOT_SCHEDULED"} and started is not None:
                raise ValueError(f"unscheduled cell has a start timestamp: {key}")
            if started is not None:
                _timestamp(started, "started_at")
        else:
            start_dt = _timestamp(started, "started_at")
            end_dt = _timestamp(completed, "completed_at")
            if end_dt < start_dt:
                raise ValueError(f"completion precedes start: {key}")
            if status in {"SUCCESS", "NON_EVALUABLE"}:
                if audit != "PASS" or not isinstance(report_path, str) or not report_path:
                    raise ValueError(f"{status} requires audit PASS and report path: {key}")
                resolved = Path(report_path)
                if not resolved.is_absolute():
                    resolved = path.parent / resolved
                resolved = resolved.resolve()
                if not resolved.is_file():
                    raise ValueError(f"independent audit report is missing: {resolved}")
                report = json.loads(resolved.read_text(encoding="utf-8"))
                if report.get("status") != "PASS" or report.get("audit_status") != "PASS":
                    raise ValueError(f"referenced audit report is not PASS: {key}")
                report_key = report.get("cell_key")
                expected_key = dict(zip(KEY_FIELDS, key))
                if report_key != expected_key:
                    raise ValueError(f"audit report is bound to a different cell: {key}")
                if status == "SUCCESS" and reason != "":
                    raise ValueError(f"SUCCESS has a reason code: {key}")
                if status == "NON_EVALUABLE" and not reason:
                    raise ValueError(f"NON_EVALUABLE lacks a reason: {key}")
                output_dir = Path(planned["output_dir"])
                if not output_dir.is_absolute():
                    output_dir = path.parent / output_dir
                if not output_dir.resolve().is_dir():
                    raise ValueError(f"audited cell output directory is missing: {key}")
                verified_reports.append({"cell_key": list(key), "path": str(resolved)})
            else:
                if audit == "PASS":
                    raise ValueError(f"FAILED cell cannot claim audit PASS: {key}")
                if report_path is not None:
                    raise ValueError(f"FAILED cell unexpectedly references an audit report: {key}")
        counts[status] += 1

    if counts["NOT_RUN"]:
        run_status = "INCOMPLETE"
    elif counts["FAILED"] or counts["NON_EVALUABLE"]:
        run_status = "COMPLETED_WITH_NON_SUCCESS"
    else:
        run_status = "PASS"
    return {
        "status": "PASS",
        "run_status": run_status,
        "plan_sha256": checkpoint["plan_sha256"],
        "planned_cells": len(plan),
        "status_counts": counts,
        "verified_audit_reports": verified_reports,
        "audit_scope": (
            "Independent structural audit of plan/checkpoint and referenced cell-audit reports. "
            "It does not rerun the numerical cell reconstruction; each referenced report "
            "must have been produced by the separate cell auditor."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    report = audit_checkpoint(args.checkpoint)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
