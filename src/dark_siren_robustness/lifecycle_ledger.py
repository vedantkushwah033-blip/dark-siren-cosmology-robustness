"""Strict cell-level lifecycle ledger primitives for confirmatory experiments.

This module is deliberately independent of the current replicate runner. It
provides a tested contract that the runner and independent auditor can adopt
without silently changing legacy artifact semantics.
"""
from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Iterable, Mapping

STATUSES = frozenset({"SUCCESS", "NON_EVALUABLE", "FAILED", "NOT_RUN"})
AUDIT_STATUSES = frozenset({"PASS", "NOT_PASS", "NOT_RUN"})
CELL_KEY_FIELDS = ("replicate", "gamma", "validation_size", "arm")

REASONS_BY_STATUS = {
    "SUCCESS": frozenset({""}),
    "NON_EVALUABLE": frozenset({
        "NONFINITE_LIKELIHOOD",
        "INVALID_POSTERIOR_NORMALIZATION",
        "SUPPORT_BOUNDARY_MODE",
        "PREDECLARED_WEIGHT_INSTABILITY",
    }),
    "FAILED": frozenset({
        "CHILD_PROCESS_ERROR",
        "CHILD_AUDIT_FAILED",
        "DATA_INTEGRITY_FAILURE",
    }),
    "NOT_RUN": frozenset({
        "RUN_NOT_STARTED",
        "RUN_INTERRUPTED",
        "UPSTREAM_FAILURE",
        "NOT_SCHEDULED",
    }),
}


def _key(row: Mapping) -> tuple:
    missing = [field for field in CELL_KEY_FIELDS if field not in row]
    if missing:
        raise ValueError(f"cell row is missing key fields: {missing}")
    return tuple(row[field] for field in CELL_KEY_FIELDS)


def initial_ledger(planned_rows: Iterable[Mapping]) -> list[dict]:
    """Create one NOT_RUN ledger row for every planned cell, before execution."""
    rows = []
    seen = set()
    for planned in planned_rows:
        row = dict(planned)
        key = _key(row)
        if key in seen:
            raise ValueError(f"duplicate planned cell: {key}")
        seen.add(key)
        row.update({
            "status": "NOT_RUN",
            "reason_code": "RUN_NOT_STARTED",
            "reason_detail": "",
            "audit_status": "NOT_RUN",
            "started_at": None,
            "completed_at": None,
            "audit_report_path": None,
        })
        rows.append(row)
    if not rows:
        raise ValueError("planned cell set cannot be empty")
    return rows


def set_cell_outcome(
    row: Mapping,
    *,
    status: str,
    audit_status: str,
    reason_code: str = "",
    reason_detail: str = "",
    started_at: str | None = None,
    completed_at: str | None = None,
    audit_report_path: str | None = None,
) -> dict:
    """Return a validated outcome row; does not mutate the input mapping."""
    if status not in STATUSES:
        raise ValueError(f"unknown lifecycle status: {status}")
    if audit_status not in AUDIT_STATUSES:
        raise ValueError(f"unknown audit status: {audit_status}")
    if reason_code not in REASONS_BY_STATUS[status]:
        raise ValueError(f"reason {reason_code!r} is not allowed for {status}")
    if status in {"SUCCESS", "NON_EVALUABLE"} and audit_status != "PASS":
        raise ValueError(f"{status} requires an independent audit PASS")
    if status == "FAILED" and audit_status == "PASS":
        raise ValueError("FAILED cannot claim audit PASS")
    if status == "NOT_RUN" and audit_status != "NOT_RUN":
        raise ValueError("NOT_RUN requires audit status NOT_RUN")
    if status == "SUCCESS" and reason_code:
        raise ValueError("SUCCESS must have an empty reason code")
    if status != "SUCCESS" and not reason_code:
        raise ValueError(f"{status} requires a reason code")
    if status in {"SUCCESS", "NON_EVALUABLE", "FAILED"} and not completed_at:
        raise ValueError(f"{status} requires completed_at")
    if status in {"SUCCESS", "NON_EVALUABLE", "FAILED"} and not started_at:
        raise ValueError(f"{status} requires started_at")
    if status in {"SUCCESS", "NON_EVALUABLE"} and (
        not isinstance(audit_report_path, str) or not audit_report_path.strip()
    ):
        raise ValueError(f"{status} requires audit_report_path")
    if started_at and completed_at and completed_at < started_at:
        raise ValueError("completed_at cannot precede started_at")
    result = dict(row)
    result.update({
        "status": status,
        "reason_code": reason_code,
        "reason_detail": reason_detail,
        "audit_status": audit_status,
        "started_at": started_at,
        "completed_at": completed_at,
        "audit_report_path": audit_report_path,
    })
    return result


def validate_ledger(planned_rows: Iterable[Mapping], ledger_rows: Iterable[Mapping]) -> dict:
    """Validate exact plan/ledger key coverage and all lifecycle invariants."""
    planned = list(planned_rows)
    ledger = list(ledger_rows)
    planned_keys = [_key(row) for row in planned]
    ledger_keys = [_key(row) for row in ledger]
    if len(set(planned_keys)) != len(planned_keys):
        raise ValueError("planned cells contain duplicate keys")
    if len(set(ledger_keys)) != len(ledger_keys):
        raise ValueError("ledger contains duplicate cell keys")
    if set(planned_keys) != set(ledger_keys):
        missing = set(planned_keys) - set(ledger_keys)
        extra = set(ledger_keys) - set(planned_keys)
        raise ValueError(f"ledger/plan key mismatch: missing={len(missing)}, extra={len(extra)}")
    counts = {status: 0 for status in sorted(STATUSES)}
    for row in ledger:
        status = row.get("status")
        if status not in STATUSES:
            raise ValueError(f"unknown lifecycle status: {status}")
        audit = row.get("audit_status")
        if audit not in AUDIT_STATUSES:
            raise ValueError(f"unknown audit status: {audit}")
        reason = row.get("reason_code", "")
        if reason not in REASONS_BY_STATUS[status]:
            raise ValueError(f"invalid reason {reason!r} for status {status}")
        if status in {"SUCCESS", "NON_EVALUABLE"} and audit != "PASS":
            raise ValueError(f"{status} cell lacks audit PASS")
        if status == "FAILED" and audit == "PASS":
            raise ValueError("FAILED cell cannot claim audit PASS")
        if status == "NOT_RUN" and audit != "NOT_RUN":
            raise ValueError("NOT_RUN cell has an audit status other than NOT_RUN")
        if status != "SUCCESS" and not reason:
            raise ValueError(f"{status} cell has no reason code")
        if status == "SUCCESS" and reason:
            raise ValueError("SUCCESS cell has a failure reason")
        if status in {"SUCCESS", "NON_EVALUABLE", "FAILED"} and not row.get("started_at"):
            raise ValueError(f"{status} cell has no start timestamp")
        if status in {"SUCCESS", "NON_EVALUABLE", "FAILED"} and not row.get("completed_at"):
            raise ValueError(f"{status} cell has no completion timestamp")
        audit_path = row.get("audit_report_path")
        if status in {"SUCCESS", "NON_EVALUABLE"} and (
            not isinstance(audit_path, str) or not audit_path.strip()
        ):
            raise ValueError(f"{status} cell has no audit report path")
        started_at = row.get("started_at")
        completed_at = row.get("completed_at")
        if started_at and completed_at and completed_at < started_at:
            raise ValueError("completed_at cannot precede started_at")
        counts[status] += 1
    return {"status": "PASS", "planned_cells": len(planned), "ledger_cells": len(ledger),
            "status_counts": counts}


def atomic_write_text(path: str | Path, content: str) -> None:
    """Durably replace a text file so interruption cannot leave a partial file."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{target.name}.", suffix=".tmp", dir=target.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, target)
        directory_fd = os.open(target.parent, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def atomic_write_json(path: str | Path, value) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")
