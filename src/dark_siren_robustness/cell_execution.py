"""Checkpoint-owning orchestration for one independently audited cell.

This generic lifecycle adapter is not yet wired into the legacy five-arm
replicate runner. Numerical evaluation and independent auditing are injected
callbacks and remain separate from lifecycle persistence.
"""
from __future__ import annotations

from datetime import datetime, timezone
from typing import Callable, Mapping

from .cell_ledger_store import mark_cell_started, read_checkpoint, record_cell_outcome


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="microseconds").replace("+00:00", "Z")


def _cell_row(checkpoint: Mapping, key: tuple) -> dict:
    matches = [
        row for row in checkpoint["cell_ledger"]
        if (row["replicate"], row["gamma"], row["validation_size"], row["arm"]) == tuple(key)
    ]
    if len(matches) != 1:
        raise KeyError(f"planned cell key must match exactly once: {key}")
    return matches[0]


def _record_failure(
    checkpoint_path, key, *, reason_code: str, audit_status: str, detail: str
) -> dict:
    checkpoint = record_cell_outcome(
        checkpoint_path, tuple(key), status="FAILED", audit_status=audit_status,
        reason_code=reason_code, reason_detail=detail[:4000], completed_at=_utc_now(),
    )
    return _cell_row(checkpoint, key)


def execute_planned_cell(
    checkpoint_path,
    key: tuple,
    *,
    evaluator: Callable[[Mapping], Mapping],
    auditor: Callable[[Mapping, Mapping], Mapping],
) -> dict:
    """Evaluate and independently audit exactly one planned cell.

    The evaluator receives the immutable planned-cell row and returns an
    artifact descriptor. The auditor receives that row and descriptor and must
    return {"status": "PASS", "audit_report_path": "..."} for success.

    Callback exceptions become FAILED outcomes. This adapter deliberately does
    not infer NON_EVALUABLE: scientific exclusion must be predeclared and
    independently audited by the experiment protocol.
    """
    checkpoint = read_checkpoint(checkpoint_path)
    plan_matches = [
        row for row in checkpoint["planned_cells"]
        if (row["replicate"], row["gamma"], row["validation_size"], row["arm"]) == tuple(key)
    ]
    if len(plan_matches) != 1:
        raise KeyError(f"planned cell key must match exactly once: {key}")
    plan_row = plan_matches[0]

    mark_cell_started(checkpoint_path, tuple(key), _utc_now())
    try:
        artifacts = evaluator(plan_row)
        if not isinstance(artifacts, Mapping):
            raise TypeError("cell evaluator must return a mapping artifact descriptor")
    except Exception as exc:
        return _record_failure(
            checkpoint_path, key, reason_code="CHILD_PROCESS_ERROR",
            audit_status="NOT_RUN", detail=f"{type(exc).__name__}: {exc}",
        )

    try:
        report = auditor(plan_row, artifacts)
        if not isinstance(report, Mapping) or report.get("status") != "PASS":
            detail = "independent cell auditor did not return status PASS"
            if isinstance(report, Mapping):
                detail = str(report.get("reason_detail") or detail)
            return _record_failure(
                checkpoint_path, key, reason_code="CHILD_AUDIT_FAILED",
                audit_status="NOT_PASS", detail=detail,
            )
        report_path = report.get("audit_report_path")
        if not isinstance(report_path, str) or not report_path.strip():
            return _record_failure(
                checkpoint_path, key, reason_code="CHILD_AUDIT_FAILED",
                audit_status="NOT_PASS",
                detail="independent audit PASS omitted audit_report_path",
            )
    except Exception as exc:
        return _record_failure(
            checkpoint_path, key, reason_code="CHILD_AUDIT_FAILED",
            audit_status="NOT_PASS", detail=f"{type(exc).__name__}: {exc}",
        )

    checkpoint = record_cell_outcome(
        checkpoint_path, tuple(key), status="SUCCESS", audit_status="PASS",
        completed_at=_utc_now(), audit_report_path=report_path,
    )
    return _cell_row(checkpoint, key)
