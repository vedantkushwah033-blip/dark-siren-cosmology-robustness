"""Execute, independently audit, and lifecycle-record exactly one arm cell.

Recovery contract: after a process interruption, stop all workers using the
checkpoint, call finalize_interrupted_cells on the checkpoint, then rerun the
cell. This command never globally finalizes other cells because that would be
unsafe while a multi-cell worker pool is active.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from dark_siren_robustness.cell_ledger_store import (
    mark_cell_started,
    read_checkpoint,
    record_cell_outcome,
)
from dark_siren_robustness.lifecycle_ledger import atomic_write_json
from scripts.audit_arm_cell_from_bundle import audit_cell
from scripts.run_arm_cell_from_bundle import run_cell


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _key(replicate: int, gamma: float, validation_size: int, arm: str) -> tuple:
    return (int(replicate), float(gamma), int(validation_size), str(arm))


def _row_key(row: dict) -> tuple:
    return (int(row["replicate"]), float(row["gamma"]),
            int(row["validation_size"]), str(row["arm"]))


def run_audited_cell(
    checkpoint_path: Path,
    bundle_dir: Path,
    output_dir: Path,
    *,
    replicate: int,
    gamma: float,
    validation_size: int,
    arm: str,
    audit_report_path: Path | None = None,
) -> dict:
    checkpoint_path = Path(checkpoint_path).resolve()
    output_dir = Path(output_dir).resolve()
    audit_path = (
        Path(audit_report_path).resolve()
        if audit_report_path is not None
        else output_dir.with_name(output_dir.name + ".independent_audit.json")
    )
    key = _key(replicate, gamma, validation_size, arm)

    checkpoint = read_checkpoint(checkpoint_path)
    matches = [row for row in checkpoint["cell_ledger"] if _row_key(row) == key]
    if len(matches) != 1:
        raise KeyError(f"checkpoint must contain exactly one planned cell: {key}")
    row = matches[0]
    if row["status"] != "NOT_RUN":
        raise ValueError(f"cell already has terminal status {row['status']}")
    if row.get("started_at") is not None and row.get("reason_code") != "RUN_INTERRUPTED":
        raise RuntimeError(
            "cell has a persisted start but is not marked RUN_INTERRUPTED; "
            "stop all workers, finalize interrupted cells, then retry"
        )

    mark_cell_started(checkpoint_path, key, _now())

    if not output_dir.exists():
        try:
            run_cell(
                bundle_dir, output_dir, replicate=replicate, gamma=gamma,
                validation_size=validation_size, arm=arm,
            )
        except Exception as exc:
            detail = f"{type(exc).__name__}: {exc}"
            record_cell_outcome(
                checkpoint_path, key, status="FAILED", audit_status="NOT_RUN",
                reason_code="CHILD_PROCESS_ERROR", reason_detail=detail,
                completed_at=_now(),
            )
            return {"status": "FAILED", "reason_code": "CHILD_PROCESS_ERROR",
                    "reason_detail": detail, "cell_key": list(key)}

    try:
        report = audit_cell(bundle_dir, output_dir)
        if report.get("status") != "PASS" or report.get("audit_status") != "PASS":
            raise ValueError("independent cell auditor did not return PASS")
        atomic_write_json(audit_path, report)
    except Exception as exc:
        detail = f"{type(exc).__name__}: {exc}"
        record_cell_outcome(
            checkpoint_path, key, status="FAILED", audit_status="NOT_PASS",
            reason_code="CHILD_AUDIT_FAILED", reason_detail=detail,
            completed_at=_now(),
        )
        return {"status": "FAILED", "reason_code": "CHILD_AUDIT_FAILED",
                "reason_detail": detail, "cell_key": list(key)}

    record_cell_outcome(
        checkpoint_path, key, status="SUCCESS", audit_status="PASS",
        reason_code="", reason_detail="", completed_at=_now(),
        audit_report_path=str(audit_path),
    )
    return {
        "status": "SUCCESS",
        "audit_status": "PASS",
        "cell_key": list(key),
        "output_dir": str(output_dir),
        "audit_report_path": str(audit_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--bundle-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--audit-report", type=Path, default=None)
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--gamma", type=float, required=True)
    parser.add_argument("--validation-size", type=int, required=True)
    parser.add_argument("--arm", required=True)
    args = parser.parse_args()
    result = run_audited_cell(
        args.checkpoint, args.bundle_dir, args.output_dir,
        replicate=args.replicate, gamma=args.gamma,
        validation_size=args.validation_size, arm=args.arm,
        audit_report_path=args.audit_report,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
