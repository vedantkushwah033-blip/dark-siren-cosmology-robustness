import json
from pathlib import Path

import pytest

from scripts.audit_arm_cell_plan import audit_checkpoint


def _checkpoint(tmp_path, *, status="SUCCESS"):
    output = tmp_path / "cell-output"
    output.mkdir()
    audit = tmp_path / "cell-audit.json"
    key = {"replicate": 1, "gamma": 0.0, "validation_size": 100,
           "arm": "A_catalogue_only"}
    audit.write_text(json.dumps({
        "status": "PASS", "audit_status": "PASS", "cell_key": key,
    }))
    plan = [{**key, "bundle_dir": str(tmp_path / "bundle"),
             "output_dir": str(output)}]
    import hashlib
    payload = json.dumps(plan, sort_keys=True, separators=(",", ":"), allow_nan=False)
    digest = hashlib.sha256(payload.encode()).hexdigest()
    row = {
        **plan[0], "status": status,
        "reason_code": "" if status == "SUCCESS" else "CHILD_PROCESS_ERROR",
        "reason_detail": "", "audit_status": "PASS" if status == "SUCCESS" else "NOT_RUN",
        "started_at": "2026-11-22T00:00:00+00:00",
        "completed_at": "2026-11-22T00:00:01+00:00",
        "audit_report_path": str(audit) if status == "SUCCESS" else None,
    }
    checkpoint = {"schema_version": 1, "plan_sha256": digest,
                  "planned_cells": plan, "cell_ledger": [row]}
    path = tmp_path / "checkpoint.json"
    path.write_text(json.dumps(checkpoint))
    return path, checkpoint


def test_independent_auditor_accepts_successful_audited_cell(tmp_path):
    path, _ = _checkpoint(tmp_path)
    result = audit_checkpoint(path)
    assert result["status"] == "PASS"
    assert result["run_status"] == "PASS"
    assert result["status_counts"]["SUCCESS"] == 1
    assert len(result["verified_audit_reports"]) == 1


def test_auditor_rejects_tampered_plan_without_refreshed_hash(tmp_path):
    path, checkpoint = _checkpoint(tmp_path)
    checkpoint["planned_cells"][0]["gamma"] = 1.0
    path.write_text(json.dumps(checkpoint))
    with pytest.raises(ValueError, match="plan SHA-256"):
        audit_checkpoint(path)


def test_auditor_rejects_audit_report_bound_to_another_cell(tmp_path):
    path, checkpoint = _checkpoint(tmp_path)
    report_path = Path(checkpoint["cell_ledger"][0]["audit_report_path"])
    report = json.loads(report_path.read_text())
    report["cell_key"]["arm"] = "B_parent_oracle"
    report_path.write_text(json.dumps(report))
    with pytest.raises(ValueError, match="different cell"):
        audit_checkpoint(path)


def test_auditor_distinguishes_incomplete_from_completed_failure(tmp_path):
    path, checkpoint = _checkpoint(tmp_path, status="FAILED")
    result = audit_checkpoint(path)
    assert result["run_status"] == "COMPLETED_WITH_NON_SUCCESS"
    checkpoint["cell_ledger"][0].update({
        "status": "NOT_RUN", "reason_code": "RUN_NOT_STARTED",
        "started_at": None, "completed_at": None, "audit_status": "NOT_RUN",
    })
    path.write_text(json.dumps(checkpoint))
    assert audit_checkpoint(path)["run_status"] == "INCOMPLETE"


def test_auditor_rejects_success_if_audit_report_missing(tmp_path):
    path, checkpoint = _checkpoint(tmp_path)
    Path(checkpoint["cell_ledger"][0]["audit_report_path"]).unlink()
    with pytest.raises(ValueError, match="report is missing"):
        audit_checkpoint(path)
