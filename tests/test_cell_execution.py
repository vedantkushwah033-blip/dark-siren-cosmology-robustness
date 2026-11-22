from pathlib import Path

import pytest

from dark_siren_robustness.cell_execution import execute_planned_cell
from dark_siren_robustness.cell_ledger_store import initialize_checkpoint, read_checkpoint


KEY = (0, 1.0, 100, "A_catalogue_only")
PLAN = [{
    "replicate": KEY[0],
    "gamma": KEY[1],
    "validation_size": KEY[2],
    "arm": KEY[3],
    "plan_token": "frozen-cell-001",
}]


def _checkpoint(tmp_path: Path) -> Path:
    path = tmp_path / "checkpoint.json"
    initialize_checkpoint(path, PLAN)
    return path


def test_cell_execution_persists_success_only_after_independent_audit(tmp_path):
    path = _checkpoint(tmp_path)
    calls = []

    def evaluate(plan):
        calls.append(("evaluate", plan["plan_token"]))
        return {"artifact_manifest": "cells/a/manifest.json"}

    def audit(plan, artifacts):
        calls.append(("audit", artifacts["artifact_manifest"]))
        return {"status": "PASS", "audit_report_path": "cells/a/audit.json"}

    row = execute_planned_cell(path, KEY, evaluator=evaluate, auditor=audit)

    assert calls == [
        ("evaluate", "frozen-cell-001"),
        ("audit", "cells/a/manifest.json"),
    ]
    assert row["status"] == "SUCCESS"
    assert row["audit_status"] == "PASS"
    assert row["audit_report_path"] == "cells/a/audit.json"
    assert row["started_at"] and row["completed_at"]
    assert read_checkpoint(path)["cell_ledger"][0]["status"] == "SUCCESS"


def test_cell_execution_records_evaluator_exception_as_failed(tmp_path):
    path = _checkpoint(tmp_path)

    def evaluate(_plan):
        raise RuntimeError("synthetic worker failure")

    def audit(_plan, _artifacts):
        pytest.fail("auditor must not run after evaluator failure")

    row = execute_planned_cell(path, KEY, evaluator=evaluate, auditor=audit)

    assert row["status"] == "FAILED"
    assert row["reason_code"] == "CHILD_PROCESS_ERROR"
    assert row["audit_status"] == "NOT_RUN"
    assert "synthetic worker failure" in row["reason_detail"]


def test_cell_execution_records_audit_rejection_as_failed(tmp_path):
    path = _checkpoint(tmp_path)

    row = execute_planned_cell(
        path, KEY,
        evaluator=lambda _plan: {"artifact_manifest": "cells/a/manifest.json"},
        auditor=lambda _plan, _artifacts: {
            "status": "NOT_PASS", "reason_detail": "reconstruction mismatch"
        },
    )

    assert row["status"] == "FAILED"
    assert row["reason_code"] == "CHILD_AUDIT_FAILED"
    assert row["audit_status"] == "NOT_PASS"
    assert row["reason_detail"] == "reconstruction mismatch"


def test_cell_execution_rejects_audit_pass_without_report_path(tmp_path):
    path = _checkpoint(tmp_path)

    row = execute_planned_cell(
        path, KEY,
        evaluator=lambda _plan: {"artifact_manifest": "cells/a/manifest.json"},
        auditor=lambda _plan, _artifacts: {"status": "PASS"},
    )

    assert row["status"] == "FAILED"
    assert row["reason_code"] == "CHILD_AUDIT_FAILED"
    assert row["audit_status"] == "NOT_PASS"


def test_cell_execution_does_not_rerun_terminal_cell(tmp_path):
    path = _checkpoint(tmp_path)
    calls = []

    def evaluate(_plan):
        calls.append("evaluate")
        return {"artifact_manifest": "cells/a/manifest.json"}

    auditor = lambda _plan, _artifacts: {
        "status": "PASS", "audit_report_path": "cells/a/audit.json"
    }
    execute_planned_cell(path, KEY, evaluator=evaluate, auditor=auditor)

    with pytest.raises(ValueError, match="terminal outcome"):
        execute_planned_cell(path, KEY, evaluator=evaluate, auditor=auditor)

    assert calls == ["evaluate"]
    assert read_checkpoint(path)["cell_ledger"][0]["status"] == "SUCCESS"


def test_cell_execution_callbacks_receive_immutable_plan_row(tmp_path):
    path = _checkpoint(tmp_path)
    observed = {}

    def evaluate(plan):
        with pytest.raises(TypeError):
            plan["gamma"] = 99.0
        return {"artifact_manifest": "cells/a/manifest.json"}

    def audit(plan, _artifacts):
        observed["gamma"] = plan["gamma"]
        with pytest.raises(TypeError):
            plan["plan_token"] = "mutated"
        return {"status": "PASS", "audit_report_path": "cells/a/audit.json"}

    row = execute_planned_cell(path, KEY, evaluator=evaluate, auditor=audit)

    assert observed["gamma"] == KEY[1]
    assert row["gamma"] == KEY[1]
    assert read_checkpoint(path)["planned_cells"][0]["plan_token"] == "frozen-cell-001"
