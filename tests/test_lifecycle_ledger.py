"""Tests for strict cell-level lifecycle ledger primitives."""
import json

import pytest

from dark_siren_robustness.lifecycle_ledger import (
    atomic_write_json,
    initial_ledger,
    set_cell_outcome,
    validate_ledger,
)


def _plan():
    return [
        {"replicate": 0, "gamma": -1.0, "validation_size": 100, "arm": "A"},
        {"replicate": 0, "gamma": -1.0, "validation_size": 100, "arm": "E"},
        {"replicate": 1, "gamma": 1.0, "validation_size": 300, "arm": "A"},
        {"replicate": 1, "gamma": 1.0, "validation_size": 300, "arm": "E"},
    ]


def test_initial_plan_is_complete_and_every_cell_starts_not_run():
    plan = _plan()
    ledger = initial_ledger(plan)
    report = validate_ledger(plan, ledger)
    assert report["status"] == "PASS"
    assert report["status_counts"] == {
        "FAILED": 0, "NON_EVALUABLE": 0, "NOT_RUN": 4, "SUCCESS": 0
    }
    assert all(row["reason_code"] == "RUN_NOT_STARTED" for row in ledger)


def test_all_four_lifecycle_states_validate_with_distinct_rules():
    plan = _plan()
    rows = initial_ledger(plan)
    rows[0] = set_cell_outcome(
        rows[0], status="SUCCESS", audit_status="PASS",
        started_at="2026-11-21T00:00:00Z",
        completed_at="2026-11-21T00:01:00Z", audit_report_path="audit/a.json",
    )
    rows[1] = set_cell_outcome(
        rows[1], status="NON_EVALUABLE", audit_status="PASS",
        reason_code="PREDECLARED_WEIGHT_INSTABILITY",
        reason_detail="Threshold frozen in protocol.",
        started_at="2026-11-21T00:00:00Z",
        completed_at="2026-11-21T00:01:00Z", audit_report_path="audit/e.json",
    )
    rows[2] = set_cell_outcome(
        rows[2], status="FAILED", audit_status="NOT_PASS",
        reason_code="CHILD_AUDIT_FAILED", reason_detail="Independent reconstruction failed.",
        started_at="2026-11-21T00:00:00Z",
    )
    report = validate_ledger(plan, rows)
    assert report["status_counts"] == {
        "FAILED": 1, "NON_EVALUABLE": 1, "NOT_RUN": 1, "SUCCESS": 1
    }


def test_non_evaluable_requires_audit_pass_and_frozen_reason():
    row = initial_ledger(_plan())[0]
    with pytest.raises(ValueError, match="requires an independent audit PASS"):
        set_cell_outcome(
            row, status="NON_EVALUABLE", audit_status="NOT_PASS",
            reason_code="PREDECLARED_WEIGHT_INSTABILITY",
            completed_at="2026-11-21T00:01:00Z",
        )
    with pytest.raises(ValueError, match="not allowed"):
        set_cell_outcome(
            row, status="NON_EVALUABLE", audit_status="PASS",
            reason_code="ESS_TOO_LOW", completed_at="2026-11-21T00:01:00Z",
        )


def test_audit_failure_cannot_be_downgraded_to_non_evaluable():
    row = initial_ledger(_plan())[0]
    with pytest.raises(ValueError, match="requires an independent audit PASS"):
        set_cell_outcome(
            row, status="NON_EVALUABLE", audit_status="NOT_PASS",
            reason_code="INVALID_POSTERIOR_NORMALIZATION",
            completed_at="2026-11-21T00:01:00Z",
        )


def test_failed_and_not_run_cannot_claim_audit_pass():
    row = initial_ledger(_plan())[0]
    with pytest.raises(ValueError, match="cannot claim audit PASS"):
        set_cell_outcome(
            row, status="FAILED", audit_status="PASS",
            reason_code="CHILD_PROCESS_ERROR",
        )
    with pytest.raises(ValueError, match="requires audit status NOT_RUN"):
        set_cell_outcome(
            row, status="NOT_RUN", audit_status="NOT_PASS",
            reason_code="RUN_INTERRUPTED",
        )


def test_plan_and_ledger_reject_missing_duplicate_and_extra_keys():
    plan = _plan()
    rows = initial_ledger(plan)
    with pytest.raises(ValueError, match="duplicate planned cell"):
        initial_ledger(plan + [plan[0]])
    with pytest.raises(ValueError, match="duplicate cell keys"):
        validate_ledger(plan, rows + [rows[0]])
    with pytest.raises(ValueError, match="key mismatch"):
        validate_ledger(plan, rows[:-1])
    extra = dict(rows[-1], arm="X")
    with pytest.raises(ValueError, match="key mismatch"):
        validate_ledger(plan, rows[:-1] + [extra])


def test_unknown_status_and_reason_are_rejected():
    plan = _plan()
    rows = initial_ledger(plan)
    rows[0]["status"] = "RUNNING"
    with pytest.raises(ValueError, match="unknown lifecycle status"):
        validate_ledger(plan, rows)
    rows = initial_ledger(plan)
    rows[0]["reason_code"] = "MADE_UP_AFTER_LOOKING"
    with pytest.raises(ValueError, match="invalid reason"):
        validate_ledger(plan, rows)


def test_success_requires_audit_and_completion_timestamp():
    row = initial_ledger(_plan())[0]
    with pytest.raises(ValueError, match="requires completed_at"):
        set_cell_outcome(row, status="SUCCESS", audit_status="PASS")
    with pytest.raises(ValueError, match="requires an independent audit PASS"):
        set_cell_outcome(
            row, status="SUCCESS", audit_status="NOT_RUN",
            completed_at="2026-11-21T00:01:00Z",
        )


def test_atomic_json_write_leaves_parseable_complete_file(tmp_path):
    target = tmp_path / "nested" / "ledger.json"
    atomic_write_json(target, {"rows": [{"status": "NOT_RUN"}]})
    assert json.loads(target.read_text()) == {"rows": [{"status": "NOT_RUN"}]}
    assert not list(target.parent.glob("*.tmp"))


def test_ledger_validator_rejects_mutated_planned_provenance():
    plan = [
        {"cell_id": "rep000_g0_n100_A", "replicate": 0, "gamma": 0.0,
         "validation_size": 100, "arm": "A", "analysis_seed": 123,
         "h0_true": 70.0},
    ]
    rows = initial_ledger(plan)
    rows[0]["analysis_seed"] = 124
    with pytest.raises(ValueError, match="planned provenance mismatch"):
        validate_ledger(plan, rows)

    rows = initial_ledger(plan)
    rows[0]["h0_true"] = 71.0
    with pytest.raises(ValueError, match="planned provenance mismatch"):
        validate_ledger(plan, rows)

    rows = initial_ledger(plan)
    rows[0]["cell_id"] = "mutated-cell-id"
    with pytest.raises(ValueError, match="planned provenance mismatch"):
        validate_ledger(plan, rows)



def test_timestamp_validation_requires_iso8601_with_explicit_offset():
    row = initial_ledger(_plan())[0]
    with pytest.raises(ValueError, match="valid ISO-8601"):
        set_cell_outcome(
            row, status="SUCCESS", audit_status="PASS",
            started_at="yesterday", completed_at="today",
            audit_report_path="audit/a.json",
        )
    with pytest.raises(ValueError, match="UTC offset"):
        set_cell_outcome(
            row, status="SUCCESS", audit_status="PASS",
            started_at="2026-11-21T00:00:00",
            completed_at="2026-11-21T00:01:00Z",
            audit_report_path="audit/a.json",
        )


def test_timestamp_order_uses_instants_not_lexicographic_strings():
    row = initial_ledger(_plan())[0]
    result = set_cell_outcome(
        row, status="SUCCESS", audit_status="PASS",
        started_at="2026-11-21T00:30:00Z",
        completed_at="2026-11-21T00:00:00-01:00",
        audit_report_path="audit/a.json",
    )
    assert result["status"] == "SUCCESS"
    plan = _plan()
    ledger = initial_ledger(plan)
    ledger[0] = result
    validate_ledger(plan, ledger)



def test_success_and_non_evaluable_require_audit_report_path():
    plan = _plan()
    for status, reason in [
        ("SUCCESS", ""),
        ("NON_EVALUABLE", "PREDECLARED_WEIGHT_INSTABILITY"),
    ]:
        row = initial_ledger(plan)[0]
        with pytest.raises(ValueError, match="requires audit_report_path"):
            set_cell_outcome(
                row, status=status, audit_status="PASS", reason_code=reason,
                started_at="2026-11-21T00:00:00Z",
                completed_at="2026-11-21T00:01:00Z",
            )

    row = initial_ledger(plan)[0]
    completed = set_cell_outcome(
        row, status="SUCCESS", audit_status="PASS",
        started_at="2026-11-21T00:00:00Z",
        completed_at="2026-11-21T00:01:00Z",
        audit_report_path="audit/report.json",
    )
    completed["audit_report_path"] = "   "
    with pytest.raises(ValueError, match="no audit report path"):
        validate_ledger(plan, [completed])
