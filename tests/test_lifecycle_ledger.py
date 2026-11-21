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
        started_at="2026-11-21T00:00:00Z", completed_at="2026-11-21T00:01:00Z",
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



def test_success_and_non_evaluable_require_start_timestamp():
    row = initial_ledger(_plan())[0]
    with pytest.raises(ValueError, match="requires started_at"):
        set_cell_outcome(
            row, status="SUCCESS", audit_status="PASS",
            completed_at="2026-11-21T00:01:00Z",
        )
    with pytest.raises(ValueError, match="requires started_at"):
        set_cell_outcome(
            row, status="NON_EVALUABLE", audit_status="PASS",
            reason_code="PREDECLARED_WEIGHT_INSTABILITY",
            completed_at="2026-11-21T00:01:00Z",
        )


def test_ledger_validator_rejects_missing_or_reversed_timestamps():
    plan = _plan()
    rows = initial_ledger(plan)
    rows[0] = set_cell_outcome(
        rows[0], status="SUCCESS", audit_status="PASS",
        started_at="2026-11-21T00:00:00Z",
        completed_at="2026-11-21T00:01:00Z", audit_report_path="audit/a.json",
    )
    missing_start = [dict(row) for row in rows]
    missing_start[0]["started_at"] = None
    with pytest.raises(ValueError, match="no start timestamp"):
        validate_ledger(plan, missing_start)

    reversed_times = [dict(row) for row in rows]
    reversed_times[0]["started_at"] = "2026-11-21T00:02:00Z"
    with pytest.raises(ValueError, match="cannot precede"):
        validate_ledger(plan, reversed_times)


def test_atomic_json_write_leaves_parseable_complete_file(tmp_path):
    target = tmp_path / "nested" / "ledger.json"
    atomic_write_json(target, {"rows": [{"status": "NOT_RUN"}]})
    assert json.loads(target.read_text()) == {"rows": [{"status": "NOT_RUN"}]}
    assert not list(target.parent.glob("*.tmp"))


from dark_siren_robustness.cell_ledger_store import (
    finalize_interrupted_cells,
    initialize_checkpoint,
    mark_cell_started,
    read_checkpoint,
    record_cell_outcome,
)


def _checkpoint_plan():
    return [
        {"replicate": 0, "gamma": 0.0, "validation_size": 100, "arm": "A",
         "analysis_seed": 123, "h0_true": 70.0},
        {"replicate": 0, "gamma": 0.0, "validation_size": 100, "arm": "E",
         "analysis_seed": 123, "h0_true": 70.0},
    ]


def test_checkpoint_atomically_starts_with_full_not_run_plan(tmp_path):
    path = tmp_path / "run" / "checkpoint.json"
    checkpoint = initialize_checkpoint(path, _checkpoint_plan())
    loaded = read_checkpoint(path)
    assert checkpoint == loaded
    assert len(loaded["planned_cells"]) == 2
    assert [row["status"] for row in loaded["cell_ledger"]] == ["NOT_RUN", "NOT_RUN"]
    assert len(loaded["plan_sha256"]) == 64
    with pytest.raises(FileExistsError):
        initialize_checkpoint(path, _checkpoint_plan())


def test_checkpoint_persists_start_and_success_transitions(tmp_path):
    path = tmp_path / "checkpoint.json"
    initialize_checkpoint(path, _checkpoint_plan())
    key = (0, 0.0, 100, "A")
    mark_cell_started(path, key, "2026-11-21T05:00:00Z")
    started = read_checkpoint(path)["cell_ledger"][0]
    assert started["status"] == "NOT_RUN"
    assert started["reason_code"] == "RUN_NOT_STARTED"
    assert started["started_at"] == "2026-11-21T05:00:00Z"
    record_cell_outcome(
        path, key, status="SUCCESS", audit_status="PASS",
        started_at="2026-11-21T05:00:00Z",
        completed_at="2026-11-21T05:01:00Z",
        audit_report_path="audit/A.json",
    )
    completed = read_checkpoint(path)["cell_ledger"][0]
    assert completed["status"] == "SUCCESS"
    assert completed["audit_status"] == "PASS"
    assert completed["started_at"] < completed["completed_at"]


def test_checkpoint_recovery_keeps_finished_and_marks_unfinished(tmp_path):
    path = tmp_path / "checkpoint.json"
    initialize_checkpoint(path, _checkpoint_plan())
    key = (0, 0.0, 100, "A")
    mark_cell_started(path, key, "2026-11-21T05:00:00Z")
    record_cell_outcome(
        path, key, status="SUCCESS", audit_status="PASS",
        started_at="2026-11-21T05:00:00Z",
        completed_at="2026-11-21T05:01:00Z",
        audit_report_path="audit/A.json",
    )
    finalize_interrupted_cells(path, "Test interruption after first cell.")
    rows = read_checkpoint(path)["cell_ledger"]
    assert [row["status"] for row in rows] == ["SUCCESS", "NOT_RUN"]
    assert rows[1]["reason_code"] == "NOT_SCHEDULED"
    assert rows[1]["reason_detail"] == "Test interruption after first cell."


def test_checkpoint_rejects_plan_tampering_and_terminal_overwrite(tmp_path):
    path = tmp_path / "checkpoint.json"
    initialize_checkpoint(path, _checkpoint_plan())
    checkpoint = json.loads(path.read_text())
    checkpoint["planned_cells"][0]["analysis_seed"] = 999
    path.write_text(json.dumps(checkpoint))
    with pytest.raises(ValueError, match="planned-cell hash mismatch"):
        read_checkpoint(path)



def test_checkpoint_refuses_to_overwrite_a_terminal_cell(tmp_path):
    path = tmp_path / "checkpoint.json"
    initialize_checkpoint(path, _checkpoint_plan())
    key = (0, 0.0, 100, "A")
    mark_cell_started(path, key, "2026-11-21T05:00:00Z")
    record_cell_outcome(
        path, key, status="SUCCESS", audit_status="PASS",
        completed_at="2026-11-21T05:01:00Z",
        audit_report_path="audit/A.json",
    )
    with pytest.raises(ValueError, match="cannot be overwritten"):
        record_cell_outcome(
            path, key, status="FAILED", audit_status="NOT_PASS",
            reason_code="CHILD_PROCESS_ERROR",
        )

def test_checkpoint_cannot_rewrite_persisted_start_timestamp(tmp_path):
    path = tmp_path / "checkpoint.json"
    initialize_checkpoint(path, _checkpoint_plan())
    key = (0, 0.0, 100, "A")
    mark_cell_started(path, key, "2026-11-21T05:00:00Z")
    with pytest.raises(ValueError, match="start timestamp cannot be changed"):
        record_cell_outcome(
            path, key, status="SUCCESS", audit_status="PASS",
            started_at="2026-11-21T05:00:30Z",
            completed_at="2026-11-21T05:01:00Z",
        )
    assert read_checkpoint(path)["cell_ledger"][0]["started_at"] == "2026-11-21T05:00:00Z"


def test_success_and_non_evaluable_require_audit_report_path():
    row = initial_ledger(_plan())[0]
    with pytest.raises(ValueError, match="requires audit_report_path"):
        set_cell_outcome(
            row, status="SUCCESS", audit_status="PASS",
            started_at="2026-11-21T00:00:00Z",
            completed_at="2026-11-21T00:01:00Z",
        )
    with pytest.raises(ValueError, match="requires audit_report_path"):
        set_cell_outcome(
            row, status="NON_EVALUABLE", audit_status="PASS",
            reason_code="PREDECLARED_WEIGHT_INSTABILITY",
            started_at="2026-11-21T00:00:00Z",
            completed_at="2026-11-21T00:01:00Z",
        )


def test_ledger_validator_rejects_removed_audit_report_path():
    plan = _plan()
    rows = initial_ledger(plan)
    rows[0] = set_cell_outcome(
        rows[0], status="SUCCESS", audit_status="PASS",
        started_at="2026-11-21T00:00:00Z",
        completed_at="2026-11-21T00:01:00Z",
        audit_report_path="audit/a.json",
    )
    rows[0]["audit_report_path"] = None
    with pytest.raises(ValueError, match="no audit report path"):
        validate_ledger(plan, rows)


def test_failed_outcomes_require_start_and_completion_timestamps():
    row = initial_ledger(_plan())[0]
    with pytest.raises(ValueError, match="requires started_at"):
        set_cell_outcome(
            row, status="FAILED", audit_status="NOT_PASS",
            reason_code="CHILD_PROCESS_ERROR",
            completed_at="2026-11-21T00:01:00Z",
        )
    with pytest.raises(ValueError, match="requires completed_at"):
        set_cell_outcome(
            row, status="FAILED", audit_status="NOT_PASS",
            reason_code="CHILD_PROCESS_ERROR",
            started_at="2026-11-21T00:00:00Z",
        )


def test_ledger_validator_rejects_failed_cell_without_terminal_timestamps():
    plan = _plan()
    rows = initial_ledger(plan)
    rows[0] = set_cell_outcome(
        rows[0], status="FAILED", audit_status="NOT_PASS",
        reason_code="CHILD_PROCESS_ERROR",
        started_at="2026-11-21T00:00:00Z",
        completed_at="2026-11-21T00:01:00Z",
    )
    rows[0]["completed_at"] = None
    with pytest.raises(ValueError, match="no completion timestamp"):
        validate_ledger(plan, rows)


def test_concurrent_checkpoint_updates_do_not_lose_completed_cells(tmp_path):
    from concurrent.futures import ThreadPoolExecutor

    from dark_siren_robustness.cell_ledger_store import (
        mark_cell_started,
        record_cell_outcome,
    )

    plan = [
        {
            "replicate": index,
            "gamma": 0.0,
            "validation_size": 100,
            "arm": f"arm_{index}",
            "analysis_seed": 1000 + index,
            "h0_true": 70.0,
        }
        for index in range(16)
    ]
    path = tmp_path / "concurrent" / "checkpoint.json"
    initialize_checkpoint(path, plan)

    def complete(row):
        key = (row["replicate"], row["gamma"], row["validation_size"], row["arm"])
        started_at = "2026-11-21T05:00:00Z"
        mark_cell_started(path, key, started_at)
        record_cell_outcome(
            path,
            key,
            status="SUCCESS",
            audit_status="PASS",
            started_at=started_at,
            completed_at="2026-11-21T05:01:00Z",
            audit_report_path=f"audit/{row['arm']}.json",
        )

    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(complete, plan))

    checkpoint = read_checkpoint(path)
    assert len(checkpoint["cell_ledger"]) == len(plan)
    assert all(row["status"] == "SUCCESS" for row in checkpoint["cell_ledger"])
    assert all(row["audit_status"] == "PASS" for row in checkpoint["cell_ledger"])


def test_lifecycle_timestamps_require_iso8601_explicit_offsets():
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


def test_timestamp_order_uses_instants_across_utc_offsets():
    row = initial_ledger(_plan())[0]
    result = set_cell_outcome(
        row, status="SUCCESS", audit_status="PASS",
        started_at="2026-11-21T00:30:00Z",
        completed_at="2026-11-21T00:00:00-01:00",
        audit_report_path="audit/a.json",
    )
    assert result["status"] == "SUCCESS"
    with pytest.raises(ValueError, match="cannot precede"):
        set_cell_outcome(
            row, status="SUCCESS", audit_status="PASS",
            started_at="2026-11-21T00:30:00Z",
            completed_at="2026-11-20T23:15:00-01:00",
            audit_report_path="audit/a.json",
        )


def test_checkpoint_rejects_terminal_failure_without_persisted_start(tmp_path):
    path = tmp_path / "checkpoint.json"
    initialize_checkpoint(path, _checkpoint_plan())
    key = (0, 0.0, 100, "A")
    with pytest.raises(ValueError, match="persisted start transition"):
        record_cell_outcome(
            path, key, status="FAILED", audit_status="NOT_PASS",
            reason_code="CHILD_PROCESS_ERROR",
            started_at="2026-11-21T05:00:00Z",
            completed_at="2026-11-21T05:01:00Z",
        )


def test_checkpoint_recovery_marks_active_and_unstarted_cells_distinctly(tmp_path):
    path = tmp_path / "checkpoint.json"
    initialize_checkpoint(path, _checkpoint_plan())
    key = (0, 0.0, 100, "A")
    mark_cell_started(path, key, "2026-11-21T05:00:00Z")
    live = read_checkpoint(path)["cell_ledger"]
    assert live[0]["reason_code"] == "RUN_NOT_STARTED"
    assert live[0]["started_at"] == "2026-11-21T05:00:00Z"
    finalize_interrupted_cells(path, "Process stopped.")
    rows = read_checkpoint(path)["cell_ledger"]
    assert rows[0]["reason_code"] == "RUN_INTERRUPTED"
    assert rows[1]["reason_code"] == "NOT_SCHEDULED"


def test_interrupted_started_cell_can_resume_without_rewriting_first_start(tmp_path):
    path = tmp_path / "resume" / "checkpoint.json"
    initialize_checkpoint(path, _checkpoint_plan())
    key = (0, 0.0, 100, "A")
    first_start = "2026-11-21T05:00:00Z"
    mark_cell_started(path, key, first_start)
    finalize_interrupted_cells(path, "Simulated process interruption.")

    interrupted = read_checkpoint(path)["cell_ledger"][0]
    assert interrupted["status"] == "NOT_RUN"
    assert interrupted["reason_code"] == "RUN_INTERRUPTED"
    assert interrupted["started_at"] == first_start

    # A retry is allowed only after the interrupted state has been persisted.
    mark_cell_started(path, key, "2026-11-21T05:05:00Z")
    resumed = read_checkpoint(path)["cell_ledger"][0]
    assert resumed["status"] == "NOT_RUN"
    assert resumed["reason_code"] == "RUN_NOT_STARTED"
    assert resumed["started_at"] == first_start

    record_cell_outcome(
        path, key, status="SUCCESS", audit_status="PASS",
        started_at=first_start,
        completed_at="2026-11-21T05:06:00Z",
        audit_report_path="audit/A.json",
    )
    completed = read_checkpoint(path)["cell_ledger"][0]
    assert completed["status"] == "SUCCESS"
    assert completed["started_at"] == first_start


def _complete_checkpoint_cell_in_process(args):
    """Top-level worker so the process-lock test exercises independent processes."""
    path, row = args
    from dark_siren_robustness.cell_ledger_store import mark_cell_started, record_cell_outcome

    key = (row["replicate"], row["gamma"], row["validation_size"], row["arm"])
    started_at = "2026-11-21T05:00:00Z"
    mark_cell_started(path, key, started_at)
    record_cell_outcome(
        path,
        key,
        status="SUCCESS",
        audit_status="PASS",
        started_at=started_at,
        completed_at="2026-11-21T05:01:00Z",
        audit_report_path=f"audit/{row['arm']}.json",
    )


def test_multiprocess_checkpoint_updates_preserve_every_terminal_cell(tmp_path):
    """Exercise sidecar flock across processes, not merely threads."""
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor

    if "fork" not in multiprocessing.get_all_start_methods():
        pytest.skip("fcntl/flock process test requires a POSIX fork context")

    plan = [
        {
            "replicate": index,
            "gamma": 0.0,
            "validation_size": 100,
            "arm": f"process_arm_{index}",
            "analysis_seed": 2000 + index,
            "h0_true": 70.0,
        }
        for index in range(20)
    ]
    path = tmp_path / "multiprocess" / "checkpoint.json"
    initialize_checkpoint(path, plan)

    context = multiprocessing.get_context("fork")
    with ProcessPoolExecutor(max_workers=8, mp_context=context) as pool:
        list(pool.map(_complete_checkpoint_cell_in_process, [(str(path), row) for row in plan]))

    checkpoint = read_checkpoint(path)
    assert len(checkpoint["cell_ledger"]) == len(plan)
    assert all(row["status"] == "SUCCESS" for row in checkpoint["cell_ledger"])
    assert all(row["audit_status"] == "PASS" for row in checkpoint["cell_ledger"])
    assert len({row["audit_report_path"] for row in checkpoint["cell_ledger"]}) == len(plan)



def test_checkpoint_rejects_ledger_provenance_tampering(tmp_path):
    """The immutable plan hash alone is insufficient; ledger rows must match it."""
    path = tmp_path / "checkpoint.json"
    initialize_checkpoint(path, _checkpoint_plan())
    checkpoint = json.loads(path.read_text())
    checkpoint["cell_ledger"][0]["analysis_seed"] = 999
    path.write_text(json.dumps(checkpoint))
    with pytest.raises(ValueError, match="planned provenance mismatch"):
        read_checkpoint(path)


def test_checkpoint_rejects_unsupported_schema_version(tmp_path):
    path = tmp_path / "checkpoint.json"
    initialize_checkpoint(path, _checkpoint_plan())
    checkpoint = json.loads(path.read_text())
    checkpoint["schema_version"] = 999
    path.write_text(json.dumps(checkpoint))
    with pytest.raises(ValueError, match="unsupported.*schema"):
        read_checkpoint(path)


def _initialize_checkpoint_in_process(args):
    path, plan = args
    from dark_siren_robustness.cell_ledger_store import initialize_checkpoint

    try:
        initialize_checkpoint(path, plan)
        return "CREATED"
    except FileExistsError:
        return "ALREADY_EXISTS"


def test_concurrent_checkpoint_initialization_has_exactly_one_winner(tmp_path):
    """Two processes racing to create one run ID must not replace each other's plan."""
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor

    if "fork" not in multiprocessing.get_all_start_methods():
        pytest.skip("fcntl/flock process test requires a POSIX fork context")

    path = tmp_path / "initialize-race" / "checkpoint.json"
    plan_a = _checkpoint_plan()
    plan_b = [dict(row, analysis_seed=999) for row in _checkpoint_plan()]
    context = multiprocessing.get_context("fork")
    with ProcessPoolExecutor(max_workers=2, mp_context=context) as pool:
        outcomes = list(pool.map(
            _initialize_checkpoint_in_process,
            [(str(path), plan_a), (str(path), plan_b)],
        ))

    assert sorted(outcomes) == ["ALREADY_EXISTS", "CREATED"]
    checkpoint = read_checkpoint(path)
    persisted_plan = checkpoint["planned_cells"]
    assert persisted_plan in (plan_a, plan_b)


def test_not_run_cannot_claim_completion_or_an_audit_report():
    row = initial_ledger(_plan())[0]
    with pytest.raises(ValueError, match="NOT_RUN cannot have completed_at"):
        set_cell_outcome(
            row, status="NOT_RUN", audit_status="NOT_RUN",
            reason_code="RUN_INTERRUPTED",
            started_at="2026-11-21T00:00:00Z",
            completed_at="2026-11-21T00:01:00Z",
        )
    with pytest.raises(ValueError, match="NOT_RUN cannot have audit_report_path"):
        set_cell_outcome(
            row, status="NOT_RUN", audit_status="NOT_RUN",
            reason_code="RUN_INTERRUPTED",
            started_at="2026-11-21T00:00:00Z",
            audit_report_path="audit/unfinished.json",
        )


def test_ledger_validator_rejects_not_run_completion_metadata():
    plan = _plan()
    rows = initial_ledger(plan)
    rows[0]["completed_at"] = "2026-11-21T00:01:00Z"
    with pytest.raises(ValueError, match="NOT_RUN cell has a completion timestamp"):
        validate_ledger(plan, rows)
    rows = initial_ledger(plan)
    rows[0]["audit_report_path"] = "audit/unfinished.json"
    with pytest.raises(ValueError, match="NOT_RUN cell has an audit report path"):
        validate_ledger(plan, rows)



def test_not_run_reason_codes_have_consistent_start_timestamp_semantics(tmp_path):
    plan = _plan()
    rows = initial_ledger(plan)

    # An unstarted cell and an active cell may both be RUN_NOT_STARTED.
    rows[0] = set_cell_outcome(
        rows[0], status="NOT_RUN", audit_status="NOT_RUN",
        reason_code="RUN_NOT_STARTED",
    )
    rows[1] = set_cell_outcome(
        rows[1], status="NOT_RUN", audit_status="NOT_RUN",
        reason_code="RUN_NOT_STARTED",
        started_at="2026-11-21T00:00:00Z",
    )
    rows[2] = set_cell_outcome(
        rows[2], status="NOT_RUN", audit_status="NOT_RUN",
        reason_code="RUN_INTERRUPTED",
        started_at="2026-11-21T00:00:00Z",
    )
    rows[3] = set_cell_outcome(
        rows[3], status="NOT_RUN", audit_status="NOT_RUN",
        reason_code="NOT_SCHEDULED",
    )
    assert validate_ledger(plan, rows)["status"] == "PASS"

    with pytest.raises(ValueError, match="RUN_INTERRUPTED requires"):
        set_cell_outcome(
            initial_ledger(plan)[0], status="NOT_RUN", audit_status="NOT_RUN",
            reason_code="RUN_INTERRUPTED",
        )
    with pytest.raises(ValueError, match="NOT_SCHEDULED cannot have"):
        set_cell_outcome(
            initial_ledger(plan)[0], status="NOT_RUN", audit_status="NOT_RUN",
            reason_code="NOT_SCHEDULED",
            started_at="2026-11-21T00:00:00Z",
        )
    with pytest.raises(ValueError, match="UPSTREAM_FAILURE cannot have"):
        set_cell_outcome(
            initial_ledger(plan)[0], status="NOT_RUN", audit_status="NOT_RUN",
            reason_code="UPSTREAM_FAILURE",
            started_at="2026-11-21T00:00:00Z",
        )


def test_checkpoint_start_claim_requires_timestamp(tmp_path):
    path = tmp_path / "checkpoint.json"
    initialize_checkpoint(path, _checkpoint_plan())
    key = (0, 0.0, 100, "A")
    with pytest.raises(ValueError, match="requires a non-empty"):
        mark_cell_started(path, key, None)
    with pytest.raises(ValueError, match="UTC offset"):
        mark_cell_started(path, key, "2026-11-21T05:00:00")
    assert read_checkpoint(path)["cell_ledger"][0]["started_at"] is None


def test_checkpoint_rejects_terminal_outcome_without_active_start_claim(tmp_path):
    path = tmp_path / "checkpoint.json"
    key = (0, 0.0, 100, "A")
    initialize_checkpoint(path, _checkpoint_plan())

    # Supplying timestamps to the terminal call must not bypass mark_cell_started.
    with pytest.raises(ValueError, match="active cell with a persisted start"):
        record_cell_outcome(
            path, key, status="SUCCESS", audit_status="PASS",
            started_at="2026-11-21T05:00:00Z",
            completed_at="2026-11-21T05:01:00Z",
            audit_report_path="audit/A.json",
        )

    # A cell finalized as not scheduled cannot later be committed as evaluated.
    finalize_interrupted_cells(path, "No worker was scheduled.")
    with pytest.raises(ValueError, match="active cell with a persisted start"):
        record_cell_outcome(
            path, key, status="SUCCESS", audit_status="PASS",
            started_at="2026-11-21T05:00:00Z",
            completed_at="2026-11-21T05:01:00Z",
            audit_report_path="audit/A.json",
        )
