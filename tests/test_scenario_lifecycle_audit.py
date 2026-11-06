from __future__ import annotations

import json
import pandas as pd
import pytest

from dark_siren_robustness.scenario_lifecycle import (
    ScenarioLifecycleError,
    run_scenario_campaign,
)
from dark_siren_robustness.scenario_lifecycle_audit import audit_scenario_lifecycle


def executor_success(scenario):
    return pd.DataFrame([
        {"draw": d, "validation_size": size, "status": "success"}
        for d in range(2) for size in (10, 20)
    ])


def execute(tmp_path, executor=executor_success):
    return run_scenario_campaign(
        n_scenarios=2, n_draws=2, validation_sizes=(10, 20),
        executor=executor, checkpoint_dir=tmp_path,
    )


def test_independent_audit_accepts_complete_success(tmp_path):
    execute(tmp_path)
    result = audit_scenario_lifecycle(tmp_path)
    assert result["status"] == "SCENARIO_LIFECYCLE_AUDIT_PASS"
    assert result["n_planned"] == 8
    assert result["n_success"] == 8
    assert result["n_not_run"] == 0


def test_independent_audit_accepts_upstream_scenario_failure(tmp_path):
    def executor(scenario):
        if scenario == 0:
            raise ScenarioLifecycleError("oracle_reference", "oracle_failed", "oracle fit failed")
        return executor_success(scenario)

    execute(tmp_path, executor)
    result = audit_scenario_lifecycle(tmp_path)
    assert result["n_non_evaluable"] == 4
    assert result["n_success"] == 4
    assert result["run_status"] == "COMPLETED_WITH_FAILURES"


def test_independent_audit_accepts_draw_level_failure(tmp_path):
    def executor(scenario):
        frame = executor_success(scenario)
        if scenario == 1:
            frame.loc[0, "status"] = "failed"
            frame.loc[0, "failure_code"] = "optimizer"
            frame.loc[0, "failure_reason"] = "failed fit"
        return frame

    execute(tmp_path, executor)
    result = audit_scenario_lifecycle(tmp_path)
    assert result["n_failed"] == 1
    assert result["n_success"] == 7


def test_independent_audit_accepts_interrupted_running_checkpoint(tmp_path):
    class StopNow(BaseException):
        pass

    def executor(scenario):
        raise StopNow()

    with pytest.raises(StopNow):
        execute(tmp_path, executor)
    result = audit_scenario_lifecycle(tmp_path)
    assert result["n_not_run"] == 8
    assert result["run_status"] == "RUNNING"


def test_independent_audit_rejects_tampered_denominator(tmp_path):
    execute(tmp_path)
    checkpoint_path = tmp_path / "checkpoint.json"
    snapshot = json.loads(checkpoint_path.read_text())
    snapshot["ledger_records"] = snapshot["ledger_records"][:-1]
    checkpoint_path.write_text(json.dumps(snapshot))
    with pytest.raises(AssertionError, match="denominator mismatch"):
        audit_scenario_lifecycle(tmp_path)


def test_independent_audit_rejects_tampered_state_counts(tmp_path):
    execute(tmp_path)
    checkpoint_path = tmp_path / "checkpoint.json"
    snapshot = json.loads(checkpoint_path.read_text())
    snapshot["run_state"]["attempt_counts"]["n_success"] = 7
    checkpoint_path.write_text(json.dumps(snapshot))
    with pytest.raises(AssertionError, match="counts disagree"):
        audit_scenario_lifecycle(tmp_path)


def test_independent_audit_rejects_invalid_status_and_transition_count(tmp_path):
    execute(tmp_path)
    checkpoint_path = tmp_path / "checkpoint.json"
    snapshot = json.loads(checkpoint_path.read_text())
    snapshot["ledger_records"][0]["status"] = "MAYBE"
    checkpoint_path.write_text(json.dumps(snapshot))
    with pytest.raises(AssertionError, match="unknown attempt status"):
        audit_scenario_lifecycle(tmp_path)

def test_atomic_snapshot_preserves_last_complete_scenario_when_next_is_interrupted(tmp_path):
    class StopNow(BaseException):
        pass

    def executor(scenario):
        if scenario == 1:
            raise StopNow()
        return executor_success(scenario)

    with pytest.raises(StopNow):
        execute(tmp_path, executor)
    snapshot = json.loads((tmp_path / "checkpoint.json").read_text())
    ledger = pd.DataFrame(snapshot["ledger_records"])
    state = snapshot["run_state"]
    assert set(ledger.loc[ledger.scenario.eq(0), "status"]) == {"SUCCESS"}
    assert set(ledger.loc[ledger.scenario.eq(1), "status"]) == {"NOT_RUN"}
    assert state["scenario_states"] == {"0": "SUCCESS", "1": "RUNNING"}
    audited = audit_scenario_lifecycle(tmp_path)
    assert audited["n_success"] == 4
    assert audited["n_not_run"] == 4
    assert audited["run_status"] == "RUNNING"
