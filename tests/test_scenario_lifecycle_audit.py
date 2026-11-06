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
    ledger_path = tmp_path / "attempt_ledger.csv"
    frame = pd.read_csv(ledger_path)
    frame = frame.iloc[:-1]
    frame.to_csv(ledger_path, index=False)
    with pytest.raises(AssertionError, match="denominator mismatch"):
        audit_scenario_lifecycle(tmp_path)


def test_independent_audit_rejects_tampered_state_counts(tmp_path):
    execute(tmp_path)
    state_path = tmp_path / "run_state.json"
    state = json.loads(state_path.read_text())
    state["attempt_counts"]["n_success"] = 7
    state_path.write_text(json.dumps(state))
    with pytest.raises(AssertionError, match="counts disagree"):
        audit_scenario_lifecycle(tmp_path)


def test_independent_audit_rejects_invalid_status_and_transition_count(tmp_path):
    execute(tmp_path)
    ledger_path = tmp_path / "attempt_ledger.csv"
    frame = pd.read_csv(ledger_path)
    frame.loc[0, "status"] = "MAYBE"
    frame.to_csv(ledger_path, index=False)
    with pytest.raises(AssertionError, match="unknown attempt status"):
        audit_scenario_lifecycle(tmp_path)
