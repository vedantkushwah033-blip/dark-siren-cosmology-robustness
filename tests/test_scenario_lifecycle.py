from __future__ import annotations

import json
import pandas as pd
import pytest

from dark_siren_robustness.scenario_lifecycle import (
    ScenarioLifecycleError,
    run_scenario_campaign,
)


def successful_records(scenario: int) -> pd.DataFrame:
    return pd.DataFrame([
        {"draw": draw, "validation_size": size, "status": "success"}
        for draw in range(2) for size in (10, 20)
    ])


def run(tmp_path, executor):
    return run_scenario_campaign(
        n_scenarios=2,
        n_draws=2,
        validation_sizes=(10, 20),
        executor=executor,
        checkpoint_dir=tmp_path,
        run_metadata={"study": "unit-test"},
    )


def test_all_planned_attempts_reconciled_and_checkpoints_written(tmp_path):
    ledger, state = run(tmp_path, successful_records)
    assert len(ledger) == 8
    assert set(ledger.status) == {"SUCCESS"}
    assert state["run_status"] == "SUCCESS"
    assert state["attempt_counts"]["n_planned"] == 8
    assert state["scenario_states"] == {"0": "SUCCESS", "1": "SUCCESS"}
    saved = pd.read_csv(tmp_path / "attempt_ledger.csv")
    assert len(saved) == 8
    assert set(saved.status) == {"SUCCESS"}
    assert json.loads((tmp_path / "run_state.json").read_text())["run_status"] == "SUCCESS"


@pytest.mark.parametrize(
    "stage",
    [
        "scenario_generation",
        "oracle_reference",
        "catalogue_reference",
        "finite_parent_reference",
    ],
)
def test_upstream_scenario_failure_marks_all_its_planned_draws_non_evaluable(tmp_path, stage):
    def executor(scenario):
        if scenario == 0:
            raise ScenarioLifecycleError(stage, "upstream_reference_failure", "reference unavailable")
        return successful_records(scenario)

    ledger, state = run(tmp_path, executor)
    failed_scenario = ledger.loc[ledger.scenario.eq(0)]
    assert len(failed_scenario) == 4
    assert set(failed_scenario.status) == {"NON_EVALUABLE"}
    assert set(failed_scenario.stage) == {stage}
    assert set(failed_scenario.failure_code) == {"upstream_reference_failure"}
    assert set(ledger.loc[ledger.scenario.eq(1), "status"]) == {"SUCCESS"}
    assert state["scenario_states"]["0"] == "FAILED"
    assert state["scenario_failures"]["0"]["stage"] == stage
    assert state["attempt_counts"]["n_non_evaluable"] == 4
    assert state["attempt_counts"]["n_planned"] == 8
    assert state["run_status"] == "COMPLETED_WITH_FAILURES"


def test_draw_level_failures_remain_distinct_from_scenario_failure(tmp_path):
    def executor(scenario):
        records = successful_records(scenario)
        if scenario == 0:
            records.loc[0, "status"] = "failed"
            records.loc[0, "failure_code"] = "optimizer_failure"
            records.loc[0, "failure_reason"] = "optimizer did not converge"
            records.loc[0, "failure_stage"] = "validation_fit"
        return records

    ledger, state = run(tmp_path, executor)
    row = ledger.loc[ledger.scenario.eq(0) & ledger.draw.eq(0) & ledger.validation_size.eq(10)].iloc[0]
    assert row.status == "FAILED"
    assert row.stage == "validation_fit"
    assert row.failure_code == "optimizer_failure"
    assert state["scenario_states"]["0"] == "PARTIAL_FAILURE"
    assert state["attempt_counts"]["n_failed"] == 1
    assert state["attempt_counts"]["n_success"] == 7


def test_missing_result_record_becomes_non_evaluable_not_silently_dropped(tmp_path):
    def executor(scenario):
        records = successful_records(scenario)
        return records.iloc[:-1].copy() if scenario == 0 else records

    ledger, state = run(tmp_path, executor)
    missing = ledger.loc[ledger.scenario.eq(0) & ledger.draw.eq(1) & ledger.validation_size.eq(20)].iloc[0]
    assert missing.status == "NON_EVALUABLE"
    assert missing.stage == "result_reconciliation"
    assert missing.failure_code == "missing_attempt_record"
    assert state["attempt_counts"]["n_planned"] == 8
    assert state["attempt_counts"]["n_non_evaluable"] == 1


def test_unexpected_result_key_fails_scenario_and_accounts_all_attempts(tmp_path):
    def executor(scenario):
        records = successful_records(scenario)
        if scenario == 0:
            records.loc[0, "draw"] = 99
        return records

    ledger, state = run(tmp_path, executor)
    group = ledger.loc[ledger.scenario.eq(0)]
    assert set(group.status) == {"NON_EVALUABLE"}
    assert set(group.stage) == {"result_reconciliation"}
    assert state["scenario_states"]["0"] == "FAILED"


def test_process_interruption_leaves_unstarted_attempts_checkpointed(tmp_path):
    class SimulatedInterruption(BaseException):
        pass

    def executor(scenario):
        if scenario == 0:
            raise SimulatedInterruption()
        return successful_records(scenario)

    with pytest.raises(SimulatedInterruption):
        run(tmp_path, executor)
    # The plan was committed before the first executor call. An abrupt process
    # interruption therefore cannot erase the planned denominator.
    saved = pd.read_csv(tmp_path / "attempt_ledger.csv")
    assert len(saved) == 8
    assert set(saved.status) == {"NOT_RUN"}
    state = json.loads((tmp_path / "run_state.json").read_text())
    assert state["scenario_states"]["0"] == "RUNNING"
    assert state["scenario_states"]["1"] == "NOT_STARTED"


def test_duplicate_attempt_keys_are_classified_as_reconciliation_failure(tmp_path):
    def executor(scenario):
        records = successful_records(scenario)
        if scenario == 0:
            return pd.concat([records, records.iloc[[0]]], ignore_index=True)
        return records

    ledger, state = run(tmp_path, executor)
    group = ledger.loc[ledger.scenario.eq(0)]
    assert set(group.status) == {"NON_EVALUABLE"}
    assert set(group.stage) == {"result_reconciliation"}
    assert state["scenario_failures"]["0"]["failure_code"] == "scenario_result_reconciliation_failure"


def test_empty_or_invalid_result_does_not_shrink_denominator(tmp_path):
    def executor(scenario):
        return pd.DataFrame(columns=["draw", "validation_size", "status"])

    ledger, state = run(tmp_path, executor)
    assert len(ledger) == 8
    assert set(ledger.status) == {"NON_EVALUABLE"}
    assert state["attempt_counts"]["n_planned"] == 8

def test_fractional_validation_size_is_not_silently_truncated(tmp_path):
    with pytest.raises(TypeError):
        run_scenario_campaign(
            n_scenarios=1, n_draws=2, validation_sizes=(10.5,),
            executor=successful_records, checkpoint_dir=tmp_path,
        )
    assert not (tmp_path / "checkpoint.json").exists()


def test_non_json_metadata_fails_before_checkpoint_or_execution(tmp_path):
    called = []
    with pytest.raises(TypeError):
        run_scenario_campaign(
            n_scenarios=1, n_draws=2, validation_sizes=(10,),
            executor=lambda scenario: called.append(scenario),
            checkpoint_dir=tmp_path, run_metadata={"invalid": object()},
        )
    assert called == []
    assert not (tmp_path / "checkpoint.json").exists()

def test_existing_checkpoint_directory_is_not_overwritten(tmp_path):
    run(tmp_path, successful_records)
    previous = (tmp_path / "checkpoint.json").read_text(encoding="utf-8")
    with pytest.raises(FileExistsError, match="fresh output directory"):
        run(tmp_path, successful_records)
    assert (tmp_path / "checkpoint.json").read_text(encoding="utf-8") == previous
