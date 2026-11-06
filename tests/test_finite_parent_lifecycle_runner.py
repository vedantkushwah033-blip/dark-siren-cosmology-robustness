from __future__ import annotations

import numpy as np
import pandas as pd

from dark_siren_robustness import finite_parent_sampling as finite_parent
from dark_siren_robustness.scenario_lifecycle_audit import audit_scenario_lifecycle


def _small_design():
    return {
        "n_scenarios": 2,
        "n_validation_draws": 2,
        "n_parent": 1200,
        "n_detected": 6,
        "validation_sizes": (100,),
        "seed": 20261110,
        "h0_points": 31,
        "z_points": 101,
    }


def test_checkpointed_finite_parent_runner_emits_full_ledger_and_auditable_snapshot(tmp_path):
    records, ledger, state = finite_parent.run_finite_parent_benchmark_checkpointed(
        checkpoint_dir=tmp_path, **_small_design()
    )
    assert len(records) == 4
    assert len(ledger) == 4
    assert state["attempt_counts"]["n_planned"] == 4
    assert state["attempt_counts"]["n_success"] + state["attempt_counts"]["n_failed"] + state["attempt_counts"]["n_non_evaluable"] == 4
    assert audit_scenario_lifecycle(tmp_path)["n_planned"] == 4
    assert (tmp_path / "checkpoint.json").exists()


def test_checkpointed_runner_records_generation_failure_without_losing_planned_draws(
    tmp_path, monkeypatch
):
    original = finite_parent.generate_synthetic_population
    design = _small_design()
    root = np.random.SeedSequence([design["seed"], 0])
    failing_sim_seed = int(root.spawn(3)[0].generate_state(1, dtype=np.uint32)[0])

    def selective_generation_failure(*, seed, **kwargs):
        if seed == failing_sim_seed:
            raise RuntimeError("injected parent-generation failure")
        return original(seed=seed, **kwargs)

    monkeypatch.setattr(
        finite_parent, "generate_synthetic_population", selective_generation_failure
    )
    records, ledger, state = finite_parent.run_finite_parent_benchmark_checkpointed(
        checkpoint_dir=tmp_path, **design
    )
    failed = ledger.loc[ledger.scenario.eq(0)]
    succeeded = ledger.loc[ledger.scenario.eq(1)]
    assert len(failed) == 2
    assert set(failed.status) == {"NON_EVALUABLE"}
    assert set(failed.stage) == {"scenario_generation"}
    assert set(failed.failure_code) == {"model_fit_failure"}
    assert set(succeeded.status).issubset({"SUCCESS", "FAILED"})
    assert len(succeeded) == 2
    assert set(records["scenario"]) == {1}
    assert state["scenario_states"]["0"] == "FAILED"
    audit = audit_scenario_lifecycle(tmp_path)
    assert audit["n_planned"] == 4
    assert audit["n_non_evaluable"] == 2


def test_checkpointed_runner_rejects_invalid_global_design_before_creating_plan(tmp_path):
    design = _small_design()
    design["validation_sizes"] = (100.5,)
    try:
        finite_parent.run_finite_parent_benchmark_checkpointed(
            checkpoint_dir=tmp_path, **design
        )
    except TypeError:
        pass
    else:
        raise AssertionError("fractional validation size should be rejected")
    assert not (tmp_path / "checkpoint.json").exists()


import pytest


@pytest.mark.parametrize(
    "target_name,fail_on_call,expected_stage",
    [
        ("estimate_host_density", 1, "oracle_reference"),
        ("estimate_host_density", 2, "catalogue_reference"),
        ("estimate_host_density_from_validation", 1, "finite_parent_reference"),
    ],
)
def test_checkpointed_runner_classifies_reference_fit_failure_by_stage(
    tmp_path, monkeypatch, target_name, fail_on_call, expected_stage
):
    original = getattr(finite_parent, target_name)
    calls = {"count": 0}

    def selective_failure(*args, **kwargs):
        calls["count"] += 1
        if calls["count"] == fail_on_call:
            raise RuntimeError("injected reference-stage failure")
        return original(*args, **kwargs)

    monkeypatch.setattr(finite_parent, target_name, selective_failure)
    design = _small_design()
    records, ledger, state = finite_parent.run_finite_parent_benchmark_checkpointed(
        checkpoint_dir=tmp_path, **design
    )
    failed = ledger.loc[ledger.scenario.eq(0)]
    assert set(failed.status) == {"NON_EVALUABLE"}
    assert set(failed.stage) == {expected_stage}
    assert state["scenario_states"]["0"] == "FAILED"
    assert set(records["scenario"]) == {1}
    assert audit_scenario_lifecycle(tmp_path)["n_non_evaluable"] == 2


def test_checkpointed_runner_treats_catalogue_size_gate_as_non_evaluable(tmp_path):
    design = _small_design()
    design.update(true_intercept=-20.0, true_slope=0.0, true_quadratic=0.0)
    records, ledger, state = finite_parent.run_finite_parent_benchmark_checkpointed(
        checkpoint_dir=tmp_path, **design
    )
    assert len(records) == 4
    assert set(records.status) == {"non_evaluable"}
    assert set(ledger.status) == {"NON_EVALUABLE"}
    assert set(ledger.stage) == {"catalogue_reference"}
    assert audit_scenario_lifecycle(tmp_path)["n_non_evaluable"] == 4
