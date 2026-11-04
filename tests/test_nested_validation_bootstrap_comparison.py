from __future__ import annotations

import pytest

from scripts.run_nested_validation_bootstrap_comparison import _seed, run_comparison


def test_nested_comparison_seed_derivation_is_reproducible_and_keyed():
    first = _seed(20261030, 2, 4, 100, 1)
    assert first == _seed(20261030, 2, 4, 100, 1)
    assert first != _seed(20261030, 2, 4, 100, 2)


def test_nested_comparison_rejects_duplicate_validation_sizes_before_simulation():
    with pytest.raises(ValueError, match="must not contain duplicates"):
        run_comparison(
            scenarios=1, draws=2, n_parent=100, n_detected=2,
            validation_sizes=(20, 20), bootstrap_replicates=2, seed=7,
        )


def test_nested_comparison_rejects_invalid_bootstrap_count_before_simulation():
    with pytest.raises(ValueError, match="bootstrap_replicates"):
        run_comparison(
            scenarios=1, draws=2, n_parent=100, n_detected=2,
            validation_sizes=(20,), bootstrap_replicates=1, seed=7,
        )


@pytest.mark.parametrize("sizes", [(20.5,), (True,), (20, 30.5)])
def test_nested_comparison_rejects_noninteger_validation_sizes(sizes):
    with pytest.raises(TypeError, match="validation sizes must be integers"):
        run_comparison(
            scenarios=1, draws=2, n_parent=100, n_detected=2,
            validation_sizes=sizes, bootstrap_replicates=2, seed=7,
        )


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"scenarios": True}, "scenarios must be an integer"),
        ({"draws": 2.5}, "draws must be an integer"),
        ({"h0_true": float("nan")}, "h0_true must be finite and positive"),
        ({"sigma_mpc": 0.0}, "sigma_mpc must be finite and positive"),
        ({"threshold_mpc": float("inf")}, "threshold_mpc must be finite and positive"),
    ],
)
def test_nested_comparison_rejects_malformed_parameters_before_simulation(kwargs, message):
    parameters = dict(
        scenarios=1, draws=2, n_parent=100, n_detected=2,
        validation_sizes=(20,), bootstrap_replicates=2, seed=7,
    )
    parameters.update(kwargs)
    with pytest.raises((TypeError, ValueError), match=message):
        run_comparison(**parameters)


def test_nested_comparison_persists_wrapper_exception_metadata(monkeypatch):
    import numpy as np
    import pandas as pd

    import scripts.run_nested_validation_bootstrap_comparison as runner

    def fake_population(**kwargs):
        return {
            "parent_truth": pd.DataFrame({
                "redshift": np.linspace(0.01, 0.40, kwargs["n_parent"])
            }),
            "gw_events_observed": pd.DataFrame({
                "distance_mpc": np.full(kwargs["n_detected"], 350.0)
            }),
        }

    def fake_posterior(*args, **kwargs):
        return {}, {"posterior_mean": 70.0}

    def fail_bootstrap(*args, **kwargs):
        raise ValueError("forced wrapper failure for serialization test")

    monkeypatch.setattr(runner, "generate_synthetic_population", fake_population)
    monkeypatch.setattr(runner, "_posterior", fake_posterior)
    monkeypatch.setattr(runner, "bootstrap_validation_h0_posterior", fail_bootstrap)

    records, attempts = runner.run_comparison(
        scenarios=1, draws=2, n_parent=100, n_detected=2,
        validation_sizes=(20,), bootstrap_replicates=2, seed=19,
    )

    failed = attempts.loc[attempts["method"].isin(
        ["paired_bootstrap", "design_aware_srswor"]
    )]
    assert len(failed) == 4
    assert failed["bootstrap_wrapper_failure_class"].eq("ValueError").all()
    assert failed["bootstrap_wrapper_failure_message"].eq(
        "forced wrapper failure for serialization test"
    ).all()
    assert records.loc[
        records["method"].isin(["paired_bootstrap", "design_aware_srswor"]),
        "bootstrap_status",
    ].eq("failed").all()
    assert failed["n_bootstrap_attempted"].isna().all()
    assert failed["n_bootstrap_successful"].isna().all()
    assert failed["n_bootstrap_failed"].isna().all()
    assert failed["bootstrap_failure_records"].map(len).eq(0).all()


def test_nested_comparison_preserves_partial_bootstrap_attempt_ledger(monkeypatch):
    import numpy as np
    import pandas as pd

    import scripts.run_nested_validation_bootstrap_comparison as runner

    def fake_population(**kwargs):
        return {
            "parent_truth": pd.DataFrame({
                "redshift": np.linspace(0.01, 0.40, kwargs["n_parent"])
            }),
            "gw_events_observed": pd.DataFrame({
                "distance_mpc": np.full(kwargs["n_detected"], 350.0)
            }),
        }

    def fake_posterior(*args, **kwargs):
        return {}, {"posterior_mean": 70.0}

    def partial_bootstrap(*args, **kwargs):
        return {
            "n_attempted": 3, "n_successful": 1, "n_failed": 2,
            "success_mask": np.array([True, False, False]),
            "posterior_mean_percentile_interval_95": (np.nan, np.nan),
            "raw_bootstrap_standard_error": np.nan,
            "failure_records": [
                {"attempt": 1, "reason": "forced failure 1"},
                {"attempt": 2, "reason": "forced failure 2"},
            ],
            "fpc_factor": 0.9, "input_fingerprints": {},
            "attempt_posterior_mean": np.array([70.1, np.nan, np.nan]),
            "reference_posterior_mean_by_attempt": np.array([70.0, np.nan, np.nan]),
            "posterior_median_draws_by_attempt": np.array([70.1, np.nan, np.nan]),
            "ci68_low_draws_by_attempt": np.array([69.0, np.nan, np.nan]),
            "ci68_high_draws_by_attempt": np.array([71.0, np.nan, np.nan]),
            "inverse_weight_ess_by_attempt": np.array([10.0, np.nan, np.nan]),
            "inverse_weight_maximum_by_attempt": np.array([0.2, np.nan, np.nan]),
            "coefficient_boundary_limited_by_attempt": np.array([False, False, False]),
            "posterior_density_draws": np.empty((1, 81)),
            "fixed_full_validation_reference_mean": 70.0,
            "sample_indices": [np.array([1]), np.array([2]), np.array([3])],
            "pseudo_population_indices": [np.array([], dtype=int)] * 3,
            "selected_pseudo_population_positions": [np.array([], dtype=int)] * 3,
        }

    monkeypatch.setattr(runner, "generate_synthetic_population", fake_population)
    monkeypatch.setattr(runner, "_posterior", fake_posterior)
    monkeypatch.setattr(runner, "bootstrap_validation_h0_posterior", partial_bootstrap)

    records, attempts = runner.run_comparison(
        scenarios=1, draws=2, n_parent=100, n_detected=2,
        validation_sizes=(20,), bootstrap_replicates=3, seed=29,
    )

    independent = attempts.loc[attempts["method"].isin(
        ["paired_bootstrap", "design_aware_srswor"]
    )]
    assert independent["bootstrap_status"].eq("non_evaluable").all()
    assert independent["n_bootstrap_attempted"].eq(3).all()
    assert independent["n_bootstrap_successful"].eq(1).all()
    assert independent["n_bootstrap_failed"].eq(2).all()
    assert independent["bootstrap_success_mask"].map(sum).eq(1).all()
    assert independent["bootstrap_failure_records"].map(len).eq(2).all()
    assert records.loc[
        records["method"].isin(["paired_bootstrap", "design_aware_srswor"]),
        "bootstrap_se",
    ].isna().all()
    derived = attempts.loc[attempts["method"].eq("paired_bootstrap_fpc")]
    assert derived["bootstrap_is_derived"].all()
    assert derived["bootstrap_source_method"].eq("paired_bootstrap").all()


def test_nested_comparison_records_scenario_target_failure_and_continues(monkeypatch):
    import numpy as np
    import pandas as pd

    import scripts.run_nested_validation_bootstrap_comparison as runner

    def fake_population(**kwargs):
        return {
            "parent_truth": pd.DataFrame({
                "redshift": np.linspace(0.01, 0.40, kwargs["n_parent"])
            }),
            "gw_events_observed": pd.DataFrame({
                "distance_mpc": np.full(kwargs["n_detected"], 350.0)
            }),
        }

    def fail_target(*args, **kwargs):
        raise RuntimeError("forced finite-parent target failure")

    monkeypatch.setattr(runner, "generate_synthetic_population", fake_population)
    monkeypatch.setattr(runner, "_posterior", fail_target)

    records, attempts = runner.run_comparison(
        scenarios=1, draws=2, n_parent=100, n_detected=2,
        validation_sizes=(20,), bootstrap_replicates=2, seed=23,
    )

    assert len(records) == 2 * 1 * 3
    assert len(attempts) == len(records)
    assert records["outer_status"].eq("failed").all()
    assert records["bootstrap_status"].eq("not_run").all()
    assert records["scenario_failure"].str.contains(
        "forced finite-parent target failure"
    ).all()
    assert attempts["scenario_failure"].str.contains(
        "forced finite-parent target failure"
    ).all()


def test_replay_auditor_rejects_partial_manifest_before_reconstruction(tmp_path):
    import json

    import pandas as pd

    from scripts.audit_nested_validation_bootstrap_replay import audit

    (tmp_path / "manifest.json").write_text(json.dumps({
        "schema": "nested_validation_bootstrap_comparison_smoke_v2",
        "status": "PARTIAL_INCOMPLETE_SCENARIOS",
        "scenarios": 1,
        "completed_scenarios": 0,
        "failed_scenarios": 1,
    }), encoding="utf-8")
    pd.DataFrame([{
        "scenario": 0,
        "scenario_status": "failed",
        "scenario_failure": "forced scenario target failure",
    }]).to_csv(tmp_path / "scenario_failure_records.csv", index=False)
    pd.DataFrame(columns=[
        "scenario", "validation_size", "outer_draw", "method",
        "bootstrap_status",
    ]).to_csv(tmp_path / "calibration_records.csv", index=False)
    (tmp_path / "attempt_audit.jsonl").write_text(
        '{"scenario":0,"validation_size":20,"outer_draw":0,"method":"paired_bootstrap"}\n',
        encoding="utf-8",
    )

    with pytest.raises(AssertionError, match="partial/incomplete artifacts"):
        audit(tmp_path)


def test_replay_auditor_rejects_failure_ledger_row_mismatch(tmp_path):
    import json

    import pandas as pd

    from scripts.audit_nested_validation_bootstrap_replay import audit

    (tmp_path / "manifest.json").write_text(json.dumps({
        "schema": "nested_validation_bootstrap_comparison_smoke_v2",
        "status": "SMOKE_ONLY_NOT_CONFIRMATORY",
        "scenarios": 1,
        "planned_scenarios": 1,
        "completed_scenarios": 1,
        "failed_scenarios": 0,
        "outer_draws": 1,
        "validation_sizes": [20],
        "planned_method_rows": 3,
        "completed_method_rows": 1,
        "summary_scenario_ids": [0],
        "parent_size": 100,
        "detected_event_count": 2,
        "seed": 7,
        "model_parameters": {},
        "bootstrap_replicates_per_method_per_outer_draw": 2,
    }), encoding="utf-8")
    pd.DataFrame([{
        "scenario": 0,
        "scenario_status": "complete",
        "scenario_failure": "",
    }]).to_csv(tmp_path / "scenario_failure_records.csv", index=False)
    pd.DataFrame([{
        "scenario": 0,
        "validation_size": 20,
        "outer_draw": 0,
        "method": "paired_bootstrap",
        "scenario_failure": "different row failure",
        "outer_status": "failed",
        "bootstrap_status": "not_run",
    }]).to_csv(tmp_path / "calibration_records.csv", index=False)
    (tmp_path / "attempt_audit.jsonl").write_text(
        '{"scenario":0,"validation_size":20,"outer_draw":0,"method":"paired_bootstrap",'
        '"scenario_failure":""}\n',
        encoding="utf-8",
    )

    with pytest.raises(AssertionError, match="calibration rows do not reconcile"):
        audit(tmp_path)


def _write_key_integrity_fixture(directory, mutation):
    import json

    import pandas as pd

    methods = ["paired_bootstrap", "paired_bootstrap_fpc", "design_aware_srswor"]
    rows = [{
        "scenario": 0, "validation_size": 20, "outer_draw": 0,
        "method": method, "scenario_failure": "",
        "outer_status": "success", "bootstrap_status": "success",
    } for method in methods]
    attempts = [{
        "scenario": 0, "validation_size": 20, "outer_draw": 0,
        "method": method, "scenario_failure": "",
    } for method in methods]
    if mutation == "missing":
        rows.pop()
        attempts.pop()
    elif mutation == "duplicate":
        rows.append(dict(rows[0]))
        attempts.append(dict(attempts[0]))
    elif mutation == "unexpected":
        rows[0]["validation_size"] = 999
        attempts[0]["validation_size"] = 999
    elif mutation == "altered_manifest":
        pass

    pd.DataFrame(rows).to_csv(directory / "calibration_records.csv", index=False)
    pd.DataFrame(attempts).to_json(
        directory / "attempt_audit.jsonl", orient="records", lines=True
    )
    pd.DataFrame([{
        "scenario": 0, "scenario_status": "complete", "scenario_failure": "",
    }]).to_csv(directory / "scenario_failure_records.csv", index=False)
    (directory / "manifest.json").write_text(json.dumps({
        "schema": "nested_validation_bootstrap_comparison_smoke_v2",
        "status": "SMOKE_ONLY_NOT_CONFIRMATORY",
        "scenarios": 1, "planned_scenarios": 1,
        "completed_scenarios": 1, "failed_scenarios": 0,
        "outer_draws": 1, "validation_sizes": [20],
        "planned_method_rows": 3 if mutation != "altered_manifest" else 4,
        "completed_method_rows": 3,
        "summary_scenario_ids": [0],
        "parent_size": 100,
        "detected_event_count": 2,
        "seed": 7,
        "model_parameters": {},
        "bootstrap_replicates_per_method_per_outer_draw": 2,
    }), encoding="utf-8")


@pytest.mark.parametrize(
    "mutation, message",
    [
        ("missing", "missing or unexpected planned keys"),
        ("duplicate", "duplicate scenario/draw/size/method keys"),
        ("unexpected", "missing or unexpected planned keys"),
        ("altered_manifest", "planned/completed method-row denominators"),
    ],
)
def test_replay_auditor_rejects_tampered_artifact_key_sets(tmp_path, mutation, message):
    from scripts.audit_nested_validation_bootstrap_replay import audit

    _write_key_integrity_fixture(tmp_path, mutation)

    with pytest.raises(AssertionError, match=message):
        audit(tmp_path)



def test_replay_auditor_rejects_missing_scenarios_before_ledger_access(tmp_path):
    import json

    from scripts.audit_nested_validation_bootstrap_replay import audit

    (tmp_path / "manifest.json").write_text(json.dumps({
        "schema": "nested_validation_bootstrap_comparison_smoke_v2",
        "status": "SMOKE_ONLY_NOT_CONFIRMATORY",
    }), encoding="utf-8")

    with pytest.raises(AssertionError, match="manifest is missing required fields.*scenarios"):
        audit(tmp_path)


def test_replay_auditor_rejects_missing_manifest_fields():
    from scripts.audit_nested_validation_bootstrap_replay import _require_manifest_fields

    with pytest.raises(AssertionError, match="manifest is missing required fields"):
        _require_manifest_fields({
            "schema": "nested_validation_bootstrap_comparison_smoke_v2",
            "status": "SMOKE_ONLY_NOT_CONFIRMATORY",
        })


@pytest.mark.parametrize(
    ("filename", "column"),
    [
        ("scenario_summary.csv", "bootstrap_se"),
        ("scenario_balanced_summary.csv", "bootstrap_se"),
    ],
)
def test_replay_auditor_rejects_modified_summary_values(tmp_path, filename, column):
    import pandas as pd

    from scripts.audit_nested_validation_bootstrap_replay import _require_summary_match

    rebuilt = pd.DataFrame({
        "scenario": [0, 0],
        "validation_size": [20, 40],
        "method": ["paired_bootstrap", "paired_bootstrap"],
        "bootstrap_se": [1.25, 2.5],
    })
    stored = rebuilt.copy()
    stored.loc[0, column] += 0.5
    stored.to_csv(tmp_path / filename, index=False)

    with pytest.raises(AssertionError, match="numeric column mismatch"):
        _require_summary_match(tmp_path, filename, rebuilt)

def test_nested_comparison_does_not_misclassify_invalid_returned_ledger_as_wrapper_failure(monkeypatch):
    import numpy as np
    import pandas as pd

    import scripts.run_nested_validation_bootstrap_comparison as runner

    def fake_population(**kwargs):
        return {
            "parent_truth": pd.DataFrame({
                "redshift": np.linspace(0.01, 0.40, kwargs["n_parent"])
            }),
            "gw_events_observed": pd.DataFrame({
                "distance_mpc": np.full(kwargs["n_detected"], 350.0)
            }),
        }

    def fake_posterior(*args, **kwargs):
        return {}, {"posterior_mean": 70.0}

    def invalid_ledger(*args, **kwargs):
        return {
            "n_attempted": 3,
            "n_successful": 1,
            "n_failed": 1,
            "success_mask": np.array([True, False, False]),
        }

    monkeypatch.setattr(runner, "generate_synthetic_population", fake_population)
    monkeypatch.setattr(runner, "_posterior", fake_posterior)
    monkeypatch.setattr(runner, "bootstrap_validation_h0_posterior", invalid_ledger)

    with pytest.raises(ValueError, match="counts do not reconcile with success mask"):
        runner.run_comparison(
            scenarios=1, draws=1, n_parent=100, n_detected=2,
            validation_sizes=(20,), bootstrap_replicates=3, seed=31,
        )
