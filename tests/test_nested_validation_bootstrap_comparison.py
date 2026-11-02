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
    }), encoding="utf-8")
    pd.DataFrame(columns=[
        "scenario", "validation_size", "outer_draw", "method",
        "bootstrap_status",
    ]).to_csv(tmp_path / "calibration_records.csv", index=False)
    (tmp_path / "attempt_audit.jsonl").write_text(
        '{"scenario":0,"validation_size":20,"outer_draw":0,"method":"paired_bootstrap"}\\n',
        encoding="utf-8",
    )

    with pytest.raises(AssertionError, match="partial/incomplete artifacts"):
        audit(tmp_path)
