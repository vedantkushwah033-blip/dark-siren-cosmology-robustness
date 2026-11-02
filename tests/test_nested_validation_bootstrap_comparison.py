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
        "status": "PARTIAL_INCOMPLETE_SCENARIOS",
        "scenarios": 1,
        "completed_scenarios": 0,
        "failed_scenarios": 1,
    }), encoding="utf-8")
    pd.DataFrame([{
        "scenario": 0,
        "scenario_status": "failed",
        "scenario_failure": "ledger failure",
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
        '"scenario_failure":"ledger failure"}\\n',
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
