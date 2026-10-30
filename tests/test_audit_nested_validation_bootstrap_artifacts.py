from __future__ import annotations

import json

import pandas as pd
import pytest

from scripts.audit_nested_validation_bootstrap_artifacts import audit
from scripts.audit_nested_validation_bootstrap_replay import audit as audit_replay
from scripts.run_nested_validation_bootstrap_comparison import run_comparison
from dark_siren_robustness.bootstrap_calibration import summarize_bootstrap_calibration


@pytest.fixture
def nested_artifact(tmp_path):
    records, attempts = run_comparison(
        scenarios=1, draws=2, n_parent=1200, n_detected=8,
        validation_sizes=(100,), bootstrap_replicates=3, seed=20261030,
    )
    scenario, aggregate = summarize_bootstrap_calibration(records)
    records.to_csv(tmp_path / "calibration_records.csv", index=False)
    attempts.to_json(tmp_path / "attempt_audit.jsonl", orient="records", lines=True)
    scenario.to_csv(tmp_path / "scenario_summary.csv", index=False)
    aggregate.to_csv(tmp_path / "scenario_balanced_summary.csv", index=False)
    manifest = {
        "schema": "nested_validation_bootstrap_comparison_smoke_v1",
        "status": "SMOKE_ONLY_NOT_CONFIRMATORY",
        "seed": 20261030,
        "scenarios": 1,
        "outer_draws": 2,
        "parent_size": 1200,
        "detected_event_count": 8,
        "model_parameters": {
            "h0_true_km_s_mpc": 70.0,
            "distance_noise_sigma_mpc": 140.0,
            "detection_threshold_mpc": 250.0,
            "redshift_grid": {"start": 0.01, "stop": 0.40, "points": 401},
            "h0_grid": {"start": 50.0, "stop": 90.0, "points": 81},
            "merger_weight": "linear from 1.0 at z=0.01 to 2.5 at z=0.40",
            "catalogue_retention_probability": "expit(0.25 - 1.2*x), x=(z-0.205)/0.195",
        },
        "validation_sizes": [100],
        "bootstrap_replicates_per_method_per_outer_draw": 3,
        "methods": ["paired_bootstrap", "paired_bootstrap_fpc", "design_aware_srswor"],
    }
    (tmp_path / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return tmp_path


def test_independent_nested_artifact_auditor_accepts_runner_output(nested_artifact):
    result = audit(nested_artifact)
    assert result["status"] == "STRUCTURAL_AND_SUMMARY_AUDIT_PASS"
    assert result["calibration_rows"] == 6
    assert result["attempt_rows"] == 6


def test_independent_nested_artifact_auditor_rejects_tampered_summary(nested_artifact):
    path = nested_artifact / "scenario_summary.csv"
    frame = pd.read_csv(path)
    frame.loc[0, "mean_interval_width"] += 0.25
    frame.to_csv(path, index=False)
    with pytest.raises(ValueError, match="independently reconstructed mean_interval_width differs"):
        audit(nested_artifact)


def test_independent_nested_artifact_auditor_rejects_non_nested_outer_sample(nested_artifact):
    path = nested_artifact / "attempt_audit.jsonl"
    frame = pd.read_json(path, lines=True)
    row = frame.index[(frame["method"] == "paired_bootstrap")][0]
    frame.at[row, "validation_indices"] = frame.at[row, "validation_indices"][:-1]
    frame.to_json(path, orient="records", lines=True)
    with pytest.raises(ValueError, match="outer validation index length mismatch"):
        audit(nested_artifact)


def test_independent_nested_artifact_auditor_rejects_wrong_bootstrap_seed(nested_artifact):
    path = nested_artifact / "attempt_audit.jsonl"
    frame = pd.read_json(path, lines=True)
    row = frame.index[(frame["method"] == "paired_bootstrap")][0]
    frame.at[row, "bootstrap_seed"] = int(frame.at[row, "bootstrap_seed"]) + 1
    frame.to_json(path, orient="records", lines=True)
    with pytest.raises(ValueError, match="bootstrap seed mismatch"):
        audit(nested_artifact)


def test_independent_nested_artifact_auditor_rejects_corrupt_pseudopopulation_mapping(nested_artifact):
    path = nested_artifact / "attempt_audit.jsonl"
    frame = pd.read_json(path, lines=True)
    row = frame.index[(frame["method"] == "design_aware_srswor")][0]
    frame.at[row, "bootstrap_sample_indices"][0][0] = (
        int(frame.at[row, "bootstrap_sample_indices"][0][0]) + 1
    ) % 100
    frame.to_json(path, orient="records", lines=True)
    with pytest.raises(ValueError, match="selected pseudo-population positions do not map"):
        audit(nested_artifact)


def test_independent_nested_replay_reconstructs_smoke_scenario_and_posteriors(nested_artifact):
    result = audit_replay(nested_artifact)
    assert result["status"] == "INDEPENDENT_SCENARIO_AND_POSTERIOR_REPLAY_PASS"
    assert result["scenarios_reconstructed"] == 1
    assert result["outer_fits_reconstructed"] == 2
    assert result["successful_bootstrap_methods_replayed"] >= 2


def test_independent_nested_replay_rejects_tampered_parent_fingerprint(nested_artifact):
    path = nested_artifact / "attempt_audit.jsonl"
    frame = pd.read_json(path, lines=True)
    frame.loc[:, "parent_frame_sha256"] = "0" * 64
    frame.to_json(path, orient="records", lines=True)
    with pytest.raises(AssertionError, match="parent-frame fingerprint mismatch"):
        audit_replay(nested_artifact)


def test_independent_nested_replay_rejects_tampered_bootstrap_posterior_mean(nested_artifact):
    path = nested_artifact / "attempt_audit.jsonl"
    frame = pd.read_json(path, lines=True)
    records = pd.read_csv(nested_artifact / "calibration_records.csv")
    successful = records.loc[
        records["method"].eq("paired_bootstrap")
        & records["bootstrap_status"].eq("success")
    ]
    assert not successful.empty
    key = successful.iloc[0][["scenario", "validation_size", "outer_draw", "method"]]
    mask = (
        frame["scenario"].eq(key["scenario"])
        & frame["validation_size"].eq(key["validation_size"])
        & frame["outer_draw"].eq(key["outer_draw"])
        & frame["method"].eq(key["method"])
    )
    row = frame.index[mask][0]
    frame.at[row, "bootstrap_attempt_posterior_mean"][0] += 0.5
    frame.to_json(path, orient="records", lines=True)
    with pytest.raises(AssertionError, match="posterior mean"):
        audit_replay(nested_artifact)


def test_independent_replay_rejects_fabricated_wrapper_failure(nested_artifact):
    records_path = nested_artifact / "calibration_records.csv"
    records = pd.read_csv(records_path)
    selected = records.index[
        records["method"].eq("paired_bootstrap")
        & records["bootstrap_status"].eq("success")
    ][0]
    records.at[selected, "bootstrap_status"] = "failed"
    records.to_csv(records_path, index=False)
    with pytest.raises(AssertionError, match="wrapper failure|did not reproduce"):
        audit_replay(nested_artifact)
