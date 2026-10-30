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
    frame.at[row, "bootstrap_sample_indices"][0] = (
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
