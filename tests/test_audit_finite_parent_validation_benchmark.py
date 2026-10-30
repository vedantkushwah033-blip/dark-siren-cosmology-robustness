from __future__ import annotations

import json

import pytest

from dark_siren_robustness.finite_parent_sampling import run_finite_parent_benchmark
from scripts.audit_finite_parent_validation_benchmark import audit


def _artifact(tmp_path):
    frame = run_finite_parent_benchmark(
        n_scenarios=1, n_validation_draws=2, n_parent=1200,
        n_detected=6, validation_sizes=(100, 200), seed=20261033,
        h0_points=41, z_points=201,
    )
    frame.to_csv(tmp_path / "records.csv", index=False)
    summary = {
        "schema": "finite_parent_srswor_smoke_v1",
        "status": "SMOKE_ONLY_NOT_CONFIRMATORY",
        "n_attempted": len(frame),
        "n_successful": int(frame["status"].eq("success").sum()),
        "n_failed": int(frame["status"].eq("failed").sum()),
        "validation_sizes": [100, 200],
        "scenarios": 1,
        "draws_per_scenario_and_size": 2,
    }
    (tmp_path / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    return frame


def test_finite_parent_auditor_replays_srswor_and_provenance(tmp_path):
    frame = _artifact(tmp_path)
    result = audit(tmp_path)
    assert result["status"] == "STRUCTURAL_AUDIT_PASS"
    assert result["n_attempted"] == 4
    for draw in range(2):
        small = frame.loc[(frame["draw"] == draw) & (frame["validation_size"] == 100)].iloc[0]
        large = frame.loc[(frame["draw"] == draw) & (frame["validation_size"] == 200)].iloc[0]
        small_indices = [int(v) for v in small["validation_indices"].split(",")]
        large_indices = [int(v) for v in large["validation_indices"].split(",")]
        assert large_indices[:100] == small_indices


def test_finite_parent_auditor_rejects_repeated_validation_index(tmp_path):
    frame = _artifact(tmp_path)
    encoded = frame.iloc[0]["validation_indices"]
    indices = encoded.split(",")
    indices[1] = indices[0]
    frame.loc[0, "validation_indices"] = ",".join(indices)
    frame.to_csv(tmp_path / "records.csv", index=False)
    with pytest.raises(AssertionError, match="exact SRSWOR"):
        audit(tmp_path)


def test_finite_parent_auditor_rejects_provenance_mismatch(tmp_path):
    frame = _artifact(tmp_path)
    frame.loc[0, "catalogue_sha256"] = "0" * 64
    frame.to_csv(tmp_path / "records.csv", index=False)
    with pytest.raises(AssertionError, match="inconsistent"):
        audit(tmp_path)
