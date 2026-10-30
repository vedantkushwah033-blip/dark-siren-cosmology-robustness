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
    successful = frame.loc[frame["status"].eq("success")]
    metrics_by_size = {}
    for size in (100, 200):
        group = successful.loc[successful["validation_size"].eq(size)]
        if group.empty:
            continue
        errors = group["posterior_mean"].to_numpy(dtype=float) - group["h0_true"].to_numpy(dtype=float)
        by_scenario = group.assign(error=errors).groupby("scenario")
        scenario_coverage = group.groupby("scenario")["covered_68"].mean()
        metrics_by_size[str(size)] = {
            "posterior_mean_bias": float(errors.mean()),
            "posterior_mean_rmse": float((errors**2).mean()**0.5),
            "finite_parent_target_bias": float(group["error_vs_finite_parent_target"].mean()),
            "finite_parent_target_rmse": float(
                (group["error_vs_finite_parent_target"]**2).mean()**0.5
            ),
            "coverage_68": float(group["covered_68"].astype(float).mean()),
            "mean_ci68_width": float(group["ci68_width"].mean()),
            "mean_density_ise": float(group["density_ise"].mean()),
            "mean_inclusion_probability_mae": float(group["inclusion_probability_mae"].mean()),
            "mean_weight_ess": float(group["weight_ess"].mean()),
            "mean_max_inverse_weight": float(group["max_inverse_weight"].mean()),
            "boundary_limited_fraction": float(group["fit_boundary_limited"].astype(float).mean()),
            "n_scenarios_with_successes": int(group["scenario"].nunique()),
            "scenario_balanced_bias_mean": float(by_scenario["error"].mean().mean()),
            "scenario_balanced_rmse_mean": float(
                group.assign(squared_error=errors**2).groupby("scenario")["squared_error"].mean().pow(0.5).mean()
            ),
            "scenario_balanced_coverage_mean": float(scenario_coverage.mean()),
            "scenario_coverage_sd": float(scenario_coverage.std(ddof=1)) if len(scenario_coverage) > 1 else None,
        }
    summary = {
        "schema": "finite_parent_srswor_smoke_v1",
        "status": "SMOKE_ONLY_NOT_CONFIRMATORY",
        "n_attempted": len(frame),
        "n_successful": int(successful.shape[0]),
        "n_failed": int(frame["status"].eq("failed").sum()),
        "validation_sizes": [100, 200],
        "scenarios": 1,
        "draws_per_scenario_and_size": 2,
        "h0_points": 41,
        "z_points": 201,
        "h0_grid_bounds": [50.0, 90.0],
        "redshift_grid_bounds": [0.01, 0.40],
        "metrics_by_validation_size": metrics_by_size,
    }
    (tmp_path / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    return frame


def test_finite_parent_auditor_replays_srswor_and_provenance(tmp_path):
    frame = _artifact(tmp_path)
    result = audit(tmp_path)
    assert result["status"] == "INDEPENDENT_NUMERICAL_AUDIT_PASS"
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


def test_finite_parent_auditor_rejects_tampered_posterior_metric(tmp_path):
    frame = _artifact(tmp_path)
    successful = frame.index[frame["status"].eq("success")]
    assert len(successful) > 0
    frame.loc[successful[0], "posterior_mean"] += 0.5
    frame.to_csv(tmp_path / "records.csv", index=False)
    with pytest.raises(AssertionError, match="posterior mean"):
        audit(tmp_path)


def test_finite_parent_auditor_requires_grid_manifest(tmp_path):
    _artifact(tmp_path)
    manifest = json.loads((tmp_path / "summary.json").read_text(encoding="utf-8"))
    del manifest["h0_points"]
    (tmp_path / "summary.json").write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(AssertionError, match="grid specification"):
        audit(tmp_path)

def test_finite_parent_auditor_rejects_reproducibly_successful_attempt_marked_failed(tmp_path):
    frame = _artifact(tmp_path)
    successful = frame.index[frame["status"].eq("success")]
    assert len(successful) > 0
    index = successful[0]
    frame.loc[index, "status"] = "failed"
    frame.loc[index, "failure_code"] = "unexpected_exception"
    frame.loc[index, "failure_reason"] = "injected false failure"
    frame.to_csv(tmp_path / "records.csv", index=False)

    summary_path = tmp_path / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["n_successful"] -= 1
    summary["n_failed"] += 1
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    with pytest.raises(AssertionError, match="independently reproducible as a successful"):
        audit(tmp_path)


def test_finite_parent_auditor_rejects_failed_attempt_without_reason(tmp_path):
    frame = _artifact(tmp_path)
    successful = frame.index[frame["status"].eq("success")]
    assert len(successful) > 0
    index = successful[0]
    frame.loc[index, "status"] = "failed"
    frame.loc[index, "failure_code"] = "unexpected_exception"
    frame.loc[index, "failure_reason"] = ""
    frame.to_csv(tmp_path / "records.csv", index=False)

    summary_path = tmp_path / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["n_successful"] -= 1
    summary["n_failed"] += 1
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    with pytest.raises(AssertionError, match="missing its failure reason"):
        audit(tmp_path)



def test_finite_parent_auditor_rejects_tampered_finite_frame_target(tmp_path):
    frame = _artifact(tmp_path)
    successful = frame.index[frame["status"].eq("success")]
    assert len(successful) > 0
    frame.loc[successful[0], "finite_parent_target_posterior_mean"] += 0.25
    frame.to_csv(tmp_path / "records.csv", index=False)
    with pytest.raises(AssertionError, match="finite-parent target posterior mean"):
        audit(tmp_path)



def test_failure_reason_audit_requires_exact_replayed_exception():
    from scripts.audit_finite_parent_validation_benchmark import (
        _assert_failure_reason_matches,
    )

    error = RuntimeError("independent weighted variance invalid")
    _assert_failure_reason_matches(
        "RuntimeError: independent weighted variance invalid", error
    )
    with pytest.raises(AssertionError, match="does not match independently reconstructed"):
        _assert_failure_reason_matches("RuntimeError: unrelated failure", error)
    with pytest.raises(AssertionError, match="does not match independently reconstructed"):
        _assert_failure_reason_matches("ValueError: independent weighted variance invalid", error)

def test_failure_code_categories_are_stable_and_ordered_by_exception_specificity():
    from scripts.audit_finite_parent_validation_benchmark import _independent_failure_code

    assert _independent_failure_code(FloatingPointError("overflow")) == "floating_point_failure"
    assert _independent_failure_code(ArithmeticError("invalid")) == "arithmetic_failure"
    assert _independent_failure_code(ValueError("one class")) == "invalid_validation_sample"
    assert _independent_failure_code(RuntimeError("fit failed")) == "model_fit_failure"
    assert _independent_failure_code(KeyError("unexpected")) == "unexpected_exception"
