import hashlib
import json

import numpy as np
import pandas as pd
import pytest

from dark_siren_robustness.finite_validation_density_recovery import (
    _paired_contrast_summary,
    _calibrate_intercept,
    _curvature_center,
    _analytic_target,
    run_density_recovery_benchmark,
    write_artifacts,
)
from dark_siren_robustness.finite_validation_density_recovery_audit import audit_artifact


def test_parent_curvature_center_uses_declared_nonuniform_parent_law():
    center = _curvature_center()
    assert 0.0 < center < 1.0
    assert not np.isclose(center, 1.0 / 3.0, atol=1e-3)


def test_intercept_calibration_is_finite_and_curvature_specific():
    center = _curvature_center()
    values = [_calibrate_intercept(g, 0.25, center) for g in (-1.2, 0.0, 1.2)]
    assert np.all(np.isfinite(values))
    assert len(set(np.round(values, 6))) == 3


def test_analytic_target_is_normalized_and_positive():
    grid = np.linspace(0.01, 0.40, 801)
    target = _analytic_target(grid)
    assert np.all(target > 0)
    assert np.trapezoid(target, grid) == pytest.approx(1.0, abs=1e-10)


def test_small_benchmark_has_comparison_arms_and_nested_splits():
    attempts, units, summary = run_density_recovery_benchmark(
        n_replicates=1, n_parent=1200, seed=101,
        gamma_values=(0.0,), validation_sizes=(100, 300),
        degrees=(1, 2), bandwidth_multipliers=(0.15,), grid_points=201,
    )
    expected = {
        "catalogue_only", "known_inclusion_ipw", "parent_sample_kde",
        "estimated_degree_1_ipw", "estimated_degree_2_ipw",
    }
    assert set(attempts.arm) == expected
    assert set(attempts.status).issubset({"SUCCESS", "FAILED", "NON_EVALUABLE"})
    assert len(units) == 2 * 1200
    small = units.loc[units.n_validation == 100]
    large = units.loc[units.n_validation == 300]
    small_ids = set(small.loc[small.role == "VALIDATION", "object_id"])
    large_ids = set(large.loc[large.role == "VALIDATION", "object_id"])
    assert small_ids.issubset(large_ids)
    assert summary["attempt_count"] == len(attempts)
    assert "SeedSequence([seed, replicate, stream_id])" in summary["rng_stream_protocol"]
    assert "no_H0_inference" in summary["status"]
    assert summary["paired_ipw_minus_catalogue_only_contrasts"]


def test_paired_contrast_reports_failed_and_non_evaluable_attempts():
    rows = []
    statuses = {
        0: ("SUCCESS", 0.10, "SUCCESS", 0.20),
        1: ("SUCCESS", 0.11, "FAILED", np.nan),
        2: ("FAILED", np.nan, "SUCCESS", 0.25),
        3: ("NON_EVALUABLE", np.nan, "NON_EVALUABLE", np.nan),
        4: ("NOT_RUN", np.nan, "NOT_RUN", np.nan),
    }
    for replicate, (cat_status, cat_iae, est_status, est_iae) in statuses.items():
        rows.append({
            "replicate": replicate, "gamma": 1.2, "target_retention": 0.25,
            "n_validation": 100, "bandwidth_multiplier": 0.15,
            "arm": "catalogue_only", "status": cat_status, "iae": cat_iae,
            "n_analysis_catalogue": 1000 - 100 * replicate,
        })
        rows.append({
            "replicate": replicate, "gamma": 1.2, "target_retention": 0.25,
            "n_validation": 100, "bandwidth_multiplier": 0.15,
            "arm": "estimated_degree_1_ipw", "status": est_status, "iae": est_iae,
            "n_analysis_catalogue": 1000 - 100 * replicate,
        })
    result = _paired_contrast_summary(
        pd.DataFrame(rows), n_replicates=5, bootstrap_resamples=100,
    )[0]
    assert result["n_pairs"] == 1
    assert result["n_replicates_with_both_attempts"] == 5
    assert result["mean_n_analysis_catalogue"] == pytest.approx(800.0)
    assert result["sd_n_analysis_catalogue"] == pytest.approx(np.std([1000, 900, 800, 700, 600], ddof=1))
    assert result["min_n_analysis_catalogue"] == 600
    assert result["max_n_analysis_catalogue"] == 1000
    assert result["n_catalogue_size_records"] == 5
    assert result["n_catalogue_only_success"] == 2
    assert result["n_catalogue_only_failed"] == 1
    assert result["n_catalogue_only_non_evaluable"] == 1
    assert result["n_catalogue_only_not_run"] == 1
    assert result["n_estimated_success"] == 2
    assert result["n_estimated_failed"] == 1
    assert result["n_estimated_non_evaluable"] == 1
    assert result["n_estimated_not_run"] == 1
    assert result["n_both_success"] == 1
    assert result["n_estimated_not_success_catalogue_success"] == 1
    assert result["n_catalogue_not_success_estimated_success"] == 1
    assert result["n_both_not_success"] == 2
    assert result["conditional_pair_fraction_of_planned"] == pytest.approx(0.20)


def test_independent_auditor_rejects_tampered_attempt_status_summary(tmp_path):
    attempts, units, summary = run_density_recovery_benchmark(
        n_replicates=1, n_parent=1200, seed=20261109,
        gamma_values=(0.0,), validation_sizes=(100,),
        degrees=(1,), bandwidth_multipliers=(0.15,), grid_points=201,
    )
    write_artifacts(attempts, units, summary, tmp_path)
    path = tmp_path / "summary.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["paired_ipw_minus_catalogue_only_contrasts"][0]["n_estimated_failed"] += 1
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    manifest_path = tmp_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    summary_bytes = path.read_bytes()
    manifest["files"]["summary"]["bytes"] = len(summary_bytes)
    manifest["files"]["summary"]["sha256"] = hashlib.sha256(summary_bytes).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    with pytest.raises(AssertionError, match="Paired contrast key mismatch"):
        audit_artifact(tmp_path)


def test_independent_auditor_rejects_tampered_catalogue_size_summary(tmp_path):
    attempts, units, summary = run_density_recovery_benchmark(
        n_replicates=1, n_parent=1200, seed=20261110,
        gamma_values=(0.0,), validation_sizes=(100,),
        degrees=(1,), bandwidth_multipliers=(0.15,), grid_points=201,
    )
    write_artifacts(attempts, units, summary, tmp_path)
    path = tmp_path / "summary.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    payload["paired_ipw_minus_catalogue_only_contrasts"][0]["mean_n_analysis_catalogue"] += 1.0
    path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    manifest_path = tmp_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    summary_bytes = path.read_bytes()
    manifest["files"]["summary"]["bytes"] = len(summary_bytes)
    manifest["files"]["summary"]["sha256"] = hashlib.sha256(summary_bytes).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    with pytest.raises(AssertionError, match="Paired contrast mismatch: mean_n_analysis_catalogue"):
        audit_artifact(tmp_path)


def test_independent_auditor_reconstructs_smoke_artifact(tmp_path):
    attempts, units, summary = run_density_recovery_benchmark(
        n_replicates=1, n_parent=1200, seed=202,
        gamma_values=(0.0,), validation_sizes=(100,),
        degrees=(1, 2), bandwidth_multipliers=(0.15,), grid_points=201,
    )
    write_artifacts(attempts, units, summary, tmp_path)
    report = audit_artifact(tmp_path)
    assert report["status"] == "PASS"
    assert report["checked_parent_partitions"] == 1
    assert report["checked_rng_streams"] == 1
    assert report["checked_density_estimates"] >= 3


def test_rejects_invalid_parent_size_and_validation_order():
    with pytest.raises(ValueError):
        run_density_recovery_benchmark(n_replicates=1, n_parent=100)
    with pytest.raises(ValueError, match="increasing"):
        run_density_recovery_benchmark(
            n_replicates=1, n_parent=1200, validation_sizes=(300, 100),
        )


def test_auditor_rejects_tampered_parent_unit_artifact(tmp_path):
    attempts, units, summary = run_density_recovery_benchmark(
        n_replicates=1, n_parent=1200, seed=303,
        gamma_values=(0.0,), retention_values=(0.50,),
        validation_sizes=(100,), degrees=(1,),
        bandwidth_multipliers=(0.15,), grid_points=201,
    )
    write_artifacts(attempts, units, summary, tmp_path)
    path = tmp_path / "parent_units.csv"
    path.write_text(path.read_text(encoding="utf-8") + "\\n", encoding="utf-8")
    with pytest.raises(AssertionError, match="byte count|checksum"):
        audit_artifact(tmp_path)



def _refresh_attempt_manifest(tmp_path):
    manifest_path = tmp_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload = (tmp_path / "attempt_metrics.csv").read_bytes()
    manifest["files"]["attempts"]["bytes"] = len(payload)
    manifest["files"]["attempts"]["sha256"] = hashlib.sha256(payload).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")



def _refresh_unit_manifest(tmp_path):
    manifest_path = tmp_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload = (tmp_path / "parent_units.csv").read_bytes()
    manifest["files"]["units"]["bytes"] = len(payload)
    manifest["files"]["units"]["sha256"] = hashlib.sha256(payload).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")


def test_auditor_allows_subtolerance_csv_float_roundtrip(tmp_path):
    attempts, units, summary = run_density_recovery_benchmark(
        n_replicates=1, n_parent=1200, seed=20261110,
        gamma_values=(0.0,), retention_values=(0.50,),
        validation_sizes=(100,), degrees=(1,),
        bandwidth_multipliers=(0.15,), grid_points=201,
    )
    write_artifacts(attempts, units, summary, tmp_path)
    path = tmp_path / "parent_units.csv"
    frame = pd.read_csv(path)
    margin = np.abs(
        frame["retention_uniform"].to_numpy(dtype=float)
        - frame["true_inclusion_probability"].to_numpy(dtype=float)
    )
    idx = int(np.flatnonzero(margin > 1e-8)[0])
    original = float(frame.loc[idx, "retention_uniform"])
    frame.loc[idx, "retention_uniform"] = np.nextafter(original, 0.0)
    frame.to_csv(path, index=False)
    _refresh_unit_manifest(tmp_path)
    report = audit_artifact(tmp_path)
    assert report["status"] == "PASS"


def test_auditor_rejects_reused_or_tampered_retention_rng_stream(tmp_path):
    attempts, units, summary = run_density_recovery_benchmark(
        n_replicates=1, n_parent=1200, seed=20261108,
        gamma_values=(0.0,), retention_values=(0.50,),
        validation_sizes=(100,), degrees=(1,),
        bandwidth_multipliers=(0.15,), grid_points=201,
    )
    write_artifacts(attempts, units, summary, tmp_path)
    path = tmp_path / "parent_units.csv"
    frame = pd.read_csv(path)
    idx = frame.index[0]
    probability = float(frame.loc[idx, "true_inclusion_probability"])
    retained = int(frame.loc[idx, "retained"])
    # Preserve the DGP label and role while deliberately changing the raw
    # uniform. The independent auditor must detect that the declared RNG stream
    # was not used, even after the manifest is refreshed.
    frame.loc[idx, "retention_uniform"] = (
        probability / 2.0 if retained else (1.0 + probability) / 2.0
    )
    frame.to_csv(path, index=False)
    _refresh_unit_manifest(tmp_path)
    with pytest.raises(AssertionError, match="Retention uniforms do not match their independent RNG stream"):
        audit_artifact(tmp_path)



def test_auditor_rejects_missing_attempt_even_with_refreshed_manifest(tmp_path):
    attempts, units, summary = run_density_recovery_benchmark(
        n_replicates=1, n_parent=1200, seed=606,
        gamma_values=(0.0,), validation_sizes=(100,),
        degrees=(1,), bandwidth_multipliers=(0.15,), grid_points=201,
    )
    write_artifacts(attempts, units, summary, tmp_path)
    path = tmp_path / "attempt_metrics.csv"
    pd.read_csv(path).iloc[1:].to_csv(path, index=False)
    _refresh_attempt_manifest(tmp_path)
    with pytest.raises(AssertionError, match="Attempt design mismatch"):
        audit_artifact(tmp_path)


def test_auditor_rejects_incorrect_logistic_objective_even_with_refreshed_manifest(tmp_path):
    attempts, units, summary = run_density_recovery_benchmark(
        n_replicates=1, n_parent=2000, seed=101,
        gamma_values=(0.0,), retention_values=(0.50,),
        validation_sizes=(1000,), degrees=(1,),
        bandwidth_multipliers=(0.15,), grid_points=201,
    )
    write_artifacts(attempts, units, summary, tmp_path)
    path = tmp_path / "attempt_metrics.csv"
    frame = pd.read_csv(path)
    successful = frame.loc[
        (frame["arm"] == "estimated_degree_1_ipw") &
        (frame["status"] == "SUCCESS")
    ]
    assert not successful.empty, (
        "Objective-tampering fixture did not produce a successful fitted-IPW arm: "
        + frame.loc[frame["arm"] == "estimated_degree_1_ipw",
                    ["status", "message", "fit_status", "fit_message",
                     "fit_max_abs_coefficient", "fit_hessian_condition_number",
                     "fit_gradient_inf_norm"]].to_string(index=False)
    )
    idx = successful.index[0]
    frame.loc[idx, "fit_objective"] += 1.0
    frame.to_csv(path, index=False)
    _refresh_attempt_manifest(tmp_path)
    with pytest.raises(AssertionError, match="Fitted logistic objective mismatch"):
        audit_artifact(tmp_path)


def test_auditor_rejects_tampered_logistic_coefficient_even_with_refreshed_manifest(tmp_path):
    attempts, units, summary = run_density_recovery_benchmark(
        n_replicates=1, n_parent=2000, seed=101,
        gamma_values=(0.0,), retention_values=(0.50,),
        validation_sizes=(1000,), degrees=(1,),
        bandwidth_multipliers=(0.15,), grid_points=201,
    )
    write_artifacts(attempts, units, summary, tmp_path)
    path = tmp_path / "attempt_metrics.csv"
    frame = pd.read_csv(path)
    successful = frame.loc[
        (frame["arm"] == "estimated_degree_1_ipw") &
        (frame["status"] == "SUCCESS")
    ]
    assert not successful.empty, (
        "Coefficient-tampering fixture did not produce a successful fitted-IPW arm."
    )
    idx = successful.index[0]
    coefficients = np.asarray(
        json.loads(frame.loc[idx, "coefficients_json"]), dtype=float
    )
    coefficients[0] += 0.05
    frame.loc[idx, "coefficients_json"] = json.dumps(
        coefficients.tolist(), separators=(",", ":")
    )
    frame.to_csv(path, index=False)
    _refresh_attempt_manifest(tmp_path)
    with pytest.raises(
        AssertionError,
        match="Fitted logistic score mismatch|Fitted logistic objective mismatch|"
              "Fitted covariance reconstruction failed",
    ):
        audit_artifact(tmp_path)
