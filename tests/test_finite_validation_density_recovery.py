import hashlib
import json

import numpy as np
import pandas as pd
import pytest

from dark_siren_robustness.finite_validation_density_recovery import (
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
    assert "no_H0_inference" in summary["status"]
    assert summary["paired_ipw_minus_catalogue_only_contrasts"]


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


def _refresh_manifest_file(tmp_path, key, filename):
    manifest_path = tmp_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    payload = (tmp_path / filename).read_bytes()
    manifest["files"][key]["bytes"] = len(payload)
    manifest["files"][key]["sha256"] = hashlib.sha256(payload).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\\n",
                             encoding="utf-8")


def test_auditor_rejects_missing_expected_attempt_even_if_manifest_is_refreshed(tmp_path):
    attempts, units, summary = run_density_recovery_benchmark(
        n_replicates=1, n_parent=1200, seed=404,
        gamma_values=(0.0,), validation_sizes=(100,),
        degrees=(1,), bandwidth_multipliers=(0.15,), grid_points=201,
    )
    write_artifacts(attempts, units, summary, tmp_path)
    path = tmp_path / "attempt_metrics.csv"
    frame = pd.read_csv(path).iloc[1:]
    frame.to_csv(path, index=False)
    _refresh_manifest_file(tmp_path, "attempts", "attempt_metrics.csv")
    with pytest.raises(AssertionError, match="Attempt design mismatch"):
        audit_artifact(tmp_path)


def test_auditor_rejects_coefficients_that_do_not_solve_logistic_score(tmp_path):
    attempts, units, summary = run_density_recovery_benchmark(
        n_replicates=1, n_parent=1200, seed=505,
        gamma_values=(0.0,), validation_sizes=(100,),
        degrees=(1,), bandwidth_multipliers=(0.15,), grid_points=201,
    )
    write_artifacts(attempts, units, summary, tmp_path)
    path = tmp_path / "attempt_metrics.csv"
    frame = pd.read_csv(path)
    idx = frame.index[
        (frame["arm"] == "estimated_degree_1_ipw") &
        (frame["status"] == "SUCCESS")
    ][0]
    beta = np.asarray(json.loads(frame.loc[idx, "coefficients_json"]), dtype=float)
    beta[0] += 0.25
    frame.loc[idx, "coefficients_json"] = json.dumps(beta.tolist(), separators=(",", ":"))
    frame.to_csv(path, index=False)
    _refresh_manifest_file(tmp_path, "attempts", "attempt_metrics.csv")
    with pytest.raises(AssertionError, match="Logistic score tolerance failed"):
        audit_artifact(tmp_path)
