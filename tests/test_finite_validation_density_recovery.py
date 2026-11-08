import numpy as np
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
