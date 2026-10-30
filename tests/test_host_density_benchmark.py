import numpy as np
import pytest
from dark_siren_robustness.host_density_benchmark import posterior_summary, run_benchmark

def test_posterior_summary_normalizes_and_recovers_symmetric_peak():
    grid = np.linspace(50, 90, 1001)
    result = posterior_summary(grid, -0.5*((grid-70)/2)**2)
    assert result["posterior_integral"] == pytest.approx(1.0, abs=1e-10)
    assert result["posterior_mean"] == pytest.approx(70.0, abs=1e-3)
    assert result["ci68_low"] < 70 < result["ci68_high"]


def test_posterior_summary_matches_independent_trapezoidal_reconstruction():
    import json
    from scipy.integrate import cumulative_trapezoid, trapezoid

    # An asymmetric log-likelihood on a nonuniform grid exercises normalization,
    # posterior mean, and equal-tailed quantiles independently of the helper.
    grid = np.array([50.0, 51.0, 53.5, 57.0, 62.0, 68.0, 73.0, 79.0, 85.0, 90.0])
    loglike = -0.5 * ((grid - 69.0) / 5.5) ** 2 + 0.012 * (grid - 69.0)

    result = posterior_summary(grid, loglike)
    density = np.exp(loglike - np.max(loglike))
    density /= trapezoid(density, grid)
    cdf = np.r_[0.0, cumulative_trapezoid(density, grid)]
    cdf /= cdf[-1]
    expected_quantiles = [np.interp(q, cdf, grid) for q in (0.16, 0.50, 0.84)]
    stored_density = np.asarray(json.loads(result["posterior_density_grid_json"]))

    assert stored_density == pytest.approx(density, rel=1e-12, abs=1e-14)
    assert trapezoid(stored_density, grid) == pytest.approx(1.0, abs=1e-12)
    assert result["posterior_mean"] == pytest.approx(trapezoid(grid * density, grid), abs=1e-12)
    assert result["ci68_low"] == pytest.approx(expected_quantiles[0], abs=1e-12)
    assert result["posterior_median"] == pytest.approx(expected_quantiles[1], abs=1e-12)
    assert result["ci68_high"] == pytest.approx(expected_quantiles[2], abs=1e-12)
    assert result["ci68_width"] == pytest.approx(expected_quantiles[2] - expected_quantiles[0], abs=1e-12)

def test_three_arms_share_each_replicate_and_return_valid_posteriors():
    rows, summary = run_benchmark(n_replicates=1, n_parent=1200, n_detected=8,
                                  seed=55123, h0_points=51, z_points=201)
    assert set(rows.arm) == {"oracle", "catalogue_only", "known_completeness"}
    assert len(rows) == 3
    assert rows.posterior_integral.between(0.999999, 1.000001).all()
    assert rows.ci68_width.gt(0).all()
    assert len(summary["aggregate_metrics"]) == 3
    assert summary["status"] == "synthetic_exploratory_benchmark_not_real_H0_measurement"
    assert summary["bandwidth_factor"] is None
    assert "weighted Scott rule" in summary["bandwidth_rule"]

def test_rejects_too_coarse_grid():
    with pytest.raises(ValueError, match="grid"):
        run_benchmark(n_replicates=1, h0_points=10)

@pytest.mark.parametrize(
    "grid, loglike",
    [
        ([50.0, 49.0, 90.0], [0.0, 1.0, 0.0]),
        ([50.0, 60.0, 60.0], [0.0, 1.0, 0.0]),
        ([50.0, np.nan, 90.0], [0.0, 1.0, 0.0]),
        ([50.0, 70.0, 90.0], [0.0, np.inf, 0.0]),
        ([50.0, 70.0, 90.0], [0.0, 1.0]),
    ],
)
def test_posterior_summary_rejects_malformed_inputs(grid, loglike):
    with pytest.raises(ValueError):
        posterior_summary(grid, loglike)
