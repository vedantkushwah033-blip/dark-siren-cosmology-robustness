from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.stats import norm
from scipy.integrate import trapezoid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_full_parent_h0_sbc_log_grid import normalize_posterior_on_h0_grid


def test_geometric_h0_grid_normalizes_uniform_h0_prior_density():
    grid = np.geomspace(20.0, 300.0, 2241)
    log_likelihood = -0.5 * ((grid - 70.0) / 8.0) ** 2

    density, cdf = normalize_posterior_on_h0_grid(grid, log_likelihood)

    assert np.isclose(trapezoid(density, grid), 1.0, atol=1e-12)
    assert np.isclose(cdf[0], 0.0)
    assert np.isclose(cdf[-1], 1.0)
    assert np.all(np.diff(cdf) >= 0.0)


def test_geometric_h0_grid_quantiles_match_truncated_gaussian_reference():
    grid = np.geomspace(20.0, 300.0, 2241)
    log_likelihood = -0.5 * ((grid - 70.0) / 8.0) ** 2
    density, cdf = normalize_posterior_on_h0_grid(grid, log_likelihood)

    probabilities = np.array([0.025, 0.16, 0.50, 0.84, 0.975])
    observed = np.interp(probabilities, cdf, grid)
    # The Gaussian tails outside [20, 300] are negligible for this test.
    expected = 70.0 + 8.0 * norm.ppf(probabilities)

    assert np.max(np.abs(observed - expected)) < 1e-3
    assert np.isclose(np.trapz(density, grid), 1.0, atol=1e-12)


def test_geometric_grid_integration_is_with_respect_to_h0_not_log_h0():
    grid = np.geomspace(20.0, 300.0, 2241)
    log_likelihood = np.zeros_like(grid)
    density, cdf = normalize_posterior_on_h0_grid(grid, log_likelihood)

    assert np.allclose(density, 1.0 / 280.0, rtol=0, atol=1e-12)
    assert np.isclose(np.interp(0.5, cdf, grid), 160.0, atol=1e-10)
