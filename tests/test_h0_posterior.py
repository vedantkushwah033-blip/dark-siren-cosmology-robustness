import numpy as np
import pytest

from dark_siren_robustness.h0_posterior import (
    compare_h0_grid_posteriors,
    normalize_h0_posterior,
)


def test_constant_likelihood_recovers_uniform_h0_prior():
    grid = np.linspace(50.0, 90.0, 81)
    posterior = normalize_h0_posterior(grid, np.full(grid.size, -17.0))
    assert posterior.grid[0] == 50.0
    assert posterior.grid[-1] == 90.0
    assert np.allclose(posterior.density, 1.0 / 40.0, atol=1e-14)
    assert np.isclose(posterior.quantile(0.5), 70.0, atol=1e-12)
    assert np.isclose(posterior.quantile(0.05), 52.0, atol=1e-12)
    assert np.isclose(posterior.cdf[-1], 1.0)


def test_h0_posterior_is_invariant_to_common_log_likelihood_offset():
    grid = np.linspace(50.0, 90.0, 101)
    log_like = -0.5 * ((grid - 69.0) / 4.0) ** 2
    first = normalize_h0_posterior(grid, log_like)
    shifted = normalize_h0_posterior(grid, log_like + 12000.0)
    assert np.allclose(first.density, shifted.density, rtol=1e-12, atol=1e-14)
    assert np.isclose(first.quantile(0.5), shifted.quantile(0.5), atol=1e-12)


def test_grid_refinement_compares_density_and_posterior_quantiles():
    coarse_grid = np.linspace(50.0, 90.0, 101)
    fine_grid = np.linspace(50.0, 90.0, 801)
    coarse = normalize_h0_posterior(
        coarse_grid, -0.5 * ((coarse_grid - 71.0) / 5.0) ** 2
    )
    fine = normalize_h0_posterior(
        fine_grid, -0.5 * ((fine_grid - 71.0) / 5.0) ** 2
    )
    check = compare_h0_grid_posteriors(
        coarse, fine, density_tolerance=2e-3, quantile_tolerance=0.02
    )
    assert check.passed
    assert check.max_abs_density_difference < 2e-3
    assert check.max_abs_quantile_difference < 0.02


@pytest.mark.parametrize(
    "grid, log_like, kwargs",
    [
        (np.linspace(55, 85, 31), np.zeros(31), {}),
        (np.linspace(50, 90, 31), np.full(31, np.nan), {}),
        (np.linspace(50, 90, 31), np.zeros(30), {}),
        (np.linspace(50, 90, 31), np.zeros(31), {"prior_min": 90, "prior_max": 50}),
    ],
)
def test_h0_posterior_rejects_invalid_support_or_likelihood(grid, log_like, kwargs):
    with pytest.raises(ValueError):
        normalize_h0_posterior(grid, log_like, **kwargs)


def test_h0_quantile_and_grid_comparison_validate_arguments():
    posterior = normalize_h0_posterior(
        np.linspace(50.0, 90.0, 41), np.zeros(41)
    )
    with pytest.raises(ValueError, match="probability"):
        posterior.quantile(1.1)
    with pytest.raises(ValueError, match="identical prior"):
        compare_h0_grid_posteriors(
            posterior,
            normalize_h0_posterior(
                np.linspace(40.0, 90.0, 51), np.zeros(51),
                prior_min=40.0, prior_max=90.0,
            ),
        )
