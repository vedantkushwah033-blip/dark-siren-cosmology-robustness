import importlib.util
from pathlib import Path

import numpy as np
import pytest

from dark_siren_robustness.hardened_inference import posterior_h0, posterior_summary

LEGACY_PATH = Path(__file__).resolve().parents[2] / "bayesian-dark-siren-cosmology" / "src" / "bayesian_inference.py"


def test_hardened_matches_pinned_legacy_on_valid_inputs():
    # Independent reproduction of the pinned legacy calculation avoids a runtime
    # dependency on a second repository checkout.
    from scipy.integrate import trapezoid
    from scipy.stats import norm

    grid = np.linspace(50.0, 90.0, 1601)
    z = np.array([0.015, 0.022, 0.031, 0.044])
    weights = np.array([1.0, 2.0, 0.5, 1.5])
    distance, sigma = 140.0, 18.0
    missing_z = np.array([0.01, 0.025, 0.05])
    missing_weights = np.array([0.4, 0.4, 0.2])
    missing_fraction = 0.15

    observed_weights = weights / weights.sum() * (1.0 - missing_fraction)
    legacy_likelihood = np.sum(
        observed_weights[None, :]
        * norm.pdf(distance, loc=299792.458 * z[None, :] / grid[:, None], scale=sigma),
        axis=1,
    )
    legacy_likelihood += missing_fraction * np.sum(
        (missing_weights / missing_weights.sum())[None, :]
        * norm.pdf(distance, loc=299792.458 * missing_z[None, :] / grid[:, None], scale=sigma),
        axis=1,
    )
    legacy_likelihood = np.where((grid >= 50.0) & (grid <= 90.0), legacy_likelihood, 0.0)
    legacy_posterior = legacy_likelihood / trapezoid(legacy_likelihood, grid)

    hardened_posterior = posterior_h0(
        grid, z, distance, sigma, weights=weights,
        missing_weight=missing_fraction, missing_redshifts=missing_z,
        missing_weights=missing_weights,
    )
    np.testing.assert_allclose(hardened_posterior, legacy_posterior, rtol=1e-12, atol=1e-14)

    legacy_summary = _legacy_summary(grid, legacy_posterior)
    hardened_summary = posterior_summary(grid, hardened_posterior)
    for key in legacy_summary:
        assert hardened_summary[key] == pytest.approx(legacy_summary[key], abs=1e-10)


def _legacy_summary(grid, posterior):
    from scipy.integrate import trapezoid

    mean = trapezoid(grid * posterior, grid)
    variance = trapezoid((grid - mean) ** 2 * posterior, grid)
    cdf = np.concatenate(([0.0], np.cumsum((posterior[1:] + posterior[:-1]) * np.diff(grid) / 2.0)))
    cdf /= cdf[-1]
    return {
        "mean": mean,
        "median": np.interp(0.5, cdf, grid),
        "posterior_sd": np.sqrt(max(variance, 0.0)),
        "lower_68": np.interp(0.16, cdf, grid),
        "upper_68": np.interp(0.84, cdf, grid),
        "lower_95": np.interp(0.025, cdf, grid),
        "upper_95": np.interp(0.975, cdf, grid),
    }
