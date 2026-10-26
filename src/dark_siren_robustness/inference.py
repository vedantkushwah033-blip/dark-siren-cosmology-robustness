"""A deliberately simplified low-redshift dark-siren inference baseline.

This module is for controlled known-answer tests. It is not a complete
hierarchical gravitational-wave cosmology likelihood and does not model GW
detection selection.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import trapezoid
from scipy.special import logsumexp
from scipy.stats import norm


def posterior_h0_low_z(
    h0_grid: ArrayLike,
    redshifts: ArrayLike,
    observed_distance_mpc: float,
    distance_sigma_mpc: float,
    host_weights: ArrayLike | None = None,
    prior_min: float = 50.0,
    prior_max: float = 90.0,
) -> NDArray[np.float64]:
    """Compute a normalized H0 posterior under d_L ~= c*z/H0.

    Candidate-host likelihoods are combined as a weighted Gaussian mixture.
    A uniform prior is applied on [prior_min, prior_max]. This is a toy
    validation model: it omits full cosmological distances, population
    inference, and gravitational-wave detection-selection normalization.
    """
    grid = np.asarray(h0_grid, dtype=float)
    z = np.asarray(redshifts, dtype=float)

    if grid.ndim != 1 or grid.size < 2:
        raise ValueError("h0_grid must be a one-dimensional array with at least two points")
    if np.any(~np.isfinite(grid)) or np.any(np.diff(grid) <= 0):
        raise ValueError("h0_grid must be finite and strictly increasing")
    if z.ndim != 1 or z.size == 0:
        raise ValueError("redshifts must be a non-empty one-dimensional array")
    if np.any(~np.isfinite(z)) or np.any(z < 0):
        raise ValueError("redshifts must be finite and non-negative")
    if not np.isfinite(observed_distance_mpc) or observed_distance_mpc <= 0:
        raise ValueError("observed_distance_mpc must be finite and positive")
    if not np.isfinite(distance_sigma_mpc) or distance_sigma_mpc <= 0:
        raise ValueError("distance_sigma_mpc must be finite and positive")
    if (
        not np.isfinite(prior_min)
        or not np.isfinite(prior_max)
        or prior_min >= prior_max
    ):
        raise ValueError("prior bounds must be finite and prior_min < prior_max")

    if host_weights is None:
        weights = np.ones(z.size, dtype=float)
    else:
        weights = np.asarray(host_weights, dtype=float)
    if (
        weights.ndim != 1
        or weights.size != z.size
        or np.any(~np.isfinite(weights))
        or np.any(weights < 0)
        or not np.any(weights > 0)
    ):
        raise ValueError("host_weights must be finite, non-negative, and match redshifts")

    # Avoid overflow when summing large but finite relative host weights.
    scaled_weights = weights / np.max(weights)
    weights = scaled_weights / scaled_weights.sum()
    valid = (grid >= prior_min) & (grid <= prior_max)
    if not np.any(valid):
        raise ValueError("h0_grid has no points inside the prior bounds")

    predicted_distance = 299_792.458 * z[None, :] / grid[:, None]
    log_terms = norm.logpdf(
        observed_distance_mpc,
        loc=predicted_distance,
        scale=distance_sigma_mpc,
    )
    log_weights = np.full(weights.shape, -np.inf)
    positive = weights > 0
    log_weights[positive] = np.log(weights[positive])
    log_likelihood = logsumexp(log_terms + log_weights[None, :], axis=1)

    log_posterior = np.where(valid, log_likelihood, -np.inf)
    finite = np.isfinite(log_posterior)
    if not np.any(finite):
        raise ValueError("Posterior has no finite support on the supplied H0 grid")

    posterior = np.zeros_like(grid)
    posterior[finite] = np.exp(log_posterior[finite] - np.max(log_posterior[finite]))
    area = trapezoid(posterior, grid)
    if not np.isfinite(area) or area <= 0:
        raise ValueError("Posterior normalization failed")
    return posterior / area
