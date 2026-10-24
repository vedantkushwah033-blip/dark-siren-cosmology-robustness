"""Cosmology-aware toy H0 inference primitives for controlled model checks.

These functions support a simplified Gaussian distance likelihood only. They do
not implement a detector likelihood or a complete hierarchical dark-siren model.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import trapezoid
from scipy.special import logsumexp
from scipy.stats import norm

from .cosmology import luminosity_distance_flat_lcdm, low_redshift_luminosity_distance


def distance_grid_for_model(
    h0_grid: ArrayLike,
    redshifts: ArrayLike,
    model: str,
    *,
    omega_m: float = 0.3,
    integration_steps: int = 2_000,
) -> NDArray[np.float64]:
    """Return predicted distances with shape (len(H0 grid), len(hosts)), in Mpc."""
    h0 = np.asarray(h0_grid, dtype=float)
    z = np.asarray(redshifts, dtype=float)
    if h0.ndim != 1 or h0.size < 1 or np.any(~np.isfinite(h0)) or np.any(h0 <= 0) or np.any(np.diff(h0) <= 0):
        raise ValueError("h0_grid must be finite, positive, non-empty, one-dimensional, and increasing")
    if z.ndim != 1 or z.size == 0 or np.any(~np.isfinite(z)) or np.any(z < 0):
        raise ValueError("redshifts must be a non-empty finite non-negative vector")
    if model == "low_z":
        reference = low_redshift_luminosity_distance(z, 70.0)
    elif model == "flat_lcdm":
        reference = luminosity_distance_flat_lcdm(z, 70.0, omega_m=omega_m, integration_steps=integration_steps)
    else:
        raise ValueError("model must be 'low_z' or 'flat_lcdm'")
    # Both implemented relations scale exactly as 1/H0 for fixed Omega_m.
    return reference[None, :] * (70.0 / h0[:, None])


def posterior_h0_from_distance_grid(
    h0_grid: ArrayLike,
    predicted_distance_mpc: ArrayLike,
    observed_distance_mpc: float,
    distance_sigma_mpc: float,
    *,
    host_weights: ArrayLike | None = None,
    prior_min: float = 40.0,
    prior_max: float = 100.0,
) -> NDArray[np.float64]:
    """Normalize a weighted Gaussian-mixture H0 posterior from precomputed distances.

    The observed value is a Gaussian noisy distance estimator and may be
    non-positive; predicted means and the noise scale remain physical/positive.
    """
    grid = np.asarray(h0_grid, dtype=float)
    distances = np.asarray(predicted_distance_mpc, dtype=float)
    if grid.ndim != 1 or grid.size < 2 or np.any(~np.isfinite(grid)) or np.any(np.diff(grid) <= 0):
        raise ValueError("h0_grid must be finite, one-dimensional, and strictly increasing")
    if distances.ndim != 2 or distances.shape[0] != grid.size or distances.shape[1] == 0:
        raise ValueError("predicted_distance_mpc must have shape (len(h0_grid), n_hosts)")
    if np.any(~np.isfinite(distances)) or np.any(distances < 0):
        raise ValueError("predicted distances must be finite and non-negative")
    if not np.isfinite(observed_distance_mpc):
        raise ValueError("observed_distance_mpc must be finite")
    if not np.isfinite(distance_sigma_mpc) or distance_sigma_mpc <= 0:
        raise ValueError("distance_sigma_mpc must be finite and positive")
    if not np.isfinite(prior_min) or not np.isfinite(prior_max) or prior_min >= prior_max:
        raise ValueError("prior bounds must be finite and prior_min < prior_max")
    weights = np.ones(distances.shape[1]) if host_weights is None else np.asarray(host_weights, dtype=float)
    if weights.ndim != 1 or weights.size != distances.shape[1] or np.any(~np.isfinite(weights)) or np.any(weights < 0) or weights.sum() <= 0:
        raise ValueError("host_weights must be finite, non-negative, match hosts, and sum positive")
    weights = weights / weights.sum()
    positive = weights > 0
    log_weights = np.full(weights.shape, -np.inf)
    log_weights[positive] = np.log(weights[positive])
    log_likelihood = logsumexp(
        norm.logpdf(observed_distance_mpc, loc=distances, scale=distance_sigma_mpc)
        + log_weights[None, :],
        axis=1,
    )
    in_prior = (grid >= prior_min) & (grid <= prior_max)
    if np.count_nonzero(in_prior) < 2:
        raise ValueError("h0_grid must contain at least two points inside the prior")
    posterior = np.zeros_like(grid)
    finite = in_prior & np.isfinite(log_likelihood)
    if not np.any(finite):
        raise ValueError("posterior has no finite support")
    posterior[finite] = np.exp(log_likelihood[finite] - np.max(log_likelihood[finite]))
    area = float(trapezoid(posterior[in_prior], grid[in_prior]))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("posterior normalization failed")
    return posterior / area
