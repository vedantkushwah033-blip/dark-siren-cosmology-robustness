"""Validated reference implementation of the simplified synthetic dark-siren H0 model.

This module is deliberately separate from the pinned legacy simulator so that the
baseline remains reproducible. It preserves the same low-redshift distance relation
and weighted Gaussian-mixture likelihood while rejecting malformed inputs explicitly.
"""
import numpy as np
from scipy.integrate import trapezoid
from scipy.stats import norm

C_KM_S = 299792.458


def _finite_scalar(value, name):
    array = np.asarray(value)
    if array.ndim != 0:
        raise ValueError(f"{name} must be a finite scalar.")
    result = float(array)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be a finite scalar.")
    return result


def _vector(value, name, *, nonnegative=False, strictly_increasing=False, min_size=1):
    array = np.asarray(value, dtype=float)
    if array.ndim != 1 or array.size < min_size:
        raise ValueError(f"{name} must be a one-dimensional array with at least {min_size} value(s).")
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must contain only finite values.")
    if nonnegative and np.any(array < 0):
        raise ValueError(f"{name} must be non-negative.")
    if strictly_increasing and np.any(np.diff(array) <= 0):
        raise ValueError(f"{name} must be strictly increasing.")
    return array


def posterior_h0(
    H0_grid, redshifts, gw_distance, gw_sigma, weights=None,
    missing_weight=0.0, missing_redshifts=None, missing_weights=None,
    prior_min=50.0, prior_max=90.0,
):
    """Compute normalized H0 posterior for the synthetic low-z mixture model."""
    grid = _vector(H0_grid, "H0_grid", strictly_increasing=True, min_size=2)
    if np.any(grid <= 0):
        raise ValueError("H0_grid must be positive.")
    z = _vector(redshifts, "redshifts", nonnegative=True)
    if z.size == 0:
        raise ValueError("At least one catalogued redshift is required.")
    distance = _finite_scalar(gw_distance, "gw_distance")
    sigma = _finite_scalar(gw_sigma, "gw_sigma")
    if distance <= 0:
        raise ValueError("gw_distance must be positive.")
    if sigma <= 0:
        raise ValueError("gw_sigma must be positive.")
    prior_lo = _finite_scalar(prior_min, "prior_min")
    prior_hi = _finite_scalar(prior_max, "prior_max")
    if prior_lo >= prior_hi:
        raise ValueError("prior_min must be less than prior_max.")

    catalog_weights = np.ones(z.size) if weights is None else _vector(weights, "weights", nonnegative=True)
    if catalog_weights.size != z.size or catalog_weights.sum() <= 0:
        raise ValueError("Invalid catalogue weights.")
    miss_fraction = _finite_scalar(missing_weight, "missing_weight")
    if not 0.0 <= miss_fraction < 1.0:
        raise ValueError("missing_weight must be in [0, 1).")
    catalog_weights = catalog_weights / catalog_weights.sum() * (1.0 - miss_fraction)

    if miss_fraction > 0:
        if missing_redshifts is None:
            raise ValueError("missing_redshifts are required when missing_weight > 0.")
        mz = _vector(missing_redshifts, "missing_redshifts", nonnegative=True)
        mw = np.ones(mz.size) if missing_weights is None else _vector(missing_weights, "missing_weights", nonnegative=True)
        if mw.size != mz.size or mw.sum() <= 0:
            raise ValueError("Invalid missing-host support or weights.")
        mw = mw / mw.sum()
    else:
        mz = np.array([], dtype=float)
        mw = np.array([], dtype=float)

    distance_grid = C_KM_S * z[None, :] / grid[:, None]
    likelihood = np.sum(catalog_weights[None, :] * norm.pdf(distance, loc=distance_grid, scale=sigma), axis=1)
    if miss_fraction > 0:
        missing_distance_grid = C_KM_S * mz[None, :] / grid[:, None]
        likelihood += miss_fraction * np.sum(
            mw[None, :] * norm.pdf(distance, loc=missing_distance_grid, scale=sigma), axis=1
        )

    posterior = np.where((grid >= prior_lo) & (grid <= prior_hi), likelihood, 0.0)
    area = float(trapezoid(posterior, grid))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("Posterior normalization failed.")
    return posterior / area


def posterior_summary(H0_grid, posterior):
    """Summarize a finite, non-negative posterior density on an increasing H0 grid."""
    grid = _vector(H0_grid, "H0_grid", strictly_increasing=True, min_size=2)
    density = _vector(posterior, "posterior", nonnegative=True, min_size=2)
    if density.size != grid.size:
        raise ValueError("posterior must have the same length as H0_grid.")
    area = float(np.trapezoid(density, grid))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("posterior must have a finite, positive integrated area.")
    density = density / area
    mean = float(np.trapezoid(grid * density, grid))
    variance = float(np.trapezoid((grid - mean) ** 2 * density, grid))
    cdf = np.concatenate(([0.0], np.cumsum((density[1:] + density[:-1]) * np.diff(grid) / 2.0)))
    cdf /= cdf[-1]
    return {
        "mean": mean,
        "median": float(np.interp(0.5, cdf, grid)),
        "posterior_sd": float(np.sqrt(max(variance, 0.0))),
        "lower_68": float(np.interp(0.16, cdf, grid)),
        "upper_68": float(np.interp(0.84, cdf, grid)),
        "lower_95": float(np.interp(0.025, cdf, grid)),
        "upper_95": float(np.interp(0.975, cdf, grid)),
    }
