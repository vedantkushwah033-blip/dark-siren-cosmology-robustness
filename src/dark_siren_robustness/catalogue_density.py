"""Bounded KDE host-density estimators for controlled selection experiments.

Inverse-inclusion weighting assumes supplied probabilities are known. It does
not estimate the selection function and is not a real-survey correction.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import trapezoid
from scipy.stats import gaussian_kde


def estimate_host_density(
    redshift_grid: ArrayLike,
    observed_redshift: ArrayLike,
    *,
    inclusion_probability: ArrayLike | None = None,
    merger_weight: ArrayLike | None = None,
    bandwidth_factor: float | None = None,
) -> NDArray[np.float64]:
    """Estimate a normalized host-redshift density on a bounded grid.

    For non-oracle analyses, observed_redshift must come from the observable
    catalogue, never hidden parent truth. Probabilities must be strictly
    positive: zero-probability regions cannot be recovered by inverse weights.
    Reflection reduces, but does not eliminate, KDE boundary bias. The final
    density is truncated and renormalized on the specified finite domain.
    """
    grid = np.asarray(redshift_grid, dtype=float)
    sample = np.asarray(observed_redshift, dtype=float)
    if grid.ndim != 1 or grid.size < 21 or not np.all(np.isfinite(grid)):
        raise ValueError("redshift_grid must be a finite 1D array with at least 21 points")
    if np.any(np.diff(grid) <= 0) or grid[0] <= 0:
        raise ValueError("redshift_grid must be strictly increasing and positive")
    if sample.ndim != 1 or sample.size < 5 or not np.all(np.isfinite(sample)):
        raise ValueError("observed_redshift must contain at least five finite values")
    if np.any((sample < grid[0]) | (sample > grid[-1])):
        raise ValueError("all sample redshifts must lie inside redshift_grid bounds")
    if np.unique(sample).size < 3:
        raise ValueError("observed_redshift must contain at least three distinct values")

    if inclusion_probability is None:
        weights = np.ones(sample.size, dtype=float)
    else:
        probability = np.asarray(inclusion_probability, dtype=float)
        if probability.shape != sample.shape or not np.all(np.isfinite(probability)):
            raise ValueError("inclusion_probability must be finite and match observed_redshift")
        if np.any((probability <= 0) | (probability > 1)):
            raise ValueError("inclusion_probability must lie in (0, 1]")
        weights = 1.0 / probability
    if bandwidth_factor is not None and (
        not np.isfinite(bandwidth_factor) or bandwidth_factor <= 0
    ):
        raise ValueError("bandwidth_factor must be finite and positive when supplied")

    left = 2.0 * grid[0] - sample
    right = 2.0 * grid[-1] - sample
    augmented_sample = np.concatenate((sample, left, right))
    augmented_weights = np.tile(weights / 3.0, 3)
    kde = gaussian_kde(
        augmented_sample, weights=augmented_weights, bw_method=bandwidth_factor
    )
    galaxy_density = np.asarray(kde(grid), dtype=float)

    if merger_weight is None:
        rate = np.ones_like(grid)
    else:
        rate = np.asarray(merger_weight, dtype=float)
        if rate.shape != grid.shape or not np.all(np.isfinite(rate)) or np.any(rate <= 0):
            raise ValueError("merger_weight must be finite, positive, and match redshift_grid")
    unnormalized = galaxy_density * rate
    area = float(trapezoid(unnormalized, grid))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("Estimated host density has invalid normalization")
    density = unnormalized / area
    if not np.all(np.isfinite(density)) or np.any(density < 0):
        raise ArithmeticError("Estimated host density contains invalid values")
    if not np.isclose(trapezoid(density, grid), 1.0, rtol=1e-10, atol=1e-10):
        raise ArithmeticError("Estimated host density failed its unit-integral check")
    return density
