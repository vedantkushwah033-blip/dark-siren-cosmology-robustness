"""Normalized H0 posterior and grid-refinement diagnostics for toy likelihoods.

These utilities normalize a supplied log likelihood under a declared uniform
H0 prior. They do not add selection corrections or validate the likelihood.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import cumulative_trapezoid, trapezoid


@dataclass(frozen=True)
class H0GridPosterior:
    """Normalized posterior density on a grid including both prior boundaries."""

    grid: np.ndarray
    density: np.ndarray
    cdf: np.ndarray
    prior_min: float
    prior_max: float
    log_likelihood_integral: float

    def quantile(self, probability: float) -> float:
        """Return a posterior quantile by monotone interpolation of the grid CDF."""
        q = float(probability)
        if not np.isfinite(q) or not 0.0 <= q <= 1.0:
            raise ValueError("probability must lie in [0, 1]")
        return float(np.interp(q, self.cdf, self.grid))


@dataclass(frozen=True)
class H0GridConvergence:
    """Density and quantile differences between coarse and fine H0 grids."""

    max_abs_density_difference: float
    max_abs_quantile_difference: float
    quantile_probabilities: tuple[float, ...]
    passed: bool


def normalize_h0_posterior(
    H0_grid: ArrayLike,
    log_likelihood: ArrayLike,
    *,
    prior_min: float = 50.0,
    prior_max: float = 90.0,
) -> H0GridPosterior:
    """Normalize a log likelihood under a uniform prior on [prior_min, prior_max].

    Piecewise-linear interpolation is used in log likelihood at prior boundaries
    when they do not coincide with grid points. The returned grid includes those
    boundaries, so trapezoidal normalization covers the full declared prior.
    Non-finite log likelihoods are rejected rather than silently imputed.
    """
    grid = np.asarray(H0_grid, dtype=float)
    log_like = np.asarray(log_likelihood, dtype=float)
    lo, hi = float(prior_min), float(prior_max)
    if grid.ndim != 1 or grid.size < 2 or not np.all(np.isfinite(grid)):
        raise ValueError("H0_grid must be a finite one-dimensional array with >= 2 points")
    if np.any(grid <= 0) or np.any(np.diff(grid) <= 0):
        raise ValueError("H0_grid must be positive and strictly increasing")
    if log_like.shape != grid.shape or not np.all(np.isfinite(log_like)):
        raise ValueError("log_likelihood must be finite and match H0_grid")
    if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
        raise ValueError("prior bounds must be finite with prior_min < prior_max")
    if grid[0] > lo or grid[-1] < hi:
        raise ValueError("H0_grid must span the complete uniform-prior interval")

    integration_grid = np.unique(np.concatenate(([lo], grid[(grid > lo) & (grid < hi)], [hi])))
    integration_log_like = np.interp(integration_grid, grid, log_like)
    shift = float(np.max(integration_log_like))
    unnormalized = np.exp(integration_log_like - shift)
    area = float(trapezoid(unnormalized, integration_grid))
    if not np.isfinite(area) or area <= 0.0:
        raise ArithmeticError("H0 posterior normalization failed")
    density = unnormalized / area
    cdf = np.concatenate(([0.0], cumulative_trapezoid(density, integration_grid)))
    cdf /= cdf[-1]
    log_integral = shift + float(np.log(area))
    return H0GridPosterior(
        grid=integration_grid,
        density=density,
        cdf=cdf,
        prior_min=lo,
        prior_max=hi,
        log_likelihood_integral=log_integral,
    )


def compare_h0_grid_posteriors(
    coarse: H0GridPosterior,
    fine: H0GridPosterior,
    *,
    density_tolerance: float = 1e-3,
    quantile_tolerance: float = 0.05,
    quantile_probabilities: tuple[float, ...] = (0.05, 0.16, 0.50, 0.84, 0.95),
) -> H0GridConvergence:
    """Compare normalized posterior shapes on the fine grid and selected quantiles."""
    if (coarse.prior_min, coarse.prior_max) != (fine.prior_min, fine.prior_max):
        raise ValueError("coarse and fine posteriors must use identical prior bounds")
    if not np.isfinite(density_tolerance) or density_tolerance <= 0:
        raise ValueError("density_tolerance must be positive and finite")
    if not np.isfinite(quantile_tolerance) or quantile_tolerance <= 0:
        raise ValueError("quantile_tolerance must be positive and finite")
    probabilities = tuple(float(q) for q in quantile_probabilities)
    if not probabilities or any(not np.isfinite(q) or not 0.0 <= q <= 1.0 for q in probabilities):
        raise ValueError("quantile_probabilities must be a non-empty sequence in [0, 1]")

    coarse_on_fine = np.interp(fine.grid, coarse.grid, coarse.density)
    density_difference = float(np.max(np.abs(coarse_on_fine - fine.density)))
    quantile_difference = float(max(
        abs(coarse.quantile(q) - fine.quantile(q)) for q in probabilities
    ))
    return H0GridConvergence(
        max_abs_density_difference=density_difference,
        max_abs_quantile_difference=quantile_difference,
        quantile_probabilities=probabilities,
        passed=(
            density_difference <= density_tolerance
            and quantile_difference <= quantile_tolerance
        ),
    )
