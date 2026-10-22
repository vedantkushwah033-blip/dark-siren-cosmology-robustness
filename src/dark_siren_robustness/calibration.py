"""Metrics and interval utilities for repeated-sampling evaluation."""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import cumulative_trapezoid, trapezoid


def equal_tail_interval(
    parameter_grid: ArrayLike, posterior_density: ArrayLike, level: float = 0.68
) -> tuple[float, float]:
    """Return an equal-tail credible interval from a grid density."""
    grid = np.asarray(parameter_grid, dtype=float)
    density = np.asarray(posterior_density, dtype=float)
    if grid.ndim != 1 or density.ndim != 1 or grid.size != density.size or grid.size < 2:
        raise ValueError("grid and density must be same-length 1D arrays of size >= 2")
    if np.any(~np.isfinite(grid)) or np.any(np.diff(grid) <= 0):
        raise ValueError("parameter_grid must be finite and strictly increasing")
    if np.any(~np.isfinite(density)) or np.any(density < 0):
        raise ValueError("posterior_density must be finite and non-negative")
    if not np.isfinite(level) or not 0 < level < 1:
        raise ValueError("level must lie strictly between 0 and 1")
    area = trapezoid(density, grid)
    if not np.isfinite(area) or area <= 0:
        raise ValueError("posterior density must have positive finite integral")
    cdf = np.concatenate(([0.0], cumulative_trapezoid(density / area, grid)))
    cdf /= cdf[-1]
    alpha = (1.0 - level) / 2.0
    return (float(np.interp(alpha, cdf, grid)), float(np.interp(1.0 - alpha, cdf, grid)))


def summarize_interval_performance(
    true_values, estimates, lower_bounds, upper_bounds
) -> dict[str, float]:
    """Summarize bias, RMSE, interval width, coverage, and MC coverage error."""
    arrays = [
        np.asarray(true_values, dtype=float),
        np.asarray(estimates, dtype=float),
        np.asarray(lower_bounds, dtype=float),
        np.asarray(upper_bounds, dtype=float),
    ]
    if any(a.ndim != 1 for a in arrays):
        raise ValueError("all inputs must be one-dimensional arrays")
    if any(a.size == 0 for a in arrays):
        raise ValueError("inputs must not be empty")
    if len({a.size for a in arrays}) != 1:
        raise ValueError("all inputs must have the same length")
    if any(np.any(~np.isfinite(a)) for a in arrays):
        raise ValueError("all inputs must contain only finite values")

    truth, estimate, lower, upper = arrays
    if np.any(lower > upper):
        raise ValueError("each lower interval bound must be <= its upper bound")
    errors = estimate - truth
    covered = (lower <= truth) & (truth <= upper)
    coverage = float(np.mean(covered))
    n = truth.size
    return {
        "n": float(n),
        "bias": float(np.mean(errors)),
        "rmse": float(np.sqrt(np.mean(errors**2))),
        "mean_interval_width": float(np.mean(upper - lower)),
        "empirical_coverage": coverage,
        "coverage_mc_se": float(np.sqrt(coverage * (1.0 - coverage) / n)),
    }
