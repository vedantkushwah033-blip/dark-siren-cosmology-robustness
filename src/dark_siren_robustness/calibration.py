"""Metrics for repeated-sampling evaluation of interval estimates.

These summaries assume each row is one independent simulation replicate. They
do not perform simulation-based calibration (SBC), and they cannot repair a
misspecified generative model.
"""
from __future__ import annotations

import numpy as np


def summarize_interval_performance(
    true_values,
    estimates,
    lower_bounds,
    upper_bounds,
) -> dict[str, float]:
    """Summarize bias, RMSE, interval width, coverage, and MC coverage error.

    All inputs must be finite one-dimensional arrays of equal nonzero length.
    Bounds must satisfy lower_bounds <= upper_bounds. Coverage Monte Carlo
    standard error is estimated as sqrt(p_hat * (1-p_hat) / n).
    """
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
