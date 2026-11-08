"""Out-of-sample inclusion scores and inverse-weight concentration diagnostics."""
from __future__ import annotations
import numpy as np
from numpy.typing import ArrayLike


def score_inclusion_predictions(retained: ArrayLike, probability: ArrayLike) -> dict[str, float]:
    """Return Brier score and binary log loss on held-out labels."""
    y = np.asarray(retained, dtype=float)
    p = np.asarray(probability, dtype=float)
    if y.ndim != 1 or p.ndim != 1 or y.shape != p.shape or y.size == 0:
        raise ValueError("retained and probability must be nonempty matching 1D arrays")
    if not np.all(np.isfinite(y)) or np.any((y != 0) & (y != 1)):
        raise ValueError("retained labels must be finite and binary")
    if not np.all(np.isfinite(p)) or np.any((p <= 0) | (p >= 1)):
        raise ValueError("probabilities must be finite and strictly between zero and one")
    return {
        "n": int(y.size),
        "brier_score": float(np.mean((p - y) ** 2)),
        "log_loss": float(-np.mean(y * np.log(p) + (1.0 - y) * np.log1p(-p))),
    }


def inverse_probability_weight_diagnostics(probability: ArrayLike) -> dict[str, float]:
    """Summarize normalized 1/p weights without modifying or stabilizing them."""
    p = np.asarray(probability, dtype=float)
    if p.ndim != 1 or p.size == 0 or not np.all(np.isfinite(p)):
        raise ValueError("probability must be a nonempty finite 1D array")
    if np.any((p <= 0) | (p > 1)):
        raise ValueError("inclusion probabilities must lie in (0, 1]")
    # Scale inverse weights before summation to avoid overflow when p is small.
    # q_i is unchanged because the common factor cancels during normalization.
    scaled_weights = float(np.min(p)) / p
    normalized = scaled_weights / np.sum(scaled_weights)
    positive = normalized[normalized > 0]
    entropy = -float(np.sum(positive * np.log(positive)))
    with np.errstate(over="ignore"):
        maximum_weight = float(1.0 / np.min(p))
    ordered = np.sort(normalized)[::-1]
    return {
        "n": int(p.size),
        "probability_min": float(np.min(p)),
        "probability_median": float(np.median(p)),
        "probability_max": float(np.max(p)),
        "weight_max": maximum_weight,
        "kish_ess": float(np.sum(scaled_weights) ** 2 / np.sum(scaled_weights**2)),
        "entropy_ess": float(np.exp(entropy)),
        "max_normalized_weight": float(np.max(normalized)),
        "top_1pct_weight_mass": float(np.sum(ordered[:max(1, int(np.ceil(0.01 * p.size)))])),
        "top_5pct_weight_mass": float(np.sum(ordered[:max(1, int(np.ceil(0.05 * p.size)))])),
    }


def score_probability_calibration(
    retained: ArrayLike,
    probability: ArrayLike,
    *,
    n_bins: int = 10,
) -> dict[str, float]:
    """Return equal-frequency reliability-bin gaps on held-out binary labels.

    Bin edges are quantiles of predicted probabilities. Tied predictions are
    kept in the same bin; empty bins are omitted. ECE is weighted by each
    nonempty bin's share of held-out records. These are descriptive diagnostics,
    not a formal test of calibration.
    """
    y = np.asarray(retained, dtype=float)
    p = np.asarray(probability, dtype=float)
    if y.ndim != 1 or p.ndim != 1 or y.shape != p.shape or y.size == 0:
        raise ValueError("retained and probability must be nonempty matching 1D arrays")
    if not np.all(np.isfinite(y)) or np.any((y != 0) & (y != 1)):
        raise ValueError("retained labels must be finite and binary")
    if not np.all(np.isfinite(p)) or np.any((p <= 0) | (p >= 1)):
        raise ValueError("probabilities must be finite and strictly between zero and one")
    if isinstance(n_bins, bool) or not isinstance(n_bins, int) or n_bins < 2:
        raise ValueError("n_bins must be an integer >= 2")

    edges = np.quantile(p, np.linspace(0.0, 1.0, n_bins + 1))
    assignments = np.searchsorted(edges[1:-1], p, side="right")
    gaps = []
    counts = []
    for index in range(n_bins):
        mask = assignments == index
        count = int(np.sum(mask))
        if count:
            gaps.append(abs(float(np.mean(y[mask]) - np.mean(p[mask]))))
            counts.append(count)
    ece = float(np.dot(np.asarray(counts, dtype=float) / y.size, gaps))
    return {
        "calibration_ece": ece,
        "calibration_max_gap": float(max(gaps)),
        "calibration_nonempty_bins": int(len(gaps)),
    }
