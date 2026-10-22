"""Catalogue-inclusion utilities for controlled simulation experiments."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def distance_dependent_completeness(
    distance_mpc: ArrayLike, scale_mpc: float, floor: float = 0.0
) -> NDArray[np.float64]:
    """Illustrative decreasing inclusion probability; not a real-survey model.

    p(include | d) = floor + (1 - floor) * exp(-d / scale).
    Parameters must be justified before this toy mechanism is used scientifically.
    """
    distance = np.asarray(distance_mpc, dtype=float)
    if np.any(~np.isfinite(distance)) or np.any(distance < 0):
        raise ValueError("distance_mpc must contain finite, non-negative values")
    if not np.isfinite(scale_mpc) or scale_mpc <= 0:
        raise ValueError("scale_mpc must be finite and positive")
    if not np.isfinite(floor) or not 0 <= floor <= 1:
        raise ValueError("floor must lie in [0, 1]")
    return floor + (1.0 - floor) * np.exp(-distance / scale_mpc)
