"""Design-aware resampling primitives for finite-parent validation studies.

These helpers do not assert that an FPC repairs nonlinear estimator variance.
They expose the exact sampling operations so runners and auditors can share a
small, explicitly tested contract.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def finite_population_correction(n: int, population_size: int) -> float:
    """Return sqrt(1 - n/N) for an SRSWOR sample of size n from N units."""
    for name, value in (("n", n), ("population_size", population_size)):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
            raise TypeError(f"{name} must be an integer")
    if population_size < 1:
        raise ValueError("population_size must be positive")
    if n < 1 or n > population_size:
        raise ValueError("n must lie in [1, population_size]")
    return float(np.sqrt(1.0 - n / population_size))


def paired_bootstrap_indices(
    rng: np.random.Generator, sample_size: int, *, population_size: int | None = None
) -> NDArray[np.int64]:
    """Draw record indices with replacement; apply the same indices to every field."""
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng must be a numpy.random.Generator")
    if isinstance(sample_size, (bool, np.bool_)) or not isinstance(sample_size, (int, np.integer)):
        raise TypeError("sample_size must be an integer")
    if population_size is None:
        population_size = int(sample_size)
    if isinstance(population_size, (bool, np.bool_)) or not isinstance(
        population_size, (int, np.integer)
    ):
        raise TypeError("population_size must be an integer")
    if sample_size < 1 or population_size < 1:
        raise ValueError("sample_size and population_size must be positive")
    return rng.integers(0, int(population_size), size=int(sample_size), dtype=np.int64)


def bootstrap_standard_error(values: np.ndarray) -> float:
    """Sample SD of finite bootstrap replicates; require at least two values."""
    array = np.asarray(values, dtype=float)
    if array.ndim != 1 or array.size < 2 or not np.all(np.isfinite(array)):
        raise ValueError("values must contain at least two finite one-dimensional replicates")
    return float(np.std(array, ddof=1))
