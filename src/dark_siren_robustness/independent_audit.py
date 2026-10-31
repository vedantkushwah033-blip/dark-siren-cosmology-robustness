"""Independent synthetic inclusion-audit cohort generator.

This module generates an audit cohort from the declared parent redshift law,
with fresh retention draws and namespaced IDs. It is separate from the source
catalogue frame. Truth probabilities are returned separately from observable
audit records to reduce accidental inference leakage.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import expit

from dark_siren_robustness.generative_selection import latent_galaxy_intensity


def generate_independent_inclusion_audit(
    *,
    seed: int,
    n_audit: int,
    z_min: float,
    z_max: float,
    intercept: float,
    slope: float,
    quadratic: float = 0.0,
    grid_size: int = 20001,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Generate observable audit records and a separate truth-only score table.

    Redshifts are sampled independently from the normalized latent galaxy
    intensity on [z_min, z_max]. Retention labels are fresh Bernoulli draws
    under the declared logistic inclusion law. No source-frame IDs, labels,
    or retention realization are accepted as inputs.

    Returns (observable, scoring_truth). The observable frame contains only
    audit_object_id, redshift, and included. The scoring frame contains the
    true inclusion probability and must never be passed to model fitting.
    """
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, (int, np.integer)) or seed < 0:
        raise ValueError("seed must be a non-negative integer")
    if isinstance(n_audit, (bool, np.bool_)) or not isinstance(n_audit, (int, np.integer)) or n_audit < 10:
        raise ValueError("n_audit must be an integer >= 10")
    if not np.isfinite(z_min) or not np.isfinite(z_max) or not 0 < z_min < z_max:
        raise ValueError("redshift bounds must be finite and satisfy 0 < z_min < z_max")
    if not all(np.isfinite(v) for v in (intercept, slope, quadratic)):
        raise ValueError("inclusion coefficients must be finite")
    if isinstance(grid_size, (bool, np.bool_)) or not isinstance(grid_size, (int, np.integer)) or grid_size < 101:
        raise ValueError("grid_size must be an integer >= 101")

    rng = np.random.default_rng(int(seed))
    grid = np.linspace(float(z_min), float(z_max), int(grid_size))
    intensity = latent_galaxy_intensity(grid)
    increments = 0.5 * (intensity[:-1] + intensity[1:]) * np.diff(grid)
    cdf = np.concatenate(([0.0], np.cumsum(increments)))
    cdf /= cdf[-1]
    redshift = np.interp(rng.random(int(n_audit)), cdf, grid)

    midpoint = (float(z_min) + float(z_max)) / 2.0
    half_range = (float(z_max) - float(z_min)) / 2.0
    x = (redshift - midpoint) / half_range
    probability = expit(float(intercept) + float(slope) * x + float(quadratic) * (x**2 - 1.0 / 3.0))
    included = (rng.random(int(n_audit)) < probability).astype(np.int8)

    observable = pd.DataFrame({
        "audit_object_id": [f"independent-audit:{i}" for i in range(int(n_audit))],
        "redshift": redshift,
        "included": included,
    })
    scoring_truth = pd.DataFrame({
        "audit_object_id": observable["audit_object_id"].copy(),
        "true_inclusion_probability": probability,
    })
    return observable, scoring_truth
