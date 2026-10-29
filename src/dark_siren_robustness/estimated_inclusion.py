"""Estimate catalogue inclusion probabilities from an external validation sample.

This module is intentionally separate from the existing oracle-weighted KDE.
The validation data must include reference redshifts and catalogue-retention
indicators sampled independently of retention, or under a documented sampling
design. The incomplete catalogue alone does not identify the inclusion law.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import minimize
from scipy.special import expit, log_expit


def fit_logistic_inclusion_model(
    validation_redshift: ArrayLike,
    included: ArrayLike,
    *,
    z_min: float,
    z_max: float,
) -> dict[str, float | str]:
    """Fit logit P(included=1 | z) = intercept + slope * standardized_z.

    Standardization is fixed by the declared domain: standardized_z is
    (z - (z_min + z_max)/2) / ((z_max - z_min)/2). It is not estimated from
    the validation sample, so fitted coefficients are comparable across
    replicates. The caller must ensure the validation sampling design supports
    this conditional model and must not treat fitted probabilities as known.

    Raises ValueError for malformed data or a validation sample containing only
    one outcome class. Raises RuntimeError if bounded maximum likelihood fails.
    """
    z = np.asarray(validation_redshift, dtype=float)
    y_raw = np.asarray(included)
    if z.ndim != 1 or z.size < 10 or not np.all(np.isfinite(z)):
        raise ValueError("validation_redshift must contain at least 10 finite values")
    if not np.isfinite(z_min) or not np.isfinite(z_max) or not z_min < z_max:
        raise ValueError("z_min and z_max must be finite with z_min < z_max")
    if np.any((z < z_min) | (z > z_max)):
        raise ValueError("validation_redshift must lie inside the declared domain")
    if y_raw.shape != z.shape or y_raw.ndim != 1:
        raise ValueError("included must be a one-dimensional array matching validation_redshift")
    if not np.all(np.isin(y_raw, [0, 1, False, True])):
        raise ValueError("included must contain only binary 0/1 values")
    y = y_raw.astype(float)
    if np.unique(y).size != 2:
        raise ValueError("validation sample must contain both included and not-included outcomes")

    midpoint = (float(z_min) + float(z_max)) / 2.0
    half_range = (float(z_max) - float(z_min)) / 2.0
    x = (z - midpoint) / half_range

    def negative_log_likelihood(theta: NDArray[np.float64]) -> float:
        eta = theta[0] + theta[1] * x
        value = -np.sum(y * log_expit(eta) + (1.0 - y) * log_expit(-eta))
        return float(value) if np.isfinite(value) else float("inf")

    result = minimize(
        negative_log_likelihood,
        x0=np.array([0.0, 0.0], dtype=float),
        method="L-BFGS-B",
        bounds=((-25.0, 25.0), (-50.0, 50.0)),
        options={"maxiter": 2000, "ftol": 1e-12, "gtol": 1e-8},
    )
    if not result.success or not np.all(np.isfinite(result.x)) or not np.isfinite(result.fun):
        raise RuntimeError(f"logistic inclusion fit failed: {result.message}")
    return {
        "model": "logistic_inclusion_v1",
        "intercept": float(result.x[0]),
        "slope": float(result.x[1]),
        "z_min": float(z_min),
        "z_max": float(z_max),
        "n_validation": int(z.size),
        "n_included": int(np.sum(y)),
        "n_not_included": int(y.size - np.sum(y)),
        "negative_log_likelihood": float(result.fun),
    }


def predict_inclusion_probability(
    redshift: ArrayLike,
    model: Mapping[str, object],
) -> NDArray[np.float64]:
    """Predict inclusion probabilities using a fitted model, without clipping.

    No implicit probability floor or weight cap is applied. If the caller
    chooses to stabilize weights, that rule must be explicit and reported as
    part of the estimator rather than silently embedded in this function.
    """
    if model.get("model") != "logistic_inclusion_v1":
        raise ValueError("model must be a logistic_inclusion_v1 result")
    try:
        intercept = float(model["intercept"])
        slope = float(model["slope"])
        z_min = float(model["z_min"])
        z_max = float(model["z_max"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("model is missing valid coefficients or domain bounds") from exc
    if not all(np.isfinite(value) for value in (intercept, slope, z_min, z_max)) or z_min >= z_max:
        raise ValueError("model coefficients and domain bounds must be finite and valid")
    z = np.asarray(redshift, dtype=float)
    if z.ndim != 1 or not np.all(np.isfinite(z)):
        raise ValueError("redshift must be a finite one-dimensional array")
    if np.any((z < z_min) | (z > z_max)):
        raise ValueError("redshift values must lie inside the fitted model domain")
    midpoint = (z_min + z_max) / 2.0
    half_range = (z_max - z_min) / 2.0
    return np.asarray(expit(intercept + slope * ((z - midpoint) / half_range)), dtype=float)



def estimate_host_density_from_validation(
    redshift_grid: ArrayLike,
    observed_redshift: ArrayLike,
    validation_redshift: ArrayLike,
    validation_included: ArrayLike,
    *,
    merger_weight: ArrayLike | None = None,
    bandwidth_factor: float | None = None,
) -> dict[str, object]:
    """Estimate host density using inclusion probabilities learned from audit data.

    The validation sample must represent the parent-frame population under a
    known sampling design and contain both catalogue-retention outcomes. The
    returned density uses only the observed catalogue and the audit records;
    the true inclusion probabilities and hidden parent objects are never inputs.

    merger_weight, when supplied, is a positive rate multiplier defined on
    redshift_grid (not an object-level weight). This wrapper intentionally
    applies no inverse-inclusion weight clipping. Its diagnostics make
    potentially unstable inverse weights visible; scientific use should
    preregister any positivity guard or cap rather than tuning it to H0 results.
    """
    from dark_siren_robustness.catalogue_density import estimate_host_density

    grid = np.asarray(redshift_grid, dtype=float)
    observed = np.asarray(observed_redshift, dtype=float)
    validation_z = np.asarray(validation_redshift, dtype=float)
    if grid.ndim != 1 or grid.size < 21 or not np.all(np.isfinite(grid)):
        raise ValueError("redshift_grid must be a finite 1D array with at least 21 points")
    if np.any(np.diff(grid) <= 0) or grid[0] <= 0:
        raise ValueError("redshift_grid must be strictly increasing and positive")
    if observed.ndim != 1 or observed.size < 5 or not np.all(np.isfinite(observed)):
        raise ValueError("observed_redshift must contain at least five finite values")
    if np.any((observed < grid[0]) | (observed > grid[-1])):
        raise ValueError("observed_redshift must lie inside redshift_grid bounds")
    model = fit_logistic_inclusion_model(
        validation_z,
        validation_included,
        z_min=float(grid[0]),
        z_max=float(grid[-1]),
    )
    probability = predict_inclusion_probability(observed, model)
    inverse_weights = 1.0 / probability
    scaled = inverse_weights / np.max(inverse_weights)
    effective_sample_size = float(np.sum(scaled) ** 2 / np.sum(scaled**2))
    density = estimate_host_density(
        grid,
        observed,
        inclusion_probability=probability,
        merger_weight=merger_weight,
        bandwidth_factor=bandwidth_factor,
    )
    return {
        "density": density,
        "model": model,
        "observed_inclusion_probability": probability,
        "inverse_weight_effective_sample_size": effective_sample_size,
        "inverse_weight_maximum": float(np.max(inverse_weights)),
        "inverse_weight_p95": float(np.quantile(inverse_weights, 0.95)),
        "inverse_weight_p99": float(np.quantile(inverse_weights, 0.99)),
        "n_observed_catalogue": int(observed.size),
        "n_validation": int(model["n_validation"]),
    }
