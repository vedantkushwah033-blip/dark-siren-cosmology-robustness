"""Observable-only diagnostics computed before cosmological posterior construction.

Features describe validation support and fitted inverse-inclusion weights.
No simulation truth or posterior outcomes are accepted by this API.
"""
from __future__ import annotations

from collections.abc import Mapping
import numpy as np
from numpy.typing import ArrayLike


def build_preposterior_diagnostics(model: Mapping[str, object], validation_redshift: ArrayLike, validation_included: ArrayLike, observed_redshift: ArrayLike, observed_inclusion_probability: ArrayLike) -> dict[str, float | int | bool | str]:
    """Build observable-only scalar features before cosmological inference.

    The support-gap feature is the fraction of catalogue redshifts outside the
    min/max range represented in the external validation sample. It is a
    warning proxy, not a formal test of covariate shift.
    """
    if model.get("model") != "logistic_inclusion_v1":
        raise ValueError("model must be a logistic_inclusion_v1 fitted result")
    validation_z = np.asarray(validation_redshift, dtype=float)
    included_raw = np.asarray(validation_included)
    catalogue_z = np.asarray(observed_redshift, dtype=float)
    probability = np.asarray(observed_inclusion_probability, dtype=float)
    if validation_z.ndim != 1 or validation_z.size < 10 or not np.all(np.isfinite(validation_z)):
        raise ValueError("validation_redshift must contain at least 10 finite values")
    if included_raw.shape != validation_z.shape or included_raw.ndim != 1:
        raise ValueError("validation_included must match validation_redshift")
    if not np.all(np.isin(included_raw, [0, 1, False, True])):
        raise ValueError("validation_included must contain only binary values")
    if catalogue_z.ndim != 1 or catalogue_z.size < 1 or not np.all(np.isfinite(catalogue_z)):
        raise ValueError("observed_redshift must contain finite catalogue values")
    if probability.shape != catalogue_z.shape or not np.all(np.isfinite(probability)):
        raise ValueError("observed inclusion probabilities must match catalogue redshifts")
    if np.any((probability <= 0.0) | (probability >= 1.0)):
        raise ValueError("observed inclusion probabilities must lie strictly inside (0, 1)")
    try:
        intercept, slope = float(model["intercept"]), float(model["slope"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("fitted model coefficients must be present and finite") from exc
    if not np.all(np.isfinite([intercept, slope])):
        raise ValueError("fitted model coefficients must be present and finite")

    weights = 1.0 / probability
    scaled = weights / np.max(weights)
    ess = float(np.sum(scaled) ** 2 / np.sum(scaled**2))
    validation_y = included_raw.astype(int)
    validation_min, validation_max = float(np.min(validation_z)), float(np.max(validation_z))
    support_gap = float(np.mean((catalogue_z < validation_min) | (catalogue_z > validation_max)))
    return {
        "diagnostic_schema": "preposterior_weight_stability_v1",
        "n_validation": int(validation_z.size),
        "n_validation_included": int(np.sum(validation_y)),
        "validation_inclusion_fraction": float(np.mean(validation_y)),
        "validation_redshift_min": validation_min,
        "validation_redshift_max": validation_max,
        "validation_redshift_span": float(validation_max - validation_min),
        "catalogue_outside_validation_range_fraction": support_gap,
        "n_catalogue": int(catalogue_z.size),
        "fit_status": str(model.get("fit_status", "unknown")),
        "fit_boundary_limited": bool(model.get("coefficient_bound_hit", False)),
        "fitted_intercept": intercept,
        "fitted_slope": slope,
        "fitted_probability_min": float(np.min(probability)),
        "fitted_probability_p01": float(np.quantile(probability, 0.01)),
        "fitted_probability_p05": float(np.quantile(probability, 0.05)),
        "fitted_probability_median": float(np.median(probability)),
        "fitted_probability_max": float(np.max(probability)),
        "inverse_weight_ess": ess,
        "inverse_weight_ess_fraction": float(ess / catalogue_z.size),
        "inverse_weight_max": float(np.max(weights)),
        "inverse_weight_p95": float(np.quantile(weights, 0.95)),
        "inverse_weight_p99": float(np.quantile(weights, 0.99)),
    }
