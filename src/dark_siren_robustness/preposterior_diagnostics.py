"""Observable-only diagnostics computed before cosmological posterior construction.

Features describe validation support and fitted inverse-inclusion weights.
No simulation truth or posterior outcomes are accepted by this API.
"""
from __future__ import annotations

from collections.abc import Mapping
import numpy as np
from scipy.special import expit
from numpy.typing import ArrayLike


def build_preposterior_diagnostics(model: Mapping[str, object], validation_redshift: ArrayLike, validation_included: ArrayLike, observed_redshift: ArrayLike, observed_inclusion_probability: ArrayLike) -> dict[str, float | int | bool | str]:
    """Build observable-only scalar features before cosmological inference.

    The support-gap feature is the fraction of catalogue redshifts outside the
    min/max range represented in the external validation sample. It is a
    warning proxy, not a formal test of covariate shift.
    """
    model_name = model.get("model")
    if model_name not in {"logistic_inclusion_v1", "quadratic_logistic_inclusion_v1"}:
        raise ValueError("model must be a supported logistic inclusion fitted result")
    degree = 1 if model_name == "logistic_inclusion_v1" else 2
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
        quadratic = float(model.get("quadratic", 0.0))
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("fitted model coefficients must be present and finite") from exc
    if not np.all(np.isfinite([intercept, slope, quadratic])):
        raise ValueError("fitted model coefficients must be present and finite")
    if abs(intercept) > 25.0 or abs(slope) > 50.0 or abs(quadratic) > 50.0:
        raise ValueError("fitted model coefficients exceed the bounded logistic-fit domain")
    try:
        z_min, z_max = float(model["z_min"]), float(model["z_max"])
        raw_model_n = model["n_validation"]
        raw_model_included = model["n_included"]
        raw_model_excluded = model["n_not_included"]
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("fitted model is missing validation counts or domain bounds") from exc
    count_values = (raw_model_n, raw_model_included, raw_model_excluded)
    if any(
        isinstance(value, (bool, np.bool_))
        or not isinstance(value, (int, np.integer))
        for value in count_values
    ):
        raise ValueError("fitted model validation counts must be integer values")
    model_n, model_included, model_excluded = (int(value) for value in count_values)
    if not np.isfinite(z_min) or not np.isfinite(z_max) or z_min >= z_max:
        raise ValueError("fitted model domain bounds must be finite and increasing")
    if np.any((validation_z < z_min) | (validation_z > z_max)) or np.any((catalogue_z < z_min) | (catalogue_z > z_max)):
        raise ValueError("validation and catalogue redshifts must lie inside fitted model domain")
    observed_included = int(np.sum(included_raw.astype(int)))
    if (
        model_n != validation_z.size
        or model_included != observed_included
        or model_excluded != validation_z.size - observed_included
        or model_n < 10
    ):
        raise ValueError("fitted model validation counts disagree with supplied validation data")
    if model.get("fit_status") not in {"converged", "boundary_limited"}:
        raise ValueError("fitted model has an unrecognized fit status")
    if bool(model.get("coefficient_bound_hit", False)) != (model.get("fit_status") == "boundary_limited"):
        raise ValueError("fitted model boundary flag disagrees with fit status")

    x_catalogue = (catalogue_z - (z_min + z_max) / 2.0) / ((z_max - z_min) / 2.0)
    eta = intercept + slope * x_catalogue
    if degree == 2:
        eta = eta + quadratic * (x_catalogue**2 - 1.0 / 3.0)
    expected_probability = expit(eta)
    if not np.allclose(probability, expected_probability, rtol=1e-10, atol=1e-12):
        raise ValueError("supplied probabilities disagree with fitted model")
    weights = 1.0 / probability
    if not np.all(np.isfinite(weights)) or np.any(weights <= 0.0):
        raise ValueError("inverse inclusion weights must be finite and strictly positive")
    weight_sum = float(np.sum(weights))
    weight_sum_squares = float(np.sum(weights**2))
    if not np.isfinite(weight_sum) or not np.isfinite(weight_sum_squares):
        raise ValueError("inverse inclusion weight aggregates overflowed")
    normalized_weights = weights / weight_sum
    scaled = weights / np.max(weights)
    ess = float(np.sum(scaled) ** 2 / np.sum(scaled**2))
    # Tail mass and entropy expose concentration beyond the maximum weight.
    descending_mass = np.sort(normalized_weights)[::-1]
    top_1pct_n = max(1, int(np.ceil(0.01 * descending_mass.size)))
    top_5pct_n = max(1, int(np.ceil(0.05 * descending_mass.size)))
    top_1pct_mass = float(np.sum(descending_mass[:top_1pct_n]))
    top_5pct_mass = float(np.sum(descending_mass[:top_5pct_n]))
    entropy = float(-np.sum(normalized_weights * np.log(normalized_weights)))
    entropy_ess_fraction = float(np.exp(entropy) / normalized_weights.size)
    validation_y = included_raw.astype(int)
    validation_min, validation_max = float(np.min(validation_z)), float(np.max(validation_z))
    validation_span = validation_max - validation_min
    support_gap = float(np.mean((catalogue_z < validation_min) | (catalogue_z > validation_max)))
    near_margin = 0.05 * validation_span
    support_edge = float(np.mean(
        (catalogue_z >= validation_min) & (catalogue_z <= validation_max)
        & ((catalogue_z - validation_min <= near_margin)
           | (validation_max - catalogue_z <= near_margin))
    ))
    return {
        "diagnostic_schema": "preposterior_weight_stability_v1",
        "n_validation": int(validation_z.size),
        "n_validation_included": int(np.sum(validation_y)),
        "validation_inclusion_fraction": float(np.mean(validation_y)),
        "validation_redshift_min": validation_min,
        "validation_redshift_max": validation_max,
        "validation_redshift_span": float(validation_max - validation_min),
        "catalogue_outside_validation_range_fraction": support_gap,
        "catalogue_near_validation_boundary_fraction": support_edge,
        "n_catalogue": int(catalogue_z.size),
        "fit_status": str(model.get("fit_status", "unknown")),
        "fit_boundary_limited": bool(model.get("coefficient_bound_hit", False)),
        "fitted_intercept": intercept,
        "fitted_slope": slope,
        "fitted_quadratic": quadratic,
        "fitted_model_degree": degree,
        "fitted_probability_min": float(np.min(probability)),
        "fitted_probability_p01": float(np.quantile(probability, 0.01)),
        "fitted_probability_p05": float(np.quantile(probability, 0.05)),
        "fitted_probability_median": float(np.median(probability)),
        "fitted_probability_max": float(np.max(probability)),
        "inverse_weight_ess": ess,
        "inverse_weight_sum": weight_sum,
        "inverse_weight_sum_squares": weight_sum_squares,
        "inverse_weight_ess_fraction": float(ess / catalogue_z.size),
        "inverse_weight_entropy_ess_fraction": entropy_ess_fraction,
        "inverse_weight_max_share": float(np.max(normalized_weights)),
        "inverse_weight_top_1pct_mass": top_1pct_mass,
        "inverse_weight_top_5pct_mass": top_5pct_mass,
        "inverse_weight_max": float(np.max(weights)),
        "inverse_weight_p95": float(np.quantile(weights, 0.95)),
        "inverse_weight_p99": float(np.quantile(weights, 0.99)),
    }
