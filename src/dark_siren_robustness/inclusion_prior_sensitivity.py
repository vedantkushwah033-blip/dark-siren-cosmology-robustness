"""Prior-scale sensitivity of the modular inclusion-posterior H0 analysis.

This utility varies the Gaussian prior scales on non-intercept inclusion
coefficients while holding the observable validation sample, catalogue,
events, H0 grid, selection rule and population base weights fixed.
"""
from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike

from .h0_posterior import H0GridPosterior, normalize_h0_posterior
from .inclusion_posterior import (
    InclusionPosteriorQuadrature,
    fit_inclusion_posterior_quadrature,
)
from .shared_inclusion_likelihood import (
    SharedInclusionLikelihoodResult,
    shared_inclusion_selected_log_likelihood,
)


@dataclass(frozen=True)
class PriorScaleSensitivityRecord:
    """One slope-prior multiplier and its downstream H0 result."""

    slope_multiplier: float
    prior_mean: np.ndarray
    prior_sd: np.ndarray
    inclusion_posterior: InclusionPosteriorQuadrature
    selected_likelihood: SharedInclusionLikelihoodResult
    h0_posterior: H0GridPosterior
    h0_quantiles: dict[str, float]
    minimum_catalogue_kish_ess: float
    median_catalogue_kish_ess: float
    maximum_catalogue_normalized_weight: float


@dataclass(frozen=True)
class PriorScaleSensitivityResult:
    """Prior-scale records and comparisons against the multiplier-one fit."""

    records: tuple[PriorScaleSensitivityRecord, ...]
    reference_multiplier: float
    max_abs_density_difference_from_reference: tuple[float, ...]
    quantile_shift_from_reference: tuple[dict[str, float], ...]


def run_inclusion_prior_scale_sensitivity(
    validation_redshift: ArrayLike,
    retained: ArrayLike,
    H0_grid: ArrayLike,
    catalogue_redshift: ArrayLike,
    events: Sequence[Mapping[str, float]],
    *,
    detection_threshold_mpc: float,
    population_base_weights: ArrayLike,
    degree: int = 1,
    prior_mean: ArrayLike | None = None,
    prior_sd: ArrayLike | None = None,
    slope_multipliers: tuple[float, ...] = (0.5, 1.0, 2.0),
    order: int = 11,
    prior_min: float | None = None,
    prior_max: float | None = None,
) -> PriorScaleSensitivityResult:
    """Refit the validation posterior and H0 posterior over slope-prior scales.

    The intercept prior scale is held fixed. Every non-intercept coefficient
    standard deviation is multiplied by each requested positive multiplier.
    All other data and model inputs are held fixed. The reference multiplier
    must be exactly 1.0 so reported shifts have an unambiguous baseline.
    """
    multipliers = tuple(float(value) for value in slope_multipliers)
    if (
        len(multipliers) < 2
        or any(not np.isfinite(value) or value <= 0 for value in multipliers)
        or len(set(multipliers)) != len(multipliers)
        or 1.0 not in multipliers
    ):
        raise ValueError(
            "slope_multipliers must be distinct positive values and include 1.0"
        )
    dimension = degree + 1
    mean = np.zeros(dimension) if prior_mean is None else np.asarray(prior_mean, dtype=float)
    sd = (
        np.asarray([5.0] + [2.5] * (dimension - 1), dtype=float)
        if prior_sd is None else np.asarray(prior_sd, dtype=float)
    )
    if mean.shape != (dimension,) or sd.shape != (dimension,):
        raise ValueError("prior vectors must match degree + 1")
    if not np.all(np.isfinite(mean)) or not np.all(np.isfinite(sd)) or np.any(sd <= 0):
        raise ValueError("prior parameters must be finite with positive scales")

    h0 = np.asarray(H0_grid, dtype=float)
    lo = float(h0[0] if prior_min is None else prior_min)
    hi = float(h0[-1] if prior_max is None else prior_max)
    records = []
    for multiplier in multipliers:
        scaled_sd = sd.copy()
        scaled_sd[1:] *= multiplier
        fit = fit_inclusion_posterior_quadrature(
            validation_redshift,
            retained,
            degree=degree,
            prior_mean=mean,
            prior_sd=scaled_sd,
            order=order,
        )
        likelihood = shared_inclusion_selected_log_likelihood(
            h0,
            catalogue_redshift,
            events,
            fit,
            detection_threshold_mpc=detection_threshold_mpc,
            population_base_weights=population_base_weights,
        )
        posterior = normalize_h0_posterior(
            h0,
            likelihood.marginal_log_likelihood,
            prior_min=lo,
            prior_max=hi,
        )
        quantiles = {
            "q05": posterior.quantile(0.05),
            "q16": posterior.quantile(0.16),
            "q50": posterior.quantile(0.50),
            "q84": posterior.quantile(0.84),
            "q95": posterior.quantile(0.95),
        }
        records.append(PriorScaleSensitivityRecord(
            slope_multiplier=multiplier,
            prior_mean=mean.copy(),
            prior_sd=scaled_sd,
            inclusion_posterior=fit,
            selected_likelihood=likelihood,
            h0_posterior=posterior,
            h0_quantiles=quantiles,
            minimum_catalogue_kish_ess=float(np.min(likelihood.node_kish_ess)),
            median_catalogue_kish_ess=float(np.median(likelihood.node_kish_ess)),
            maximum_catalogue_normalized_weight=float(
                np.max(likelihood.node_maximum_normalized_weight)
            ),
        ))

    reference = next(record for record in records if record.slope_multiplier == 1.0)
    density_differences = []
    quantile_shifts = []
    for record in records:
        reference_density = np.interp(
            record.h0_posterior.grid,
            reference.h0_posterior.grid,
            reference.h0_posterior.density,
        )
        density_differences.append(float(np.max(np.abs(
            record.h0_posterior.density - reference_density
        ))))
        quantile_shifts.append({
            key: float(record.h0_quantiles[key] - reference.h0_quantiles[key])
            for key in reference.h0_quantiles
        })
    return PriorScaleSensitivityResult(
        records=tuple(records),
        reference_multiplier=1.0,
        max_abs_density_difference_from_reference=tuple(density_differences),
        quantile_shift_from_reference=tuple(quantile_shifts),
    )
