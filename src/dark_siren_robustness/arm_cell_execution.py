"""Independent execution of one synthetic dark-siren inference arm.

This module deliberately evaluates exactly one (validation-size, arm) unit.
It does not construct the other arms and never combines pure IPW with a
missing-host mixture correction. The caller owns lifecycle persistence.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike
from scipy.special import expit

from .catalogue_density import reflected_weighted_kde_bandwidth
from .five_arm_comparison import weight_concentration_diagnostics
from .five_arm_kde_likelihood import arm_kde_selected_log_likelihood
from .h0_posterior import normalize_h0_posterior
from .inclusion_posterior import fit_inclusion_posterior_quadrature


ARM_NAMES = (
    "A_catalogue_only",
    "B_parent_oracle",
    "C_known_probability_ipw",
    "D_map_plugin_ipw",
    "E_uncertainty_integrated_ipw",
)


@dataclass(frozen=True)
class ArmCellResult:
    """Auditable numerical outputs from one arm/validation-size calculation."""

    arm: str
    status: str
    raw_weights: np.ndarray
    component_weights: np.ndarray
    mixture_weights: np.ndarray
    component_density: np.ndarray
    component_log_likelihood: np.ndarray
    marginal_log_likelihood: np.ndarray
    component_log_detection_efficiency: np.ndarray
    posterior_summary: Mapping[str, float]
    weight_diagnostics: Mapping[str, float]
    realized_bandwidths: np.ndarray
    inclusion_mode: np.ndarray | None
    inclusion_nodes: np.ndarray | None
    inclusion_quadrature_weights: np.ndarray | None


def _positive_vector(values: ArrayLike, name: str, size: int | None = None) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if (arr.ndim != 1 or (size is not None and arr.size != size)
            or arr.size < 1 or not np.all(np.isfinite(arr)) or np.any(arr <= 0)):
        raise ValueError(f"{name} must be a finite, strictly positive vector")
    return arr.copy()


def _posterior_quantiles(grid: np.ndarray, log_likelihood: np.ndarray,
                         prior_min: float, prior_max: float) -> dict[str, float]:
    posterior = normalize_h0_posterior(
        grid, log_likelihood, prior_min=prior_min, prior_max=prior_max
    )
    return {
        "posterior_median": float(posterior.quantile(0.50)),
        "q025": float(posterior.quantile(0.025)),
        "q05": float(posterior.quantile(0.05)),
        "q16": float(posterior.quantile(0.16)),
        "q84": float(posterior.quantile(0.84)),
        "q95": float(posterior.quantile(0.95)),
        "q975": float(posterior.quantile(0.975)),
    }


def evaluate_arm_cell(
    arm: str,
    *,
    parent_redshift: ArrayLike,
    catalogue_redshift: ArrayLike,
    validation_redshift: ArrayLike,
    validation_retained: ArrayLike,
    true_catalogue_inclusion_probability: ArrayLike,
    events: Sequence[Mapping[str, float]],
    H0_grid: ArrayLike,
    redshift_grid: ArrayLike,
    merger_weight: ArrayLike,
    h0_true: float,
    h0_prior_min: float,
    h0_prior_max: float,
    detection_threshold_mpc: float,
    scott_multiplier: float,
    quadrature_order: int = 5,
    inclusion_degree: int = 1,
    coordinate_center: float = 0.205,
    coordinate_scale: float = 0.195,
    quadratic_center: float = 0.0,
    parent_base_weights: ArrayLike | None = None,
    catalogue_base_weights: ArrayLike | None = None,
) -> ArmCellResult:
    """Evaluate one arm independently from shared immutable synthetic inputs.

    A/B/C do not fit or depend on the validation inclusion model. D/E each
    perform their own validation-posterior fit. For E, the complete selected
    multi-event likelihood is marginalized over coefficient-posterior nodes.
    """
    if arm not in ARM_NAMES:
        raise ValueError(f"unknown arm: {arm}")
    parent = _positive_vector(parent_redshift, "parent_redshift")
    catalogue = _positive_vector(catalogue_redshift, "catalogue_redshift")
    validation = _positive_vector(validation_redshift, "validation_redshift")
    retained = np.asarray(validation_retained)
    if (retained.shape != validation.shape or not np.all(np.isfinite(retained))
            or np.any((retained != 0) & (retained != 1))):
        raise ValueError("validation_retained must be binary and match validation_redshift")
    if inclusion_degree not in (1, 2):
        raise ValueError("inclusion_degree must be 1 or 2")
    if not np.isfinite(h0_true) or not h0_prior_min <= h0_true <= h0_prior_max:
        raise ValueError("h0_true must lie within finite prior support")
    if not np.isfinite(h0_prior_min) or not np.isfinite(h0_prior_max) or h0_prior_min >= h0_prior_max:
        raise ValueError("H0 prior support must be finite and increasing")

    parent_base = (np.ones(parent.size) if parent_base_weights is None
                   else _positive_vector(parent_base_weights, "parent_base_weights", parent.size))
    catalogue_base = (np.ones(catalogue.size) if catalogue_base_weights is None
                      else _positive_vector(catalogue_base_weights, "catalogue_base_weights", catalogue.size))

    mode = nodes = quadrature = None
    if arm == "A_catalogue_only":
        sample, components = catalogue, catalogue_base[None, :]
        mixture = np.ones(1)
    elif arm == "B_parent_oracle":
        sample, components = parent, parent_base[None, :]
        mixture = np.ones(1)
    elif arm == "C_known_probability_ipw":
        probability = np.asarray(true_catalogue_inclusion_probability, dtype=float)
        if (probability.shape != catalogue.shape or not np.all(np.isfinite(probability))
                or np.any(probability <= 0) or np.any(probability > 1)):
            raise ValueError("true inclusion probabilities must be finite and in (0, 1]")
        with np.errstate(over="ignore", divide="ignore"):
            weights = catalogue_base / probability
        if not np.all(np.isfinite(weights)):
            raise ArithmeticError("known-probability IPW weights overflow")
        sample, components, mixture = catalogue, weights[None, :], np.ones(1)
    else:
        if not np.isfinite(coordinate_center) or not np.isfinite(coordinate_scale) or coordinate_scale <= 0:
            raise ValueError("inclusion coordinate center/scale are invalid")
        posterior = fit_inclusion_posterior_quadrature(
            validation, retained.astype(np.int8), degree=inclusion_degree,
            coordinate_center=coordinate_center, coordinate_scale=coordinate_scale,
            quadratic_center=quadratic_center, order=quadrature_order,
        )
        mode = np.asarray(posterior.mode, dtype=float).copy()
        nodes = np.asarray(posterior.nodes, dtype=float).copy()
        quadrature = np.asarray(posterior.weights, dtype=float).copy()
        x = (catalogue - coordinate_center) / coordinate_scale
        features = [np.ones_like(x), x]
        if inclusion_degree == 2:
            features.append(x * x - quadratic_center)
        design = np.column_stack(features)
        if arm == "D_map_plugin_ipw":
            probability = expit(design @ mode)
            with np.errstate(over="ignore", divide="ignore"):
                weights = catalogue_base / probability
            if not np.all(np.isfinite(weights)) or np.any(weights <= 0):
                raise ArithmeticError("MAP plug-in IPW weights are numerically unusable")
            components, mixture = weights[None, :], np.ones(1)
        else:
            probabilities = expit(design @ nodes.T).T
            with np.errstate(over="ignore", divide="ignore"):
                components = catalogue_base[None, :] / probabilities
            if (not np.all(np.isfinite(components)) or np.any(components <= 0)
                    or quadrature.shape != (nodes.shape[0],)
                    or not np.isclose(np.sum(quadrature), 1.0)):
                raise ArithmeticError("posterior-integrated IPW components are numerically unusable")
            mixture = quadrature
        sample = catalogue

    result = arm_kde_selected_log_likelihood(
        H0_grid, sample, components if components.shape[0] > 1 else components[0],
        events, detection_threshold_mpc=detection_threshold_mpc,
        redshift_grid=redshift_grid, merger_weight=merger_weight,
        scott_multiplier=scott_multiplier,
        mixture_weights=(mixture if components.shape[0] > 1 else None),
    )
    summary = _posterior_quantiles(
        np.asarray(H0_grid, dtype=float), result.marginal_log_likelihood,
        h0_prior_min, h0_prior_max,
    )
    diagnostics = (weight_concentration_diagnostics(components[0])
                   if components.shape[0] == 1 else {
        "quadrature_ess": float(1.0 / np.sum(mixture ** 2)),
        "minimum_node_kish_ess": min(weight_concentration_diagnostics(w)["kish_ess"] for w in components),
        "maximum_node_top_5_percent_mass": max(weight_concentration_diagnostics(w)["top_5_percent_mass"] for w in components),
        "quadrature_nodes": int(components.shape[0]),
    })
    bandwidths = np.asarray([
        reflected_weighted_kde_bandwidth(
            sample, lower_bound=float(np.min(redshift_grid)),
            upper_bound=float(np.max(redshift_grid)), sample_weights=w,
            scott_multiplier=scott_multiplier,
        )
        for w in components
    ], dtype=float)
    return ArmCellResult(
        arm=arm, status="SUCCESS",
        raw_weights=components[0].copy(),
        component_weights=components.copy(),
        mixture_weights=mixture.copy(),
        component_density=result.component_density.copy(),
        component_log_likelihood=result.component_log_likelihood.copy(),
        marginal_log_likelihood=result.marginal_log_likelihood.copy(),
        component_log_detection_efficiency=result.component_log_detection_efficiency.copy(),
        posterior_summary=summary, weight_diagnostics=diagnostics,
        realized_bandwidths=bandwidths,
        inclusion_mode=mode, inclusion_nodes=nodes,
        inclusion_quadrature_weights=quadrature,
    )
