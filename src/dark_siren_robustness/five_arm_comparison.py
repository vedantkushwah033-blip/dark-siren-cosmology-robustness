"""Matched five-arm inclusion-weight construction for synthetic comparisons.

This module constructs host-level probabilities and raw weights only. It does
not run the GW likelihood or claim that any arm is calibrated.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike
from scipy.special import expit

from .generative_selection import catalogue_inclusion_probability
from .inclusion_posterior import InclusionPosteriorQuadrature, fit_inclusion_posterior_quadrature


@dataclass(frozen=True)
class FiveArmWeights:
    """Raw arm-specific host weights and the fitted validation posterior."""

    parent_redshift: np.ndarray
    catalogue_redshift: np.ndarray
    arm_a_catalogue_only: np.ndarray
    arm_b_parent_oracle: np.ndarray
    arm_c_known_probability: np.ndarray
    arm_d_map_plugin: np.ndarray
    arm_e_node_weights: np.ndarray
    arm_e_quadrature_weights: np.ndarray
    arm_c_true_inclusion_probability: np.ndarray
    arm_d_map_inclusion_probability: np.ndarray
    arm_e_node_inclusion_probability: np.ndarray
    inclusion_posterior: InclusionPosteriorQuadrature


def _validate_redshifts(values: ArrayLike, name: str, *, minimum: int = 5) -> np.ndarray:
    z = np.asarray(values, dtype=float)
    if z.ndim != 1 or z.size < minimum or not np.all(np.isfinite(z)):
        raise ValueError(f"{name} must be a finite one-dimensional array with at least {minimum} entries")
    if np.any(z <= 0):
        raise ValueError(f"{name} must contain strictly positive redshifts")
    return z.copy()


def _validate_base_weights(
    values: ArrayLike | None, size: int, name: str
) -> np.ndarray:
    if values is None:
        return np.ones(size, dtype=float)
    weights = np.asarray(values, dtype=float)
    if (
        weights.shape != (size,)
        or not np.all(np.isfinite(weights))
        or np.any(weights <= 0)
    ):
        raise ValueError(
            f"{name} must be finite, strictly positive, and match its redshift array"
        )
    return weights.copy()


def _probabilities_for_nodes(
    redshift: np.ndarray,
    nodes: np.ndarray,
    *,
    degree: int,
    center: float,
    scale: float,
) -> np.ndarray:
    x = (redshift - center) / scale
    design = np.column_stack([x ** power for power in range(degree + 1)])
    probability = expit(design @ nodes.T).T
    if probability.ndim != 2 or probability.shape[1] != redshift.size:
        raise ArithmeticError("node-specific inclusion probabilities have invalid shape")
    if not np.all(np.isfinite(probability)) or np.any(probability <= 0) or np.any(probability > 1):
        raise ArithmeticError("node-specific inclusion probabilities are numerically unusable")
    return probability


def build_five_arm_weights(
    parent_redshift: ArrayLike,
    catalogue_redshift: ArrayLike,
    validation_redshift: ArrayLike,
    validation_retained: ArrayLike,
    *,
    catalogue_floor: float,
    catalogue_midpoint: float = 0.22,
    catalogue_width: float = 0.05,
    degree: int = 1,
    coordinate_center: float = 0.205,
    coordinate_scale: float = 0.195,
    prior_mean: ArrayLike | None = None,
    prior_sd: ArrayLike | None = None,
    quadrature_order: int = 15,
    parent_base_weights: ArrayLike | None = None,
    catalogue_base_weights: ArrayLike | None = None,
) -> FiveArmWeights:
    """Construct A-E raw weights from one matched synthetic replicate.

    A: declared base-population weights on the observed catalogue.
    B: the same base-population weighting rule on the full parent population.
    C: catalogue base weights divided by true inclusion probabilities.
    D: catalogue base weights divided by probabilities at the validation MAP.
    E: one catalogue-base-weighted IPW vector per posterior quadrature node.

    Supply parent_base_weights and catalogue_base_weights evaluated from the
    same declared population/rate function. If omitted, equal base weights are
    used in both populations. No clipping or normalization is performed. The
    caller must use the same fixed KDE and downstream selection-conditioned
    likelihood for every arm. E's coefficient nodes are integrated at the
    complete joint-likelihood level, not by averaging host weights first.
    """
    parent = _validate_redshifts(parent_redshift, "parent_redshift")
    catalogue = _validate_redshifts(catalogue_redshift, "catalogue_redshift")
    validation = _validate_redshifts(validation_redshift, "validation_redshift")
    parent_base = _validate_base_weights(
        parent_base_weights, parent.size, "parent_base_weights"
    )
    catalogue_base = _validate_base_weights(
        catalogue_base_weights, catalogue.size, "catalogue_base_weights"
    )
    retained = np.asarray(validation_retained)
    if retained.shape != validation.shape or np.any((retained != 0) & (retained != 1)):
        raise ValueError("validation_retained must be binary and match validation_redshift")
    if not np.isfinite(catalogue_floor) or not 0 <= catalogue_floor <= 1:
        raise ValueError("catalogue_floor must lie in [0, 1]")
    if not np.isfinite(catalogue_midpoint) or catalogue_midpoint <= 0:
        raise ValueError("catalogue_midpoint must be finite and positive")
    if not np.isfinite(catalogue_width) or catalogue_width <= 0:
        raise ValueError("catalogue_width must be finite and positive")
    if not np.isfinite(coordinate_center) or not np.isfinite(coordinate_scale) or coordinate_scale <= 0:
        raise ValueError("coefficient coordinates must have finite center and positive scale")

    posterior = fit_inclusion_posterior_quadrature(
        validation,
        retained,
        degree=degree,
        prior_mean=prior_mean,
        prior_sd=prior_sd,
        coordinate_center=coordinate_center,
        coordinate_scale=coordinate_scale,
        order=quadrature_order,
    )

    true_probability = catalogue_inclusion_probability(
        catalogue,
        floor=float(catalogue_floor),
        midpoint=float(catalogue_midpoint),
        width=float(catalogue_width),
    )
    map_probability = _probabilities_for_nodes(
        catalogue,
        posterior.mode[None, :],
        degree=posterior.degree,
        center=posterior.center,
        scale=posterior.scale,
    )[0]
    node_probability = _probabilities_for_nodes(
        catalogue,
        posterior.nodes,
        degree=posterior.degree,
        center=posterior.center,
        scale=posterior.scale,
    )

    with np.errstate(over="ignore", divide="ignore"):
        known_weights = catalogue_base / true_probability
        map_weights = catalogue_base / map_probability
        node_weights = catalogue_base[None, :] / node_probability
    if not all(np.all(np.isfinite(value)) for value in (known_weights, map_weights, node_weights)):
        raise ArithmeticError("inverse inclusion weights overflow; refusing silent clipping")

    return FiveArmWeights(
        parent_redshift=parent,
        catalogue_redshift=catalogue,
        arm_a_catalogue_only=catalogue_base,
        arm_b_parent_oracle=parent_base,
        arm_c_known_probability=known_weights,
        arm_d_map_plugin=map_weights,
        arm_e_node_weights=node_weights,
        arm_e_quadrature_weights=np.asarray(posterior.weights, dtype=float).copy(),
        arm_c_true_inclusion_probability=true_probability,
        arm_d_map_inclusion_probability=map_probability,
        arm_e_node_inclusion_probability=node_probability,
        inclusion_posterior=posterior,
    )


def weight_concentration_diagnostics(raw_weights: ArrayLike) -> dict[str, float]:
    """Return untrimmed Kish ESS and normalized upper-tail mass diagnostics."""
    weights = np.asarray(raw_weights, dtype=float)
    if weights.ndim != 1 or weights.size < 1 or not np.all(np.isfinite(weights)) or np.any(weights <= 0):
        raise ValueError("raw_weights must be a finite, strictly positive one-dimensional array")
    scaled = weights / np.max(weights)
    if np.any(scaled == 0.0):
        raise ArithmeticError("weight scaling underflowed positive mass; refusing to hide support loss")
    total = float(np.sum(scaled))
    normalized = scaled / total
    ess = float(total * total / np.sum(scaled * scaled))
    ordered = np.sort(normalized)[::-1]
    n1 = max(1, int(np.ceil(0.01 * weights.size)))
    n5 = max(1, int(np.ceil(0.05 * weights.size)))
    return {
        "kish_ess": ess,
        "kish_ess_fraction": ess / weights.size,
        "maximum_normalized_weight": float(ordered[0]),
        "top_1_percent_mass": float(np.sum(ordered[:n1])),
        "top_5_percent_mass": float(np.sum(ordered[:n5])),
    }
