"""Matched five-arm selected H0 likelihoods for the synthetic Euclidean model.

The same detected events and threshold are used for every arm. E is
marginalized over one shared inclusion-coefficient node for the entire
multi-event likelihood.
"""
from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike
from scipy.special import log_ndtr, logsumexp
from scipy.stats import norm

from .dataset_likelihood import C_KM_S
from .five_arm_comparison import FiveArmWeights


@dataclass(frozen=True)
class FiveArmLikelihoodResult:
    """Selected log likelihoods on the common H0 grid."""

    H0_grid: np.ndarray
    arm_a_log_likelihood: np.ndarray
    arm_b_log_likelihood: np.ndarray
    arm_c_log_likelihood: np.ndarray
    arm_d_log_likelihood: np.ndarray
    arm_e_node_log_likelihood: np.ndarray
    arm_e_marginal_log_likelihood: np.ndarray
    arm_e_node_log_detection_efficiency: np.ndarray


def _validate_grid(H0_grid: ArrayLike) -> np.ndarray:
    grid = np.asarray(H0_grid, dtype=float)
    if (
        grid.ndim != 1 or grid.size < 2 or not np.all(np.isfinite(grid))
        or np.any(grid <= 0) or np.any(np.diff(grid) <= 0)
    ):
        raise ValueError("H0_grid must be finite, positive, and strictly increasing")
    return grid


def _validate_events(events: Sequence[Mapping[str, float]], threshold: float) -> tuple[np.ndarray, np.ndarray]:
    if not events:
        raise ValueError("at least one detected event is required")
    distances = np.asarray([float(event["gw_distance"]) for event in events], dtype=float)
    sigmas = np.asarray([float(event["gw_sigma"]) for event in events], dtype=float)
    if not np.all(np.isfinite(distances)) or np.any(distances < threshold):
        raise ValueError("every event distance must be finite and meet the detection threshold")
    if not np.all(np.isfinite(sigmas)) or np.any(sigmas <= 0):
        raise ValueError("every event gw_sigma must be finite and positive")
    return distances, sigmas


def _single_population_log_likelihood(
    H0_grid: np.ndarray,
    redshift: np.ndarray,
    raw_weights: np.ndarray,
    distances: np.ndarray,
    sigmas: np.ndarray,
    threshold: float,
) -> tuple[np.ndarray, np.ndarray]:
    if (
        redshift.ndim != 1 or redshift.size < 2
        or raw_weights.shape != redshift.shape
        or not np.all(np.isfinite(redshift)) or np.any(redshift <= 0)
        or not np.all(np.isfinite(raw_weights)) or np.any(raw_weights <= 0)
    ):
        raise ValueError("host redshifts and raw weights must be finite, positive, and shape-matched")
    log_weights = np.log(raw_weights)
    log_weights -= logsumexp(log_weights)
    means = C_KM_S * redshift[None, :] / H0_grid[:, None]
    event_log_likelihood = np.zeros((distances.size, H0_grid.size), dtype=float)
    event_log_efficiency = np.zeros_like(event_log_likelihood)

    for event_index, (distance, sigma) in enumerate(zip(distances, sigmas)):
        event_log_likelihood[event_index] = logsumexp(
            log_weights[None, :]
            + norm.logpdf(distance, loc=means, scale=sigma),
            axis=1,
        )
        # P(measured distance >= threshold | z,H0) = Phi((mean-threshold)/sigma).
        event_log_efficiency[event_index] = logsumexp(
            log_weights[None, :]
            + log_ndtr((means - threshold) / sigma),
            axis=1,
        )

    selected = np.sum(event_log_likelihood - event_log_efficiency, axis=0)
    return selected, event_log_efficiency


def matched_five_arm_selected_log_likelihood(
    H0_grid: ArrayLike,
    arms: FiveArmWeights,
    events: Sequence[Mapping[str, float]],
    *,
    detection_threshold_mpc: float,
) -> FiveArmLikelihoodResult:
    """Evaluate A-D plug-in selected likelihoods and E's shared-node integral.

    The complete event likelihood and its event-specific selection denominator
    are evaluated for each arm. E uses a single shared coefficient node across
    all events, then marginalizes that joint likelihood once with the
    validation-posterior quadrature weights. No H0 posterior normalization is
    performed here.
    """
    grid = _validate_grid(H0_grid)
    threshold = float(detection_threshold_mpc)
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError("detection_threshold_mpc must be finite and positive")
    distances, sigmas = _validate_events(events, threshold)

    arm_a, _ = _single_population_log_likelihood(
        grid, arms.catalogue_redshift, arms.arm_a_catalogue_only,
        distances, sigmas, threshold,
    )
    arm_b, _ = _single_population_log_likelihood(
        grid, arms.parent_redshift, arms.arm_b_parent_oracle,
        distances, sigmas, threshold,
    )
    arm_c, _ = _single_population_log_likelihood(
        grid, arms.catalogue_redshift, arms.arm_c_known_probability,
        distances, sigmas, threshold,
    )
    arm_d, _ = _single_population_log_likelihood(
        grid, arms.catalogue_redshift, arms.arm_d_map_plugin,
        distances, sigmas, threshold,
    )

    node_weights = np.asarray(arms.arm_e_node_weights, dtype=float)
    quadrature_weights = np.asarray(arms.arm_e_quadrature_weights, dtype=float)
    if (
        node_weights.ndim != 2
        or node_weights.shape[1] != arms.catalogue_redshift.size
        or node_weights.shape[0] != quadrature_weights.size
        or not np.all(np.isfinite(quadrature_weights))
        or np.any(quadrature_weights < 0)
        or not np.isclose(np.sum(quadrature_weights), 1.0)
        or not np.any(quadrature_weights > 0)
    ):
        raise ValueError("Arm E node weights and quadrature weights are inconsistent")

    node_log_likelihood = np.empty((node_weights.shape[0], grid.size), dtype=float)
    node_log_efficiency = np.empty((node_weights.shape[0], distances.size, grid.size), dtype=float)
    for node_index, weights in enumerate(node_weights):
        node_log_likelihood[node_index], node_log_efficiency[node_index] = (
            _single_population_log_likelihood(
                grid, arms.catalogue_redshift, weights,
                distances, sigmas, threshold,
            )
        )
    log_quadrature = np.full(quadrature_weights.shape, -np.inf, dtype=float)
    positive = quadrature_weights > 0
    log_quadrature[positive] = np.log(quadrature_weights[positive])
    arm_e = logsumexp(node_log_likelihood + log_quadrature[:, None], axis=0)

    # Keep the per-event efficiency contributions: their sum is the node's
    # log(beta_k) total and allows an independent audit of the selected factor.
    return FiveArmLikelihoodResult(
        H0_grid=grid.copy(),
        arm_a_log_likelihood=arm_a,
        arm_b_log_likelihood=arm_b,
        arm_c_log_likelihood=arm_c,
        arm_d_log_likelihood=arm_d,
        arm_e_node_log_likelihood=node_log_likelihood,
        arm_e_marginal_log_likelihood=arm_e,
        arm_e_node_log_detection_efficiency=node_log_efficiency,
    )
