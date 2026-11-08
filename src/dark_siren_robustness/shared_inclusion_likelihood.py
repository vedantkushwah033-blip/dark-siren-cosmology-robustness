"""Shared inclusion-posterior propagation through a discrete toy host mixture.

This is a synthetic low-redshift Euclidean likelihood primitive. One inclusion
coefficient vector is shared by every event at a quadrature node. The same
node-specific IPW population weights enter the event likelihood and the
one-sided GW detection efficiency. No event-rate model is included.
"""
from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike
from scipy.special import log_expit, logsumexp, log_ndtr

from .dataset_likelihood import C_KM_S, event_log_likelihood
from .inclusion_posterior import InclusionPosteriorQuadrature, log_marginalize_shared_parameter


@dataclass(frozen=True)
class SharedInclusionLikelihoodResult:
    """Node-conditional and validation-posterior-marginalized likelihoods."""

    node_log_selected_likelihood: np.ndarray
    node_log_detection_efficiency: np.ndarray
    marginal_log_likelihood: np.ndarray
    node_catalogue_weights: np.ndarray
    node_kish_ess: np.ndarray
    node_maximum_normalized_weight: np.ndarray
    node_minimum_log_inclusion_probability: np.ndarray


def shared_inclusion_selected_log_likelihood(
    H0_grid: ArrayLike,
    catalogue_redshift: ArrayLike,
    events: Sequence[Mapping[str, float]],
    inclusion_posterior: InclusionPosteriorQuadrature,
    *,
    detection_threshold_mpc: float,
    population_base_weights: ArrayLike | None = None,
    positive_truncated: bool = False,
) -> SharedInclusionLikelihoodResult:
    """Evaluate the selected joint likelihood at each shared theta node.

    Each event mapping must contain gw_distance and gw_sigma. The catalogue
    is discrete and fixed across theta nodes. Optional population_base_weights
    can encode a declared merger-rate weight; they are multiplied by 1/pi(z;theta)
    before normalization. The default corresponds to equal base population mass.

    The generator-selection convention is one-sided: measured distance >=
    detection_threshold_mpc. The returned marginal likelihood integrates the
    *joint* event likelihood over one shared validation-posterior theta node.
    It does not apply a second missing-host mixture correction.
    """
    h0 = np.asarray(H0_grid, dtype=float)
    z = np.asarray(catalogue_redshift, dtype=float)
    threshold = float(detection_threshold_mpc)
    if h0.ndim != 1 or h0.size < 2 or not np.all(np.isfinite(h0)) or np.any(h0 <= 0) or np.any(np.diff(h0) <= 0):
        raise ValueError("H0_grid must be finite, positive, and strictly increasing")
    if z.ndim != 1 or z.size < 2 or not np.all(np.isfinite(z)) or np.any(z < 0):
        raise ValueError("catalogue_redshift must be a finite non-negative 1D array with at least two hosts")
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError("detection_threshold_mpc must be finite and positive")
    if not events:
        raise ValueError("at least one detected event is required")
    if not isinstance(positive_truncated, (bool, np.bool_)):
        raise ValueError("positive_truncated must be boolean")

    base = np.ones(z.size, dtype=float) if population_base_weights is None else np.asarray(population_base_weights, dtype=float)
    if base.shape != z.shape or not np.all(np.isfinite(base)) or np.any(base < 0) or not np.any(base > 0):
        raise ValueError("population_base_weights must be finite, non-negative, match catalogue, and have positive mass")
    log_base = np.full(base.shape, -np.inf)
    positive_base = base > 0
    log_base[positive_base] = np.log(base[positive_base])

    nodes = np.asarray(inclusion_posterior.nodes, dtype=float)
    posterior_weights = np.asarray(inclusion_posterior.weights, dtype=float)
    if nodes.ndim != 2 or nodes.shape[0] != posterior_weights.size or nodes.shape[1] != inclusion_posterior.degree + 1:
        raise ValueError("inclusion posterior nodes have inconsistent dimensions")
    if not np.all(np.isfinite(nodes)) or posterior_weights.ndim != 1 or not np.all(np.isfinite(posterior_weights)) or np.any(posterior_weights < 0) or not np.isclose(posterior_weights.sum(), 1.0):
        raise ValueError("inclusion posterior nodes/weights must be finite with normalized non-negative weights")

    x = (z - inclusion_posterior.center) / inclusion_posterior.scale
    design = np.column_stack([x ** power for power in range(inclusion_posterior.degree + 1)])
    eta = design @ nodes.T
    if not np.all(np.isfinite(eta)):
        raise ArithmeticError("inclusion linear predictors are non-finite")
    log_probability = log_expit(eta)
    if not np.all(np.isfinite(log_probability)):
        raise ArithmeticError("log inclusion probabilities are non-finite; IPW weights are not evaluable")
    log_raw_weights = log_base[:, None] - log_probability
    log_normalizer = logsumexp(log_raw_weights, axis=0)
    log_host_weights = log_raw_weights - log_normalizer[None, :]
    host_weights = np.exp(log_host_weights)

    event_sigmas = np.asarray([float(event["gw_sigma"]) for event in events], dtype=float)
    if not np.all(np.isfinite(event_sigmas)) or np.any(event_sigmas <= 0):
        raise ValueError("every event gw_sigma must be finite and positive")
    if not np.allclose(event_sigmas, event_sigmas[0], rtol=1e-12, atol=0.0):
        raise ValueError("this toy selection implementation requires a common gw_sigma across events")

    n_nodes = nodes.shape[0]
    node_log_likelihood = np.zeros((n_nodes, h0.size), dtype=float)
    for event_index, event in enumerate(events):
        if not isinstance(event, Mapping) or set(event) != {"gw_distance", "gw_sigma"}:
            raise ValueError(f"events[{event_index}] must contain exactly gw_distance and gw_sigma")
        distance = float(event["gw_distance"])
        sigma = float(event["gw_sigma"])
        if not np.isfinite(distance) or distance < threshold:
            raise ValueError(f"events[{event_index}] must be finite and satisfy the one-sided detection threshold")
        if not np.isfinite(sigma) or sigma <= 0:
            raise ValueError(f"events[{event_index}].gw_sigma must be finite and positive")
        for node_index in range(n_nodes):
            node_log_likelihood[node_index] += event_log_likelihood(
                h0, redshifts=z, gw_distance=distance, gw_sigma=sigma,
                weights=host_weights[:, node_index],
                positive_truncated=positive_truncated,
            )

    # The generator accepts a noisy distance iff d_obs >= threshold. Use the
    # same weighted host population as the event mixture at each theta node.
    means = C_KM_S * z[None, :] / h0[:, None]
    if not np.all(np.isfinite(means)):
        raise ArithmeticError("predicted host distances are non-finite on the H0 grid")
    log_detection_by_host = log_ndtr((means - threshold) / event_sigmas[0])
    node_log_beta = np.empty((n_nodes, h0.size), dtype=float)
    for node_index in range(n_nodes):
        node_log_beta[node_index] = logsumexp(
            log_detection_by_host + log_host_weights[:, node_index][None, :], axis=1
        )

    if not np.all(np.isfinite(node_log_beta)):
        raise ArithmeticError("one-sided detection efficiency is non-finite")
    node_log_selected = node_log_likelihood - len(events) * node_log_beta
    marginal = log_marginalize_shared_parameter(node_log_selected, posterior_weights)
    kish_ess = 1.0 / np.sum(host_weights * host_weights, axis=0)
    maximum_weight = np.max(host_weights, axis=0)
    minimum_log_probability = np.min(log_probability, axis=0)
    return SharedInclusionLikelihoodResult(
        node_log_selected_likelihood=node_log_selected,
        node_log_detection_efficiency=node_log_beta,
        marginal_log_likelihood=marginal,
        node_catalogue_weights=host_weights,
        node_kish_ess=kish_ess,
        node_maximum_normalized_weight=maximum_weight,
        node_minimum_log_inclusion_probability=minimum_log_probability,
    )
