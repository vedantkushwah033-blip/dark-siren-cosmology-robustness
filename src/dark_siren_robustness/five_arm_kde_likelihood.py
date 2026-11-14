"""Five-arm selected H0 likelihoods using a common reflected weighted KDE.

All arms use the same finite redshift grid and weighted-Scott multiplier. Arm E integrates the complete selected
multi-event likelihood over shared inclusion-coefficient nodes.
"""
from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike
from scipy.integrate import trapezoid
from scipy.special import log_ndtr, logsumexp
from scipy.stats import norm

from .catalogue_density import estimate_host_density
from .dataset_likelihood import C_KM_S
from .five_arm_comparison import FiveArmWeights


@dataclass(frozen=True)
class FiveArmKDELikelihoodResult:
    H0_grid: np.ndarray
    redshift_grid: np.ndarray
    arm_a_log_likelihood: np.ndarray
    arm_b_log_likelihood: np.ndarray
    arm_c_log_likelihood: np.ndarray
    arm_d_log_likelihood: np.ndarray
    arm_e_node_log_likelihood: np.ndarray
    arm_e_marginal_log_likelihood: np.ndarray
    arm_e_node_log_detection_efficiency: np.ndarray
    arm_density: dict[str, np.ndarray]
    arm_e_node_density: np.ndarray
    scott_multiplier: float | None
    absolute_bandwidth_sd: float | None


def _validate_grid(values: ArrayLike, name: str, minimum: int = 21) -> np.ndarray:
    grid = np.asarray(values, dtype=float)
    if (grid.ndim != 1 or grid.size < minimum or not np.all(np.isfinite(grid))
            or np.any(np.diff(grid) <= 0) or grid[0] <= 0):
        raise ValueError(f"{name} must be finite, positive, strictly increasing, and have at least {minimum} points")
    return grid


def _trapezoid_weights(grid: np.ndarray) -> np.ndarray:
    """Return linear weights whose dot product reproduces composite trapezoids."""
    delta = np.diff(grid)
    weights = np.empty_like(grid)
    weights[0] = 0.5 * delta[0]
    weights[-1] = 0.5 * delta[-1]
    weights[1:-1] = 0.5 * (delta[:-1] + delta[1:])
    return weights


def _batched_density_likelihood(
    H0_grid: np.ndarray,
    z_grid: np.ndarray,
    densities: np.ndarray,
    distances: np.ndarray,
    sigmas: np.ndarray,
    threshold: float,
    *,
    h0_chunk_size: int = 128,
) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate many host densities while reusing each event/H0 redshift kernel.

    Returns log likelihoods with shape (density, H0) and log detection
    efficiencies with shape (density, event, H0). H0 chunking bounds the
    temporary event/H0/redshift kernel memory; the returned likelihood and
    efficiency arrays still scale with the full H0 grid.
    """
    if (densities.ndim != 2 or densities.shape[1] != z_grid.size
            or not np.all(np.isfinite(densities)) or np.any(densities < 0)):
        raise ValueError("host densities must be finite, nonnegative, and match redshift_grid")
    areas = np.asarray([trapezoid(row, z_grid) for row in densities])
    if not np.allclose(areas, 1.0, rtol=1e-9, atol=1e-9):
        raise ValueError("every host density must integrate to one on redshift_grid")
    if h0_chunk_size < 1:
        raise ValueError("h0_chunk_size must be positive")

    quadrature = _trapezoid_weights(z_grid)
    n_density, n_h0 = densities.shape[0], H0_grid.size
    log_numerator = np.empty((n_density, distances.size, n_h0))
    log_efficiency = np.empty_like(log_numerator)

    for start in range(0, n_h0, h0_chunk_size):
        stop = min(start + h0_chunk_size, n_h0)
        means = C_KM_S * z_grid[None, :] / H0_grid[start:stop, None]
        for event_index, (distance, sigma) in enumerate(zip(distances, sigmas)):
            # The expensive event/H0/redshift kernels are formed once and
            # contracted against every arm/node density by matrix multiply.
            numerator_kernel = norm.pdf(distance, loc=means, scale=sigma) * quadrature[None, :]
            efficiency_kernel = np.exp(log_ndtr((means - threshold) / sigma)) * quadrature[None, :]
            numerator = numerator_kernel @ densities.T
            efficiency = efficiency_kernel @ densities.T
            log_numerator[:, event_index, start:stop] = np.log(numerator.T)
            log_efficiency[:, event_index, start:stop] = np.log(efficiency.T)

    if not np.all(np.isfinite(log_numerator)) or not np.all(np.isfinite(log_efficiency)):
        raise ArithmeticError("KDE-integrated event likelihood or selection efficiency is non-finite")
    log_likelihood = np.sum(log_numerator - log_efficiency, axis=1)
    return log_likelihood, log_efficiency


def matched_five_arm_kde_selected_log_likelihood(
    H0_grid: ArrayLike,
    arms: FiveArmWeights,
    events: Sequence[Mapping[str, float]],
    *,
    detection_threshold_mpc: float,
    redshift_grid: ArrayLike,
    merger_weight: ArrayLike,
    scott_multiplier: float | None = None,
    absolute_bandwidth_sd: float | None = None,
) -> FiveArmKDELikelihoodResult:
    """Evaluate all five arms through one shared reflected-KDE bandwidth rule.

    Supply either a common Scott multiplier or a common absolute redshift-space
    kernel SD. FiveArmWeights carries the complete
    base-population and inverse-inclusion weights. Use a unit merger_weight
    when those arrays already include the base measure.
    """
    h0 = np.asarray(H0_grid, dtype=float)
    if h0.ndim != 1 or h0.size < 2 or not np.all(np.isfinite(h0)) or np.any(h0 <= 0) or np.any(np.diff(h0) <= 0):
        raise ValueError("H0_grid must be finite, positive, and strictly increasing")
    z_grid = _validate_grid(redshift_grid, "redshift_grid")
    rate = np.asarray(merger_weight, dtype=float)
    if rate.shape != z_grid.shape or not np.all(np.isfinite(rate)) or np.any(rate <= 0):
        raise ValueError("merger_weight must be finite, positive, and match redshift_grid")
    if (scott_multiplier is None) == (absolute_bandwidth_sd is None):
        raise ValueError("supply exactly one of scott_multiplier or absolute_bandwidth_sd")
    if scott_multiplier is not None and (not np.isfinite(scott_multiplier) or scott_multiplier <= 0):
        raise ValueError("scott_multiplier must be finite and positive")
    if absolute_bandwidth_sd is not None and (not np.isfinite(absolute_bandwidth_sd) or absolute_bandwidth_sd <= 0):
        raise ValueError("absolute_bandwidth_sd must be finite and positive")
    threshold = float(detection_threshold_mpc)
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError("detection_threshold_mpc must be finite and positive")
    if not events:
        raise ValueError("at least one detected event is required")
    distances = np.asarray([float(e["gw_distance"]) for e in events])
    sigmas = np.asarray([float(e["gw_sigma"]) for e in events])
    if (not np.all(np.isfinite(distances)) or np.any(distances < threshold)
            or not np.all(np.isfinite(sigmas)) or np.any(sigmas <= 0)):
        raise ValueError("events must have finite detected distances and positive finite uncertainties")

    # FiveArmWeights already contains b(z), b(z)/pi_true(z),
    # b(z)/pi_MAP(z), and b(z)/pi(z;theta_k). Pass those complete raw
    # galaxy weights directly; using only 1/pi here would silently drop the
    # common base-population/merger-rate factor and make the arms incomparable.
    def density(sample: np.ndarray, raw_weights: np.ndarray) -> np.ndarray:
        return estimate_host_density(
            z_grid, sample, sample_weights=raw_weights,
            merger_weight=rate, scott_multiplier=scott_multiplier,
            absolute_bandwidth_sd=absolute_bandwidth_sd,
        )

    p_a = density(arms.catalogue_redshift, arms.arm_a_catalogue_only)
    p_b = density(arms.parent_redshift, arms.arm_b_parent_oracle)
    p_c = density(arms.catalogue_redshift, arms.arm_c_known_probability)
    p_d = density(arms.catalogue_redshift, arms.arm_d_map_plugin)
    node_weights = np.asarray(arms.arm_e_node_weights, dtype=float)
    if node_weights.ndim != 2 or node_weights.shape[1] != arms.catalogue_redshift.size:
        raise ValueError("Arm E node weights have invalid shape")
    p_e = np.stack([density(arms.catalogue_redshift, w) for w in node_weights])
    q = np.asarray(arms.arm_e_quadrature_weights, dtype=float)
    if (q.shape != (p_e.shape[0],) or not np.all(np.isfinite(q))
            or np.any(q < 0) or not np.isclose(np.sum(q), 1.0) or not np.any(q > 0)):
        raise ValueError("Arm E quadrature weights are invalid")

    densities = {
        "A_catalogue_only": p_a, "B_parent_oracle": p_b,
        "C_known_probability_ipw": p_c, "D_map_plugin_ipw": p_d,
    }
    # Stack A-D and every E quadrature-node density so each event/H0 kernel
    # is evaluated once and reused across all arms and coefficient nodes.
    density_names = list(densities)
    all_densities = np.vstack([*(densities[name] for name in density_names), p_e])
    all_log_likelihood, all_log_efficiency = _batched_density_likelihood(
        h0, z_grid, all_densities, distances, sigmas, threshold,
    )
    curves = {
        name: all_log_likelihood[index]
        for index, name in enumerate(density_names)
    }
    node_ll = all_log_likelihood[len(density_names):]
    node_eff = all_log_efficiency[len(density_names):]
    positive = q > 0
    marginal = logsumexp(node_ll[positive] + np.log(q[positive, None]), axis=0)
    return FiveArmKDELikelihoodResult(
        H0_grid=h0.copy(), redshift_grid=z_grid.copy(),
        arm_a_log_likelihood=curves["A_catalogue_only"],
        arm_b_log_likelihood=curves["B_parent_oracle"],
        arm_c_log_likelihood=curves["C_known_probability_ipw"],
        arm_d_log_likelihood=curves["D_map_plugin_ipw"],
        arm_e_node_log_likelihood=node_ll,
        arm_e_marginal_log_likelihood=marginal,
        arm_e_node_log_detection_efficiency=node_eff,
        arm_density=densities, arm_e_node_density=p_e,
        scott_multiplier=None if scott_multiplier is None else float(scott_multiplier),
        absolute_bandwidth_sd=None if absolute_bandwidth_sd is None else float(absolute_bandwidth_sd),
    )
