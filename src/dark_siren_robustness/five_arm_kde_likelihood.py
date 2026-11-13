"""Five-arm selected H0 likelihoods using a common reflected weighted KDE.

All arms use the same finite redshift grid, merger-rate function, and absolute
SciPy covariance bandwidth factor. Arm E integrates the complete selected
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


def _validate_grid(values: ArrayLike, name: str, minimum: int = 21) -> np.ndarray:
    grid = np.asarray(values, dtype=float)
    if (grid.ndim != 1 or grid.size < minimum or not np.all(np.isfinite(grid))
            or np.any(np.diff(grid) <= 0) or grid[0] <= 0):
        raise ValueError(f"{name} must be finite, positive, strictly increasing, and have at least {minimum} points")
    return grid


def _density_likelihood(
    H0_grid: np.ndarray, z_grid: np.ndarray, density: np.ndarray,
    distances: np.ndarray, sigmas: np.ndarray, threshold: float,
) -> tuple[np.ndarray, np.ndarray]:
    if density.shape != z_grid.shape or not np.all(np.isfinite(density)) or np.any(density < 0):
        raise ValueError("host density must be finite, nonnegative, and match redshift_grid")
    area = float(trapezoid(density, z_grid))
    if not np.isclose(area, 1.0, rtol=1e-9, atol=1e-9):
        raise ValueError("host density must integrate to one on redshift_grid")
    means = C_KM_S * z_grid[None, :] / H0_grid[:, None]
    event_log_likelihood = np.empty((distances.size, H0_grid.size))
    event_log_efficiency = np.empty_like(event_log_likelihood)
    for i, (distance, sigma) in enumerate(zip(distances, sigmas)):
        numerator = density[None, :] * norm.pdf(distance, loc=means, scale=sigma)
        efficiency = density[None, :] * np.exp(log_ndtr((means - threshold) / sigma))
        event_log_likelihood[i] = np.log(trapezoid(numerator, z_grid, axis=1))
        event_log_efficiency[i] = np.log(trapezoid(efficiency, z_grid, axis=1))
    if not np.all(np.isfinite(event_log_likelihood)) or not np.all(np.isfinite(event_log_efficiency)):
        raise ArithmeticError("KDE-integrated event likelihood or selection efficiency is non-finite")
    return np.sum(event_log_likelihood - event_log_efficiency, axis=0), event_log_efficiency


def matched_five_arm_kde_selected_log_likelihood(
    H0_grid: ArrayLike,
    arms: FiveArmWeights,
    events: Sequence[Mapping[str, float]],
    *,
    detection_threshold_mpc: float,
    redshift_grid: ArrayLike,
    merger_weight: ArrayLike,
    bandwidth_factor: float,
) -> FiveArmKDELikelihoodResult:
    """Evaluate all five arms through one fixed-bandwidth reflected weighted KDE."""
    h0 = np.asarray(H0_grid, dtype=float)
    if h0.ndim != 1 or h0.size < 2 or not np.all(np.isfinite(h0)) or np.any(h0 <= 0) or np.any(np.diff(h0) <= 0):
        raise ValueError("H0_grid must be finite, positive, and strictly increasing")
    z_grid = _validate_grid(redshift_grid, "redshift_grid")
    rate = np.asarray(merger_weight, dtype=float)
    if rate.shape != z_grid.shape or not np.all(np.isfinite(rate)) or np.any(rate <= 0):
        raise ValueError("merger_weight must be finite, positive, and match redshift_grid")
    if not np.isfinite(bandwidth_factor) or bandwidth_factor <= 0:
        raise ValueError("bandwidth_factor must be finite and positive")
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
            merger_weight=rate, bandwidth_factor=float(bandwidth_factor),
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
    curves = {}
    for name, p in densities.items():
        curves[name], _ = _density_likelihood(h0, z_grid, p, distances, sigmas, threshold)
    node_ll = np.empty((p_e.shape[0], h0.size))
    node_eff = np.empty((p_e.shape[0], distances.size, h0.size))
    for j, p in enumerate(p_e):
        node_ll[j], node_eff[j] = _density_likelihood(h0, z_grid, p, distances, sigmas, threshold)
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
    )
