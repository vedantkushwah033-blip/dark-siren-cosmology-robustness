"""Downstream H0-grid convergence for the matched five-arm KDE likelihood."""
import numpy as np
import pytest

from dark_siren_robustness.five_arm_comparison import build_five_arm_weights
from dark_siren_robustness.five_arm_kde_likelihood import (
    matched_five_arm_kde_selected_log_likelihood,
)
from dark_siren_robustness.h0_posterior import (
    compare_h0_grid_posteriors,
    normalize_h0_posterior,
)


def _fixed_replicate():
    rng = np.random.default_rng(20261113)
    parent = np.sort(rng.uniform(0.0101, 0.3999, 900))
    inclusion = 0.10 + 0.90 / (1.0 + np.exp((parent - 0.22) / 0.05))
    catalogue = parent[rng.random(parent.size) < inclusion]
    validation = np.sort(rng.uniform(0.0101, 0.3999, 300))
    p_validation = 0.10 + 0.90 / (1.0 + np.exp((validation - 0.22) / 0.05))
    retained = (rng.random(validation.size) < p_validation).astype(int)
    parent_base = 1.0 + 1.5 * (parent - 0.01) / 0.39
    catalogue_base = 1.0 + 1.5 * (catalogue - 0.01) / 0.39
    arms = build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.10, quadrature_order=5,
        parent_base_weights=parent_base,
        catalogue_base_weights=catalogue_base,
    )
    events = [
        {"gw_distance": 275.0, "gw_sigma": 140.0},
        {"gw_distance": 340.0, "gw_sigma": 140.0},
        {"gw_distance": 420.0, "gw_sigma": 140.0},
    ]
    return arms, events


def test_five_arm_kde_posteriors_stabilize_under_h0_grid_refinement():
    arms, events = _fixed_replicate()
    z_grid = np.linspace(0.01, 0.40, 501)
    outputs = {}
    for points in (81, 161, 321):
        h0_grid = np.linspace(50.0, 90.0, points)
        result = matched_five_arm_kde_selected_log_likelihood(
            h0_grid, arms, events,
            detection_threshold_mpc=250.0,
            redshift_grid=z_grid,
            merger_weight=np.ones_like(z_grid),
            scott_multiplier=0.7,
        )
        curves = {
            "A": result.arm_a_log_likelihood,
            "B": result.arm_b_log_likelihood,
            "C": result.arm_c_log_likelihood,
            "D": result.arm_d_log_likelihood,
            "E": result.arm_e_marginal_log_likelihood,
        }
        outputs[points] = {
            arm: normalize_h0_posterior(
                h0_grid, curve, prior_min=50.0, prior_max=90.0
            )
            for arm, curve in curves.items()
        }

    for arm in ("A", "B", "C", "D", "E"):
        coarse_to_fine = compare_h0_grid_posteriors(
            outputs[161][arm], outputs[321][arm],
            density_tolerance=0.003,
            quantile_tolerance=0.20,
        )
        assert coarse_to_fine.passed, (
            f"Arm {arm} did not stabilize from 161 to 321 H0 points: "
            f"density={coarse_to_fine.max_abs_density_difference:.6g}, "
            f"quantile={coarse_to_fine.max_abs_quantile_difference:.6g} km/s"
        )
