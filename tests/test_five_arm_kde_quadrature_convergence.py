"""Coefficient quadrature convergence after propagation through the full H0 likelihood."""
import numpy as np

from dark_siren_robustness.five_arm_comparison import build_five_arm_weights
from dark_siren_robustness.five_arm_kde_likelihood import (
    matched_five_arm_kde_selected_log_likelihood,
)
from dark_siren_robustness.h0_posterior import (
    compare_h0_grid_posteriors,
    normalize_h0_posterior,
)


def test_shared_coefficient_quadrature_stabilizes_full_arm_e_posterior():
    rng = np.random.default_rng(20261114)
    parent = np.sort(rng.uniform(0.0101, 0.3999, 900))
    inclusion = 0.10 + 0.90 / (1.0 + np.exp((parent - 0.22) / 0.05))
    catalogue = parent[rng.random(parent.size) < inclusion]
    validation = np.sort(rng.uniform(0.0101, 0.3999, 300))
    p_validation = 0.10 + 0.90 / (1.0 + np.exp((validation - 0.22) / 0.05))
    retained = (rng.random(validation.size) < p_validation).astype(int)
    parent_base = 1.0 + 1.5 * (parent - 0.01) / 0.39
    catalogue_base = 1.0 + 1.5 * (catalogue - 0.01) / 0.39
    h0_grid = np.linspace(50.0, 90.0, 121)
    z_grid = np.linspace(0.01, 0.40, 501)
    events = [
        {"gw_distance": 275.0, "gw_sigma": 140.0},
        {"gw_distance": 340.0, "gw_sigma": 140.0},
        {"gw_distance": 420.0, "gw_sigma": 140.0},
    ]

    posteriors = {}
    for order in (7, 9, 11):
        arms = build_five_arm_weights(
            parent, catalogue, validation, retained,
            catalogue_floor=0.10, quadrature_order=order,
            parent_base_weights=parent_base,
            catalogue_base_weights=catalogue_base,
        )
        result = matched_five_arm_kde_selected_log_likelihood(
            h0_grid, arms, events,
            detection_threshold_mpc=250.0,
            redshift_grid=z_grid,
            merger_weight=np.ones_like(z_grid),
            scott_multiplier=0.7,
        )
        posteriors[order] = normalize_h0_posterior(
            h0_grid, result.arm_e_marginal_log_likelihood,
            prior_min=50.0, prior_max=90.0,
        )

    comparison = compare_h0_grid_posteriors(
        posteriors[9], posteriors[11],
        density_tolerance=0.005,
        quantile_tolerance=0.25,
    )
    assert comparison.passed, (
        "Arm E changed under coefficient quadrature refinement: "
        f"density={comparison.max_abs_density_difference:.6g}, "
        f"quantile={comparison.max_abs_quantile_difference:.6g} km/s"
    )
