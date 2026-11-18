"""Curved-logit DGP regression for shared-node selected H0 quadrature."""
import numpy as np
from scipy.integrate import trapezoid
from scipy.special import expit

from dark_siren_robustness.five_arm_comparison import build_five_arm_weights
from dark_siren_robustness.five_arm_kde_likelihood import (
    matched_five_arm_kde_selected_log_likelihood,
)
from dark_siren_robustness.h0_posterior import (
    compare_h0_grid_posteriors,
    normalize_h0_posterior,
)


def test_curved_inclusion_dgp_arm_e_h0_posterior_stabilizes_with_quadrature_order():
    rng = np.random.default_rng(20261118)
    parent = np.sort(rng.uniform(0.0101, 0.3999, 1200))
    x_parent = (parent - 0.205) / 0.195
    gamma = 1.2
    alpha = -0.15
    beta = -3.9
    center_square = 0.2845005579191241
    true_probability = expit(
        alpha + beta * x_parent + gamma * (x_parent**2 - center_square)
    )
    catalogue_mask = rng.random(parent.size) < true_probability
    catalogue = parent[catalogue_mask]
    x_catalogue = (catalogue - 0.205) / 0.195
    catalogue_probability = expit(
        alpha + beta * x_catalogue
        + gamma * (x_catalogue**2 - center_square)
    )

    validation = np.sort(rng.uniform(0.0101, 0.3999, 300))
    x_validation = (validation - 0.205) / 0.195
    validation_probability = expit(
        alpha + beta * x_validation
        + gamma * (x_validation**2 - center_square)
    )
    retained = (rng.random(validation.size) < validation_probability).astype(int)

    parent_base = 1.0 + 1.5 * (parent - 0.01) / 0.39
    catalogue_base = 1.0 + 1.5 * (catalogue - 0.01) / 0.39
    events = [
        {"gw_distance": 275.0, "gw_sigma": 140.0},
        {"gw_distance": 340.0, "gw_sigma": 140.0},
        {"gw_distance": 420.0, "gw_sigma": 140.0},
    ]
    h0_grid = np.linspace(50.0, 90.0, 121)
    z_grid = np.linspace(0.01, 0.40, 501)
    posteriors = {}
    posterior_mean_densities = {}

    for order in (7, 9, 11):
        arms = build_five_arm_weights(
            parent,
            catalogue,
            validation,
            retained,
            degree=1,
            coordinate_center=0.205,
            coordinate_scale=0.195,
            quadrature_order=order,
            parent_base_weights=parent_base,
            catalogue_base_weights=catalogue_base,
            true_catalogue_inclusion_probability=catalogue_probability,
        )
        result = matched_five_arm_kde_selected_log_likelihood(
            h0_grid,
            arms,
            events,
            detection_threshold_mpc=250.0,
            redshift_grid=z_grid,
            merger_weight=np.ones_like(z_grid),
            scott_multiplier=0.7,
        )
        posterior_mean_densities[order] = np.average(
            result.arm_e_node_density,
            axis=0,
            weights=arms.arm_e_quadrature_weights,
        )
        assert np.isclose(
            trapezoid(posterior_mean_densities[order], z_grid),
            1.0,
            rtol=1e-8,
            atol=1e-8,
        )
        posteriors[order] = normalize_h0_posterior(
            h0_grid,
            result.arm_e_marginal_log_likelihood,
            prior_min=50.0,
            prior_max=90.0,
        )

    comparison = compare_h0_grid_posteriors(
        posteriors[9],
        posteriors[11],
        density_tolerance=0.01,
        quantile_tolerance=0.5,
    )
    assert comparison.passed, (
        "Curved-DGP Arm E changed under coefficient quadrature refinement: "
        f"density={comparison.max_abs_density_difference:.6g}, "
        f"quantile={comparison.max_abs_quantile_difference:.6g} km/s"
    )

    density_total_variation = 0.5 * trapezoid(
        np.abs(posterior_mean_densities[9] - posterior_mean_densities[11]),
        z_grid,
    )
    assert density_total_variation <= 0.01, (
        "Curved-DGP posterior-mean Arm-E host density changed under "
        f"quadrature refinement: TV={density_total_variation:.6g}"
    )
