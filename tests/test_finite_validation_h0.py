import numpy as np
import pytest
from scipy.integrate import trapezoid
from scipy.special import logsumexp

from dark_siren_robustness.bayesian_inclusion import weighted_log_marginal_likelihood
from dark_siren_robustness.finite_validation_h0 import (
    marginalize_selection_conditioned_likelihood,
    normalize_h0_posterior,
    reconstruct_host_density_draws,
    selected_log_likelihood_by_draw,
)
from dark_siren_robustness.selection_toy import selection_aware_toy_likelihood


def _density_pair():
    z = np.linspace(0.01, 0.40, 151)
    base = 1.0 + 2.2 * np.exp(-0.5 * ((z - 0.105) / 0.025) ** 2)
    merger = 1.0 + 1.5 * (z - z[0]) / (z[-1] - z[0])
    first = base * merger
    first /= trapezoid(first, z)
    second = base * np.exp(0.8 * (z - z[0]) / (z[-1] - z[0])) * merger
    second /= trapezoid(second, z)
    return z, np.stack([first, second])


def test_vectorized_per_draw_likelihood_matches_existing_reference():
    z, densities = _density_pair()
    h0 = np.linspace(55.0, 85.0, 21)
    distances = np.array([270.0, 330.0, 410.0])
    result = selected_log_likelihood_by_draw(
        z, densities, h0, distances, sigma=140.0, detection_threshold=250.0,
        h0_batch_size=5,
    )
    assert result["selected_log_likelihoods"].shape == (2, 21)
    for draw in range(2):
        reference = selection_aware_toy_likelihood(
            z, densities[draw], h0, distances, sigma=140.0,
            detection_threshold=250.0,
        )
        np.testing.assert_allclose(
            result["unconditional_log_likelihoods"][draw],
            reference["unconditional_log_likelihood"], rtol=1e-11, atol=1e-11,
        )
        np.testing.assert_allclose(
            result["detection_efficiencies"][draw],
            reference["detection_efficiency"], rtol=1e-11, atol=1e-12,
        )
        np.testing.assert_allclose(
            result["selected_log_likelihoods"][draw],
            reference["selected_log_likelihood"], rtol=1e-11, atol=1e-11,
        )


def test_selection_correction_is_applied_inside_each_draw_before_marginalization():
    z, densities = _density_pair()
    h0 = np.linspace(55.0, 85.0, 15)
    distances = np.array([270.0, 330.0, 410.0])
    result = selected_log_likelihood_by_draw(
        z, densities, h0, distances, sigma=140.0, detection_threshold=250.0
    )
    log_weights = np.log(np.array([0.3, 0.7]))
    correct = marginalize_selection_conditioned_likelihood(
        result["selected_log_likelihoods"], log_weights
    )
    wrong_unconditional = weighted_log_marginal_likelihood(
        result["unconditional_log_likelihoods"], log_weights
    )
    average_beta = np.average(
        result["detection_efficiencies"], axis=0, weights=np.exp(log_weights)
    )
    wrong_shortcut = wrong_unconditional - distances.size * np.log(average_beta)
    assert not np.allclose(correct, wrong_shortcut, rtol=1e-6, atol=1e-6)


def test_weighted_h0_marginalization_matches_explicit_logsumexp():
    log_likelihoods = np.array([[-3.0, -1.0], [-0.5, -2.0], [-1.5, -0.2]])
    log_weights = np.log(np.array([0.1, 0.6, 0.3]))
    result = marginalize_selection_conditioned_likelihood(log_likelihoods, log_weights)
    expected = logsumexp(log_likelihoods + log_weights[:, None], axis=0) - logsumexp(log_weights)
    np.testing.assert_allclose(result, expected, atol=1e-12)


def test_h0_posterior_normalization_and_boundary_diagnostic():
    h0 = np.linspace(50.0, 90.0, 801)
    loglike = -0.5 * ((h0 - 70.0) / 4.0) ** 2
    prior = np.ones_like(h0)
    result = normalize_h0_posterior(h0, loglike, prior_density=prior)
    assert abs(result["posterior_integral"] - 1.0) < 1e-12
    assert abs(result["posterior_mean"] - 70.0) < 1e-8
    assert abs(result["posterior_median"] - 70.0) < 1e-8
    assert result["ci68_low"] < 70.0 < result["ci68_high"]
    assert result["ci95_low"] < 70.0 < result["ci95_high"]
    assert result["ci95_width"] > result["ci68_width"]
    assert result["boundary_mass_fraction"] < 1e-5


def test_coefficient_proposals_reconstruct_normalized_host_density_draws():
    rng = np.random.default_rng(101)
    analysis_z = rng.uniform(0.01, 0.40, 80)
    z_grid = np.linspace(0.01, 0.40, 201)
    merger = 1.0 + 1.5 * (z_grid - z_grid[0]) / (z_grid[-1] - z_grid[0])
    coefficients = np.array([[-0.7, -0.8], [-0.5, -0.6]])
    densities = reconstruct_host_density_draws(
        z_grid, analysis_z, coefficients, degree=1, merger_weight=merger,
        scott_multiplier=0.15,
    )
    assert densities.shape == (2, 201)
    np.testing.assert_allclose(trapezoid(densities, z_grid, axis=1), 1.0, atol=1e-10)


def test_invalid_selection_inputs_are_rejected():
    z, densities = _density_pair()
    h0 = np.linspace(55.0, 85.0, 10)
    with pytest.raises(ValueError):
        selected_log_likelihood_by_draw(
            z, densities, h0, np.array([100.0]), sigma=140.0, detection_threshold=250.0
        )
    with pytest.raises(ValueError):
        normalize_h0_posterior(h0, np.zeros(9), prior_density=np.ones(10))



def test_density_and_likelihood_draws_have_full_domain_and_memory_guards():
    rng = np.random.default_rng(101)
    analysis_z = rng.uniform(0.01, 0.40, 80)
    z_grid = np.linspace(0.01, 0.40, 201)
    merger = 1.0 + 1.5 * (z_grid - z_grid[0]) / (z_grid[-1] - z_grid[0])
    coefficients = np.array([[-0.7, -0.8], [-0.5, -0.6]])
    with pytest.raises(ValueError, match="full fixed inclusion-prior domain"):
        reconstruct_host_density_draws(
            z_grid[1:], analysis_z, coefficients, degree=1,
            merger_weight=merger[1:], scott_multiplier=0.15,
        )
    with pytest.raises(MemoryError, match="reduce draws or grid size"):
        reconstruct_host_density_draws(
            z_grid, analysis_z, coefficients, degree=1,
            merger_weight=merger, scott_multiplier=0.15, max_output_bytes=100,
        )

    densities = reconstruct_host_density_draws(
        z_grid, analysis_z, coefficients, degree=1,
        merger_weight=merger, scott_multiplier=0.15,
    )
    with pytest.raises(MemoryError, match="reduce draws, grid size"):
        selected_log_likelihood_by_draw(
            z_grid, densities, np.linspace(55.0, 85.0, 21),
            np.array([270.0, 330.0, 410.0]), sigma=140.0,
            detection_threshold=250.0, max_output_bytes=100,
        )


def test_real_weighted_inclusion_posterior_reaches_normalized_h0_posterior():
    from scipy.special import expit
    from dark_siren_robustness.bayesian_inclusion import (
        importance_sample_inclusion_posterior,
    )

    rng = np.random.default_rng(910)
    validation_z = rng.uniform(0.01, 0.40, 160)
    validation_u = (validation_z - 0.01) / 0.39
    validation_labels = rng.binomial(1, expit(-0.8 - 0.7 * validation_u))
    inclusion_posterior = importance_sample_inclusion_posterior(
        validation_z,
        validation_labels,
        degree=1,
        prior_sd=1.0,
        n_draws=1000,
        seed=12,
        min_ess=30,
        min_ess_fraction=0.01,
    )
    assert inclusion_posterior.status == "SUCCESS", inclusion_posterior.message

    analysis_parent_z = rng.uniform(0.01, 0.40, 400)
    analysis_u = (analysis_parent_z - 0.01) / 0.39
    analysis_labels = rng.binomial(1, expit(-0.8 - 0.7 * analysis_u))
    analysis_catalogue = analysis_parent_z[analysis_labels == 1]
    assert analysis_catalogue.size >= 20

    z_grid = np.linspace(0.01, 0.40, 101)
    h0_grid = np.linspace(55.0, 85.0, 41)
    merger = 1.0 + 1.5 * (z_grid - 0.01) / 0.39
    density_draws = reconstruct_host_density_draws(
        z_grid,
        analysis_catalogue,
        inclusion_posterior.coefficient_proposals,
        degree=1,
        merger_weight=merger,
        scott_multiplier=0.3,
    )
    likelihood = selected_log_likelihood_by_draw(
        z_grid,
        density_draws,
        h0_grid,
        np.array([330.0, 410.0, 520.0]),
        sigma=120.0,
        detection_threshold=250.0,
        h0_batch_size=8,
    )
    log_marginal = marginalize_selection_conditioned_likelihood(
        likelihood["selected_log_likelihoods"],
        inclusion_posterior.log_importance_weights,
    )
    posterior = normalize_h0_posterior(
        h0_grid, log_marginal, prior_density=np.ones_like(h0_grid)
    )
    assert density_draws.shape == (1000, 101)
    assert likelihood["selected_log_likelihoods"].shape == (1000, 41)
    assert abs(posterior["posterior_integral"] - 1.0) < 1e-10
    assert posterior["ci95_low"] < posterior["ci68_low"]
    assert posterior["ci95_high"] > posterior["ci68_high"]
    assert 0.0 <= posterior["boundary_mass_fraction"] <= 1.0
