import numpy as np
import pytest
from scipy.integrate import trapezoid
from scipy.special import expit, logsumexp

from dark_siren_robustness.inclusion_posterior import sample_inclusion_posterior
from dark_siren_robustness.bayesian_inclusion import (
    fixed_redshift_design,
    importance_sample_inclusion_posterior,
    weighted_log_marginal_likelihood,
)


def test_fixed_basis_is_domain_anchored_not_sample_standardized():
    z = np.array([0.01, 0.10, 0.40])
    design = fixed_redshift_design(z, degree=2)
    np.testing.assert_allclose(design, [[1, 0, 0], [1, 0.09 / 0.39, (0.09 / 0.39) ** 2], [1, 1, 1]])
    assert fixed_redshift_design(z[::-1], degree=2)[0, 1] == 1.0


def test_importance_sampler_returns_normalized_weights_and_posterior_predictions():
    rng = np.random.default_rng(18)
    z = np.linspace(0.01, 0.40, 100)
    u = (z - 0.01) / 0.39
    y = rng.binomial(1, expit(-0.7 - 0.8 * u))
    result = importance_sample_inclusion_posterior(
        z, y, degree=1, n_draws=5000, seed=31, min_ess=50, min_ess_fraction=0.01
    )
    assert result.status == "SUCCESS"
    assert result.importance_ess > 50
    assert result.max_normalized_weight < 0.02
    np.testing.assert_allclose(result.normalized_importance_weights.sum(), 1.0, atol=1e-12)
    predicted = result.predict(np.array([0.01, 0.20, 0.40]))
    assert predicted.shape == (3,)
    assert np.all((predicted > 0) & (predicted < 1))


def test_one_class_validation_has_a_proper_prior_posterior():
    z = np.linspace(0.01, 0.40, 40)
    for outcome in (0, 1):
        y = np.full(z.size, outcome, dtype=int)
        result = importance_sample_inclusion_posterior(
            z, y, degree=1, n_draws=5000, seed=12 + outcome,
            min_ess=50, min_ess_fraction=0.01
        )
        assert result.status == "SUCCESS", result.message
        assert result.coefficient_proposals.shape == (5000, 2)
        assert np.all(np.isfinite(result.posterior_mean))
        np.testing.assert_allclose(result.normalized_importance_weights.sum(), 1.0, atol=1e-12)


def test_importance_sampler_matches_existing_multichain_mcmc_posterior():
    rng = np.random.default_rng(81)
    z = np.linspace(0.01, 0.40, 160)
    u = (z - 0.01) / 0.39
    y = rng.binomial(1, expit(-0.7 - 0.8 * u))
    mcmc = sample_inclusion_posterior(
        z, y, degree=1, warmup=1000, draws_per_chain=1000, seed=82
    )
    importance = importance_sample_inclusion_posterior(
        z, y, degree=1, n_draws=20000, seed=83, min_ess=200,
        min_ess_fraction=0.05
    )
    assert mcmc.status == "SUCCESS", mcmc.message
    assert importance.status == "SUCCESS", importance.message
    mcmc_mean = mcmc.draws.mean(axis=(0, 1))
    np.testing.assert_allclose(importance.posterior_mean, mcmc_mean, atol=0.12, rtol=0)

    def toy_log_likelihood(beta):
        return np.column_stack([
            -0.5 * (beta[:, 0] - 0.2) ** 2 - 0.2 * (beta[:, 1] + 0.6) ** 2,
            -0.5 * (beta[:, 0] + 0.5) ** 2 - 0.15 * (beta[:, 1] - 0.1) ** 2,
        ])

    mcmc_beta = mcmc.draws.reshape(-1, 2)
    mcmc_log_marginal = logsumexp(toy_log_likelihood(mcmc_beta), axis=0) - np.log(len(mcmc_beta))
    importance_log_marginal = weighted_log_marginal_likelihood(
        toy_log_likelihood(importance.coefficient_proposals),
        importance.log_importance_weights,
    )
    np.testing.assert_allclose(importance_log_marginal, mcmc_log_marginal, atol=0.10, rtol=0)


def test_rank_deficient_validation_design_is_non_evaluable():
    z = np.full(30, 0.2)
    y = np.tile([0, 1], 15)
    result = importance_sample_inclusion_posterior(
        z, y, degree=1, n_draws=3000, seed=13
    )
    assert result.status == "NON_EVALUABLE"
    assert result.coefficient_proposals is None


def test_importance_sampler_matches_two_dimensional_grid_quadrature():
    rng = np.random.default_rng(202)
    z = np.linspace(0.01, 0.40, 36)
    u = (z - 0.01) / 0.39
    y = rng.binomial(1, expit(-0.6 - 0.9 * u))
    result = importance_sample_inclusion_posterior(
        z, y, degree=1, n_draws=30000, seed=77, min_ess=100, min_ess_fraction=0.02
    )
    assert result.status == "SUCCESS"

    b0 = np.linspace(-3.5, 2.0, 501)
    b1 = np.linspace(-3.5, 2.0, 501)
    g0, g1 = np.meshgrid(b0, b1, indexing="ij")
    eta = g0[..., None] + g1[..., None] * u
    log_like = np.sum(y * eta - np.logaddexp(0.0, eta), axis=-1)
    log_prior = -0.5 * (g0**2 + g1**2) - np.log(2 * np.pi)
    log_mass = log_like + log_prior
    mass = np.exp(log_mass - np.max(log_mass))
    norm = trapezoid(trapezoid(mass, b1, axis=1), b0)
    mean0 = trapezoid(trapezoid(g0 * mass, b1, axis=1), b0) / norm
    mean1 = trapezoid(trapezoid(g1 * mass, b1, axis=1), b0) / norm
    np.testing.assert_allclose(result.posterior_mean, [mean0, mean1], atol=0.035, rtol=0)


def test_weighted_log_marginal_likelihood_uses_self_normalized_importance_weights():
    log_likelihood = np.array([[-2.0, -1.0], [-0.5, -3.0], [-1.0, -2.0]])
    log_weights = np.log(np.array([0.2, 0.3, 0.5]))
    result = weighted_log_marginal_likelihood(log_likelihood, log_weights)
    expected = np.log(np.sum(np.exp(log_likelihood) * np.array([0.2, 0.3, 0.5])[:, None], axis=0))
    np.testing.assert_allclose(result, expected, atol=1e-12)


def test_weighted_log_marginal_likelihood_is_invariant_to_common_weight_scale():
    ll = np.array([-2.0, -1.0, -3.0])
    lw = np.log(np.array([1.0, 4.0, 2.0]))
    first = weighted_log_marginal_likelihood(ll, lw)
    second = weighted_log_marginal_likelihood(ll, lw + 50.0)
    np.testing.assert_allclose(first, second, atol=1e-12)


def test_invalid_validation_and_proposal_inputs_are_rejected():
    z = np.linspace(0.01, 0.40, 10)
    y = np.zeros(10)
    with pytest.raises(ValueError):
        fixed_redshift_design(z, degree=3)
    with pytest.raises(ValueError):
        importance_sample_inclusion_posterior(z, y, degree=1, n_draws=500)
    with pytest.raises(ValueError):
        weighted_log_marginal_likelihood(np.array([1.0, np.nan]), np.array([0.0, 0.0]))
