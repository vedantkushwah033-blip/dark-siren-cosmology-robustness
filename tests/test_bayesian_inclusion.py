import numpy as np
import pytest
from scipy.integrate import trapezoid
from scipy.special import expit

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


def test_proper_prior_keeps_one_class_validation_posterior_evaluable():
    z = np.linspace(0.01, 0.40, 40)
    y = np.ones(z.size, dtype=int)
    result = importance_sample_inclusion_posterior(
        z, y, degree=1, n_draws=3000, seed=12, min_ess=50, min_ess_fraction=0.01
    )
    assert result.status == "SUCCESS"
    assert np.all(np.isfinite(result.posterior_mean))
    assert result.importance_ess > 50


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
