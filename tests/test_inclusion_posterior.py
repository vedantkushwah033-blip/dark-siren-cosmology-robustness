import numpy as np
import pytest
from scipy.special import logsumexp

from dark_siren_robustness.inclusion_posterior import (
    fit_inclusion_posterior_quadrature,
    log_marginalize_shared_parameter,
    posterior_predictive_inclusion,
)


def test_proper_posterior_quadrature_normalizes_and_predicts():
    rng = np.random.default_rng(4182)
    z = rng.uniform(0.01, 0.40, 500)
    x = (z - z.mean()) / z.std()
    probability = 1.0 / (1.0 + np.exp(-(-0.35 + 0.85 * x + 0.25 * x * x)))
    retained = rng.binomial(1, probability)

    posterior = fit_inclusion_posterior_quadrature(
        z, retained, degree=2, order=13
    )
    assert posterior.nodes.shape == (13**3, 3)
    assert posterior.weights.shape == (13**3,)
    assert np.all(np.isfinite(posterior.nodes))
    assert np.all(np.isfinite(posterior.weights))
    assert np.all(posterior.weights >= 0)
    assert np.isclose(posterior.weights.sum(), 1.0)
    assert 1.0 <= posterior.quadrature_ess <= posterior.weights.size
    assert 0 < posterior.maximum_normalized_weight <= 1
    assert np.all(np.linalg.eigvalsh(posterior.posterior_covariance) > 0)

    predictive = posterior_predictive_inclusion(posterior, np.linspace(0.01, 0.40, 51))
    assert predictive.shape == (51,)
    assert np.all(np.isfinite(predictive))
    assert np.all((predictive > 0) & (predictive < 1))


def test_quadrature_order_refinement_is_stable_for_regular_validation_sample():
    rng = np.random.default_rng(9102)
    z = rng.uniform(0.01, 0.40, 800)
    x = (z - z.mean()) / z.std()
    p = 1.0 / (1.0 + np.exp(-(-0.2 + 0.65*x - 0.18*x*x)))
    retained = rng.binomial(1, p)
    low = fit_inclusion_posterior_quadrature(z, retained, degree=2, order=9)
    high = fit_inclusion_posterior_quadrature(z, retained, degree=2, order=15)
    assert np.max(np.abs(low.posterior_mean - high.posterior_mean)) < 0.04
    grid = np.linspace(0.01, 0.40, 101)
    assert np.max(np.abs(
        posterior_predictive_inclusion(low, grid)
        - posterior_predictive_inclusion(high, grid)
    )) < 0.015


def test_proper_prior_keeps_one_class_validation_posterior_finite():
    z = np.linspace(0.01, 0.40, 80)
    retained = np.ones(z.size, dtype=int)
    posterior = fit_inclusion_posterior_quadrature(
        z, retained, degree=1, order=11
    )
    assert np.all(np.isfinite(posterior.nodes))
    assert np.isclose(posterior.weights.sum(), 1.0)
    assert np.all(np.isfinite(posterior.posterior_mean))


def test_rejects_rank_deficient_or_invalid_validation_data():
    with pytest.raises(ValueError, match="positive finite spread"):
        fit_inclusion_posterior_quadrature(np.ones(20), np.zeros(20), degree=1)
    z = np.linspace(0.01, 0.40, 20)
    bad = np.zeros(20)
    bad[0] = 2
    with pytest.raises(ValueError, match="binary"):
        fit_inclusion_posterior_quadrature(z, bad, degree=1)
    with pytest.raises(ValueError, match="prior_sd"):
        fit_inclusion_posterior_quadrature(
            z, np.arange(20) % 2, degree=1, prior_sd=[5.0, 0.0]
        )


def test_shared_parameter_marginalization_is_not_eventwise_factorization():
    # Two events prefer opposite values of the same shared inclusion parameter.
    # E_theta[L1(theta)L2(theta)] = 0.5*(1*9 + 9*1) = 9,
    # while E[L1]E[L2] = 5*5 = 25.
    event_log_likelihoods = np.log(np.asarray([[1.0, 9.0], [9.0, 1.0]]))
    joint_by_node = event_log_likelihoods.sum(axis=1)
    marginalized_joint = log_marginalize_shared_parameter(
        joint_by_node, np.asarray([0.5, 0.5])
    )
    eventwise = log_marginalize_shared_parameter(
        event_log_likelihoods, np.asarray([0.5, 0.5])
    )
    assert np.isclose(np.exp(marginalized_joint), 9.0)
    assert np.isclose(np.prod(np.exp(eventwise)), 25.0)
    assert not np.isclose(marginalized_joint, np.sum(eventwise))


def test_log_marginalization_is_stable_to_common_likelihood_offset():
    log_likelihood = np.asarray([[-10000.0, -10002.0], [-10001.0, -9999.0]])
    weights = np.asarray([0.3, 0.7])
    result = log_marginalize_shared_parameter(log_likelihood, weights)
    shifted = log_marginalize_shared_parameter(log_likelihood + 12000.0, weights)
    assert np.allclose(shifted - 12000.0, result)


def test_log_marginalization_validates_weights_and_log_likelihoods():
    with pytest.raises(ValueError, match="sum to one"):
        log_marginalize_shared_parameter(np.zeros(2), np.asarray([0.2, 0.2]))
    with pytest.raises(ValueError, match="NaN"):
        log_marginalize_shared_parameter(np.asarray([0.0, np.nan]), np.asarray([0.5, 0.5]))
    with pytest.raises(ValueError, match="match posterior weights"):
        log_marginalize_shared_parameter(np.zeros((3, 2)), np.asarray([0.5, 0.5]))


def test_quadrature_ladder_reports_convergence_and_nonconvergence_explicitly():
    from dark_siren_robustness.inclusion_posterior import (
        fit_inclusion_posterior_quadrature_ladder,
    )

    rng = np.random.default_rng(22061)
    z = rng.uniform(0.01, 0.40, 350)
    x = (z - z.mean()) / z.std()
    p = 1.0 / (1.0 + np.exp(-(-0.25 + 0.7 * x)))
    retained = rng.binomial(1, p)

    permissive = fit_inclusion_posterior_quadrature_ladder(
        z, retained, degree=1, orders=(5, 7, 9),
        mean_tolerance=1.0, covariance_tolerance=1.0,
        predictive_tolerance=1.0, consecutive_passes=1,
    )
    assert permissive.status == "CONVERGED"
    assert len(permissive.fits) >= 2
    assert len(permissive.mean_max_abs_change) == len(permissive.fits) - 1

    strict = fit_inclusion_posterior_quadrature_ladder(
        z, retained, degree=1, orders=(5, 7, 9),
        mean_tolerance=1e-14, covariance_tolerance=1e-14,
        predictive_tolerance=1e-14, consecutive_passes=2,
    )
    assert strict.status == "NON_CONVERGED"
    assert len(strict.fits) == 3


def test_quadrature_ladder_rejects_invalid_order_schedule():
    from dark_siren_robustness.inclusion_posterior import (
        fit_inclusion_posterior_quadrature_ladder,
    )

    with pytest.raises(ValueError, match="strictly increasing"):
        fit_inclusion_posterior_quadrature_ladder(
            np.linspace(0.01, 0.4, 30), np.arange(30) % 2,
            orders=(7, 7, 11),
        )


def test_log_target_gradient_and_hessian_match_finite_differences():
    from dark_siren_robustness.inclusion_posterior import (
        _log_target, _prior_vectors,
    )

    rng = np.random.default_rng(7712)
    z = rng.uniform(0.01, 0.40, 70)
    x = (z - z.mean()) / z.std()
    design = np.column_stack([np.ones(z.size), x, x*x])
    y = rng.binomial(1, 1.0 / (1.0 + np.exp(-(-0.3 + 0.5*x - 0.2*x*x))))
    mu, sd = _prior_vectors(3, None, None)
    precision = 1.0 / sd**2
    theta = np.asarray([0.2, -0.4, 0.15])
    eps = 2e-5

    def objective(b):
        return -_log_target(b, design, y, mu, sd)

    def gradient(b):
        from scipy.special import expit
        return design.T @ (expit(design @ b) - y) + (b-mu)*precision

    p = 3
    numerical_gradient = np.empty(p)
    numerical_hessian = np.empty((p, p))
    for j in range(p):
        step = np.zeros(p)
        step[j] = eps
        numerical_gradient[j] = (objective(theta+step)-objective(theta-step))/(2*eps)
        numerical_hessian[:, j] = (gradient(theta+step)-gradient(theta-step))/(2*eps)

    from scipy.special import expit
    probability = expit(design @ theta)
    analytic_hessian = design.T @ (
        (probability*(1-probability))[:, None] * design
    ) + np.diag(precision)
    assert np.allclose(numerical_gradient, gradient(theta), rtol=2e-6, atol=2e-6)
    assert np.allclose(numerical_hessian, analytic_hessian, rtol=2e-6, atol=2e-6)


def test_shared_node_marginalization_constant_likelihood_is_exact():
    log_likelihood = np.full((7, 4), -3.75)
    weights = np.asarray([0.02, 0.08, 0.1, 0.2, 0.25, 0.15, 0.2])
    result = log_marginalize_shared_parameter(log_likelihood, weights)
    assert np.allclose(result, -3.75, rtol=0, atol=1e-14)
