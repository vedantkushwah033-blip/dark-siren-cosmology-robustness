import numpy as np
import pytest

from dark_siren_robustness.inclusion_posterior import (
    fit_inclusion_posterior_quadrature,
)
from dark_siren_robustness.independent_posterior_reference import (
    reconstruct_logistic_posterior,
)


@pytest.mark.parametrize(
    ("degree", "truth", "prior_sd", "all_retained"),
    [
        (1, np.array([-0.4, 0.9]), np.array([5.0, 2.5]), False),
        (2, np.array([-0.2, 0.7, -0.45]), np.array([5.0, 2.5, 2.5]), False),
        (1, np.array([1.0, 0.0]), np.array([5.0, 2.5]), True),
    ],
)
def test_reference_posterior_reconstructs_production_rule(
    degree, truth, prior_sd, all_retained
):
    rng = np.random.default_rng(78123 + degree + int(all_retained))
    redshift = rng.uniform(0.01, 0.40, 180)
    center, scale = 0.205, 0.195
    x = (redshift - center) / scale
    design = np.column_stack([x**power for power in range(degree + 1)])
    probability = 1.0 / (1.0 + np.exp(-(design @ truth)))
    retained = (
        np.ones(redshift.size, dtype=int)
        if all_retained
        else rng.binomial(1, probability)
    )
    prior_mean = np.zeros(degree + 1)

    production = fit_inclusion_posterior_quadrature(
        redshift,
        retained,
        degree=degree,
        prior_mean=prior_mean,
        prior_sd=prior_sd,
        coordinate_center=center,
        coordinate_scale=scale,
        order=9,
    )
    reference = reconstruct_logistic_posterior(
        redshift,
        retained,
        degree=degree,
        prior_mean=prior_mean,
        prior_sd=prior_sd,
        center=center,
        scale=scale,
        order=9,
    )

    assert reference.mode_score_inf < 2e-7
    np.testing.assert_allclose(reference.mode, production.mode, rtol=2e-6, atol=2e-7)
    np.testing.assert_allclose(reference.nodes, production.nodes, rtol=2e-6, atol=2e-7)
    np.testing.assert_allclose(reference.weights, production.weights, rtol=2e-5, atol=2e-8)
    np.testing.assert_allclose(
        reference.posterior_mean, production.posterior_mean, rtol=2e-5, atol=2e-7
    )
    np.testing.assert_allclose(
        reference.posterior_covariance,
        production.posterior_covariance,
        rtol=3e-5,
        atol=3e-7,
    )


def test_reference_posterior_matches_independent_tensor_legendre_integral():
    """Compare the Laplace-Hermite rule with a separate prior-coordinate integral."""
    from scipy.special import expit, logsumexp, roots_legendre

    rng = np.random.default_rng(914207)
    redshift = rng.uniform(0.01, 0.40, 95)
    center, scale = 0.205, 0.195
    x = (redshift - center) / scale
    design = np.column_stack((np.ones_like(x), x))
    truth = np.array([-0.35, 0.85])
    retained = rng.binomial(1, expit(design @ truth))
    prior_mean = np.array([0.0, 0.0])
    prior_sd = np.array([3.0, 1.8])

    reference = reconstruct_logistic_posterior(
        redshift,
        retained,
        degree=1,
        prior_mean=prior_mean,
        prior_sd=prior_sd,
        center=center,
        scale=scale,
        order=31,
    )

    # Independent tensor Gauss-Legendre integration in standardized prior
    # coordinates. This does not use the Laplace proposal or Hermite nodes.
    abscissa, weights_1d = roots_legendre(260)
    standardized = 8.0 * abscissa
    quadrature_weights = 8.0 * weights_1d
    u0, u1 = np.meshgrid(standardized, standardized, indexing="ij")
    beta0 = prior_mean[0] + prior_sd[0] * u0
    beta1 = prior_mean[1] + prior_sd[1] * u1
    eta = (
        design[:, 0, None, None] * beta0[None, :, :]
        + design[:, 1, None, None] * beta1[None, :, :]
    )
    log_likelihood = np.sum(
        retained[:, None, None] * eta - np.logaddexp(0.0, eta),
        axis=0,
    )
    # Standard-normal prior density in u coordinates; the Jacobian cancels
    # the coefficient-scale factors exactly.
    log_integrand = (
        log_likelihood
        - 0.5 * (u0 * u0 + u1 * u1)
        - np.log(2.0 * np.pi)
    )
    log_measure = (
        np.log(quadrature_weights[:, None])
        + np.log(quadrature_weights[None, :])
    )
    log_mass = log_integrand + log_measure
    log_evidence = float(logsumexp(log_mass))
    normalized_mass = np.exp(log_mass - log_evidence)

    direct_mean = np.array([
        np.sum(normalized_mass * beta0),
        np.sum(normalized_mass * beta1),
    ])
    centered0 = beta0 - direct_mean[0]
    centered1 = beta1 - direct_mean[1]
    direct_covariance = np.array([
        [
            np.sum(normalized_mass * centered0 * centered0),
            np.sum(normalized_mass * centered0 * centered1),
        ],
        [
            np.sum(normalized_mass * centered1 * centered0),
            np.sum(normalized_mass * centered1 * centered1),
        ],
    ])

    np.testing.assert_allclose(reference.log_normalizer, log_evidence, rtol=0, atol=3e-5)
    np.testing.assert_allclose(reference.posterior_mean, direct_mean, rtol=0, atol=3e-4)
    np.testing.assert_allclose(
        reference.posterior_covariance, direct_covariance, rtol=0, atol=5e-4
    )
