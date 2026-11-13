"""Direct two-dimensional integration reference for the proper logistic posterior."""
import numpy as np
from scipy.integrate import trapezoid
from scipy.special import expit

from dark_siren_robustness.inclusion_posterior import (
    fit_inclusion_posterior_quadrature,
    posterior_predictive_inclusion,
)


def test_laplace_hermite_posterior_moments_match_direct_2d_integration():
    rng = np.random.default_rng(20261113)
    z = np.linspace(0.012, 0.398, 90)
    center, scale = 0.205, 0.195
    x = (z - center) / scale
    true_probability = expit(-0.35 - 1.1 * x)
    retained = (rng.random(z.size) < true_probability).astype(int)

    fit = fit_inclusion_posterior_quadrature(
        z, retained, degree=1, coordinate_center=center,
        coordinate_scale=scale, order=19,
    )

    # Integrate the declared likelihood x Gaussian prior directly on a dense
    # coefficient grid. This does not use the Laplace proposal or GH nodes.
    sd = np.array([5.0, 2.5])
    half_width = 8.0 * np.sqrt(np.diag(fit.laplace_covariance))
    beta0 = np.linspace(fit.mode[0] - half_width[0], fit.mode[0] + half_width[0], 801)
    beta1 = np.linspace(fit.mode[1] - half_width[1], fit.mode[1] + half_width[1], 801)
    b0, b1 = np.meshgrid(beta0, beta1, indexing="ij")
    log_likelihood = np.zeros_like(b0)
    for xi, yi in zip(x, retained):
        eta_i = b0 + b1 * xi
        log_likelihood += yi * eta_i - np.logaddexp(0.0, eta_i)
    log_prior = -0.5 * ((b0 / sd[0]) ** 2 + (b1 / sd[1]) ** 2)
    log_target = log_likelihood + log_prior
    mass = np.exp(log_target - np.max(log_target))
    normalization = trapezoid(trapezoid(mass, beta1, axis=1), beta0)
    density = mass / normalization

    marginal0 = trapezoid(density, beta1, axis=1)
    marginal1 = trapezoid(density, beta0, axis=0)
    mean0 = trapezoid(beta0 * marginal0, beta0)
    mean1 = trapezoid(beta1 * marginal1, beta1)
    second00 = trapezoid(beta0**2 * marginal0, beta0)
    second11 = trapezoid(beta1**2 * marginal1, beta1)
    second01 = trapezoid(trapezoid(b0 * b1 * density, beta1, axis=1), beta0)
    reference_mean = np.array([mean0, mean1])
    reference_cov = np.array([
        [second00 - mean0**2, second01 - mean0 * mean1],
        [second01 - mean0 * mean1, second11 - mean1**2],
    ])

    np.testing.assert_allclose(fit.posterior_mean, reference_mean, rtol=2e-3, atol=2e-3)
    np.testing.assert_allclose(
        fit.posterior_covariance, reference_cov, rtol=5e-3, atol=2e-3
    )

    prediction_z = np.array([0.05, 0.20, 0.35])
    prediction_x = (prediction_z - center) / scale
    direct_predictive = []
    for xp in prediction_x:
        predictive_integrand = expit(b0 + b1 * xp) * density
        direct_predictive.append(
            trapezoid(trapezoid(predictive_integrand, beta1, axis=1), beta0)
        )
    np.testing.assert_allclose(
        posterior_predictive_inclusion(fit, prediction_z),
        np.asarray(direct_predictive),
        rtol=3e-3, atol=2e-3,
    )
