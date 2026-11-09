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
