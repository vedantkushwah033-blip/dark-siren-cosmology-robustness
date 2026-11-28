import numpy as np
import pytest
from scipy.special import expit

from dark_siren_robustness.design_based_logistic import (
    srswor_logistic_linearization_covariance,
)


def test_srswor_logistic_covariance_matches_explicit_linearization():
    x1 = np.linspace(-1.0, 1.0, 80)
    design = np.column_stack([np.ones_like(x1), x1])
    beta = np.array([-0.3, 0.8])
    probability = expit(design @ beta)
    # Deterministic, non-degenerate finite validation sample.
    retained = (np.arange(x1.size) % 3 != 0).astype(int)

    result = srswor_logistic_linearization_covariance(
        design, retained, beta, population_size=400
    )

    score = design * (retained - probability)[:, None]
    score_centered = score - score.mean(axis=0)
    score_cov = score_centered.T @ score_centered / (len(retained) - 1)
    sensitivity = (design.T * (probability * (1 - probability))) @ design / len(retained)
    fpc = 1 - len(retained) / 400
    expected = np.linalg.inv(sensitivity) @ (fpc / len(retained) * score_cov) @ np.linalg.inv(sensitivity).T

    np.testing.assert_allclose(result.sensitivity, sensitivity, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(result.score_covariance, score_cov, rtol=1e-13, atol=1e-13)
    np.testing.assert_allclose(result.covariance, expected, rtol=1e-12, atol=1e-12)
    assert result.fpc == pytest.approx(fpc)
    assert np.all(np.linalg.eigvalsh(result.covariance) >= -1e-12)


def test_srswor_logistic_covariance_is_zero_for_census():
    x1 = np.linspace(-1.0, 1.0, 40)
    design = np.column_stack([np.ones_like(x1), x1])
    beta = np.array([0.2, -0.4])
    retained = (np.arange(x1.size) % 2).astype(int)

    result = srswor_logistic_linearization_covariance(
        design, retained, beta, population_size=len(retained)
    )
    assert result.fpc == 0.0
    np.testing.assert_allclose(result.covariance, np.zeros((2, 2)), atol=1e-14)


@pytest.mark.parametrize(
    "design,retained,beta,population_size",
    [
        (np.ones((5, 2)), np.array([0, 1, 0, 1, 0]), np.zeros(2), 10),
        (np.ones((5, 1)), np.array([0, 1, 0, 1, 0]), np.zeros(1), 4),
        (np.ones((5, 1)), np.array([0, 1, 0, 1, 0]), np.zeros(2), 10),
    ],
)
def test_srswor_logistic_covariance_rejects_invalid_inputs(
    design, retained, beta, population_size
):
    with pytest.raises(ValueError):
        srswor_logistic_linearization_covariance(
            design, retained, beta, population_size=population_size
        )


def test_srswor_linearization_tracks_repeated_fixed_frame_sampling():
    from scipy.optimize import minimize

    rng = np.random.default_rng(20261128)
    population_size = 1200
    sample_size = 300
    replicates = 500
    coordinate = np.linspace(-1.0, 1.0, population_size)
    full_design = np.column_stack([np.ones(population_size), coordinate])
    generating_beta = np.array([-0.2, 0.8])
    retention_probability = expit(full_design @ generating_beta)
    finite_retention = (
        rng.random(population_size) < retention_probability
    ).astype(int)

    def fit_logistic(indices):
        design = full_design[indices]
        outcome = finite_retention[indices]

        def objective(beta):
            eta = design @ beta
            return float(np.sum(np.logaddexp(0.0, eta) - outcome * eta))

        def gradient(beta):
            return design.T @ (expit(design @ beta) - outcome)

        fit = minimize(
            objective, np.zeros(2), jac=gradient, method="BFGS",
            options={"gtol": 1e-9, "maxiter": 500},
        )
        assert np.all(np.isfinite(fit.x))
        assert np.linalg.norm(gradient(fit.x), ord=np.inf) < 1e-6
        return np.asarray(fit.x)

    # The estimand is the logistic score root for this one realized finite
    # frame, not the coefficient used to generate the Bernoulli labels.
    finite_frame_beta = fit_logistic(np.arange(population_size))
    estimates = []
    estimated_design_variances = []
    for _ in range(replicates):
        selected = rng.choice(population_size, size=sample_size, replace=False)
        estimate = fit_logistic(selected)
        estimates.append(estimate)
        estimated_design_variances.append(
            srswor_logistic_linearization_covariance(
                full_design[selected], finite_retention[selected], estimate,
                population_size=population_size,
            ).covariance
        )

    empirical = np.cov(np.asarray(estimates).T, ddof=1)
    mean_linearized = np.mean(estimated_design_variances, axis=0)
    # A moderate tolerance is appropriate for a first-order approximation
    # at n=300; this is a fixed-frame repeated-sampling check, not a theorem.
    for index in range(2):
        assert mean_linearized[index, index] == pytest.approx(
            empirical[index, index], rel=0.20
        )
    assert np.all(np.isfinite(finite_frame_beta))
