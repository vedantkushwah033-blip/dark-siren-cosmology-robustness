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
