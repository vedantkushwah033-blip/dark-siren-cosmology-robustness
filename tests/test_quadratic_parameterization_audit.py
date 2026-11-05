"""Deterministic algebra checks for quadratic inclusion parameterizations.

These tests validate identities and feasible-set mappings only. They do not
claim optimizer equivalence or empirical downstream robustness.
"""
import numpy as np
import pytest
from scipy.special import expit


@pytest.mark.parametrize(
    "alpha,beta,gamma,center",
    [
        (0.2, -1.1, 0.0, 0.2845),
        (1.7, 0.4, -2.0, 0.2845),
        (-3.0, 2.2, 1.5, 0.0),
        (0.0, 0.0, 4.0, 0.731),
    ],
)
def test_centered_and_uncentered_quadratic_predictors_are_identical(
    alpha, beta, gamma, center
):
    x = np.linspace(-1.0, 1.0, 1001)

    centered_eta = alpha + beta * x + gamma * (x**2 - center)
    uncentered_intercept = alpha - gamma * center
    uncentered_eta = (
        uncentered_intercept + beta * x + gamma * x**2
    )

    np.testing.assert_allclose(
        centered_eta, uncentered_eta, rtol=0.0, atol=2e-15
    )
    np.testing.assert_allclose(
        expit(centered_eta), expit(uncentered_eta),
        rtol=0.0, atol=2e-16,
    )


@pytest.mark.parametrize(
    "quadratic,center",
    [(-50.0, 0.2845), (-2.0, 0.2845), (0.0, 0.2845),
     (1.5, 0.2845), (50.0, 0.2845), (3.0, 0.731)],
)
def test_transformed_centered_intercept_bounds_map_to_uncentered_box(
    quadratic, center
):
    uncentered_lower, uncentered_upper = -25.0, 25.0
    centered_lower = uncentered_lower + quadratic * center
    centered_upper = uncentered_upper + quadratic * center

    # The coefficient transformation is a = alpha - q*c.
    mapped_lower = centered_lower - quadratic * center
    mapped_upper = centered_upper - quadratic * center

    assert mapped_lower == pytest.approx(uncentered_lower, abs=1e-14)
    assert mapped_upper == pytest.approx(uncentered_upper, abs=1e-14)


def test_same_numeric_centered_and_uncentered_intercept_boxes_are_not_equivalent():
    center = 0.28450056
    quadratic = 2.0

    # A centered intercept alpha in [-25, 25] maps to a in
    # [-25 - q*c, 25 - q*c], not generally [-25, 25].
    mapped_lower = -25.0 - quadratic * center
    mapped_upper = 25.0 - quadratic * center

    assert mapped_lower != pytest.approx(-25.0)
    assert mapped_upper != pytest.approx(25.0)


def test_zero_center_is_the_special_case_where_intercept_boxes_match():
    center = 0.0
    quadratic = 7.5
    assert -25.0 - quadratic * center == -25.0
    assert 25.0 - quadratic * center == 25.0
