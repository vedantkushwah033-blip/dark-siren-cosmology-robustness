"""Deterministic algebra and constrained-fit checks for quadratic inclusion bases.

These tests validate basis invariance under equivalent feasible bounds. They
do not claim robustness to a wider coefficient box or downstream H0 validity.
"""
import numpy as np
import pytest
from scipy.optimize import minimize
from scipy.special import expit, log_expit


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
    uncentered_eta = uncentered_intercept + beta * x + gamma * x**2
    np.testing.assert_allclose(centered_eta, uncentered_eta, rtol=0.0, atol=2e-15)
    np.testing.assert_allclose(
        expit(centered_eta), expit(uncentered_eta), rtol=0.0, atol=2e-16
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
    mapped_lower = centered_lower - quadratic * center
    mapped_upper = centered_upper - quadratic * center
    assert mapped_lower == pytest.approx(uncentered_lower, abs=1e-14)
    assert mapped_upper == pytest.approx(uncentered_upper, abs=1e-14)


def test_same_numeric_centered_and_uncentered_intercept_boxes_are_not_equivalent():
    center = 0.28450056
    quadratic = 2.0
    mapped_lower = -25.0 - quadratic * center
    mapped_upper = 25.0 - quadratic * center
    assert mapped_lower != pytest.approx(-25.0)
    assert mapped_upper != pytest.approx(25.0)


def test_zero_center_is_the_special_case_where_intercept_boxes_match():
    center = 0.0
    quadratic = 7.5
    assert -25.0 - quadratic * center == -25.0
    assert 25.0 - quadratic * center == 25.0


def test_centered_and_uncentered_constrained_mles_agree_under_transformed_bounds():
    """Compare distinct optimizer coordinates over the same feasible predictor set."""
    rng = np.random.default_rng(7319)
    x = np.linspace(-1.0, 1.0, 201)
    center = 0.28450056
    probability = expit(0.3 - 0.8 * x + 1.1 * x**2)
    included = rng.binomial(1, probability)
    assert np.unique(included).size == 2

    uncentered_design = np.column_stack((np.ones_like(x), x, x**2))
    centered_design = np.column_stack((np.ones_like(x), x, x**2 - center))

    def objective(design, theta):
        eta = design @ theta
        return float(-np.sum(
            included * log_expit(eta) + (1.0 - included) * log_expit(-eta)
        ))

    def gradient(design, theta):
        return design.T @ (expit(design @ theta) - included)

    uncentered = minimize(
        lambda theta: objective(uncentered_design, theta),
        np.zeros(3),
        jac=lambda theta: gradient(uncentered_design, theta),
        method="L-BFGS-B",
        bounds=((-25.0, 25.0), (-50.0, 50.0), (-50.0, 50.0)),
        options={"maxiter": 10000, "ftol": 1e-15, "gtol": 1e-9},
    )

    # In centered coordinates a = alpha - gamma*c must remain in [-25, 25].
    # This is a pair of linear inequalities, not a fixed alpha box.
    constraints = [
        {
            "type": "ineq",
            "fun": lambda theta: 25.0 - (theta[0] - center * theta[2]),
            "jac": lambda theta: np.array([-1.0, 0.0, center]),
        },
        {
            "type": "ineq",
            "fun": lambda theta: 25.0 + (theta[0] - center * theta[2]),
            "jac": lambda theta: np.array([1.0, 0.0, -center]),
        },
    ]
    centered = minimize(
        lambda theta: objective(centered_design, theta),
        np.zeros(3),
        jac=lambda theta: gradient(centered_design, theta),
        method="SLSQP",
        bounds=((-100.0, 100.0), (-50.0, 50.0), (-50.0, 50.0)),
        constraints=constraints,
        options={"maxiter": 10000, "ftol": 1e-12},
    )

    assert uncentered.success, uncentered.message
    assert centered.success, centered.message
    assert uncentered.fun == pytest.approx(centered.fun, rel=0.0, abs=1e-8)

    # Transform centered coefficients back to the uncentered basis.
    alpha, beta, gamma = centered.x
    transformed = np.array([alpha - gamma * center, beta, gamma])
    np.testing.assert_allclose(transformed, uncentered.x, rtol=0.0, atol=2e-5)

    eta_u = uncentered_design @ uncentered.x
    eta_c = centered_design @ centered.x
    np.testing.assert_allclose(eta_c, eta_u, rtol=0.0, atol=2e-8)
    np.testing.assert_allclose(expit(eta_c), expit(eta_u), rtol=0.0, atol=2e-8)
    assert np.max(np.abs(gradient(uncentered_design, uncentered.x))) < 1e-6
    assert np.max(np.abs(gradient(centered_design, centered.x))) < 1e-5


def test_centered_and_uncentered_fit_agree_when_transformed_intercept_bound_is_active():
    """Exercise the q-dependent intercept constraint at an actual optimum."""
    x = np.linspace(-1.0, 1.0, 201)
    center = 0.28450056
    fixed_quadratic = 1.2
    included = np.ones(x.size, dtype=float)

    uncentered_design = np.column_stack((np.ones_like(x), x))
    centered_design = uncentered_design.copy()
    uncentered_offset = fixed_quadratic * x**2
    centered_offset = fixed_quadratic * (x**2 - center)

    def objective(design, offset, theta):
        eta = design @ theta + offset
        return float(-np.sum(log_expit(eta)))

    def gradient(design, offset, theta):
        eta = design @ theta + offset
        residual = expit(eta) - included
        return design.T @ residual

    uncentered = minimize(
        lambda theta: objective(uncentered_design, uncentered_offset, theta),
        np.zeros(2),
        jac=lambda theta: gradient(uncentered_design, uncentered_offset, theta),
        method="L-BFGS-B",
        bounds=((-25.0, 25.0), (-50.0, 50.0)),
        options={"maxiter": 10000, "ftol": 1e-15, "gtol": 1e-9},
    )

    # a = alpha - q*c, so the equivalent centered intercept interval
    # is [-25 + q*c, 25 + q*c], not the same fixed numeric box.
    centered_intercept_bounds = (
        -25.0 + fixed_quadratic * center,
        25.0 + fixed_quadratic * center,
    )
    centered = minimize(
        lambda theta: objective(centered_design, centered_offset, theta),
        np.zeros(2),
        jac=lambda theta: gradient(centered_design, centered_offset, theta),
        method="L-BFGS-B",
        bounds=(centered_intercept_bounds, (-50.0, 50.0)),
        options={"maxiter": 10000, "ftol": 1e-15, "gtol": 1e-9},
    )

    assert uncentered.success, uncentered.message
    assert centered.success, centered.message
    assert uncentered.x[0] == pytest.approx(25.0, abs=1e-7)
    assert centered.x[0] == pytest.approx(centered_intercept_bounds[1], abs=1e-7)
    assert uncentered.fun == pytest.approx(centered.fun, rel=0.0, abs=1e-10)

    eta_uncentered = uncentered_design @ uncentered.x + uncentered_offset
    eta_centered = centered_design @ centered.x + centered_offset
    np.testing.assert_allclose(eta_centered, eta_uncentered, rtol=0.0, atol=1e-9)
