import numpy as np
import pytest
from scipy.special import expit

from dark_siren_robustness.estimated_inclusion import (
    fit_quadratic_logistic_inclusion_model,
    predict_inclusion_probability,
)
from dark_siren_robustness.estimated_inclusion_v2 import (
    fit_logistic_inclusion_model_v2,
    fit_quadratic_logistic_inclusion_model_v2,
    predict_inclusion_probability_v2,
    _logistic_nll_gradient_hessian,
)


def test_v2_linear_solver_recovers_well_conditioned_probability_curve():
    rng = np.random.default_rng(202611051)
    z_min, z_max = 0.01, 0.40
    z = rng.uniform(z_min, z_max, 12000)
    x = (z - (z_min + z_max) / 2) / ((z_max - z_min) / 2)
    y = rng.binomial(1, expit(0.25 - 1.15 * x))
    fit = fit_logistic_inclusion_model_v2(z, y, z_min=z_min, z_max=z_max)
    grid = np.linspace(z_min, z_max, 101)
    expected = expit(0.25 - 1.15 * ((grid - 0.205) / 0.195))
    predicted = predict_inclusion_probability_v2(grid, fit)
    assert fit["model"] == "logistic_inclusion_v2"
    assert fit["fit_status"] == "converged"
    assert fit["projected_gradient_inf_norm"] <= 1e-7
    assert np.mean(np.abs(predicted - expected)) < 0.015


def test_v2_quadratic_solver_recovers_nonlinear_probability_curve():
    rng = np.random.default_rng(202611052)
    z_min, z_max = 0.01, 0.40
    z = rng.uniform(z_min, z_max, 16000)
    x = (z - 0.205) / 0.195
    y = rng.binomial(1, expit(0.15 - 0.8 * x + 1.4 * x**2))
    fit = fit_quadratic_logistic_inclusion_model_v2(z, y, z_min=z_min, z_max=z_max)
    grid = np.linspace(z_min, z_max, 101)
    x_grid = (grid - 0.205) / 0.195
    predicted = predict_inclusion_probability_v2(grid, fit)
    assert fit["model"] == "quadratic_logistic_inclusion_v2"
    assert fit["fit_status"] == "converged"
    assert fit["projected_gradient_inf_norm"] <= 1e-7
    assert np.mean(np.abs(predicted - expit(0.15 - 0.8*x_grid + 1.4*x_grid**2))) < 0.015
    # V1 production quadratic fitter uses the uncentered x**2 basis.
    assert fit["intercept"] == pytest.approx(0.15, abs=0.10)
    assert fit["quadratic"] == pytest.approx(1.4, abs=0.10)


def test_v2_solver_reports_separated_data_as_boundary_limited():
    z_min, z_max = 0.01, 0.40
    z = np.linspace(z_min, z_max, 200)
    y = (z > 0.205).astype(int)
    fit = fit_logistic_inclusion_model_v2(z, y, z_min=z_min, z_max=z_max)
    assert fit["fit_status"] == "boundary_limited"
    assert fit["coefficient_bound_hit"]
    assert abs(fit["slope"]) == pytest.approx(50.0, abs=1e-5)
    assert fit["projected_gradient_inf_norm"] <= 1e-7


@pytest.mark.parametrize(
    "fitter",
    [fit_logistic_inclusion_model_v2, fit_quadratic_logistic_inclusion_model_v2],
)
def test_v2_solver_rejects_rank_deficient_design(fitter):
    z_min, z_max = 0.01, 0.40
    z = np.full(100, 0.205)
    y = np.tile([0, 1], 50)
    with pytest.raises(ValueError, match="rank deficient.*coefficients are not identifiable"):
        fitter(z, y, z_min=z_min, z_max=z_max)


def test_v2_solver_rejects_degenerate_outcome_and_predictor_domain():
    with pytest.raises(ValueError, match="both"):
        fit_logistic_inclusion_model_v2(
            np.linspace(0.01, 0.40, 20), np.ones(20), z_min=0.01, z_max=0.40
        )
    model = {
        "model": "logistic_inclusion_v2", "intercept": 0.0, "slope": 0.0,
        "z_min": 0.01, "z_max": 0.40,
    }
    with pytest.raises(ValueError, match="inside"):
        predict_inclusion_probability_v2([0.5], model)


def test_v2_analytic_gradient_and_hessian_match_finite_differences():
    rng = np.random.default_rng(202611053)
    design = np.column_stack([np.ones(80), rng.normal(size=80), rng.normal(size=80)])
    y = rng.binomial(1, expit(design @ np.array([0.2, -0.4, 0.3])))
    theta = np.array([0.15, -0.25, 0.35])
    objective, gradient, hessian = _logistic_nll_gradient_hessian(design, y, theta)
    eps = 1e-5
    numerical_gradient = np.empty(3)
    numerical_hessian = np.empty((3, 3))
    for j in range(3):
        step = np.zeros(3); step[j] = eps
        plus = _logistic_nll_gradient_hessian(design, y, theta + step)[0]
        minus = _logistic_nll_gradient_hessian(design, y, theta - step)[0]
        numerical_gradient[j] = (plus - minus) / (2 * eps)
        numerical_hessian[:, j] = (
            _logistic_nll_gradient_hessian(design, y, theta + step)[1]
            - _logistic_nll_gradient_hessian(design, y, theta - step)[1]
        ) / (2 * eps)
    assert np.isfinite(objective)
    assert np.allclose(gradient, numerical_gradient, rtol=2e-6, atol=2e-7)
    assert np.allclose(hessian, numerical_hessian, rtol=2e-6, atol=2e-7)
    assert np.allclose(hessian, hessian.T, atol=1e-12)
    assert np.linalg.eigvalsh(hessian).min() >= -1e-10


def test_v2_reports_hessian_rank_and_conditioning_diagnostics():
    z_min, z_max = 0.01, 0.40
    z = np.linspace(z_min, z_max, 300)
    x = (z - 0.205) / 0.195
    y = np.tile([0, 1], 150)
    fit = fit_quadratic_logistic_inclusion_model_v2(z, y, z_min=z_min, z_max=z_max)
    assert fit["hessian_rank"] == 3
    assert fit["parameter_dimension"] == 3
    assert np.isfinite(fit["hessian_condition_number"])
    assert fit["hessian_min_eigenvalue"] > 0


def test_v2_rank_deficiency_is_exposed_in_hessian_diagnostics():
    z_min, z_max = 0.01, 0.40
    z = np.full(100, 0.205)
    y = np.tile([0, 1], 50)
    fit = fit_quadratic_logistic_inclusion_model_v2(z, y, z_min=z_min, z_max=z_max)
    assert fit["hessian_rank"] == 1
    assert fit["parameter_dimension"] == 3
    assert np.isinf(fit["hessian_condition_number"])



def test_v2_quadratic_matches_v1_on_identical_data_and_parameterization():
    rng = np.random.default_rng(202611054)
    z_min, z_max = 0.01, 0.40
    z = rng.uniform(z_min, z_max, 9000)
    x = (z - 0.205) / 0.195
    y = rng.binomial(1, expit(0.2 - 0.9 * x + 1.1 * x**2))
    v1 = fit_quadratic_logistic_inclusion_model(z, y, z_min=z_min, z_max=z_max)
    v2 = fit_quadratic_logistic_inclusion_model_v2(z, y, z_min=z_min, z_max=z_max)
    grid = np.linspace(z_min, z_max, 201)
    p1 = predict_inclusion_probability(grid, v1)
    p2 = predict_inclusion_probability_v2(grid, v2)
    assert v1["model"] == "quadratic_logistic_inclusion_v1"
    assert v2["model"] == "quadratic_logistic_inclusion_v2"
    assert np.allclose(p1, p2, rtol=2e-4, atol=2e-5)
    assert v1["negative_log_likelihood"] == pytest.approx(
        v2["negative_log_likelihood"], rel=1e-9, abs=1e-7
    )


def test_v2_linear_matches_v1_on_identical_data_and_parameterization():
    """V1 L-BFGS-B and V2 damped Newton solve the same bounded linear MLE."""
    from dark_siren_robustness.estimated_inclusion import (
        fit_logistic_inclusion_model,
    )

    rng = np.random.default_rng(202611055)
    z_min, z_max = 0.01, 0.40
    z = rng.uniform(z_min, z_max, 10000)
    x = (z - (z_min + z_max) / 2.0) / ((z_max - z_min) / 2.0)
    y = rng.binomial(1, expit(0.35 - 1.05 * x))

    v1 = fit_logistic_inclusion_model(z, y, z_min=z_min, z_max=z_max)
    v2 = fit_logistic_inclusion_model_v2(z, y, z_min=z_min, z_max=z_max)

    grid = np.linspace(z_min, z_max, 301)
    p1 = predict_inclusion_probability(grid, v1)
    p2 = predict_inclusion_probability_v2(grid, v2)

    assert v1["fit_status"] == "converged"
    assert v2["fit_status"] == "converged"
    assert np.allclose(
        [v1["intercept"], v1["slope"]],
        [v2["intercept"], v2["slope"]],
        rtol=2e-4,
        atol=2e-5,
    )
    assert np.allclose(p1, p2, rtol=2e-4, atol=2e-5)
    assert v1["negative_log_likelihood"] == pytest.approx(
        v2["negative_log_likelihood"], rel=1e-9, abs=1e-7
    )


@pytest.mark.parametrize(
    ("included", "expected_slope_sign"),
    [
        ("increasing", 1),
        ("decreasing", -1),
    ],
)
def test_v2_boundary_solution_satisfies_explicit_box_kkt_gradient_sign(
    included, expected_slope_sign
):
    """At an upper/lower slope face, verify the minimization KKT sign directly."""
    z_min, z_max = 0.01, 0.40
    z = np.linspace(z_min, z_max, 300)
    y = (z > 0.205).astype(int)
    if included == "decreasing":
        y = 1 - y
    fit = fit_logistic_inclusion_model_v2(
        z, y, z_min=z_min, z_max=z_max
    )
    x = (z - (z_min + z_max) / 2.0) / ((z_max - z_min) / 2.0)
    design = np.column_stack([np.ones(z.size), x])
    theta = np.array([fit["intercept"], fit["slope"]])
    _objective, gradient, _hessian = _logistic_nll_gradient_hessian(
        design, y.astype(float), theta
    )
    assert abs(fit["slope"]) == pytest.approx(50.0, abs=1e-5)
    assert np.sign(fit["slope"]) == expected_slope_sign
    if expected_slope_sign > 0:
        # Upper bound: gradient must be nonpositive for box-constrained
        # minimization under the documented KKT convention.
        assert gradient[1] <= 1e-6
    else:
        # Lower bound: gradient must be nonnegative.
        assert gradient[1] >= -1e-6
    assert fit["projected_gradient_inf_norm"] <= 1e-7


def test_v2_extreme_finite_coefficients_expose_floating_point_saturation():
    model = {
        "model": "logistic_inclusion_v2",
        "intercept": 25.0,
        "slope": 50.0,
        "z_min": 0.01,
        "z_max": 0.40,
    }
    grid = np.array([0.01, 0.205, 0.40])
    probability = predict_inclusion_probability_v2(grid, model)
    assert np.all(np.isfinite(probability))
    assert probability[0] > 0.0
    assert probability[-1] == 1.0  # expit(75) rounds to exactly one in float64
    # At x=-1 the linear predictor is 25 - 50 = -25, not -75;
    # the coefficient box therefore gives expit(-25) at this endpoint.
    assert probability[0] == pytest.approx(expit(-25.0), rel=1e-12)
