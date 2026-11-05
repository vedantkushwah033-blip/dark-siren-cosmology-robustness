import numpy as np
import pytest
from scipy.special import expit

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


def test_v2_solver_handles_rank_deficient_quadratic_design_deterministically():
    z_min, z_max = 0.01, 0.40
    z = np.full(100, 0.205)
    y = np.tile([0, 1], 50)
    first = fit_quadratic_logistic_inclusion_model_v2(z, y, z_min=z_min, z_max=z_max)
    second = fit_quadratic_logistic_inclusion_model_v2(z, y, z_min=z_min, z_max=z_max)
    assert first["fit_status"] == "converged"
    assert np.array_equal(
        [first["intercept"], first["slope"], first["quadratic"]],
        [second["intercept"], second["slope"], second["quadratic"]],
    )
    assert first["projected_gradient_inf_norm"] <= 1e-7


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
