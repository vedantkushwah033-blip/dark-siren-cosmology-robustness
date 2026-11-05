"""Cross-check V2 against an algorithmically independent optimizer.

This float64 SLSQP comparison is supplementary coverage, not the
high-precision active-set oracle required to close historical D-01.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.optimize import minimize
from scipy.special import expit

from dark_siren_robustness.estimated_inclusion_v2 import (
    fit_logistic_inclusion_model_v2,
    fit_quadratic_logistic_inclusion_model_v2,
)


def _slsqp_reference(design: np.ndarray, outcome: np.ndarray, bounds):
    """Independent bounded MLE using SLSQP and a logaddexp objective."""
    def objective(theta):
        eta = design @ theta
        return float(np.sum(np.logaddexp(0.0, eta) - outcome * eta))

    def gradient(theta):
        return design.T @ (expit(design @ theta) - outcome)

    result = minimize(
        objective, np.zeros(design.shape[1]), jac=gradient,
        method="SLSQP", bounds=bounds,
        options={"maxiter": 3000, "ftol": 1e-12, "disp": False},
    )
    assert np.all(np.isfinite(result.x)), result.message
    assert np.isfinite(result.fun), result.message
    raw = gradient(result.x)
    projected = raw.copy()
    for j, (lower, upper) in enumerate(bounds):
        if result.x[j] <= lower + 2e-6 and raw[j] > 0:
            projected[j] = 0.0
        elif result.x[j] >= upper - 2e-6 and raw[j] < 0:
            projected[j] = 0.0
    assert np.max(np.abs(projected)) <= 2e-4, (
        result.message, np.max(np.abs(projected))
    )
    return result.x, result.fun


@pytest.mark.parametrize("degree", [1, 2])
@pytest.mark.parametrize("seed", [61101, 61102, 61103])
def test_v2_matches_independent_slsqp_reference(degree: int, seed: int):
    rng = np.random.default_rng(seed)
    z_min, z_max = 0.01, 0.40
    midpoint = (z_min + z_max) / 2.0
    half_range = (z_max - z_min) / 2.0
    z = rng.uniform(z_min, z_max, 5000)
    x = (z - midpoint) / half_range

    if degree == 1:
        truth = np.array([0.35, -1.10])
        design = np.column_stack([np.ones(z.size), x])
        bounds = ((-25.0, 25.0), (-50.0, 50.0))
    else:
        truth = np.array([0.20, -0.85, 1.05])
        design = np.column_stack([np.ones(z.size), x, x * x])
        bounds = ((-25.0, 25.0), (-50.0, 50.0), (-50.0, 50.0))

    outcome = rng.binomial(1, expit(design @ truth)).astype(float)
    if degree == 1:
        fit = fit_logistic_inclusion_model_v2(
            z, outcome, z_min=z_min, z_max=z_max
        )
        theta_v2 = np.array([fit["intercept"], fit["slope"]])
    else:
        fit = fit_quadratic_logistic_inclusion_model_v2(
            z, outcome, z_min=z_min, z_max=z_max
        )
        theta_v2 = np.array([fit["intercept"], fit["slope"], fit["quadratic"]])

    theta_ref, objective_ref = _slsqp_reference(design, outcome, bounds)
    eta_v2 = design @ theta_v2
    objective_v2 = float(np.sum(np.logaddexp(0.0, eta_v2) - outcome * eta_v2))
    assert fit["fit_status"] == "converged"
    assert fit["projected_gradient_inf_norm"] <= 1e-7
    assert objective_v2 <= objective_ref + 2e-7
    assert np.allclose(theta_v2, theta_ref, rtol=2e-5, atol=2e-5)


@pytest.mark.parametrize("degree", [1, 2])
@pytest.mark.parametrize("direction", ["increasing", "decreasing"])
def test_v2_separated_fit_satisfies_independent_box_kkt_conditions(
    degree: int, direction: str
):
    """Stress the constrained optimum and independently check every KKT face."""
    z_min, z_max = 0.01, 0.40
    z = np.linspace(z_min, z_max, 401)
    x = (z - (z_min + z_max) / 2.0) / ((z_max - z_min) / 2.0)
    outcome = (x > 0.12).astype(float)
    if direction == "decreasing":
        outcome = 1.0 - outcome

    if degree == 1:
        design = np.column_stack([np.ones(z.size), x])
        bounds = ((-25.0, 25.0), (-50.0, 50.0))
        fit = fit_logistic_inclusion_model_v2(
            z, outcome, z_min=z_min, z_max=z_max
        )
        theta = np.array([fit["intercept"], fit["slope"]])
    else:
        design = np.column_stack([np.ones(z.size), x, x * x])
        bounds = ((-25.0, 25.0), (-50.0, 50.0), (-50.0, 50.0))
        fit = fit_quadratic_logistic_inclusion_model_v2(
            z, outcome, z_min=z_min, z_max=z_max
        )
        theta = np.array([fit["intercept"], fit["slope"], fit["quadratic"]])

    # Recompute objective derivatives without calling the production solver's
    # private likelihood/gradient/Hessian implementation.
    eta = design @ theta
    probability = expit(eta)
    independent_gradient = design.T @ (probability - outcome)
    independent_hessian = design.T @ (
        design * (probability * (1.0 - probability))[:, None]
    )
    assert np.linalg.eigvalsh(independent_hessian).min() >= -1e-8

    # Box-constrained minimization KKT signs: gradient >= 0 at a lower
    # face, <= 0 at an upper face, and zero in the interior.
    kkt_components = []
    for coefficient, gradient, (lower, upper) in zip(
        theta, independent_gradient, bounds
    ):
        if coefficient <= lower + 1e-5:
            kkt_components.append(max(0.0, -float(gradient)))
        elif coefficient >= upper - 1e-5:
            kkt_components.append(max(0.0, float(gradient)))
        else:
            kkt_components.append(abs(float(gradient)))

    assert fit["coefficient_bound_hit"]
    assert fit["fit_status"] == "boundary_limited"
    assert max(kkt_components) <= 2e-7
    assert fit["projected_gradient_inf_norm"] <= 1e-7
