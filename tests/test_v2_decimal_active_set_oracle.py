"""High-precision active-set oracle tests for bounded logistic V2.

The Decimal reference is deliberately test-owned: it does not import the
production objective, derivatives, or optimizer. It enumerates every box face,
solves score equations on free coordinates, and independently checks KKT signs.
"""
from __future__ import annotations

from decimal import Decimal, localcontext
from itertools import product

import numpy as np
import pytest
from scipy.special import expit

from dark_siren_robustness.estimated_inclusion_v2 import (
    fit_logistic_inclusion_model_v2,
    fit_quadratic_logistic_inclusion_model_v2,
)

D = Decimal


def _softplus(value: Decimal) -> Decimal:
    if value > 0:
        return value + (D(1) + (-value).exp()).ln()
    return (D(1) + value.exp()).ln()


def _probability(value: Decimal) -> Decimal:
    if value >= 0:
        return D(1) / (D(1) + (-value).exp())
    exp_value = value.exp()
    return exp_value / (D(1) + exp_value)


def _decimal_objective_gradient_hessian(design, outcome, theta):
    dimension = len(theta)
    objective = D(0)
    gradient = [D(0) for _ in range(dimension)]
    hessian = [[D(0) for _ in range(dimension)] for _ in range(dimension)]
    for row, observed in zip(design, outcome):
        eta = sum(value * coefficient for value, coefficient in zip(row, theta))
        probability = _probability(eta)
        residual = probability - observed
        variance = probability * (D(1) - probability)
        objective += _softplus(eta) - observed * eta
        for j in range(dimension):
            gradient[j] += row[j] * residual
            for k in range(dimension):
                hessian[j][k] += row[j] * row[k] * variance
    return objective, gradient, hessian


def _decimal_solve(matrix, rhs):
    """Pivoted Gaussian elimination in the active-set reference arithmetic."""
    n = len(rhs)
    augmented = [list(matrix[i]) + [rhs[i]] for i in range(n)]
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < D("1e-65"):
            raise ArithmeticError("singular Decimal Hessian")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [value / divisor for value in augmented[column]]
        for row in range(n):
            if row == column:
                continue
            multiplier = augmented[row][column]
            if multiplier:
                augmented[row] = [
                    augmented[row][j] - multiplier * augmented[column][j]
                    for j in range(n + 1)
                ]
    return [augmented[i][-1] for i in range(n)]


def _decimal_active_set_oracle(design_float, outcome_float, bounds_float):
    """Enumerate lower/free/upper faces and retain feasible KKT candidates."""
    with localcontext() as context:
        context.prec = 70
        design = [[D(str(float(value))) for value in row] for row in design_float]
        outcome = [D(str(float(value))) for value in outcome_float]
        bounds = [(D(str(lo)), D(str(hi))) for lo, hi in bounds_float]
        dimension = len(bounds)
        tolerance = D("1e-35")
        candidates = []

        for face in product((-1, 0, 1), repeat=dimension):
            theta = [D(0) for _ in range(dimension)]
            free = []
            for j, state in enumerate(face):
                if state == -1:
                    theta[j] = bounds[j][0]
                elif state == 1:
                    theta[j] = bounds[j][1]
                else:
                    free.append(j)

            converged = not free
            for _ in range(120):
                objective, gradient, hessian = _decimal_objective_gradient_hessian(
                    design, outcome, theta
                )
                if not free or max(abs(gradient[j]) for j in free) <= tolerance:
                    converged = True
                    break
                reduced_hessian = [[hessian[j][k] for k in free] for j in free]
                reduced_gradient = [gradient[j] for j in free]
                try:
                    step = _decimal_solve(reduced_hessian, reduced_gradient)
                except ArithmeticError:
                    break

                accepted = False
                scale = D(1)
                while scale >= D("1e-35"):
                    candidate = theta.copy()
                    for j, delta in zip(free, step):
                        candidate[j] = theta[j] - scale * delta
                    inside = all(
                        bounds[j][0] <= candidate[j] <= bounds[j][1]
                        for j in free
                    )
                    if inside:
                        candidate_objective = _decimal_objective_gradient_hessian(
                            design, outcome, candidate
                        )[0]
                        if candidate_objective < objective:
                            theta = candidate
                            accepted = True
                            break
                    scale /= 2
                if not accepted:
                    break

            if not converged:
                continue
            objective, gradient, _ = _decimal_objective_gradient_hessian(
                design, outcome, theta
            )
            kkt_ok = True
            for j, state in enumerate(face):
                if state == -1 and gradient[j] < -tolerance:
                    kkt_ok = False
                elif state == 1 and gradient[j] > tolerance:
                    kkt_ok = False
                elif state == 0 and abs(gradient[j]) > tolerance:
                    kkt_ok = False
            if kkt_ok:
                candidates.append((objective, tuple(theta)))

        if not candidates:
            raise AssertionError("Decimal active-set oracle found no KKT-feasible face")
        return min(candidates, key=lambda item: item[0])


def _fit_and_compare(x, outcome, *, degree):
    z_min, z_max = 0.01, 0.40
    z = 0.205 + 0.195 * np.asarray(x, dtype=float)
    x_used = (z - 0.205) / 0.195
    columns = [np.ones(z.size), x_used]
    if degree == 2:
        columns.append(x_used**2)
        fit = fit_quadratic_logistic_inclusion_model_v2(
            z, outcome, z_min=z_min, z_max=z_max
        )
        coefficients = np.array([fit["intercept"], fit["slope"], fit["quadratic"]])
        bounds = ((-25.0, 25.0), (-50.0, 50.0), (-50.0, 50.0))
    else:
        fit = fit_logistic_inclusion_model_v2(
            z, outcome, z_min=z_min, z_max=z_max
        )
        coefficients = np.array([fit["intercept"], fit["slope"]])
        bounds = ((-25.0, 25.0), (-50.0, 50.0))
    design = np.column_stack(columns)
    oracle_objective, oracle_theta = _decimal_active_set_oracle(
        design, np.asarray(outcome, dtype=float), bounds
    )
    production_objective = _decimal_objective_gradient_hessian(
        [[D(str(float(v))) for v in row] for row in design],
        [D(str(float(v))) for v in outcome],
        [D(str(float(v))) for v in coefficients],
    )[0]
    assert fit["projected_gradient_inf_norm"] <= 1e-7
    assert np.allclose(
        coefficients, np.asarray([float(value) for value in oracle_theta]),
        rtol=0.0, atol=3e-5,
    )
    assert abs(production_objective - oracle_objective) <= D("1e-10")
    return fit, coefficients, oracle_theta


def test_decimal_active_set_oracle_matches_well_conditioned_interior_linear_fit():
    rng = np.random.default_rng(202611061)
    x = rng.uniform(-1.0, 1.0, 350)
    outcome = rng.binomial(1, expit(0.35 - 0.9 * x)).astype(float)
    fit, _, _ = _fit_and_compare(x, outcome, degree=1)
    assert fit["fit_status"] == "converged"
    assert not fit["coefficient_bound_hit"]


@pytest.mark.parametrize(
    "x,outcome,expected_intercept,expected_slope",
    [
        (
            np.concatenate([np.zeros(100), [1.0]]),
            np.concatenate([np.zeros(1000), [1.0]]),
            -25.0, 50.0,
        ),
        (
            np.concatenate([np.zeros(1000), [-1.0]]),
            np.concatenate([np.ones(100), [0.0]]),
            25.0, 50.0,
        ),
    ],
)
def test_decimal_active_set_oracle_matches_lower_and_upper_box_faces(
    x, outcome, expected_intercept, expected_slope
):
    fit, coefficients, oracle_theta = _fit_and_compare(x, outcome, degree=1)
    assert fit["fit_status"] == "boundary_limited"
    assert fit["coefficient_bound_hit"]
    assert coefficients[0] == pytest.approx(expected_intercept, abs=3e-5)
    assert coefficients[1] == pytest.approx(expected_slope, abs=3e-5)
    assert float(oracle_theta[0]) == pytest.approx(expected_intercept, abs=3e-5)


@pytest.mark.parametrize(
    "direction,expected_slope",
    [("increasing", 50.0), ("decreasing", -50.0)],
)
def test_decimal_active_set_oracle_matches_complete_separation_faces(
    direction, expected_slope
):
    x = np.linspace(-1.0, 1.0, 201)
    outcome = (x > 0.13).astype(float)
    if direction == "decreasing":
        outcome = 1.0 - outcome
    fit, coefficients, oracle_theta = _fit_and_compare(x, outcome, degree=1)
    assert fit["fit_status"] == "boundary_limited"
    assert coefficients[1] == pytest.approx(expected_slope, abs=3e-5)
    assert float(oracle_theta[1]) == pytest.approx(expected_slope, abs=3e-5)


def test_decimal_active_set_oracle_matches_quadratic_interior_fit():
    rng = np.random.default_rng(202611062)
    x = rng.uniform(-1.0, 1.0, 500)
    outcome = rng.binomial(1, expit(0.2 - 0.7 * x + 0.55 * x**2)).astype(float)
    fit, _, _ = _fit_and_compare(x, outcome, degree=2)
    assert fit["fit_status"] == "converged"
    assert not fit["coefficient_bound_hit"]
