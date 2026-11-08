import numpy as np
import pytest

from dark_siren_robustness.inclusion_quadrature_convergence import (
    fit_inclusion_posterior_with_order_ladder,
)


def _validation_sample():
    rng = np.random.default_rng(55201)
    z = rng.uniform(0.01, 0.40, 240)
    x = (z - z.mean()) / z.std()
    probability = 1.0 / (1.0 + np.exp(-(-0.25 + 0.7 * x)))
    return z, rng.binomial(1, probability)


def test_order_ladder_returns_converged_with_full_transition_record():
    z, retained = _validation_sample()
    result = fit_inclusion_posterior_with_order_ladder(
        z, retained, degree=1, orders=(7, 9),
        mean_tolerance=10.0, sd_tolerance=10.0,
        predictive_tolerance=10.0, log_normalizer_tolerance=10.0,
    )
    assert result.status == "CONVERGED"
    assert result.attempted_orders == (7, 9)
    assert len(result.checks) == 1
    assert result.checks[0].passed
    assert result.posterior.order == 9
    assert result.checks[0].quadrature_ess >= 1.0
    assert 0.0 < result.checks[0].maximum_normalized_weight <= 1.0


def test_order_ladder_reports_non_convergence_at_maximum_order():
    z, retained = _validation_sample()
    result = fit_inclusion_posterior_with_order_ladder(
        z, retained, degree=1, orders=(5, 7),
        mean_tolerance=1e-15, sd_tolerance=1e-15,
        predictive_tolerance=1e-15, log_normalizer_tolerance=1e-15,
    )
    assert result.status == "NON_CONVERGED"
    assert result.attempted_orders == (5, 7)
    assert len(result.checks) == 1
    assert not result.checks[0].passed
    assert result.posterior.order == 7


@pytest.mark.parametrize("orders", [(7,), (7, 7), (9, 7), (True, 9), (2, 7)])
def test_order_ladder_rejects_invalid_order_sequences(orders):
    z, retained = _validation_sample()
    with pytest.raises(ValueError):
        fit_inclusion_posterior_with_order_ladder(z, retained, orders=orders)


def test_order_ladder_rejects_nonpositive_tolerance():
    z, retained = _validation_sample()
    with pytest.raises(ValueError, match="tolerances"):
        fit_inclusion_posterior_with_order_ladder(
            z, retained, orders=(7, 9), mean_tolerance=0.0
        )
