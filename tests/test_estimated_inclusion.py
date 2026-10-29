import numpy as np
import pytest
from scipy.special import expit

from dark_siren_robustness.estimated_inclusion import (
    fit_logistic_inclusion_model,
    predict_inclusion_probability,
)


def test_logistic_inclusion_fit_recovers_declared_synthetic_curve():
    rng = np.random.default_rng(20261029)
    z_min, z_max = 0.01, 0.40
    z = rng.uniform(z_min, z_max, 6000)
    x = (z - (z_min + z_max) / 2) / ((z_max - z_min) / 2)
    expected = expit(0.35 - 1.8 * x)
    included = rng.binomial(1, expected)

    model = fit_logistic_inclusion_model(z, included, z_min=z_min, z_max=z_max)
    grid = np.linspace(z_min, z_max, 101)
    predicted = predict_inclusion_probability(grid, model)
    truth = expit(0.35 - 1.8 * ((grid - (z_min + z_max) / 2) / ((z_max - z_min) / 2)))

    assert model["model"] == "logistic_inclusion_v1"
    assert model["n_validation"] == len(z)
    assert model["n_included"] + model["n_not_included"] == len(z)
    assert np.all(np.isfinite(predicted))
    assert np.all((predicted > 0) & (predicted < 1))
    assert np.mean(np.abs(predicted - truth)) < 0.025


def test_logistic_inclusion_fit_is_deterministic_with_separation():
    z = np.linspace(0.01, 0.40, 100)
    included = (z < 0.23).astype(int)
    first = fit_logistic_inclusion_model(z, included, z_min=0.01, z_max=0.40)
    second = fit_logistic_inclusion_model(z, included, z_min=0.01, z_max=0.40)
    assert first == second
    assert np.isfinite(first["intercept"])
    assert np.isfinite(first["slope"])
    assert abs(first["slope"]) <= 50.0


@pytest.mark.parametrize(
    "z, included, kwargs, message",
    [
        (np.ones(9), np.ones(9), {"z_min": 0.01, "z_max": 0.40}, "at least 10"),
        (np.linspace(0.01, 0.40, 20), np.ones(20), {"z_min": 0.01, "z_max": 0.40}, "both"),
        (np.linspace(0.01, 0.40, 20), np.arange(20), {"z_min": 0.01, "z_max": 0.40}, "binary"),
        (np.linspace(0.00, 0.40, 20), np.tile([0, 1], 10), {"z_min": 0.01, "z_max": 0.40}, "inside"),
        (np.linspace(0.01, 0.40, 20), np.tile([0, 1], 10), {"z_min": 0.40, "z_max": 0.01}, "z_min"),
    ],
)
def test_logistic_inclusion_fit_rejects_invalid_inputs(z, included, kwargs, message):
    with pytest.raises(ValueError, match=message):
        fit_logistic_inclusion_model(z, included, **kwargs)


def test_predict_rejects_unknown_model_and_out_of_domain_values():
    model = {
        "model": "logistic_inclusion_v1",
        "intercept": 0.0,
        "slope": -1.0,
        "z_min": 0.01,
        "z_max": 0.40,
    }
    with pytest.raises(ValueError, match="model"):
        predict_inclusion_probability([0.1], {**model, "model": "other"})
    with pytest.raises(ValueError, match="inside"):
        predict_inclusion_probability([0.5], model)
