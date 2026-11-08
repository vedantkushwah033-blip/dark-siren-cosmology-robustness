import numpy as np
import pytest
from scipy.special import expit
from dark_siren_robustness.estimated_inclusion import fit_inclusion_logistic


def test_linear_logistic_fit_recovers_probabilities_on_large_validation_sample():
    rng = np.random.default_rng(20261108)
    z = rng.uniform(0.01, 0.40, 5000)
    x = (z - 0.205) / 0.11
    truth = expit(0.25 - 0.9 * x)
    labels = rng.binomial(1, truth)
    fit = fit_inclusion_logistic(z, labels, degree=1)
    assert fit.status == "SUCCESS"
    assert fit.coefficients is not None and fit.covariance is not None
    assert np.all(np.linalg.eigvalsh(fit.covariance) > 0)
    prediction = fit.predict(z)
    assert prediction.shape == z.shape
    assert np.all((prediction > 0) & (prediction < 1))
    assert np.mean((prediction - truth) ** 2) < 0.002


def test_quadratic_model_is_supported_and_predicts_finite_probabilities():
    rng = np.random.default_rng(11)
    z = rng.uniform(0.01, 0.40, 1800)
    x = (z - 0.205) / 0.11
    labels = rng.binomial(1, expit(0.1 - 0.5 * x + 0.35 * (x**2 - 1)))
    fit = fit_inclusion_logistic(z, labels, degree=2)
    assert fit.status == "SUCCESS"
    assert np.all(np.isfinite(fit.predict(np.linspace(0.01, 0.40, 100))))


def test_single_class_validation_is_non_evaluable():
    z = np.linspace(0.01, 0.40, 100)
    fit = fit_inclusion_logistic(z, np.ones(100, dtype=int))
    assert fit.status == "NON_EVALUABLE"
    with pytest.raises(RuntimeError):
        fit.predict(z)


def test_rank_deficient_design_is_failed():
    z = np.repeat(0.2, 20)
    fit = fit_inclusion_logistic(z, np.tile([0, 1], 10))
    assert fit.status == "FAILED"
    assert "rank deficient" in fit.message


def test_invalid_inputs_are_rejected():
    with pytest.raises(ValueError):
        fit_inclusion_logistic([0.1, 0.2, 0.3], [0, 1, 0])
    with pytest.raises(ValueError):
        fit_inclusion_logistic(np.linspace(0.1, 0.3, 20), np.tile([0, 1], 10), degree=3)


def test_complete_separation_is_not_reported_as_success():
    z = np.linspace(0.01, 0.40, 200)
    labels = (z > 0.205).astype(int)
    fit = fit_inclusion_logistic(z, labels)
    assert fit.status == "FAILED"
    assert fit.coefficients is None


def test_precision_loss_is_recovered_by_damped_newton_polish(monkeypatch):
    from types import SimpleNamespace
    import dark_siren_robustness.estimated_inclusion as inclusion_module

    rng = np.random.default_rng(731)
    z = rng.uniform(0.01, 0.40, 600)
    x = (z - z.mean()) / z.std()
    labels = rng.binomial(1, expit(-0.3 + 0.65 * x))

    # Simulate BFGS returning a finite but not fully converged point with a
    # precision-loss status. The analytic Newton fallback should finish it.
    monkeypatch.setattr(
        inclusion_module,
        "minimize",
        lambda *args, **kwargs: SimpleNamespace(
            x=np.zeros(2), success=False,
            message="Desired error not necessarily achieved due to precision loss.",
            nit=3,
        ),
    )
    fit = fit_inclusion_logistic(z, labels, degree=1)
    assert fit.status == "SUCCESS"
    assert fit.gradient_inf_norm <= 1e-10
    assert fit.iterations > 3
    assert fit.coefficients is not None
    assert fit.covariance is not None
