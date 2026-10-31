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
    assert model["fit_status"] == "converged"
    assert model["coefficient_bound_hit"] is False
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
    assert first["fit_status"] == "boundary_limited"
    assert first["coefficient_bound_hit"] is True


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
    with pytest.raises(ValueError, match="coefficients exceed"):
        predict_inclusion_probability([0.1], {**model, "slope": 50.01})



def test_validation_estimator_integrates_with_bounded_host_density():
    from scipy.integrate import trapezoid

    rng = np.random.default_rng(20261030)
    z_min, z_max = 0.01, 0.40
    grid = np.linspace(z_min, z_max, 401)
    parent_z = rng.uniform(z_min, z_max, 10000)
    x = (parent_z - (z_min + z_max) / 2) / ((z_max - z_min) / 2)
    true_probability = expit(0.25 - 1.2 * x)
    retained = rng.random(parent_z.size) < true_probability
    observed = parent_z[retained]
    validation_idx = rng.choice(parent_z.size, size=3000, replace=False)

    from dark_siren_robustness.estimated_inclusion import estimate_host_density_from_validation

    result = estimate_host_density_from_validation(
        grid,
        observed,
        parent_z[validation_idx],
        retained[validation_idx].astype(int),
        bandwidth_factor=0.8,
    )
    density = result["density"]
    assert np.all(np.isfinite(density))
    assert np.all(density >= 0)
    assert trapezoid(density, grid) == pytest.approx(1.0, abs=1e-10)
    assert result["n_validation"] == len(validation_idx)
    assert 0 < result["inverse_weight_effective_sample_size"] <= len(observed)
    assert result["inverse_weight_maximum"] >= 1.0
    assert len(result["observed_inclusion_probability"]) == len(observed)


def test_validation_estimator_rejects_audit_sample_with_single_class():
    from dark_siren_robustness.estimated_inclusion import estimate_host_density_from_validation

    grid = np.linspace(0.01, 0.40, 101)
    observed = np.linspace(0.02, 0.39, 20)
    validation = np.linspace(0.01, 0.40, 50)
    with pytest.raises(ValueError, match="both"):
        estimate_host_density_from_validation(
            grid, observed, validation, np.ones(validation.size, dtype=int)
        )


def test_validation_estimator_passes_grid_level_merger_weight():
    from scipy.integrate import trapezoid
    from dark_siren_robustness.estimated_inclusion import estimate_host_density_from_validation

    rng = np.random.default_rng(20261031)
    grid = np.linspace(0.01, 0.40, 301)
    parent = rng.uniform(grid[0], grid[-1], 5000)
    x = (parent - (grid[0] + grid[-1]) / 2) / ((grid[-1] - grid[0]) / 2)
    p = expit(0.1 - 0.8 * x)
    retained = rng.binomial(1, p)
    idx = rng.choice(parent.size, 2500, replace=False)
    observed = parent[retained.astype(bool)]
    merger_weight = 1.0 + 2.0 * grid

    result = estimate_host_density_from_validation(
        grid, observed, parent[idx], retained[idx], merger_weight=merger_weight
    )
    baseline = estimate_host_density_from_validation(
        grid, observed, parent[idx], retained[idx]
    )
    assert trapezoid(result["density"], grid) == pytest.approx(1.0, abs=1e-10)
    assert trapezoid(grid * result["density"], grid) > trapezoid(
        grid * baseline["density"], grid
    )


def test_preposterior_diagnostics_are_observable_only_and_reconstruct_weight_metrics():
    from dark_siren_robustness.preposterior_diagnostics import build_preposterior_diagnostics

    validation_z = np.linspace(0.02, 0.39, 200)
    validation_y = np.tile([0, 1], 100)
    catalogue_z = np.linspace(0.03, 0.38, 80)
    probability = expit(0.2 - 1.1 * ((catalogue_z - 0.205) / 0.195))
    model = {
        "model": "logistic_inclusion_v1", "fit_status": "converged",
        "coefficient_bound_hit": False, "intercept": 0.2, "slope": -1.1,
        "z_min": 0.01, "z_max": 0.40, "n_validation": 200, "n_included": 100, "n_not_included": 100,
    }
    result = build_preposterior_diagnostics(
        model, validation_z, validation_y, catalogue_z, probability
    )
    weights = 1.0 / probability
    expected_ess = np.sum(weights) ** 2 / np.sum(weights**2)
    assert result["diagnostic_schema"] == "preposterior_weight_stability_v1"
    assert result["n_validation"] == 200
    assert result["n_validation_included"] == 100
    assert result["inverse_weight_ess"] == pytest.approx(expected_ess)
    assert result["inverse_weight_ess_fraction"] == pytest.approx(expected_ess / 80)
    normalized = weights / weights.sum()
    entropy_fraction = np.exp(-np.sum(normalized * np.log(normalized))) / len(weights)
    assert result["inverse_weight_entropy_ess_fraction"] == pytest.approx(entropy_fraction)
    assert result["inverse_weight_max_share"] == pytest.approx(normalized.max())
    assert result["inverse_weight_top_1pct_mass"] >= result["inverse_weight_max_share"]
    assert result["inverse_weight_top_5pct_mass"] >= result["inverse_weight_top_1pct_mass"]
    assert 0 <= result["catalogue_near_validation_boundary_fraction"] <= 1
    assert result["catalogue_outside_validation_range_fraction"] == 0.0
    assert not any("true_" in key or "posterior" in key or "coverage" in key for key in result)


def test_preposterior_diagnostics_flag_catalogue_outside_validation_support():
    from dark_siren_robustness.preposterior_diagnostics import build_preposterior_diagnostics

    validation_z = np.linspace(0.10, 0.30, 100)
    catalogue_z = np.array([0.02, 0.15, 0.25, 0.39])
    # A zero-intercept, zero-slope logistic model predicts p=0.5 at every z.
    # Keep the fixture internally consistent so this test isolates support-gap logic.
    probability = np.full(catalogue_z.size, 0.5)
    model = {"model": "logistic_inclusion_v1", "fit_status": "converged", "coefficient_bound_hit": False, "intercept": 0.0, "slope": 0.0, "z_min": 0.01, "z_max": 0.40, "n_validation": 100, "n_included": 50, "n_not_included": 50}
    result = build_preposterior_diagnostics(
        model, validation_z, np.tile([0, 1], 50), catalogue_z, probability
    )
    assert result["catalogue_outside_validation_range_fraction"] == pytest.approx(0.5)


def test_preposterior_diagnostics_reject_invalid_probabilities_and_mismatched_arrays():
    from dark_siren_robustness.preposterior_diagnostics import build_preposterior_diagnostics

    validation_z = np.linspace(0.01, 0.40, 20)
    validation_y = np.tile([0, 1], 10)
    catalogue_z = np.linspace(0.02, 0.39, 10)
    model = {"model": "logistic_inclusion_v1", "intercept": 0.0, "slope": 0.0}
    with pytest.raises(ValueError, match="strictly inside"):
        build_preposterior_diagnostics(model, validation_z, validation_y, catalogue_z, np.zeros(10))
    with pytest.raises(ValueError, match="match catalogue"):
        build_preposterior_diagnostics(model, validation_z, validation_y, catalogue_z, np.full(9, 0.5))


def test_preposterior_diagnostics_reject_coefficients_outside_fit_bounds():
    from dark_siren_robustness.preposterior_diagnostics import build_preposterior_diagnostics

    validation_z = np.linspace(0.02, 0.39, 200)
    validation_y = np.tile([0, 1], 100)
    catalogue_z = np.linspace(0.03, 0.38, 80)
    model = {
        "model": "logistic_inclusion_v1", "fit_status": "boundary_limited",
        "coefficient_bound_hit": True, "intercept": 0.0, "slope": 50.01,
        "z_min": 0.01, "z_max": 0.40, "n_validation": 200,
        "n_included": 100, "n_not_included": 100,
    }
    with pytest.raises(ValueError, match="coefficients exceed"):
        build_preposterior_diagnostics(
            model, validation_z, validation_y, catalogue_z, np.full(80, 0.5)
        )

@pytest.mark.parametrize(
    "field, value",
    [
        ("n_validation", 200.5),
        ("n_included", True),
        ("n_not_included", 99),
    ],
)
def test_preposterior_diagnostics_reject_noninteger_or_inconsistent_model_counts(field, value):
    from dark_siren_robustness.preposterior_diagnostics import build_preposterior_diagnostics

    validation_z = np.linspace(0.02, 0.39, 200)
    validation_y = np.tile([0, 1], 100)
    catalogue_z = np.linspace(0.03, 0.38, 80)
    probability = np.full(80, 0.5)
    model = {
        "model": "logistic_inclusion_v1", "fit_status": "converged",
        "coefficient_bound_hit": False, "intercept": 0.0, "slope": 0.0,
        "z_min": 0.01, "z_max": 0.40, "n_validation": 200,
        "n_included": 100, "n_not_included": 100,
    }
    model[field] = value
    with pytest.raises(ValueError, match="validation counts"):
        build_preposterior_diagnostics(
            model, validation_z, validation_y, catalogue_z, probability
        )

def test_preposterior_diagnostics_reject_probabilities_not_generated_by_fitted_model():
    from dark_siren_robustness.preposterior_diagnostics import build_preposterior_diagnostics

    validation_z = np.linspace(0.02, 0.39, 200)
    validation_y = np.tile([0, 1], 100)
    catalogue_z = np.linspace(0.03, 0.38, 80)
    model = {
        "model": "logistic_inclusion_v1", "fit_status": "converged",
        "coefficient_bound_hit": False, "intercept": 0.2, "slope": -1.1,
        "z_min": 0.01, "z_max": 0.40, "n_validation": 200,
        "n_included": 100, "n_not_included": 100,
    }
    with pytest.raises(ValueError, match="disagree with fitted model"):
        build_preposterior_diagnostics(
            model, validation_z, validation_y, catalogue_z, np.full(80, 0.5)
        )


def test_quadratic_logistic_inclusion_fit_recovers_curvature_and_predicts():
    rng = np.random.default_rng(20261031)
    z_min, z_max = 0.01, 0.40
    z = rng.uniform(z_min, z_max, 30000)
    x = (z - (z_min + z_max) / 2) / ((z_max - z_min) / 2)
    truth = expit(0.25 - 1.2 * x + 1.8 * (x**2 - 1.0 / 3.0))
    included = rng.binomial(1, truth)
    model = fit_logistic_inclusion_model(z, included, z_min=z_min, z_max=z_max, degree=2)
    grid = np.linspace(z_min, z_max, 301)
    predicted = predict_inclusion_probability(grid, model)
    gx = (grid - (z_min + z_max) / 2) / ((z_max - z_min) / 2)
    expected = expit(0.25 - 1.2 * gx + 1.8 * (gx**2 - 1.0 / 3.0))
    assert model["model"] == "quadratic_logistic_inclusion_v1"
    assert model["degree"] == 2
    assert model["fit_status"] == "converged"
    assert model["quadratic"] == pytest.approx(1.8, abs=0.10)
    assert np.mean(np.abs(predicted - expected)) < 0.02


def test_inclusion_fit_rejects_unsupported_polynomial_degree():
    with pytest.raises(ValueError, match="degree"):
        fit_logistic_inclusion_model(
            np.linspace(0.01, 0.40, 20), np.tile([0, 1], 10),
            z_min=0.01, z_max=0.40, degree=3,
        )


def test_quadratic_preposterior_diagnostics_reconstruct_fitted_probabilities():
    from dark_siren_robustness.preposterior_diagnostics import build_preposterior_diagnostics
    rng = np.random.default_rng(20261032)
    z_min, z_max = 0.01, 0.40
    validation_z = rng.uniform(z_min, z_max, 12000)
    x = (validation_z - (z_min + z_max) / 2) / ((z_max - z_min) / 2)
    included = rng.binomial(1, expit(0.2 - 0.8*x + 1.2*(x*x - 1/3)))
    model = fit_logistic_inclusion_model(
        validation_z, included, z_min=z_min, z_max=z_max, degree=2
    )
    catalogue_z = np.linspace(0.02, 0.39, 80)
    probability = predict_inclusion_probability(catalogue_z, model)
    result = build_preposterior_diagnostics(
        model, validation_z, included, catalogue_z, probability
    )
    assert result["fitted_model_degree"] == 2
    assert result["fitted_quadratic"] == pytest.approx(model["quadratic"])
    assert 0 < result["inverse_weight_ess"] <= len(catalogue_z)
