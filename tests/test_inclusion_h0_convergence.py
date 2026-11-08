from __future__ import annotations

import numpy as np

from dark_siren_robustness.inclusion_h0_convergence import (
    run_inclusion_h0_quadrature_convergence,
)


def _case():
    rng = np.random.default_rng(20261108)
    validation_z = rng.uniform(0.01, 0.40, 120)
    x = (validation_z - validation_z.mean()) / validation_z.std()
    retained = rng.binomial(1, 1.0 / (1.0 + np.exp(-(-0.25 + 0.55 * x))))
    catalogue_z = np.linspace(0.015, 0.395, 35)
    events = [
        {"gw_distance": 520.0, "gw_sigma": 140.0},
        {"gw_distance": 760.0, "gw_sigma": 140.0},
        {"gw_distance": 910.0, "gw_sigma": 140.0},
    ]
    base = 1.0 + 1.5 * (catalogue_z - 0.01) / 0.39
    h0 = np.linspace(50.0, 90.0, 41)
    return validation_z, retained, h0, catalogue_z, events, base


def test_downstream_h0_quadrature_convergence_recomputes_full_likelihood():
    validation_z, retained, h0, catalogue_z, events, base = _case()
    result = run_inclusion_h0_quadrature_convergence(
        validation_z, retained, h0, catalogue_z, events,
        detection_threshold_mpc=250.0,
        population_base_weights=base,
        degree=1,
        orders=(5, 7, 9),
        centered_log_likelihood_tolerance=1e6,
        posterior_density_tolerance=1e6,
        posterior_quantile_tolerance=1e6,
        consecutive_passes=1,
    )
    assert result.status == "CONVERGED"
    assert result.orders == (5, 7)
    assert len(result.fits) == len(result.likelihoods) == len(result.posteriors) == 2
    assert result.centered_log_likelihood_max_abs_change[0] >= 0.0
    assert result.posterior_density_max_abs_change[0] >= 0.0
    assert result.posterior_quantile_max_abs_change[0] >= 0.0
    for likelihood, posterior in zip(result.likelihoods, result.posteriors):
        assert likelihood.marginal_log_likelihood.shape == (41,)
        assert np.isclose(np.trapezoid(posterior.density, posterior.grid), 1.0)


def test_downstream_h0_quadrature_convergence_can_remain_nonconverged():
    validation_z, retained, h0, catalogue_z, events, base = _case()
    result = run_inclusion_h0_quadrature_convergence(
        validation_z, retained, h0, catalogue_z, events,
        detection_threshold_mpc=250.0,
        population_base_weights=base,
        degree=1,
        orders=(5, 7, 9),
        centered_log_likelihood_tolerance=1e-14,
        posterior_density_tolerance=1e-14,
        posterior_quantile_tolerance=1e-14,
        consecutive_passes=2,
    )
    assert result.status == "NON_CONVERGED"
    assert result.orders == (5, 7, 9)
    assert len(result.centered_log_likelihood_max_abs_change) == 2


def test_downstream_convergence_rejects_invalid_order_schedule():
    validation_z, retained, h0, catalogue_z, events, base = _case()
    try:
        run_inclusion_h0_quadrature_convergence(
            validation_z, retained, h0, catalogue_z, events,
            detection_threshold_mpc=250.0,
            population_base_weights=base,
            orders=(7, 7, 11),
        )
    except ValueError as error:
        assert "strictly increasing" in str(error)
    else:
        raise AssertionError("invalid quadrature order schedule was accepted")


def test_downstream_convergence_preserves_fixed_domain_coordinates():
    validation_z, retained, h0, catalogue_z, events, base = _case()
    result = run_inclusion_h0_quadrature_convergence(
        validation_z, retained, h0, catalogue_z, events,
        detection_threshold_mpc=250.0,
        population_base_weights=base,
        degree=1,
        coordinate_center=0.205,
        coordinate_scale=0.195,
        orders=(5, 7),
        centered_log_likelihood_tolerance=1e6,
        posterior_density_tolerance=1e6,
        posterior_quantile_tolerance=1e6,
        consecutive_passes=1,
    )
    assert result.orders == (5, 7)
    assert all(np.isclose(fit.center, 0.205) for fit in result.fits)
    assert all(np.isclose(fit.scale, 0.195) for fit in result.fits)

