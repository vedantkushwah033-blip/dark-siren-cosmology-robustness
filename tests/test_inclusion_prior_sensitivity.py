from __future__ import annotations

import numpy as np
import pytest

from dark_siren_robustness.inclusion_prior_sensitivity import (
    run_inclusion_prior_scale_sensitivity,
)


def _case():
    rng = np.random.default_rng(20261108)
    validation_z = rng.uniform(0.01, 0.40, 100)
    x = (validation_z - validation_z.mean()) / validation_z.std()
    retained = rng.binomial(1, 1.0 / (1.0 + np.exp(-(-0.3 + 0.6 * x))))
    catalogue_z = np.linspace(0.02, 0.39, 28)
    events = [
        {"gw_distance": 640.0, "gw_sigma": 140.0},
        {"gw_distance": 880.0, "gw_sigma": 140.0},
    ]
    base = 1.0 + 1.5 * (catalogue_z - 0.01) / 0.39
    h0 = np.linspace(50.0, 90.0, 41)
    return validation_z, retained, catalogue_z, events, base, h0


def test_prior_scale_sensitivity_recomputes_shared_selected_h0_posterior():
    validation_z, retained, catalogue_z, events, base, h0 = _case()
    result = run_inclusion_prior_scale_sensitivity(
        validation_z,
        retained,
        h0,
        catalogue_z,
        events,
        detection_threshold_mpc=250.0,
        population_base_weights=base,
        degree=1,
        slope_multipliers=(0.5, 1.0, 2.0),
        order=5,
    )
    assert result.reference_multiplier == 1.0
    assert tuple(record.slope_multiplier for record in result.records) == (
        0.5, 1.0, 2.0
    )
    for record in result.records:
        assert record.prior_sd[0] == 5.0
        assert np.isclose(record.prior_sd[1], 2.5 * record.slope_multiplier)
        assert np.isclose(
            np.trapezoid(
                record.h0_posterior.density, record.h0_posterior.grid
            ),
            1.0,
        )
        assert np.all(np.isfinite(record.selected_likelihood.marginal_log_likelihood))
        assert record.minimum_catalogue_kish_ess >= 1.0
    reference_index = 1
    assert result.max_abs_density_difference_from_reference[reference_index] < 1e-12
    assert all(
        abs(value) < 1e-12
        for value in result.quantile_shift_from_reference[reference_index].values()
    )


def test_prior_scale_sensitivity_requires_reference_multiplier_one():
    validation_z, retained, catalogue_z, events, base, h0 = _case()
    with pytest.raises(ValueError, match="include 1.0"):
        run_inclusion_prior_scale_sensitivity(
            validation_z,
            retained,
            h0,
            catalogue_z,
            events,
            detection_threshold_mpc=250.0,
            population_base_weights=base,
            slope_multipliers=(0.5, 2.0),
            order=5,
        )


def test_prior_scale_sensitivity_rejects_nonpositive_multiplier():
    validation_z, retained, catalogue_z, events, base, h0 = _case()
    with pytest.raises(ValueError, match="distinct positive"):
        run_inclusion_prior_scale_sensitivity(
            validation_z,
            retained,
            h0,
            catalogue_z,
            events,
            detection_threshold_mpc=250.0,
            population_base_weights=base,
            slope_multipliers=(0.0, 1.0, 2.0),
            order=5,
        )


def test_prior_scale_sensitivity_preserves_fixed_domain_coordinates():
    validation_z, retained, catalogue_z, events, base, h0 = _case()
    result = run_inclusion_prior_scale_sensitivity(
        validation_z,
        retained,
        h0,
        catalogue_z,
        events,
        detection_threshold_mpc=250.0,
        population_base_weights=base,
        degree=1,
        coordinate_center=0.205,
        coordinate_scale=0.195,
        slope_multipliers=(0.5, 1.0, 2.0),
        order=5,
    )
    for record in result.records:
        assert record.inclusion_posterior.center == 0.205
        assert record.inclusion_posterior.scale == 0.195
        assert record.prior_sd[0] == 5.0
        assert np.isclose(
            record.prior_sd[1], 2.5 * record.slope_multiplier
        )
