import numpy as np
import pytest

from dark_siren_robustness.cosmology import (
    luminosity_distance_flat_lcdm,
    low_redshift_luminosity_distance,
)


def test_zero_redshift_has_zero_luminosity_distance():
    assert luminosity_distance_flat_lcdm(0.0, 70.0) == pytest.approx(0.0)


def test_luminosity_distance_increases_with_redshift():
    distances = luminosity_distance_flat_lcdm([0.01, 0.05, 0.1], 70.0)
    assert np.all(np.diff(distances) > 0)


def test_low_redshift_approximation_scales_linearly():
    distances = low_redshift_luminosity_distance([0.01, 0.02], 70.0)
    assert distances[1] == pytest.approx(2 * distances[0])


@pytest.mark.parametrize("kwargs", [
    {"redshift": -0.1, "h0_km_s_mpc": 70.0},
    {"redshift": 0.1, "h0_km_s_mpc": 0.0},
    {"redshift": 0.1, "h0_km_s_mpc": 70.0, "omega_m": 1.0},
])
def test_invalid_cosmology_inputs_raise(kwargs):
    with pytest.raises(ValueError):
        luminosity_distance_flat_lcdm(**kwargs)
