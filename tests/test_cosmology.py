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


def test_flat_lcdm_scales_inversely_with_h0():
    z = np.array([0.03, 0.1, 0.4])
    d70 = luminosity_distance_flat_lcdm(z, 70.0)
    d80 = luminosity_distance_flat_lcdm(z, 80.0)
    assert d80 == pytest.approx(d70 * (70.0 / 80.0), rel=1e-12)


def test_flat_lcdm_approaches_low_redshift_limit():
    z = 1e-5
    exact = luminosity_distance_flat_lcdm(z, 70.0)
    approximation = low_redshift_luminosity_distance(z, 70.0)
    assert exact == pytest.approx(approximation, rel=1e-4)


def test_flat_lcdm_agrees_with_independent_quadrature():
    from scipy.integrate import quad

    z = 0.7
    h0 = 70.0
    omega_m = 0.3
    integral, _ = quad(
        lambda zp: 1.0 / np.sqrt(omega_m * (1.0 + zp) ** 3 + 1.0 - omega_m),
        0.0,
        z,
        epsabs=1e-12,
        epsrel=1e-12,
    )
    reference = (1.0 + z) * (299_792.458 / h0) * integral
    computed = luminosity_distance_flat_lcdm(
        z, h0, omega_m=omega_m, integration_steps=20_000
    )
    assert computed == pytest.approx(reference, rel=2e-8)


def test_flat_lcdm_preserves_array_shape_and_finite_values():
    z = np.array([[0.0, 0.05], [0.1, 0.3]])
    distances = luminosity_distance_flat_lcdm(z, 70.0)
    assert distances.shape == z.shape
    assert np.all(np.isfinite(distances))
    assert np.all(distances >= 0)
