import numpy as np
import pytest
from scipy.integrate import quad

from dark_siren_robustness.cosmology import luminosity_distance_flat_lcdm


def _quad_reference(z: float, h0: float = 70.0, omega_m: float = 0.3) -> float:
    integral, _ = quad(
        lambda zp: 1.0 / np.sqrt(omega_m * (1.0 + zp) ** 3 + 1.0 - omega_m),
        0.0, z, epsabs=1e-12, epsrel=1e-12,
    )
    return (1.0 + z) * (299_792.458 / h0) * integral


@pytest.mark.parametrize("z", [0.01, 0.1, 0.7, 1.5])
def test_trapezoid_resolution_converges_toward_quad_reference(z):
    reference = _quad_reference(z)
    estimates = [
        float(luminosity_distance_flat_lcdm(z, 70.0, integration_steps=steps))
        for steps in (500, 2_000, 10_000)
    ]
    errors = [abs(value - reference) for value in estimates]
    assert errors[2] < errors[1] < errors[0]
    assert errors[1] / reference < 1e-6
    assert errors[2] / reference < 1e-7


def test_default_resolution_is_accurate_over_reference_redshift_grid():
    for z in [0.01, 0.1, 0.3, 0.7, 1.0, 1.5]:
        reference = _quad_reference(z)
        computed = float(luminosity_distance_flat_lcdm(z, 70.0))
        assert computed == pytest.approx(reference, rel=1e-6)