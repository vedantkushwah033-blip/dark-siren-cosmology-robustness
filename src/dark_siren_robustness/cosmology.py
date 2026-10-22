"""Small, explicit cosmology utilities used by the simulation study."""

from __future__ import annotations

import numpy as np
from scipy.integrate import trapezoid
from numpy.typing import ArrayLike, NDArray

C_KM_S = 299_792.458


def luminosity_distance_flat_lcdm(
    redshift: ArrayLike,
    h0_km_s_mpc: float,
    omega_m: float = 0.3,
    integration_steps: int = 2_000,
) -> NDArray[np.float64]:
    """Compute flat Lambda-CDM luminosity distance in Mpc by numerical integration.

    This transparent reference implementation should be independently validated
    over the redshift range used in any scientific analysis.
    """
    z = np.asarray(redshift, dtype=float)
    if np.any(~np.isfinite(z)) or np.any(z < 0):
        raise ValueError("redshift must contain only finite, non-negative values")
    if not np.isfinite(h0_km_s_mpc) or h0_km_s_mpc <= 0:
        raise ValueError("h0_km_s_mpc must be finite and positive")
    if not np.isfinite(omega_m) or not 0 < omega_m < 1:
        raise ValueError("omega_m must lie strictly between 0 and 1")
    if not isinstance(integration_steps, int) or integration_steps < 10:
        raise ValueError("integration_steps must be an integer >= 10")

    flat_z = z.ravel()
    result = np.zeros_like(flat_z)
    for index, zi in enumerate(flat_z):
        if zi == 0:
            continue
        grid = np.linspace(0.0, zi, integration_steps + 1)
        e_z = np.sqrt(omega_m * (1.0 + grid) ** 3 + (1.0 - omega_m))
        integral = trapezoid(1.0 / e_z, grid)
        result[index] = (1.0 + zi) * (C_KM_S / h0_km_s_mpc) * integral
    return result.reshape(z.shape)


def low_redshift_luminosity_distance(
    redshift: ArrayLike, h0_km_s_mpc: float
) -> NDArray[np.float64]:
    """Low-redshift approximation d_L ~= c z / H0, in Mpc.

    This approximation omits higher-order cosmological terms and is not a
    substitute for a full cosmological likelihood.
    """
    z = np.asarray(redshift, dtype=float)
    if np.any(~np.isfinite(z)) or np.any(z < 0):
        raise ValueError("redshift must contain only finite, non-negative values")
    if not np.isfinite(h0_km_s_mpc) or h0_km_s_mpc <= 0:
        raise ValueError("h0_km_s_mpc must be finite and positive")
    return C_KM_S * z / h0_km_s_mpc
