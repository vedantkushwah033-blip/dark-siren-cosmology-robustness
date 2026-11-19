import numpy as np
from scipy.stats import norm

from audit_full_parent_h0_sbc_prior_support import _posterior_mass_below


def test_mass_below_cutoff_matches_normal_cdf():
    mu, sigma = 22.0, 0.7
    log_density = lambda x: float(norm.logpdf(x, loc=mu, scale=sigma))
    scale = float(norm.logpdf(mu, loc=mu, scale=sigma))
    total = float(np.sqrt(2.0 * np.pi) * sigma)
    reference = {
        "adaptive_log_scale": scale,
        "adaptive_normalization_scaled": total,
    }
    observed = _posterior_mass_below(log_density, reference, 10.0, 20.0)
    expected = float(norm.cdf((20.0 - mu) / sigma))
    assert np.isclose(observed, expected, atol=1e-10)


def test_mass_below_cutoff_is_zero_at_lower_boundary():
    log_density = lambda x: -0.5 * (x - 22.0) ** 2
    reference = {"adaptive_log_scale": 0.0, "adaptive_normalization_scaled": 1.0}
    assert _posterior_mass_below(log_density, reference, 10.0, 10.0) == 0.0
