import math
import numpy as np
import pytest
from src.dark_siren_robustness.uncertainty_likelihoods import (
    C_KM_S, candidate_distance_mpc, gaussian_log_likelihood,
    log_likelihood_fixed_sigma, log_likelihood_fractional_sigma,
    log_likelihood_fixed_reference_sigma,
)

def test_distance_relation():
    assert candidate_distance_mpc(70.0, 0.01) == pytest.approx(C_KM_S * 0.01 / 70.0)

def test_normalized_gaussian():
    got = gaussian_log_likelihood(100.0, 95.0, 10.0)
    want = -0.5 * (5.0 / 10.0) ** 2 - math.log(10.0) - 0.5 * math.log(2.0 * math.pi)
    assert got == pytest.approx(want)

def test_fixed_sigma_matches_formula():
    h = np.array([60.0, 70.0, 80.0])
    z, obs, sigma = 0.02, 90.0, 12.0
    mu = C_KM_S * z / h
    want = -0.5 * ((obs - mu) / sigma) ** 2 - np.log(sigma) - 0.5 * np.log(2 * np.pi)
    np.testing.assert_allclose(log_likelihood_fixed_sigma(obs, h, z, sigma), want)

def test_fractional_sigma_contains_candidate_normalization():
    h = np.array([60.0, 70.0, 80.0])
    z, obs, f = 0.015, 70.0, 0.1
    mu = C_KM_S * z / h
    sigma = f * mu
    want = -0.5 * ((obs - mu) / sigma) ** 2 - np.log(sigma) - 0.5 * np.log(2 * np.pi)
    got = log_likelihood_fractional_sigma(obs, h, z, f)
    np.testing.assert_allclose(got, want)
    assert not np.allclose(got + 0.5 * ((obs - mu) / sigma) ** 2, (got + 0.5 * ((obs - mu) / sigma) ** 2)[0])

def test_reference_sigma_is_fixed_over_candidate_grid():
    h = np.array([55.0, 70.0, 85.0])
    z, obs, f, href = 0.02, 100.0, 0.1, 70.0
    sigma = f * C_KM_S * z / href
    want = log_likelihood_fixed_sigma(obs, h, z, sigma)
    got = log_likelihood_fixed_reference_sigma(obs, h, z, f, reference_h0=href)
    np.testing.assert_allclose(got, want)

@pytest.mark.parametrize("sigma", [0.0, -1.0, float("nan"), float("inf")])
def test_invalid_sigma_rejected(sigma):
    with pytest.raises(ValueError):
        gaussian_log_likelihood(1.0, 1.0, sigma)
