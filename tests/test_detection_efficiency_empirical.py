"""Empirical checks for the toy detection-efficiency calculation.

These tests validate internal consistency of the declared Gaussian threshold
toy model only; they do not approximate a gravitational-wave detector.
"""
import numpy as np

from dark_siren_robustness.dataset_likelihood import event_detection_efficiency


def test_empirical_detection_rate_matches_untruncated_gaussian_efficiency():
    rng = np.random.default_rng(20261028)
    h0 = 70.0
    z = 0.03
    sigma = 35.0
    threshold = 130.0
    mean = 299_792.458 * z / h0

    expected = event_detection_efficiency(
        [h0, 71.0], redshifts=[z], gw_sigma=sigma,
        detection_distance_min=threshold,
    )[0]
    observed_distance = rng.normal(loc=mean, scale=sigma, size=200_000)
    empirical = np.mean(observed_distance >= threshold)

    # Fixed seed and 200k draws make this a high-signal stochastic regression
    # test while allowing several Monte Carlo standard errors of tolerance.
    mc_se = np.sqrt(expected * (1.0 - expected) / observed_distance.size)
    assert abs(empirical - expected) < 5.0 * mc_se + 1e-4


def test_empirical_detection_rate_matches_positive_truncated_efficiency():
    rng = np.random.default_rng(20261029)
    h0 = 70.0
    z = 0.005
    sigma = 35.0
    threshold = 25.0
    mean = 299_792.458 * z / h0

    expected = event_detection_efficiency(
        [h0, 71.0], redshifts=[z], gw_sigma=sigma,
        detection_distance_min=threshold, positive_truncated=True,
    )[0]

    # Rejection sampling explicitly represents the model d | d > 0 used by
    # positive_truncated=True. Generate in batches to keep runtime bounded.
    accepted = []
    total = 0
    target = 200_000
    while total < target:
        draw = rng.normal(loc=mean, scale=sigma, size=250_000)
        draw = draw[draw > 0.0]
        take = min(target - total, draw.size)
        if take:
            accepted.append(draw[:take])
            total += take
    positive_distance = np.concatenate(accepted)
    empirical = np.mean(positive_distance >= threshold)

    mc_se = np.sqrt(expected * (1.0 - expected) / positive_distance.size)
    assert abs(empirical - expected) < 5.0 * mc_se + 1e-4


def test_empirical_detection_rate_matches_weighted_catalogue_and_missing_host_mixture():
    rng = np.random.default_rng(20261030)
    h0 = 70.0
    sigma = 30.0
    threshold = 115.0
    catalog_z = np.array([0.018, 0.032, 0.045])
    catalog_weights = np.array([0.2, 0.5, 0.3])
    missing_z = np.array([0.055, 0.075])
    missing_weights = np.array([0.35, 0.65])
    missing_mass = 0.4
    n = 300_000

    expected = event_detection_efficiency(
        [h0, 71.0],
        redshifts=catalog_z,
        weights=catalog_weights,
        gw_sigma=sigma,
        detection_distance_min=threshold,
        missing_weight=missing_mass,
        missing_redshifts=missing_z,
        missing_weights=missing_weights,
    )[0]

    # Draw hosts from exactly the same two-component mixture as the analytic
    # beta(H0): first choose population component, then choose a host by weight.
    is_missing = rng.random(n) < missing_mass
    z = np.empty(n)
    n_missing = int(is_missing.sum())
    n_catalog = n - n_missing
    z[~is_missing] = rng.choice(catalog_z, size=n_catalog, p=catalog_weights)
    z[is_missing] = rng.choice(missing_z, size=n_missing, p=missing_weights)
    observed_distance = 299_792.458 * z / h0 + rng.normal(0.0, sigma, size=n)
    empirical = np.mean(observed_distance >= threshold)

    mc_se = np.sqrt(expected * (1.0 - expected) / n)
    assert abs(empirical - expected) < 5.0 * mc_se + 1e-4
