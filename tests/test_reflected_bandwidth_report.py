"""Tests for reporting realized physical reflected-KDE bandwidths."""
import numpy as np
import pytest
from dark_siren_robustness.catalogue_density import reflected_weighted_kde_bandwidth


def test_reported_bandwidth_matches_reflected_weighted_kde_covariance():
    z = np.array([0.03, 0.06, 0.09, 0.14, 0.19, 0.23, 0.28, 0.33, 0.37])
    weights = np.array([1.0, 1.2, 1.4, 1.8, 2.0, 2.4, 3.0, 3.5, 4.0])
    multiplier = 0.73
    scaled = weights / weights.max()
    probability = scaled / scaled.sum()
    mean = float(probability @ z)
    variance = float(
        np.sum(probability * (z - mean) ** 2)
        / (1.0 - float(probability @ probability))
    )
    neff = scaled.sum() ** 2 / np.sum(scaled**2)
    expected = multiplier * neff ** (-1.0 / 5.0) * np.sqrt(variance)
    observed = reflected_weighted_kde_bandwidth(
        z, lower_bound=0.01, upper_bound=0.40,
        sample_weights=weights, scott_multiplier=multiplier,
    )
    assert observed == pytest.approx(expected, rel=1e-12)


def test_bandwidth_report_rejects_invalid_weights():
    with pytest.raises(ValueError, match="sample_weights"):
        reflected_weighted_kde_bandwidth(
            np.linspace(0.02, 0.38, 9), lower_bound=0.01,
            upper_bound=0.40, sample_weights=np.zeros(9),
            scott_multiplier=1.0,
        )
