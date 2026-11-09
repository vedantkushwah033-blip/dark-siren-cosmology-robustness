from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import trapezoid
from scipy.stats import gaussian_kde

from dark_siren_robustness.catalogue_density import estimate_host_density


def test_density_is_normalized_and_finite_on_bounded_grid():
    rng = np.random.default_rng(4101)
    grid = np.linspace(0.01, 0.40, 2001)
    sample = rng.uniform(0.01, 0.40, 600)
    density = estimate_host_density(grid, sample)
    assert np.all(np.isfinite(density))
    assert np.all(density >= 0)
    assert trapezoid(density, grid) == pytest.approx(1.0, abs=1e-10)


def test_inverse_inclusion_weights_reduce_selection_induced_mean_shift():
    grid = np.linspace(0.01, 0.40, 4001)
    parent = np.linspace(0.0101, 0.3999, 5000)
    probability = 0.08 + 0.92 / (1.0 + np.exp((parent - 0.22) / 0.045))
    rng = np.random.default_rng(4102)
    keep = rng.random(parent.size) < probability
    observed, observed_p = parent[keep], probability[keep]
    assert observed.size > 1000

    unweighted = estimate_host_density(grid, observed, bandwidth_factor=0.8)
    corrected = estimate_host_density(
        grid, observed, inclusion_probability=observed_p, bandwidth_factor=0.8
    )
    parent_mean = (grid[-1] + grid[0]) / 2.0
    unweighted_mean = trapezoid(grid * unweighted, grid)
    corrected_mean = trapezoid(grid * corrected, grid)
    assert unweighted_mean < parent_mean - 0.005
    assert abs(corrected_mean - parent_mean) < abs(unweighted_mean - parent_mean)


def test_merger_weight_is_applied_before_final_normalization():
    grid = np.linspace(0.01, 0.40, 2001)
    sample = np.linspace(0.011, 0.399, 500)
    base = estimate_host_density(grid, sample)
    weighted = estimate_host_density(grid, sample, merger_weight=1.0 + 2.0 * grid)
    assert trapezoid(weighted, grid) == pytest.approx(1.0, abs=1e-10)
    assert trapezoid(grid * weighted, grid) > trapezoid(grid * base, grid)


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"inclusion_probability": np.array([1.0, 0.0, 1.0, 1.0, 1.0])}, "inclusion_probability"),
        ({"bandwidth_factor": 0.0}, "bandwidth_factor"),
        ({"merger_weight": np.ones(21)}, "merger_weight"),
    ],
)
def test_invalid_inputs_are_rejected(kwargs, message):
    grid = np.linspace(0.01, 0.40, 101)
    sample = np.linspace(0.02, 0.39, 20)
    if "inclusion_probability" in kwargs:
        kwargs = {"inclusion_probability": np.resize(kwargs["inclusion_probability"], sample.size)}
    with pytest.raises(ValueError, match=message):
        estimate_host_density(grid, sample, **kwargs)


def test_sample_outside_declared_domain_is_rejected():
    grid = np.linspace(0.01, 0.40, 101)
    with pytest.raises(ValueError, match="inside"):
        estimate_host_density(grid, np.linspace(0.0, 0.39, 20))



def test_reflection_does_not_artificially_triple_effective_sample_size():
    grid = np.linspace(0.01, 0.40, 301)
    sample = np.linspace(0.012, 0.398, 80)
    probability = 0.15 + 0.8 / (1.0 + np.exp((sample - 0.22) / 0.05))
    weights = 1.0 / probability
    augmented = np.concatenate((sample, 2 * grid[0] - sample, 2 * grid[-1] - sample))
    augmented_weights = np.tile(weights / 3.0, 3)
    original_neff = np.sum(weights) ** 2 / np.sum(weights ** 2)
    reference_kde = gaussian_kde(
        augmented, weights=augmented_weights,
        bw_method=original_neff ** (-1.0 / 5.0),
    )
    reference = reference_kde(grid)
    reference /= trapezoid(reference, grid)
    actual = estimate_host_density(
        grid, sample, inclusion_probability=probability
    )
    np.testing.assert_allclose(actual, reference, rtol=1e-12, atol=1e-12)


def test_inverse_probability_overflow_fails_explicitly():
    grid = np.linspace(0.01, 0.40, 101)
    sample = np.linspace(0.02, 0.39, 20)
    probability = np.full(sample.size, 0.5)
    probability[0] = 1e-320
    with pytest.raises(ValueError, match="overflow"):
        estimate_host_density(
            grid, sample, inclusion_probability=probability
        )


def test_explicit_bandwidth_factor_is_absolute_scipy_covariance_factor():
    grid = np.linspace(0.01, 0.40, 301)
    sample = np.linspace(0.012, 0.398, 80)
    probability = 0.15 + 0.8 / (1.0 + np.exp((sample - 0.22) / 0.05))
    weights = 1.0 / probability
    augmented = np.concatenate((
        sample,
        2.0 * grid[0] - sample,
        2.0 * grid[-1] - sample,
    ))
    augmented_weights = np.tile(weights / 3.0, 3)
    reference_kde = gaussian_kde(
        augmented, weights=augmented_weights, bw_method=0.8
    )
    reference = reference_kde(grid)
    reference /= trapezoid(reference, grid)
    actual = estimate_host_density(
        grid, sample, inclusion_probability=probability,
        bandwidth_factor=0.8,
    )
    np.testing.assert_allclose(actual, reference, rtol=1e-12, atol=1e-12)



def test_scott_multiplier_one_matches_weighted_scott_default():
    grid = np.linspace(0.01, 0.40, 501)
    sample = np.linspace(0.012, 0.398, 120) ** 1.0
    probability = 0.12 + 0.82 / (1.0 + np.exp((sample - 0.22) / 0.05))
    default = estimate_host_density(
        grid, sample, inclusion_probability=probability
    )
    explicit = estimate_host_density(
        grid, sample, inclusion_probability=probability,
        scott_multiplier=1.0,
    )
    np.testing.assert_allclose(explicit, default, rtol=1e-12, atol=1e-12)


def test_scott_multiplier_uses_original_weighted_sample_effective_size():
    grid = np.linspace(0.01, 0.40, 501)
    sample = np.linspace(0.012, 0.398, 120)
    probability = 0.12 + 0.82 / (1.0 + np.exp((sample - 0.22) / 0.05))
    weights = 1.0 / probability
    scaled = weights / np.max(weights)
    neff = np.sum(scaled) ** 2 / np.sum(scaled ** 2)
    augmented = np.concatenate((
        sample, 2.0 * grid[0] - sample, 2.0 * grid[-1] - sample,
    ))
    augmented_weights = np.tile(scaled / 3.0, 3)
    original_probability = scaled / np.sum(scaled)
    original_mean = float(original_probability @ sample)
    original_denom = 1.0 - float(original_probability @ original_probability)
    original_variance = float(
        np.sum(original_probability * (sample - original_mean) ** 2)
        / original_denom
    )
    augmented_probability = augmented_weights / np.sum(augmented_weights)
    augmented_mean = float(augmented_probability @ augmented)
    augmented_denom = 1.0 - float(augmented_probability @ augmented_probability)
    augmented_variance = float(
        np.sum(augmented_probability * (augmented - augmented_mean) ** 2)
        / augmented_denom
    )
    expected_factor = (
        0.5 * neff ** (-1.0 / 5.0)
        * np.sqrt(original_variance / augmented_variance)
    )
    expected_kde = gaussian_kde(
        augmented, weights=augmented_weights, bw_method=expected_factor,
    )
    expected = expected_kde(grid)
    expected /= trapezoid(expected, grid)
    actual = estimate_host_density(
        grid, sample, inclusion_probability=probability,
        scott_multiplier=0.5,
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)

    # Verify the physical kernel SD, not merely scipy's covariance factor.
    assert np.sqrt(expected_kde.covariance[0, 0]) == pytest.approx(
        0.5 * neff ** (-1.0 / 5.0) * np.sqrt(original_variance),
        rel=1e-12,
    )


def test_absolute_bandwidth_and_scott_multiplier_are_mutually_exclusive():
    grid = np.linspace(0.01, 0.40, 101)
    sample = np.linspace(0.02, 0.39, 20)
    with pytest.raises(ValueError, match="mutually exclusive"):
        estimate_host_density(
            grid, sample, bandwidth_factor=0.8, scott_multiplier=0.8
        )


@pytest.mark.parametrize("multiplier", [0.0, -0.1, np.nan, np.inf])
def test_invalid_scott_multiplier_is_rejected(multiplier):
    grid = np.linspace(0.01, 0.40, 101)
    sample = np.linspace(0.02, 0.39, 20)
    with pytest.raises(ValueError, match="scott_multiplier"):
        estimate_host_density(grid, sample, scott_multiplier=multiplier)
