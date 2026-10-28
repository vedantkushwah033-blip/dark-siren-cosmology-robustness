from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import trapezoid

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
