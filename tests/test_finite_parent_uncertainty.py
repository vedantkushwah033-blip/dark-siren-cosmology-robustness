from __future__ import annotations

import numpy as np
import pytest

from dark_siren_robustness.finite_parent_uncertainty import (
    bootstrap_standard_error,
    finite_population_correction,
    paired_bootstrap_indices,
)


@pytest.mark.parametrize(
    ("n", "population", "expected"),
    [(1, 4, np.sqrt(0.75)), (2, 4, np.sqrt(0.5)), (4, 4, 0.0)],
)
def test_fpc_matches_srswor_formula(n, population, expected):
    assert finite_population_correction(n, population) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("n", "population"),
    [(0, 4), (5, 4), (1, 0), (True, 4), (1.5, 4)],
)
def test_fpc_rejects_invalid_design(n, population):
    with pytest.raises((TypeError, ValueError)):
        finite_population_correction(n, population)


def test_paired_bootstrap_indices_are_reproducible_and_allow_duplicates():
    first = paired_bootstrap_indices(np.random.default_rng(812), 100, population_size=20)
    second = paired_bootstrap_indices(np.random.default_rng(812), 100, population_size=20)
    np.testing.assert_array_equal(first, second)
    assert first.shape == (100,)
    assert np.all((first >= 0) & (first < 20))
    assert np.unique(first).size < first.size


def test_same_bootstrap_indices_preserve_pairing_between_validation_fields():
    redshift = np.linspace(0.01, 0.4, 20)
    retained = (np.arange(20) % 3 == 0).astype(int)
    indices = paired_bootstrap_indices(np.random.default_rng(19), 50, population_size=20)
    sampled_z = redshift[indices]
    sampled_retained = retained[indices]
    np.testing.assert_array_equal(sampled_retained, (np.rint((sampled_z - 0.01) / 0.39 * 19).astype(int) % 3 == 0).astype(int))


def test_bootstrap_standard_error_uses_sample_standard_deviation():
    values = np.array([1.0, 2.0, 3.0, 4.0])
    assert bootstrap_standard_error(values) == pytest.approx(np.std(values, ddof=1))


def test_bootstrap_standard_error_rejects_invalid_replicates():
    with pytest.raises(ValueError):
        bootstrap_standard_error(np.array([1.0]))
    with pytest.raises(ValueError):
        bootstrap_standard_error(np.array([1.0, np.nan]))
