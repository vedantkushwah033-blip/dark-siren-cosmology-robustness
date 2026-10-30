from __future__ import annotations

import numpy as np
import pytest

from dark_siren_robustness.design_aware_srswor_bootstrap import (
    design_aware_srswor_bootstrap_indices,
)


def test_design_bootstrap_replays_and_resamples_without_replacement_from_pseudopop():
    a = design_aware_srswor_bootstrap_indices(np.random.default_rng(441), 40, 400)
    b = design_aware_srswor_bootstrap_indices(np.random.default_rng(441), 40, 400)
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_array_equal(a[1], b[1])
    pseudo, sample = a
    assert pseudo.shape == (400,)
    assert sample.shape == (40,)
    assert np.all((pseudo >= 0) & (pseudo < 40))
    assert np.all((sample >= 0) & (sample < 40))
    # Sample positions are unique in the pseudo-population; mapped records may
    # repeat because a pseudo-population can contain replicated source units.
    assert np.all(np.isin(sample, pseudo))


def test_design_bootstrap_preserves_paired_record_fields():
    z = np.linspace(0.01, 0.4, 40)
    retained = (np.arange(40) % 4 == 0).astype(int)
    _, idx = design_aware_srswor_bootstrap_indices(np.random.default_rng(12), 40, 400)
    sampled_z, sampled_retained = z[idx], retained[idx]
    recovered = np.rint((sampled_z - 0.01) / 0.39 * 39).astype(int)
    np.testing.assert_array_equal(sampled_retained, (recovered % 4 == 0).astype(int))


@pytest.mark.parametrize(
    "args",
    [(1, 10, None), (10, 9, None), (10, 20, 21), (10, 20, 0), (True, 20, None)],
)
def test_design_bootstrap_rejects_invalid_sizes(args):
    with pytest.raises((TypeError, ValueError)):
        design_aware_srswor_bootstrap_indices(np.random.default_rng(1), *args)


def test_design_bootstrap_accepts_sample_size_equal_to_population():
    pseudo, sample = design_aware_srswor_bootstrap_indices(
        np.random.default_rng(8), 10, 10, 10
    )
    assert pseudo.shape == sample.shape == (10,)
