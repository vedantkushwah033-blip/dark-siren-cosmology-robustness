"""Deterministic known-answer and invariance tests for the ASA hidden-host model.

These tests target mathematical identities of the finite positive-truncated
Gaussian mixture. They do not establish astrophysical validity.
"""
import numpy as np
import pytest

from dark_siren_robustness.hardened_inference import posterior_h0


GRID = np.linspace(25.0, 125.0, 4001)
PRIOR = (50.0, 90.0)
DISTANCE = 120.0
SIGMA = 18.0


def fit(redshifts, weights=None, **kwargs):
    return posterior_h0(
        GRID,
        redshifts,
        DISTANCE,
        SIGMA,
        weights=weights,
        prior_min=PRIOR[0],
        prior_max=PRIOR[1],
        positive_truncated=True,
        **kwargs,
    )


def test_single_host_equals_repeated_identical_host_components():
    """Splitting a single host into identical labelled components changes nothing."""
    one = fit([0.03], [1.0])
    repeated = fit([0.03, 0.03, 0.03], [0.17, 0.31, 0.52])
    np.testing.assert_allclose(repeated, one, rtol=2e-13, atol=2e-13)


def test_permuting_host_components_and_paired_weights_changes_nothing():
    z = np.array([0.018, 0.027, 0.043, 0.061])
    w = np.array([0.11, 0.29, 0.37, 0.23])
    baseline = fit(z, w)
    permutation = np.array([2, 0, 3, 1])
    permuted = fit(z[permutation], w[permutation])
    np.testing.assert_allclose(permuted, baseline, rtol=2e-13, atol=2e-13)


def test_splitting_catalogue_weight_across_duplicate_redshift_is_invariant():
    z = np.array([0.02, 0.04])
    w = np.array([0.7, 0.3])
    baseline = fit(z, w)
    split_z = np.array([0.02, 0.02, 0.04])
    split_w = np.array([0.2, 0.5, 0.3])
    split = fit(split_z, split_w)
    np.testing.assert_allclose(split, baseline, rtol=2e-13, atol=2e-13)


def test_zero_hidden_mass_is_exact_catalogue_only_reduction():
    z = np.array([0.02, 0.04])
    w = np.array([0.7, 0.3])
    baseline = fit(z, w)
    zero_hidden = fit(
        z,
        w,
        missing_weight=0.0,
        # At zero mass, support is mathematically irrelevant and need not be supplied.
    )
    np.testing.assert_allclose(zero_hidden, baseline, rtol=1e-13, atol=1e-13)


def test_identical_hidden_and_catalogue_components_collapse_to_same_density():
    z = np.array([0.02, 0.04])
    w = np.array([0.7, 0.3])
    baseline = fit(z, w)
    with_redundant_hidden = fit(
        z,
        w,
        missing_weight=0.37,
        missing_redshifts=[0.02, 0.04],
        missing_weights=[0.7, 0.3],
    )
    np.testing.assert_allclose(
        with_redundant_hidden, baseline, rtol=2e-13, atol=2e-13
    )


def test_positive_truncation_becomes_negligible_when_all_means_are_many_sigma_positive():
    """For mu/sigma very large, Phi(mu/sigma) is numerically one."""
    grid = np.linspace(50.0, 90.0, 2001)
    z = np.array([0.20, 0.25])
    weights = np.array([0.4, 0.6])
    distance = 900.0
    sigma = 2.0
    truncated = posterior_h0(
        grid, z, distance, sigma, weights=weights,
        prior_min=50.0, prior_max=90.0, positive_truncated=True,
    )
    ordinary = posterior_h0(
        grid, z, distance, sigma, weights=weights,
        prior_min=50.0, prior_max=90.0, positive_truncated=False,
    )
    np.testing.assert_allclose(truncated, ordinary, rtol=2e-12, atol=2e-12)


def test_hidden_mixture_matches_explicit_expanded_host_list():
    grid = np.linspace(40.0, 100.0, 3001)
    z_catalogue = np.array([0.02, 0.04])
    q = np.array([0.7, 0.3])
    f = 0.28
    z_hidden = np.array([0.027, 0.052])
    r_hidden = np.array([0.35, 0.65])
    compact = posterior_h0(
        grid, z_catalogue, 145.0, 17.0,
        weights=q, missing_weight=f,
        missing_redshifts=z_hidden, missing_weights=r_hidden,
        prior_min=40.0, prior_max=100.0, positive_truncated=True,
    )
    expanded = posterior_h0(
        grid,
        np.concatenate([z_catalogue, z_hidden]),
        145.0,
        17.0,
        weights=np.concatenate([(1.0 - f) * q, f * r_hidden]),
        prior_min=40.0,
        prior_max=100.0,
        positive_truncated=True,
    )
    np.testing.assert_allclose(compact, expanded, rtol=2e-12, atol=2e-12)


@pytest.mark.parametrize("bad_fraction", [-1e-8, 1.0, 1.01, np.nan, np.inf])
def test_invalid_hidden_fraction_is_rejected(bad_fraction):
    with pytest.raises(ValueError):
        fit([0.02, 0.04], [0.7, 0.3], missing_weight=bad_fraction)
