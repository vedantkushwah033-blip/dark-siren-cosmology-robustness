import numpy as np
from scipy.integrate import trapezoid
from scipy.special import log_ndtr
from scipy.stats import norm

from dark_siren_robustness.five_arm_comparison import build_five_arm_weights
from dark_siren_robustness.five_arm_kde_likelihood import (
    _batched_density_likelihood,
    matched_five_arm_kde_selected_log_likelihood,
)
from dark_siren_robustness.dataset_likelihood import C_KM_S


def _reference(h0, z, density, events, threshold):
    means = C_KM_S * z[None, :] / h0[:, None]
    log_num = []
    log_eff = []
    for event in events:
        numerator = density[None, :] * norm.pdf(
            float(event["gw_distance"]), loc=means, scale=float(event["gw_sigma"])
        )
        efficiency = density[None, :] * np.exp(
            log_ndtr((means - threshold) / float(event["gw_sigma"]))
        )
        log_num.append(np.log(trapezoid(numerator, z, axis=1)))
        log_eff.append(np.log(trapezoid(efficiency, z, axis=1)))
    log_num = np.asarray(log_num)
    log_eff = np.asarray(log_eff)
    return np.sum(log_num - log_eff, axis=0), log_eff


def test_batched_kernel_reuse_matches_scalar_trapezoid_reference_for_all_arms():
    parent = np.linspace(0.0101, 0.3999, 420)
    catalogue = parent[::2]
    validation = np.linspace(0.0101, 0.3999, 150)
    retained = (np.arange(validation.size) % 4 != 0).astype(int)
    base_parent = 1.0 + 1.5 * (parent - 0.01) / 0.39
    base_catalogue = 1.0 + 1.5 * (catalogue - 0.01) / 0.39
    arms = build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.1, quadrature_order=5,
        parent_base_weights=base_parent,
        catalogue_base_weights=base_catalogue,
    )
    events = [
        {"gw_distance": 275.0, "gw_sigma": 140.0},
        {"gw_distance": 340.0, "gw_sigma": 125.0},
        {"gw_distance": 410.0, "gw_sigma": 155.0},
    ]
    z = np.linspace(0.01, 0.40, 241)
    h0 = np.linspace(55.0, 85.0, 37)
    threshold = 250.0
    result = matched_five_arm_kde_selected_log_likelihood(
        h0, arms, events, detection_threshold_mpc=threshold,
        redshift_grid=z, merger_weight=np.ones_like(z), scott_multiplier=0.7,
    )

    densities = [
        result.arm_density["A_catalogue_only"],
        result.arm_density["B_parent_oracle"],
        result.arm_density["C_known_probability_ipw"],
        result.arm_density["D_map_plugin_ipw"],
        *list(result.arm_e_node_density),
    ]
    actual = [
        result.arm_a_log_likelihood,
        result.arm_b_log_likelihood,
        result.arm_c_log_likelihood,
        result.arm_d_log_likelihood,
        *list(result.arm_e_node_log_likelihood),
    ]
    actual_eff = [
        *[None] * 4,
        *list(result.arm_e_node_log_detection_efficiency),
    ]
    for index, density in enumerate(densities):
        expected_ll, expected_eff = _reference(h0, z, density, events, threshold)
        np.testing.assert_allclose(actual[index], expected_ll, rtol=2e-12, atol=2e-12)
        if index >= 4:
            np.testing.assert_allclose(actual_eff[index], expected_eff, rtol=2e-12, atol=2e-12)



def test_five_arm_likelihood_accepts_one_common_absolute_physical_bandwidth():
    from dark_siren_robustness.catalogue_density import reflected_weighted_kde_bandwidth

    parent = np.linspace(0.0101, 0.3999, 420)
    catalogue = parent[::2]
    validation = np.linspace(0.0101, 0.3999, 150)
    retained = (np.arange(validation.size) % 4 != 0).astype(int)
    base_parent = 1.0 + 1.5 * (parent - 0.01) / 0.39
    base_catalogue = 1.0 + 1.5 * (catalogue - 0.01) / 0.39
    arms = build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.1, quadrature_order=5,
        parent_base_weights=base_parent,
        catalogue_base_weights=base_catalogue,
    )
    requested = 0.012
    result = matched_five_arm_kde_selected_log_likelihood(
        np.linspace(55.0, 85.0, 37), arms,
        [{"gw_distance": 275.0, "gw_sigma": 140.0}],
        detection_threshold_mpc=250.0,
        redshift_grid=np.linspace(0.01, 0.40, 241),
        merger_weight=np.ones(241),
        absolute_bandwidth_sd=requested,
    )
    samples_weights = [
        (catalogue, arms.arm_a_catalogue_only),
        (parent, arms.arm_b_parent_oracle),
        (catalogue, arms.arm_c_known_probability),
        (catalogue, arms.arm_d_map_plugin),
        *[(catalogue, weights) for weights in arms.arm_e_node_weights],
    ]
    for sample, weights in samples_weights:
        assert reflected_weighted_kde_bandwidth(
            sample, lower_bound=0.01, upper_bound=0.40,
            sample_weights=weights, absolute_bandwidth_sd=requested,
        ) == requested
    assert result.arm_density["A_catalogue_only"].shape == (241,)
    assert result.arm_e_node_density.shape[1] == 241



def test_batched_kernel_h0_chunk_sizes_match_scalar_reference_including_tail_chunk():
    z = np.linspace(0.01, 0.40, 241)
    h0 = np.linspace(55.0, 85.0, 37)  # not divisible by the smaller chunk sizes
    raw = np.stack([
        0.2 + np.exp(-0.5 * ((z - center) / width) ** 2)
        for center, width in ((0.08, 0.025), (0.19, 0.04), (0.31, 0.03))
    ])
    densities = raw / np.asarray([trapezoid(row, z) for row in raw])[:, None]
    events = [
        (275.0, 140.0),
        (410.0, 155.0),
        (720.0, 190.0),
    ]
    distances = np.asarray([item[0] for item in events])
    sigmas = np.asarray([item[1] for item in events])
    threshold = 250.0

    expected_ll = []
    expected_eff = []
    for density in densities:
        ll, eff = _reference(
            h0, z, density,
            [{"gw_distance": d, "gw_sigma": s} for d, s in events],
            threshold,
        )
        expected_ll.append(ll)
        expected_eff.append(eff)
    expected_ll = np.asarray(expected_ll)
    expected_eff = np.asarray(expected_eff)

    for chunk_size in (1, 7, 16, 128):
        actual_ll, actual_eff = _batched_density_likelihood(
            h0, z, densities, distances, sigmas, threshold,
            h0_chunk_size=chunk_size,
        )
        np.testing.assert_allclose(actual_ll, expected_ll, rtol=2e-12, atol=2e-12)
        np.testing.assert_allclose(actual_eff, expected_eff, rtol=2e-12, atol=2e-12)


def test_batched_kernel_rejects_nonpositive_h0_chunk_size():
    z = np.linspace(0.01, 0.40, 101)
    density = np.ones((1, z.size))
    density /= trapezoid(density[0], z)
    with pytest.raises(ValueError, match="h0_chunk_size"):
        _batched_density_likelihood(
            np.linspace(55.0, 85.0, 21), z, density,
            np.asarray([275.0]), np.asarray([140.0]), 250.0,
            h0_chunk_size=0,
        )
