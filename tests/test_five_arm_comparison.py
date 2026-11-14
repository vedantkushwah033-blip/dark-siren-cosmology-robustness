import numpy as np
import pytest

from dark_siren_robustness.five_arm_comparison import (
    build_five_arm_weights,
    weight_concentration_diagnostics,
)
from dark_siren_robustness.generative_selection import catalogue_inclusion_probability


def _fixture():
    rng = np.random.default_rng(20261113)
    parent = np.sort(rng.uniform(0.01, 0.40, 600))
    probability = catalogue_inclusion_probability(parent, floor=0.10)
    catalogue = parent[rng.random(parent.size) < probability]
    validation = np.sort(rng.uniform(0.01, 0.40, 240))
    retained = rng.binomial(
        1, catalogue_inclusion_probability(validation, floor=0.10)
    )
    return parent, catalogue, validation, retained


def test_five_arm_weights_have_matched_shapes_and_normalized_quadrature():
    parent, catalogue, validation, retained = _fixture()
    result = build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.10, quadrature_order=7,
    )
    assert result.arm_a_catalogue_only.shape == catalogue.shape
    assert result.arm_b_parent_oracle.shape == parent.shape
    assert result.arm_c_known_probability.shape == catalogue.shape
    assert result.arm_d_map_plugin.shape == catalogue.shape
    assert result.arm_e_node_weights.shape == (
        result.inclusion_posterior.nodes.shape[0], catalogue.size
    )
    assert result.arm_e_node_inclusion_probability.shape == result.arm_e_node_weights.shape
    assert np.isclose(result.arm_e_quadrature_weights.sum(), 1.0)
    assert np.all(result.arm_c_known_probability >= 1.0)
    assert np.all(result.arm_d_map_plugin >= 1.0)
    assert np.all(result.arm_e_node_weights >= 1.0)


def test_arm_a_and_b_are_equal_weight_oracles_on_their_respective_samples():
    parent, catalogue, validation, retained = _fixture()
    result = build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.10, quadrature_order=5,
    )
    assert np.array_equal(result.arm_a_catalogue_only, np.ones(catalogue.size))
    assert np.array_equal(result.arm_b_parent_oracle, np.ones(parent.size))


def test_weight_diagnostics_are_finite_and_report_concentration():
    result = weight_concentration_diagnostics([1.0, 1.0, 2.0, 6.0])
    assert result["kish_ess"] == pytest.approx(100.0 / 42.0)
    assert result["kish_ess_fraction"] == pytest.approx((100.0 / 42.0) / 4.0)
    assert result["maximum_normalized_weight"] == pytest.approx(0.6)
    assert result["top_1_percent_mass"] == pytest.approx(0.6)
    assert result["top_5_percent_mass"] == pytest.approx(0.6)


def test_weight_diagnostics_reject_nonpositive_or_nonfinite_weights():
    with pytest.raises(ValueError):
        weight_concentration_diagnostics([1.0, 0.0])
    with pytest.raises(ValueError):
        weight_concentration_diagnostics([1.0, np.inf])


def test_five_arm_builder_rejects_mismatched_validation_labels():
    parent, catalogue, validation, retained = _fixture()
    with pytest.raises(ValueError, match="match validation_redshift"):
        build_five_arm_weights(
            parent, catalogue, validation, retained[:-1],
            catalogue_floor=0.10,
        )

def test_weight_diagnostics_reject_positive_mass_underflow():
    with pytest.raises(ArithmeticError, match="underflowed positive mass"):
        weight_concentration_diagnostics([1e-320, 1e308])


def test_five_arm_weights_apply_common_base_population_mass_to_all_arms():
    parent, catalogue, validation, retained = _fixture()
    parent_base = 1.0 + 1.5 * (parent - 0.01) / 0.39
    catalogue_base = 1.0 + 1.5 * (catalogue - 0.01) / 0.39
    result = build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.10,
        parent_base_weights=parent_base,
        catalogue_base_weights=catalogue_base,
        quadrature_order=5,
    )
    true_probability = catalogue_inclusion_probability(catalogue, floor=0.10)
    assert np.allclose(result.arm_a_catalogue_only, catalogue_base)
    assert np.allclose(result.arm_b_parent_oracle, parent_base)
    assert np.allclose(result.arm_c_known_probability, catalogue_base / true_probability)
    assert np.allclose(
        result.arm_d_map_plugin,
        catalogue_base / result.arm_d_map_inclusion_probability,
    )
    assert np.allclose(
        result.arm_e_node_weights,
        catalogue_base[None, :] / result.arm_e_node_inclusion_probability,
    )


@pytest.mark.parametrize(
    "parent_base,catalogue_base,match",
    [
        (np.ones(599), None, "parent_base_weights"),
        (None, np.ones(2), "catalogue_base_weights"),
        (np.r_[np.ones(599), np.nan], None, "parent_base_weights"),
        (None, np.r_[np.ones(2), 0.0], "catalogue_base_weights"),
    ],
)
def test_five_arm_builder_rejects_invalid_base_population_weights(
    parent_base, catalogue_base, match
):
    parent, catalogue, validation, retained = _fixture()
    with pytest.raises(ValueError, match=match):
        build_five_arm_weights(
            parent, catalogue, validation, retained,
            catalogue_floor=0.10,
            parent_base_weights=parent_base,
            catalogue_base_weights=catalogue_base,
            quadrature_order=5,
        )


def test_five_arm_builder_accepts_external_known_inclusion_probabilities():
    parent, catalogue, validation, retained = _fixture()
    catalogue_base = 1.0 + 1.5 * (catalogue - 0.01) / 0.39
    supplied = np.linspace(0.2, 0.8, catalogue.size)
    result = build_five_arm_weights(
        parent, catalogue, validation, retained,
        catalogue_floor=0.10,
        catalogue_base_weights=catalogue_base,
        known_inclusion_probability=supplied,
        quadrature_order=5,
    )
    assert np.array_equal(result.arm_c_true_inclusion_probability, supplied)
    assert np.allclose(result.arm_c_known_probability, catalogue_base / supplied)


@pytest.mark.parametrize(
    "supplied,match",
    [
        (np.ones(2), "match catalogue_redshift"),
        (np.zeros(600), "in \(0, 1\]"),
        (np.full(600, 1.01), "in \(0, 1\]"),
        (np.full(600, np.nan), "finite"),
    ],
)
def test_five_arm_builder_rejects_invalid_external_known_probabilities(supplied, match):
    parent, catalogue, validation, retained = _fixture()
    if supplied.size == 600 and supplied.size != catalogue.size:
        # Force the mismatch case to remain a mismatch even if the fixture changes.
        supplied = np.ones(catalogue.size + 1)
    with pytest.raises(ValueError, match=match):
        build_five_arm_weights(
            parent, catalogue, validation, retained,
            catalogue_floor=0.10,
            known_inclusion_probability=supplied,
            quadrature_order=5,
        )
