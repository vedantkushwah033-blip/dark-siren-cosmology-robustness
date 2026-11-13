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
