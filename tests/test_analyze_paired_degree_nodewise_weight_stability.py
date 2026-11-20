import numpy as np
import pytest

from scripts.analyze_paired_degree_nodewise_weight_stability import (
    compute_nodewise_diagnostics,
    weighted_quantile,
)


def test_weighted_quantile_uses_posterior_mass():
    values = np.asarray([1.0, 2.0, 3.0])
    weights = np.asarray([0.1, 0.7, 0.2])
    assert weighted_quantile(values, weights, 0.05) == 1.0
    assert weighted_quantile(values, weights, 0.50) == 2.0
    assert weighted_quantile(values, weights, 0.95) == 3.0


def test_nodewise_diagnostics_reconstructs_extrema_and_weighted_summaries():
    z = np.linspace(0.01, 0.40, 40)
    nodes = np.asarray([[0.0, 0.0], [-0.2, 0.8], [0.1, -0.5]])
    posterior_weights = np.asarray([0.7, 0.2, 0.1])
    result = compute_nodewise_diagnostics(
        z,
        nodes,
        posterior_weights,
        degree=1,
        coordinate_center=0.205,
        coordinate_scale=0.195,
    )

    x = (z - 0.205) / 0.195
    design = np.column_stack([np.ones_like(x), x])
    probability = 1.0 / (1.0 + np.exp(-(design @ nodes.T)))
    base = 1.0 + 1.5 * (z - 0.01) / 0.39
    raw = base[:, None] / probability
    normalized = raw / raw.sum(axis=0, keepdims=True)
    ess = 1.0 / np.sum(normalized**2, axis=0)
    top5 = np.sort(normalized, axis=0)[-2:, :].sum(axis=0)

    assert result["minimum_node_kish_ess"] == pytest.approx(ess.min())
    assert result["maximum_node_top_5_percent_mass"] == pytest.approx(top5.max())
    assert result["posterior_weighted_q50_node_kish_ess"] == pytest.approx(
        weighted_quantile(ess, posterior_weights, 0.50)
    )
    assert result["posterior_weighted_q90_node_top_5_percent_mass"] == pytest.approx(
        weighted_quantile(top5, posterior_weights, 0.90)
    )
    assert result["quadrature_ess"] == pytest.approx(
        1.0 / np.sum((posterior_weights / posterior_weights.sum()) ** 2)
    )


def test_nodewise_diagnostics_rejects_invalid_posterior_weights():
    with pytest.raises(ValueError, match="non-negative"):
        compute_nodewise_diagnostics(
            np.linspace(0.01, 0.40, 10),
            np.asarray([[0.0, 0.0], [0.1, 0.2]]),
            np.asarray([1.1, -0.1]),
            degree=1,
            coordinate_center=0.205,
            coordinate_scale=0.195,
        )


def test_weighted_quantile_rejects_invalid_probability():
    with pytest.raises(ValueError, match="q must lie"):
        weighted_quantile(np.asarray([1.0, 2.0]), np.asarray([0.5, 0.5]), 1.1)
