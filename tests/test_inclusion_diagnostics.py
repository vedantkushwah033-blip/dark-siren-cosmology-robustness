import numpy as np
import pytest
from dark_siren_robustness.inclusion_diagnostics import (
    inverse_probability_weight_diagnostics,
    score_inclusion_predictions,
)


def test_equal_probabilities_give_equal_weights_and_full_ess():
    result = inverse_probability_weight_diagnostics(np.full(100, 0.5))
    assert result["kish_ess"] == pytest.approx(100.0)
    assert result["entropy_ess"] == pytest.approx(100.0)
    assert result["max_normalized_weight"] == pytest.approx(0.01)
    assert result["top_5pct_weight_mass"] == pytest.approx(0.05)


def test_concentrated_weights_reduce_effective_sample_size():
    p = np.r_[np.full(99, 0.9), 0.001]
    result = inverse_probability_weight_diagnostics(p)
    assert result["kish_ess"] < 10
    assert result["top_1pct_weight_mass"] > 0.5


def test_probability_scores_match_direct_calculation():
    result = score_inclusion_predictions([0, 1, 1], [0.2, 0.8, 0.6])
    assert result["brier_score"] == pytest.approx((0.04 + 0.04 + 0.16) / 3)
    assert result["log_loss"] > 0
    assert result["n"] == 3


def test_invalid_probability_inputs_are_rejected():
    with pytest.raises(ValueError):
        inverse_probability_weight_diagnostics([0.0, 0.5])
    with pytest.raises(ValueError):
        score_inclusion_predictions([0, 1], [0.0, 1.0])



def test_effective_sample_size_is_stable_for_extreme_but_representable_probabilities():
    probability = np.array([1e-250, 2e-250, 0.5, 0.9])
    result = inverse_probability_weight_diagnostics(probability)
    assert np.isfinite(result["kish_ess"])
    assert np.isfinite(result["entropy_ess"])
    assert result["kish_ess"] >= 1.0
    assert result["kish_ess"] <= probability.size
    assert result["max_normalized_weight"] == pytest.approx(
        1.0 / np.sum(np.min(probability) / probability)
    )
