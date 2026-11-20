"""Tests for deterministic prior-predictive H0 truth scheduling."""
from __future__ import annotations

import numpy as np
import pytest

from scripts.run_quadratic_five_arm_replicates import _h0_truth


def test_h0_truth_schedule_is_deterministic_and_in_support():
    first = [_h0_truth(20261121, rep, 30.0, 110.0) for rep in range(20)]
    second = [_h0_truth(20261121, rep, 30.0, 110.0) for rep in range(20)]
    assert first == second
    assert all(30.0 <= truth < 110.0 for _, truth in first)
    assert len({seed for seed, _ in first}) == 20
    assert len({truth for _, truth in first}) == 20


def test_h0_truth_changes_with_master_seed():
    assert _h0_truth(1, 4, 30.0, 110.0) != _h0_truth(2, 4, 30.0, 110.0)


@pytest.mark.parametrize("lower,upper", [(30.0, 30.0), (110.0, 30.0), (float("nan"), 90.0)])
def test_uniform_truth_schedule_rejects_invalid_bounds(lower, upper):
    with pytest.raises(ValueError):
        _h0_truth(1, 0, lower, upper)
