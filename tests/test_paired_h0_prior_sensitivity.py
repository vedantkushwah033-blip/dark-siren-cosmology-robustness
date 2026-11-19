from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_paired_h0_prior_sensitivity import interval


def test_paired_interval_matches_student_t_mean_interval():
    values = np.array([-1.0, 0.0, 1.0, 2.0, 3.0])
    mean, se, low, high = interval(values)

    expected_mean = float(np.mean(values))
    expected_se = float(np.std(values, ddof=1) / np.sqrt(len(values)))
    half_width = float(t.ppf(0.975, len(values) - 1) * expected_se)

    assert np.isclose(mean, expected_mean, rtol=0, atol=1e-14)
    assert np.isclose(se, expected_se, rtol=0, atol=1e-14)
    assert np.isclose(low, expected_mean - half_width, rtol=0, atol=1e-14)
    assert np.isclose(high, expected_mean + half_width, rtol=0, atol=1e-14)


def test_paired_interval_handles_zero_monte_carlo_variance():
    mean, se, low, high = interval(np.full(100, 2.5))

    assert mean == 2.5
    assert se == 0.0
    assert low == 2.5
    assert high == 2.5


def test_paired_interval_single_pair_has_undefined_uncertainty():
    mean, se, low, high = interval(np.array([4.0]))

    assert mean == 4.0
    assert np.isnan(se)
    assert np.isnan(low)
    assert np.isnan(high)


import pandas as pd
import pytest

from analyze_paired_h0_prior_sensitivity import validate_paired_summary_rows


def test_pairing_requires_complete_identical_key_sets():
    keys = ["replicate", "gamma", "validation_size", "arm"]
    narrow = pd.DataFrame([
        [0, -1.0, 100, "A", 70.0],
        [0, -1.0, 100, "B", 70.0],
    ], columns=keys + ["h0_true"])
    wide = narrow.copy()
    paired = validate_paired_summary_rows(narrow, wide, keys, expected=2)
    assert len(paired) == 2


def test_pairing_rejects_missing_rows():
    keys = ["replicate", "gamma", "validation_size", "arm"]
    narrow = pd.DataFrame([[0, -1.0, 100, "A"]], columns=keys)
    wide = pd.DataFrame([[0, -1.0, 100, "A"], [0, -1.0, 100, "B"]], columns=keys)
    with pytest.raises(AssertionError, match="row count"):
        validate_paired_summary_rows(narrow, wide, keys, expected=2)


def test_pairing_rejects_duplicate_keys():
    keys = ["replicate", "gamma", "validation_size", "arm"]
    narrow = pd.DataFrame([[0, -1.0, 100, "A"], [0, -1.0, 100, "A"]], columns=keys)
    with pytest.raises(AssertionError, match="duplicate"):
        validate_paired_summary_rows(narrow, narrow.copy(), keys, expected=2)


def test_pairing_rejects_truth_mismatch():
    keys = ["replicate", "gamma", "validation_size", "arm"]
    narrow = pd.DataFrame([[0, -1.0, 100, "A", 70.0]], columns=keys + ["h0_true"])
    wide = pd.DataFrame([[0, -1.0, 100, "A", 71.0]], columns=keys + ["h0_true"])
    with pytest.raises(AssertionError, match="H0 truth"):
        validate_paired_summary_rows(narrow, wide, keys, expected=1)
