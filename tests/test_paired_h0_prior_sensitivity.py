from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_paired_h0_prior_sensitivity import interval, validate_paired_summary_rows


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


KEYS = ["replicate", "gamma", "validation_size", "arm"]


def summary_frame(arm="A", truth=70.0, median=71.0):
    return pd.DataFrame([{
        "replicate": 0, "gamma": -1.0, "validation_size": 100, "arm": arm,
        "h0_true": truth, "posterior_median": median, "median_error": median-truth,
        "q025": 60.0, "q05": 62.0, "q16": 67.0, "q84": 75.0,
        "q95": 78.0, "q975": 80.0,
    }])


def test_pairing_requires_complete_identical_key_sets():
    narrow = pd.concat([summary_frame(arm="A"), summary_frame(arm="B")], ignore_index=True)
    wide = narrow.copy()
    paired = validate_paired_summary_rows(narrow, wide, KEYS, expected=2)
    assert len(paired) == 2


def test_pairing_rejects_missing_rows():
    narrow = summary_frame(arm="A")
    wide = pd.concat([narrow, summary_frame(arm="B")], ignore_index=True)
    with pytest.raises(AssertionError, match="row count"):
        validate_paired_summary_rows(narrow, wide, KEYS, expected=2)


def test_pairing_rejects_duplicate_keys():
    narrow = pd.concat([summary_frame(), summary_frame()], ignore_index=True)
    with pytest.raises(AssertionError, match="duplicate"):
        validate_paired_summary_rows(narrow, narrow.copy(), KEYS, expected=2)


def test_pairing_rejects_truth_mismatch():
    narrow = summary_frame(truth=70.0)
    wide = summary_frame(truth=71.0)
    with pytest.raises(AssertionError, match="H0 truth"):
        validate_paired_summary_rows(narrow, wide, KEYS, expected=1)


def test_pairing_rejects_inconsistent_median_error():
    narrow = summary_frame()
    wide = summary_frame()
    wide.loc[0, "median_error"] = 9.0
    with pytest.raises(AssertionError, match="median_error"):
        validate_paired_summary_rows(narrow, wide, KEYS, expected=1)


def test_pairing_rejects_nonmonotone_quantiles():
    narrow = summary_frame()
    wide = summary_frame()
    wide.loc[0, "q16"] = 76.0
    with pytest.raises(AssertionError, match="monotonically ordered"):
        validate_paired_summary_rows(narrow, wide, KEYS, expected=1)
