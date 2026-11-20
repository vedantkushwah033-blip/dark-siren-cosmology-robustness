"""Tests for the paired degree-one/degree-two replicate pilot."""
from pathlib import Path
import sys

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_paired_degree_inclusion_multireplicate import _degree_contrasts, _paired_child_inputs


def test_cross_degree_input_pairing_accepts_exact_match_and_rejects_mismatch(tmp_path):
    left, right = tmp_path / "d1", tmp_path / "d2"
    for root in (left, right):
        for rel in (
            "observable/catalogue_observed.csv",
            "observable/gw_events_observed.csv",
            "truth/parent_truth.csv",
            "truth/gw_event_truth.csv",
            "truth/validation_frame_truth.csv",
            "observable/validation_n100.csv",
        ):
            path = root / rel
            path.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame({"id": [1, 2], "z": [0.1, 0.2]}).to_csv(path, index=False)
    assert _paired_child_inputs(left, right, [100]) == 6
    pd.DataFrame({"id": [1, 2], "z": [0.1, 0.25]}).to_csv(
        right / "observable/validation_n100.csv", index=False
    )
    with pytest.raises(AssertionError):
        _paired_child_inputs(left, right, [100])


def test_degree_contrasts_are_paired_by_replicate_gamma_size_and_arm():
    base = {
        "replicate": [0, 1], "gamma": [1.0, 1.0],
        "validation_size": [100, 100], "arm": ["D_map_plugin_ipw", "D_map_plugin_ipw"],
        "status": ["SUCCESS", "SUCCESS"], "h0_true": [70.0, 70.0],
        "median_error": [2.0, -1.0], "posterior_median": [72.0, 69.0],
        "q025": [60.0, 58.0], "q16": [66.0, 65.0],
        "q84": [78.0, 73.0], "q975": [84.0, 80.0],
    }
    one = pd.DataFrame(base)
    two = one.copy()
    two["posterior_median"] += [1.5, -0.5]
    two["median_error"] += [1.5, -0.5]
    result = _degree_contrasts(one, two)
    assert len(result) == 2
    assert result.degree2_minus_degree1_median.tolist() == pytest.approx([1.5, -0.5])
    assert result.coverage68_change_d2_minus_d1.tolist() == [0, 0]


def test_degree_contrasts_excludes_unpaired_failed_cells():
    base = {
        "replicate": [0], "gamma": [0.0], "validation_size": [100],
        "arm": ["E_uncertainty_integrated_ipw"], "status": ["SUCCESS"],
        "h0_true": [70.0], "median_error": [1.0], "posterior_median": [71.0],
        "q025": [60.0], "q16": [66.0], "q84": [76.0], "q975": [82.0],
    }
    one = pd.DataFrame(base)
    two = one.iloc[0:0].copy()
    assert _degree_contrasts(one, two).empty

# This module is included in the degree-two pilot and full test workflows.
