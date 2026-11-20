from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_paired_degree_multireplicate_pilot import (
    _aggregate_paired, _paired_summary, _verify_child_pair,
)


def _summaries():
    rows = []
    for rep, e1, e2 in [(0, 2.0, 1.0), (1, 4.0, 2.0), (2, 0.0, -1.0)]:
        for degree, error in [(1, e1), (2, e2)]:
            rows.append({
                "replicate": rep, "gamma": 1.0, "validation_size": 100,
                "arm": "E_uncertainty_integrated_ipw", "median_error": error,
                "q16": 60.0 + error, "q84": 80.0 + error, "h0_true": 70.0,
            })
    return pd.DataFrame(rows)


def test_degree_contrasts_are_paired_and_use_d2_minus_d1():
    data = _summaries()
    left = data[data.replicate.isin([0, 1, 2])].copy()
    # Construct the two input frames in the same shape as the replicate runner.
    d1 = left[left.index % 2 == 0].copy()
    d2 = left[left.index % 2 == 1].copy()
    paired = _paired_summary(d1, d2)
    assert len(paired) == 3
    assert paired.median_error_difference_d2_minus_d1.tolist() == [-1.0, -2.0, -1.0]
    aggregate = _aggregate_paired(paired)
    row = aggregate.iloc[0]
    assert row.paired_replicates == 3
    assert row.degree1_bias == 2.0
    assert row.degree2_bias == pytest.approx(2.0 / 3.0)
    assert row.median_error_difference_d2_minus_d1_mean == pytest.approx(-4.0 / 3.0)
    assert row.median_error_difference_d2_minus_d1_mcse > 0


def test_degree_pair_rejects_changed_observable_inputs(tmp_path: Path):
    left, right = tmp_path / "d1", tmp_path / "d2"
    for root in (left, right):
        (root / "observable").mkdir(parents=True)
        (root / "truth").mkdir()
        (root / "results").mkdir()
        for rel in (
            "observable/catalogue_observed.csv",
            "observable/gw_events_observed.csv",
            "truth/parent_truth.csv",
            "truth/gw_event_truth.csv",
            "truth/validation_frame_truth.csv",
        ):
            (root / rel).write_text("id,value\n1,2\n")
        (root / "observable/validation_n100.csv").write_text("id,value\n1,2\n")
        np.savez(root / "results/log_likelihoods.npz", H0_grid=np.array([30., 31.]),
                 A_catalogue_only_n100=np.array([-1., -2.]),
                 B_parent_oracle_n100=np.array([-1., -2.]),
                 C_known_probability_ipw_n100=np.array([-1., -2.]))
    assert _verify_child_pair(left, right)["common_arm_likelihoods_exact"]
    (right / "observable/gw_events_observed.csv").write_text("id,value\n1,3\n")
    with pytest.raises(AssertionError):
        _verify_child_pair(left, right)


def test_degree_pair_rejects_common_arm_likelihood_mismatch(tmp_path: Path):
    left, right = tmp_path / "d1", tmp_path / "d2"
    for root in (left, right):
        (root / "observable").mkdir(parents=True)
        (root / "truth").mkdir()
        (root / "results").mkdir()
        for rel in (
            "observable/catalogue_observed.csv",
            "observable/gw_events_observed.csv",
            "truth/parent_truth.csv",
            "truth/gw_event_truth.csv",
            "truth/validation_frame_truth.csv",
        ):
            (root / rel).write_text("id,value\n1,2\n")
        (root / "observable/validation_n100.csv").write_text("id,value\n1,2\n")
    for root, value in ((left, -1.0), (right, -1.01)):
        np.savez(root / "results/log_likelihoods.npz", H0_grid=np.array([30., 31.]),
                 A_catalogue_only_n100=np.array([value, -2.]),
                 B_parent_oracle_n100=np.array([-1., -2.]),
                 C_known_probability_ipw_n100=np.array([-1., -2.]))
    with pytest.raises(AssertionError, match="common arm"):
        _verify_child_pair(left, right)
