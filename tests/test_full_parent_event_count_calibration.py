import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from dark_siren_robustness.generative_selection import generate_synthetic_population
from run_full_parent_event_count_calibration import run


def test_four_event_generation_is_prefix_of_thirty_event_generation():
    short = generate_synthetic_population(
        seed=123456, n_parent=300, n_detected=4, h0_true=70.0
    )
    long = generate_synthetic_population(
        seed=123456, n_parent=300, n_detected=30, h0_true=70.0
    )
    np.testing.assert_array_equal(
        short["parent_truth"].redshift.to_numpy(),
        long["parent_truth"].redshift.to_numpy(),
    )
    np.testing.assert_array_equal(
        short["gw_event_truth"].host_object_id.to_numpy(),
        long["gw_event_truth"].host_object_id.to_numpy()[:4],
    )
    np.testing.assert_array_equal(
        short["gw_events_observed"].distance_mpc.to_numpy(),
        long["gw_events_observed"].distance_mpc.to_numpy()[:4],
    )


def test_paired_event_count_runner_returns_two_arms_per_replicate():
    raw, summary, paired = run(
        replicates=2, master_seed=2026111928, n_parent=300,
        h0_true=70.0, prior_min=20.0, prior_max=300.0, h0_points=41,
    )
    assert len(raw) == 4
    assert set(raw.event_count) == {4, 30}
    assert raw.groupby("replicate").event_count.nunique().eq(2).all()
    assert len(summary) == 2
    assert set(summary.event_count) == {4, 30}
    assert paired.loc[0, "paired_replicates"] == 2
    assert np.isfinite(summary.mean_median_bias).all()
    assert (summary.mean_width95 >= summary.mean_width68).all()
