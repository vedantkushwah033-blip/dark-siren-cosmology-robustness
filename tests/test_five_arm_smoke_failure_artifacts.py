"""Ensure raw replicate inputs survive downstream inference failures."""
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_five_arm_smoke import run


def test_raw_inputs_are_persisted_before_inclusion_model_fit_failure(tmp_path):
    output = tmp_path / "failed_replicate"
    args = SimpleNamespace(
        seed=20261116,
        validation_seed=20261117,
        n_parent=300,
        n_detected=1,
        validation_frame=50,
        validation_sizes=[20],
        catalogue_floor=0.10,
        h0_true=70.0,
        h0_points=21,
        scott_multiplier=0.7,
        quadrature_order=2,  # Deliberately invalid; fail after raw input generation.
        sigma_mpc=140.0,
        detection_threshold_mpc=250.0,
        output_dir=output,
    )
    with pytest.raises(ValueError, match="order must be an integer"):
        run(args)
    required = [
        "observable/catalogue_observed.csv",
        "observable/gw_events_observed.csv",
        "observable/validation_n20.csv",
        "truth/parent_truth.csv",
        "truth/gw_event_truth.csv",
        "truth/gw_generation_trace_truth.csv",
        "truth/validation_frame_truth.csv",
    ]
    for relative in required:
        assert (output / relative).is_file(), relative
    assert not (output / "results/posterior_summaries.csv").exists()
