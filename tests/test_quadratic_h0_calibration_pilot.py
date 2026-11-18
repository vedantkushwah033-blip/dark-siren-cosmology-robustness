import json

import pandas as pd

from scripts.run_quadratic_h0_calibration_pilot import run_pilot
from argparse import Namespace


def test_quadratic_h0_calibration_pilot_writes_paired_auditable_outputs(tmp_path):
    args = Namespace(
        replicates=1,
        seed=901,
        gamma_values=(0.0,),
        validation_sizes=(50,),
        n_parent=300,
        n_validation_frame=100,
        n_events=2,
        beta=-3.9,
        h0_true=70.0,
        sigma_mpc=140.0,
        threshold=250.0,
        z_min=0.01,
        z_max=0.40,
        h0_min=55.0,
        h0_max=85.0,
        h0_points=31,
        output_dir=tmp_path / "pilot",
    )
    result = run_pilot(args)
    assert result["planned_metric_rows"] == 4
    assert result["actual_metric_rows"] == 4
    assert result["success_rows"] == 4
    assert result["failed_rows"] == 0

    metrics = pd.read_csv(args.output_dir / "replicate_metrics.csv")
    assert set(metrics["arm"]) == {
        "degree_1", "degree_2", "catalogue_only", "oracle_true_inclusion"
    }
    assert set(metrics["status"]) == {"SUCCESS"}
    assert metrics["pit"].between(0.0, 1.0).all()
    assert metrics["area_before_normalization"].between(0.99, 1.01).all()

    config = json.loads((args.output_dir / "config.json").read_text())
    assert config["gamma_values"] == [0.0]
    assert config["validation_sizes"] == [50]
    manifest = json.loads((args.output_dir / "manifest.json").read_text())
    assert "replicate_metrics.csv" in manifest
    assert "posterior_curves.csv" in manifest
    assert "failure_ledger.csv" in manifest
