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

    from scripts.audit_quadratic_h0_calibration_pilot import audit
    audit_result = audit(args.output_dir)
    assert audit_result["status"] == "PASS"
    assert audit_result["planned_arm_rows"] == 4
    assert audit_result["raw_cell_manifests_verified"] == 1

def test_partial_arm_failure_records_one_failed_row_per_planned_arm(tmp_path, monkeypatch):
    import numpy as np
    import scripts.run_quadratic_h0_calibration_pilot as pilot

    def fake_single(args):
        from pathlib import Path
        Path(args.output_dir).mkdir(parents=True, exist_ok=True)
        grid = np.array([55.0, 70.0, 85.0])
        curves = pd.DataFrame({
            "H0": grid,
            "density_d1_o7": [0.01, 0.08, 0.01],
            "density_d2_o7": [0.01, 0.08, 0.01],
            "density_catalogue_only": [0.01, 0.08, 0.01],
            "density_oracle_true_inclusion": [0.01, 0.08, 0.01],
        })
        curves.to_csv(str(args.output_dir) + "/h0_posterior_curves.csv", index=False)
        diagnostics = {
            "quadrature_ess": 10.0,
            "maximum_node_mass": 0.2,
            "minimum_log_inclusion_probability": -2.0,
        }
        return {
            "results": {
                "degree_1": diagnostics.copy(),
                "degree_2": diagnostics.copy(),
                "catalogue_only": {},
                "oracle_true_inclusion": {},
            }
        }

    original_metrics = pilot._posterior_metrics
    calls = {"n": 0}

    def fail_on_third_arm(grid, density, truth):
        calls["n"] += 1
        if calls["n"] == 3:
            raise ValueError("injected third-arm posterior failure")
        return original_metrics(grid, density, truth)

    monkeypatch.setattr(pilot, "run_single", fake_single)
    monkeypatch.setattr(pilot, "_posterior_metrics", fail_on_third_arm)
    args = Namespace(
        replicates=1, seed=902, gamma_values=(0.0,), validation_sizes=(50,),
        n_parent=300, n_validation_frame=100, n_events=2, beta=-3.9,
        h0_true=70.0, sigma_mpc=140.0, threshold=250.0,
        z_min=0.01, z_max=0.40, h0_min=55.0, h0_max=85.0,
        h0_points=3, output_dir=tmp_path / "atomic_failure_pilot",
    )

    result = pilot.run_pilot(args)
    metrics = pd.read_csv(args.output_dir / "replicate_metrics.csv")
    assert result["planned_metric_rows"] == 4
    assert result["actual_metric_rows"] == 4
    assert result["success_rows"] == 0
    assert result["failed_rows"] == 4
    assert len(metrics) == 4
    assert metrics["arm"].nunique() == 4
    assert set(metrics["status"]) == {"FAILED"}
    assert metrics["error"].str.contains("injected third-arm posterior failure").all()
    assert len(pd.read_csv(args.output_dir / "failure_ledger.csv")) == 4



def test_incomplete_pilot_status_fails_closed():
    import pytest
    from scripts.run_quadratic_h0_calibration_pilot import _require_complete

    with pytest.raises(SystemExit, match="3/4 succeeded"):
        _require_complete({
            "status": "INCOMPLETE",
            "planned_metric_rows": 4,
            "success_rows": 3,
            "failed_rows": 1,
        })


def test_complete_pilot_status_is_accepted():
    from scripts.run_quadratic_h0_calibration_pilot import _require_complete

    _require_complete({
        "status": "PASS",
        "planned_metric_rows": 4,
        "success_rows": 4,
        "failed_rows": 0,
    })
