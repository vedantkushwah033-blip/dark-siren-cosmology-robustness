from __future__ import annotations

import json
import subprocess
import sys

import numpy as np
import pandas as pd

from dark_siren_robustness.estimated_inclusion_study import run_estimated_inclusion_pilot


def test_persisted_validation_bootstrap_artifact_round_trips_and_audits(tmp_path):
    frame, _ = run_estimated_inclusion_pilot(
        n_replicates=1, n_parent=1200, n_detected=6,
        validation_sizes=(200,), seed=20261081,
        h0_points=41, z_points=201,
    )
    row = frame.loc[
        frame["arm"].eq("estimated_200") & frame["status"].eq("success")
    ]
    assert len(row) == 1, frame[["arm", "status", "error"]].to_dict("records")
    pilot_csv = tmp_path / "pilot.csv"
    frame.to_csv(pilot_csv, index=False)
    artifact_dir = tmp_path / "bootstrap"

    run = subprocess.run(
        [
            sys.executable, "scripts/run_validation_bootstrap_sensitivity.py",
            "--pilot-csv", str(pilot_csv), "--replicate", "0",
            "--validation-size", "200", "--n-bootstrap", "12",
            "--seed", "20261082", "--output-dir", str(artifact_dir),
        ],
        check=False, capture_output=True, text=True,
    )
    assert run.returncode == 0, run.stdout + run.stderr
    summary = json.loads((artifact_dir / "summary.json").read_text())
    assert summary["n_attempted"] == 12
    assert summary["n_successful"] + summary["n_failed"] == 12
    assert summary["finite_population_correction_applied"] is False
    with np.load(artifact_dir / "bootstrap_draws.npz", allow_pickle=False) as draws:
        assert draws["density_draws"].shape[0] == summary["n_successful"]

    audit = subprocess.run(
        [
            sys.executable, "scripts/audit_validation_bootstrap_sensitivity.py",
            "--artifact-dir", str(artifact_dir), "--pilot-csv", str(pilot_csv),
        ],
        check=False, capture_output=True, text=True,
    )
    assert audit.returncode == 0, audit.stdout + audit.stderr
    assert json.loads(audit.stdout)["status"] == "PASS"
