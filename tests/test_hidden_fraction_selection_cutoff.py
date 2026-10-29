from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

from scripts.run_hidden_fraction_selection_cutoff import (
    D_MAX_LEVELS,
    F_MODELS,
    generate_selected_dataset,
    normalized_gaussian_density,
    stable_seed,
)


def test_stable_seed_does_not_depend_on_mapping_order():
    a = stable_seed(123, {"cutoff": 500.0, "replicate": 2})
    b = stable_seed(123, {"replicate": 2, "cutoff": 500.0})
    c = stable_seed(123, {"cutoff": 1000.0, "replicate": 2})
    assert a == b
    assert a != c


def test_selected_data_generator_respects_selection_and_has_audit_labels():
    z = np.linspace(0.01, 0.40, 151)
    gc = normalized_gaussian_density(z, 0.12, 0.05)
    gm = normalized_gaussian_density(z, 0.25, 0.07)
    rows = generate_selected_dataset(
        z=z, g_catalogue=gc, g_missing=gm, h0_true=70.0, f_true=0.3,
        sigma=140.0, d_min=0.0, d_max=500.0, n_events=12, seed=77881,
    )
    assert len(rows) == 12
    assert all(0.0 <= row["distance_mpc"] <= 500.0 for row in rows)
    assert all(0.01 <= row["z_true_for_audit_only"] <= 0.40 for row in rows)
    assert {row["component_for_audit_only"] for row in rows} <= {"catalogue", "missing"}


def test_small_pilot_and_independent_artifact_audit(tmp_path):
    root = Path(__file__).resolve().parents[1]
    output = tmp_path / "pilot"
    command = [
        sys.executable, str(root / "scripts" / "run_hidden_fraction_selection_cutoff.py"),
        "--replicates", "1", "--n-events", "3", "--h0-points", "31",
        "--z-points", "101", "--seed", "20261029", "--output-dir", str(output),
    ]
    subprocess.run(command, check=True, capture_output=True, text=True, cwd=root)
    audit = subprocess.run(
        [sys.executable, str(root / "scripts" / "audit_hidden_fraction_selection_cutoff.py"), str(output)],
        check=True, capture_output=True, text=True, cwd=root,
    )
    result = json.loads(audit.stdout)
    assert result["status"] == "PASS"
    assert result["events_reconstructed_structure"] == len(D_MAX_LEVELS) * 3
    assert result["posterior_grids_independently_reconstructed"] == len(D_MAX_LEVELS) * len(F_MODELS)
    assert result["expected_posterior_grids"] == len(D_MAX_LEVELS) * len(F_MODELS)


def test_smooth_selection_misspecification_pilot_and_independent_audit(tmp_path):
    root = Path(__file__).resolve().parents[1]
    output = tmp_path / "smooth-pilot"
    subprocess.run(
        [
            sys.executable, str(root / "scripts" / "run_smooth_selection_misspecification.py"),
            "--replicates", "1", "--n-events", "2", "--h0-points", "31",
            "--z-points", "101", "--seed", "20261029", "--output-dir", str(output),
        ],
        check=True, capture_output=True, text=True, cwd=root,
    )
    audit = subprocess.run(
        [sys.executable, str(root / "scripts" / "audit_smooth_selection_misspecification.py"), str(output)],
        check=True, capture_output=True, text=True, cwd=root,
    )
    result = json.loads(audit.stdout)
    assert result["status"] == "PASS"
    assert result["events_structurally_checked"] == 3 * 2
    assert result["posterior_fits_independently_reconstructed"] == 3 * 3 * 3
    assert result["expected_posterior_fits"] == 3 * 3 * 3
