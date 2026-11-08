from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest

from dark_siren_robustness.generative_selection import (
    generate_independent_validation_sample,
)


def test_independent_validation_sampler_is_reproducible_and_separates_truth():
    kwargs = {
        "seed": 8912,
        "n_validation_frame": 900,
        "n_validation": 120,
    }
    first = generate_independent_validation_sample(**kwargs)
    second = generate_independent_validation_sample(**kwargs)
    pd.testing.assert_frame_equal(
        first["validation_observed"], second["validation_observed"]
    )
    pd.testing.assert_frame_equal(
        first["validation_frame_truth"], second["validation_frame_truth"]
    )
    assert list(first["validation_observed"].columns) == [
        "validation_id", "frame_id", "redshift", "retained"
    ]
    assert "inclusion_probability" not in first["validation_observed"].columns
    assert len(first["validation_observed"]) == 120
    assert len(first["validation_frame_truth"]) == 900
    assert first["config"]["independent_frame"] is True


def test_independent_validation_sampler_rejects_invalid_sample_sizes():
    with pytest.raises(ValueError):
        generate_independent_validation_sample(
            seed=1, n_validation_frame=100, n_validation=101
        )


def test_end_to_end_smoke_runner_writes_separate_and_checksummed_artifacts(
    tmp_path: Path,
):
    root = Path(__file__).resolve().parents[1]
    output_dir = tmp_path / "smoke"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src") + os.pathsep + env.get("PYTHONPATH", "")
    completed = subprocess.run(
        [
            sys.executable,
            str(root / "scripts" / "run_inclusion_posterior_smoke.py"),
            "--seed", "712",
            "--n-parent", "1200",
            "--n-events", "12",
            "--n-validation-frame", "1000",
            "--n-validation", "120",
            "--output-dir", str(output_dir),
        ],
        cwd=root,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    summary = json.loads(
        (output_dir / "analysis" / "summary.json").read_text()
    )
    assert summary["status"] == "SMOKE_COMPLETE_NOT_CALIBRATION"
    assert summary["analysis_seed"] != summary["validation_seed"]
    assert summary["seeds_are_distinct"] is True
    assert len(summary["quadrature_runs"]) == 3
    assert summary["quadrature_h0_shape_comparisons"]
    assert summary["h0_grid_comparison"]["fine_points"] == 81

    catalogue = pd.read_csv(
        output_dir / "observable" / "catalogue_observed.csv"
    )
    events = pd.read_csv(
        output_dir / "observable" / "gw_events_observed.csv"
    )
    validation = pd.read_csv(
        output_dir / "observable" / "validation_observed.csv"
    )
    assert len(catalogue) >= 2
    assert len(events) == 12
    assert len(validation) == 120
    assert "host_redshift" not in events.columns
    assert "inclusion_probability" not in validation.columns
    assert (events["distance_mpc"] >= 250.0).all()

    assert (output_dir / "truth" / "analysis_parent_truth.csv").exists()
    assert (output_dir / "truth" / "validation_frame_truth.csv").exists()
    assert not (output_dir / "observable" / "analysis_parent_truth.csv").exists()
    assert not (output_dir / "observable" / "validation_frame_truth.csv").exists()

    manifest = json.loads(
        (output_dir / "sha256_manifest.json").read_text()
    )
    assert "analysis/summary.json" in manifest
    for relative, expected in manifest.items():
        digest = hashlib.sha256(
            (output_dir / relative).read_bytes()
        ).hexdigest()
        assert digest == expected
