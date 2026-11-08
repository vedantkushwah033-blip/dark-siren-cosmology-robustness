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
        "validation_id", "redshift", "retained"
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
    audit = json.loads(
        (output_dir / "analysis" / "independent_audit.json").read_text()
    )
    assert audit["status"] == "PASS"
    assert audit["observable_only_reconstruction"] is True
    assert audit["truth_used_only_after_reconstruction"] is True
    assert audit["generation_trace_attempts_crosschecked"] > 0
    assert audit["selected_events_crosschecked"] == 12
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


def test_fixed_seed_validation_samples_are_nested_across_sizes():
    small = generate_independent_validation_sample(
        seed=5511, n_validation_frame=1200, n_validation=100
    )["validation_observed"]
    large = generate_independent_validation_sample(
        seed=5511, n_validation_frame=1200, n_validation=300
    )["validation_observed"]
    assert set(small["validation_id"]).issubset(set(large["validation_id"]))
    indexed = large.set_index("validation_id")
    for row in small.itertuples(index=False):
        assert indexed.loc[row.validation_id, "redshift"] == row.redshift
        assert indexed.loc[row.validation_id, "retained"] == row.retained


def test_independent_validation_h0_runner_records_prior_and_passes_independent_audit(
    tmp_path: Path,
):
    root = Path(__file__).resolve().parents[1]
    output_dir = tmp_path / "independent_h0"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src") + os.pathsep + env.get("PYTHONPATH", "")
    completed = subprocess.run(
        [
            sys.executable,
            str(root / "scripts" / "run_independent_validation_inclusion_h0_smoke.py"),
            "--seed", "812",
            "--validation-seed", "913",
            "--n-parent", "800",
            "--n-detected", "3",
            "--n-validation-frame", "500",
            "--n-validation", "80",
            "--orders", "3,5",
            "--h0-points", "21",
            "--prior-slope-multipliers", "0.5,1.0,2.0",
            "--output-dir", str(output_dir),
        ],
        cwd=root,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr

    config = json.loads((output_dir / "config.json").read_text())
    prior = config["inference"]["coefficient_prior"]
    assert prior["family"] == "independent_normal"
    assert prior["mean"] == [0.0, 0.0]
    assert prior["sd"] == [5.0, 2.5]
    coordinates = config["inference"]["coefficient_coordinates"]
    assert coordinates["basis"] == "fixed_domain_standardization"
    assert coordinates["center"] == 0.205
    assert coordinates["scale"] == 0.195

    audited = subprocess.run(
        [
            sys.executable,
            str(root / "scripts" / "audit_independent_validation_inclusion_h0_smoke.py"),
            str(output_dir),
        ],
        cwd=root,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )
    assert audited.returncode == 0, audited.stdout + audited.stderr
    report = json.loads(audited.stdout)
    assert report["status"] == "PASS"
    assert report["manifest_files_verified"] > 0
    assert report["h0_posterior_normalization_checks"] == 2
    assert report["quadrature_orders"] == [3, 5]
