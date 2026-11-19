from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from run_full_parent_h0_fixed_count_baseline import (
    draw_detected_distances,
    prepare_full_parent_likelihood,
    selected_posterior_quantiles,
    run,
)


def test_detected_distance_sampler_is_reproducible_and_respects_threshold():
    redshift = np.linspace(0.03, 0.38, 80)
    host_probability = np.linspace(1.0, 2.0, len(redshift))
    host_probability /= host_probability.sum()
    a = draw_detected_distances(
        np.random.default_rng(123), redshift, host_probability,
        70.0, 140.0, 250.0, 8,
    )
    b = draw_detected_distances(
        np.random.default_rng(123), redshift, host_probability,
        70.0, 140.0, 250.0, 8,
    )
    assert np.array_equal(a, b)
    assert a.shape == (8,)
    assert np.all(a >= 250.0)


def test_full_parent_selected_posterior_quantiles_are_finite_and_ordered():
    redshift = np.linspace(0.03, 0.38, 80)
    host_weights = 1.0 + 1.5 * (redshift - 0.01) / (0.40 - 0.01)
    grid = np.linspace(20.0, 300.0, 281)
    prepared = prepare_full_parent_likelihood(
        redshift, host_weights, grid, 140.0, 250.0
    )
    distances = draw_detected_distances(
        np.random.default_rng(456), redshift,
        host_weights / host_weights.sum(), 70.0, 140.0, 250.0, 4,
    )
    quantiles = selected_posterior_quantiles(distances, prepared, 140.0, 250.0)
    values = np.array([
        quantiles[k] for k in ("q025", "q05", "q16", "median", "q84", "q95", "q975")
    ])
    assert np.isfinite(values).all()
    assert np.all(np.diff(values) > 0)
    assert 20.0 <= values[0] < values[-1] <= 300.0


def _write_input_artifact(root: Path, status: str = "SUCCESS") -> None:
    cell = root / "replicates" / "gm1_n100_r000"
    (cell / "truth").mkdir(parents=True)
    config = {
        "analysis_generator_config": {
            "redshift_domain": [0.01, 0.40],
            "sigma_mpc": 140.0,
            "detection_threshold_mpc": 250.0,
            "n_detected_events": 4,
        },
        "inference": {"uniform_h0_prior": [20.0, 300.0]},
    }
    (cell / "config.json").write_text(json.dumps(config))
    parent = pd.DataFrame({
        "object_id": np.arange(100, dtype=np.int64),
        "redshift": np.linspace(0.01, 0.40, 100),
    })
    parent.to_csv(cell / "truth" / "analysis_parent_truth.csv", index=False)
    ledger = pd.DataFrame([{
        "replicate": 0, "gamma": -1.0, "n_validation": 100,
        "cell_id": "gm1_n100_r000", "status": status, "h0_true": 70.0,
    }])
    ledger.to_csv(root / "replicate_ledger.csv", index=False)
    manifest = {}
    for path in (root / "replicate_ledger.csv", cell / "config.json",
                 cell / "truth" / "analysis_parent_truth.csv"):
        manifest[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    (root / "manifest.json").write_text(json.dumps(manifest, sort_keys=True))


def test_small_full_parent_baseline_run_writes_verifiable_outputs(tmp_path):
    input_root = tmp_path / "input"
    input_root.mkdir()
    _write_input_artifact(input_root)
    output_root = tmp_path / "output"
    frame, summary = run(
        input_root, output_root, replicates_per_parent=3, seed=990000, h0_points=41
    )
    assert len(frame) == 3
    assert summary["status"] == "DIAGNOSTIC_COMPLETE"
    assert summary["parent_populations"] == 1
    assert summary["total_replicates"] == 3
    assert len(pd.read_csv(output_root / "parent_summary.csv")) == 1
    saved = pd.read_csv(output_root / "replicates.csv")
    assert len(saved) == 3
    manifest = json.loads((output_root / "manifest.json").read_text())
    for name, digest in manifest.items():
        assert hashlib.sha256((output_root / name).read_bytes()).hexdigest() == digest


def test_full_parent_preparation_rejects_nonpositive_weights():
    with pytest.raises(ValueError, match="invalid parent"):
        prepare_full_parent_likelihood(
            np.array([0.1, 0.2]), np.array([1.0, 0.0]),
            np.linspace(20.0, 300.0, 41), 140.0, 250.0,
        )



def test_small_full_parent_baseline_accepts_non_evaluable_fitted_cell(tmp_path):
    input_root = tmp_path / "input"
    input_root.mkdir()
    _write_input_artifact(input_root, status="NON_EVALUABLE")
    output_root = tmp_path / "output"
    frame, summary = run(
        input_root, output_root, replicates_per_parent=2, seed=990000, h0_points=41
    )
    assert len(frame) == 2
    assert summary["status"] == "DIAGNOSTIC_COMPLETE"
    assert summary["parent_populations"] == 1


def test_full_parent_run_can_override_detected_event_count(tmp_path):
    input_root = tmp_path / "input"
    input_root.mkdir()
    _write_input_artifact(input_root)
    output_root = tmp_path / "output"
    frame, summary = run(
        input_root, output_root, replicates_per_parent=2, seed=990000,
        h0_points=41, n_events_override=7,
    )
    assert len(frame) == 2
    assert set(frame.n_events) == {7}
    assert summary["n_events"] == 7
