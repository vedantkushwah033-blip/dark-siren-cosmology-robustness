from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

import pandas as pd
import pytest

from dark_siren_robustness.estimated_inclusion_benchmark import (
    run_probability_validation_benchmark,
)
from scripts.audit_finite_validation_probability_benchmark import audit


def test_independent_probability_auditor_reconstructs_raw_records(tmp_path: Path):
    raw_path = tmp_path / "raw_unit_records.csv"
    summary_path = tmp_path / "replicates.csv"
    summary = run_probability_validation_benchmark(
        n_parent=500,
        validation_sizes=(100,),
        curvature_values=(0.0, 1.2),
        degrees=(1, 2),
        replicates=2,
        seed=91,
        raw_output=raw_path,
    )
    summary.to_csv(summary_path, index=False)

    report = audit(summary_path, raw_path)
    assert report["status"] == "PASS"
    assert report["experiments_reconstructed"] == len(summary)
    raw = pd.read_csv(raw_path)
    assert set(raw["split_role"]) == {"validation", "analysis_catalogue", "nonretained_holdout"}
    assert not raw.duplicated([
        "scenario_seed", "replicate", "gamma", "retention_target",
        "n_validation", "model_degree", "object_id",
    ]).any()
    assert (raw.groupby([
        "scenario_seed", "replicate", "gamma", "retention_target",
        "n_validation", "model_degree",
    ]).size() == 500).all()


def test_independent_probability_auditor_rejects_tampered_summary(tmp_path: Path):
    raw_path = tmp_path / "raw_unit_records.csv"
    summary_path = tmp_path / "replicates.csv"
    summary = run_probability_validation_benchmark(
        n_parent=500,
        validation_sizes=(100,),
        curvature_values=(0.0,),
        retention_fractions=(0.5,),
        degrees=(1,),
        replicates=1,
        seed=92,
        raw_output=raw_path,
    )
    summary.loc[0, "heldout_brier_score"] += 0.01
    summary.to_csv(summary_path, index=False)
    with pytest.raises(ValueError, match="mismatch"):
        audit(summary_path, raw_path)



def test_independent_probability_auditor_rejects_incomplete_parent_partition(tmp_path: Path):
    raw_path = tmp_path / "raw_unit_records.csv"
    summary_path = tmp_path / "replicates.csv"
    summary = run_probability_validation_benchmark(
        n_parent=500, validation_sizes=(100,), curvature_values=(0.0,),
        degrees=(1,), replicates=1, seed=93, raw_output=raw_path,
    )
    summary.to_csv(summary_path, index=False)
    raw = pd.read_csv(raw_path).iloc[1:]
    raw.to_csv(raw_path, index=False)
    with pytest.raises(ValueError, match="expected 500 parent records"):
        audit(summary_path, raw_path)


def test_independent_probability_auditor_rejects_prediction_not_from_coefficients(tmp_path: Path):
    raw_path = tmp_path / "raw_unit_records.csv"
    summary_path = tmp_path / "replicates.csv"
    summary = run_probability_validation_benchmark(
        n_parent=500, validation_sizes=(100,), curvature_values=(0.0,),
        retention_fractions=(0.5,), degrees=(1,), replicates=1, seed=94,
        raw_output=raw_path,
    )
    summary.to_csv(summary_path, index=False)
    raw = pd.read_csv(raw_path)
    raw.loc[0, "p_hat"] = min(0.999, raw.loc[0, "p_hat"] + 0.01)
    raw.to_csv(raw_path, index=False)
    with pytest.raises(ValueError, match="predicted probabilities do not match"):
        audit(summary_path, raw_path)



def test_curvature_scenarios_share_parent_frame_and_validation_ids(tmp_path: Path):
    raw_path = tmp_path / "raw_unit_records.csv"
    run_probability_validation_benchmark(
        n_parent=500, validation_sizes=(100,), curvature_values=(-1.2, 0.0, 1.2),
        degrees=(1,), retention_fractions=(0.1,), replicates=1, seed=95, raw_output=raw_path,
    )
    raw = pd.read_csv(raw_path)
    frames = []
    for gamma in (-1.2, 0.0, 1.2):
        frame = raw.loc[raw["gamma"] == gamma].sort_values("object_id")
        frames.append(frame)
    for frame in frames[1:]:
        assert frame["object_id"].tolist() == frames[0]["object_id"].tolist()
        assert frame["redshift"].tolist() == frames[0]["redshift"].tolist()
        validation_ids = frame.loc[
            frame["split_role"] == "validation", "object_id"
        ].tolist()
        reference_validation_ids = frames[0].loc[
            frames[0]["split_role"] == "validation", "object_id"
        ].tolist()
        assert validation_ids == reference_validation_ids



def test_configured_runner_writes_provenance_manifest_and_auditable_artifacts(tmp_path: Path):
    root = Path(__file__).resolve().parents[1]
    config_path = tmp_path / "pilot.json"
    config_path.write_text(json.dumps({
        "schema_version": 1,
        "n_parent": 500,
        "validation_sizes": [100],
        "curvature_values": [0.0],
        "retention_fractions": [0.5],
        "model_degrees": [1, 2],
        "replicates": 1,
        "seed": 96,
    }), encoding="utf-8")
    summary_path = tmp_path / "output" / "replicates.csv"
    raw_path = tmp_path / "output" / "raw_unit_records.csv"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src") + os.pathsep + env.get("PYTHONPATH", "")
    completed = subprocess.run(
        [
            sys.executable,
            str(root / "scripts" / "run_finite_validation_probability_benchmark.py"),
            "--config", str(config_path),
            "--output", str(summary_path),
            "--raw-output", str(raw_path),
        ],
        cwd=root, env=env, text=True, capture_output=True, check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    manifest_path = summary_path.parent / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["summary_rows"] == 2
    assert manifest["raw_unit_records"] == 1000
    assert manifest["summary_sha256"] == hashlib.sha256(summary_path.read_bytes()).hexdigest()
    assert manifest["raw_sha256"] == hashlib.sha256(raw_path.read_bytes()).hexdigest()
    assert audit(summary_path, raw_path)["status"] == "PASS"
