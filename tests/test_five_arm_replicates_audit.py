"""Regression tests for independent five-arm lifecycle audit integrity."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

MODULE_PATH = Path(__file__).parents[1] / "scripts" / "audit_five_arm_replicates.py"
SPEC = importlib.util.spec_from_file_location("audit_five_arm_replicates", MODULE_PATH)
AUDITOR = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(AUDITOR)


def _valid_artifact(root: Path) -> None:
    root.mkdir()
    master_seed = 731
    children = np.random.SeedSequence(master_seed).spawn(2)
    seeds = [int(child.generate_state(1, dtype=np.uint32)[0]) for child in children]
    sizes = [10]
    config = {
        "status": "engineering_pilot_not_confirmatory_not_observational",
        "master_seed": master_seed,
        "replicates": 1,
        "validation_frame": 20,
        "validation_sizes": sizes,
        "h0_true": 70.0,
        "catalogue_floor": 0.10,
        "n_parent": 800,
        "n_detected": 3,
        "h0_points": 31,
        "sigma_mpc": 140.0,
        "detection_threshold_mpc": 250.0,
        "quadrature_order": 3,
        "scott_multiplier": 0.7,
        "failure_bound_assumptions": {
            "h0_prior_support": [50.0, 90.0],
            "absolute_median_error_max": 20.0,
            "squared_median_error_max": 400.0,
            "coverage_failure_indicator_bounds": [0, 1],
            "interval_width_bounds": [0.0, 40.0],
        },
        "seed_schedule": [{
            "replicate": 1,
            "analysis_seed": seeds[0],
            "validation_seed": seeds[1],
        }],
    }
    (root / "config.json").write_text(json.dumps(config), encoding="utf-8")

    ledger_rows = []
    summary_rows = []
    for arm_index, arm in enumerate(AUDITOR.ARMS):
        ledger_rows.append({
            "replicate": 1,
            "analysis_seed": seeds[0],
            "validation_seed": seeds[1],
            "validation_size": 10,
            "arm": arm,
            "status": "SUCCESS",
            "audit_status": "PASS",
            "error": "",
        })
        median = 68.0 + arm_index
        summary_rows.append({
            "replicate": 1,
            "analysis_seed": seeds[0],
            "validation_seed": seeds[1],
            "validation_size": 10,
            "arm": arm,
            "status": "SUCCESS",
            "audit_status": "PASS",
            "h0_true": 70.0,
            "posterior_median": median,
            "median_error": median - 70.0,
            "q16": 60.0,
            "q84": 80.0,
        })
    ledger = pd.DataFrame(ledger_rows)
    summary = pd.DataFrame(summary_rows)
    audits = pd.DataFrame([{
        "replicate": 1,
        "status": "PASS",
        "manifest_files_verified": 4,
        "posterior_summary_rows_reconstructed": 5,
        "maximum_log_likelihood_absolute_difference": 0.0,
    }])
    ledger.to_csv(root / "replicate_ledger.csv", index=False)
    summary.to_csv(root / "replicate_posterior_summaries.csv", index=False)
    audits.to_csv(root / "replicate_audits.csv", index=False)
    aggregate = AUDITOR._expected_aggregate(summary, ledger, 1)
    contrasts = AUDITOR._expected_contrasts(summary)
    aggregate.to_csv(root / "aggregate_metrics.csv", index=False)
    contrasts.to_csv(root / "paired_contrasts.csv", index=False)
    replicate_dir = root / "replicates" / "replicate_0001"
    replicate_dir.mkdir(parents=True)
    (replicate_dir / "config.json").write_text(json.dumps({
        "seed": seeds[0],
        "validation_seed": seeds[1],
        "n_parent": 800,
        "n_detected": 3,
        "validation_frame": 20,
        "validation_sizes": sizes,
        "h0_true": 70.0,
        "catalogue_floor": 0.10,
        "h0_grid_points": 31,
        "sigma_mpc": 140.0,
        "detection_threshold_mpc": 250.0,
        "quadrature_order": 3,
        "scott_multiplier": 0.7,
    }), encoding="utf-8")
    (replicate_dir / "audit.json").write_text(json.dumps({
        "status": "PASS",
        "manifest_files_verified": 4,
        "posterior_summary_rows_reconstructed": 5,
        "maximum_log_likelihood_absolute_difference": 0.0,
    }), encoding="utf-8")

    manifest_files = sorted(
        path for path in root.rglob("*")
        if path.is_file() and path.name != "manifest.json"
    )
    import hashlib
    (root / "manifest.json").write_text(json.dumps({
        "status": "PASS",
        "file_count": len(manifest_files),
        "files": {
            str(path.relative_to(root)): {
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for path in manifest_files
        },
    }), encoding="utf-8")




def _refresh_manifest(root: Path) -> None:
    import hashlib
    paths = sorted(
        path for path in root.rglob("*")
        if path.is_file() and path.name not in {"manifest.json", "independent_audit.json"}
    )
    (root / "manifest.json").write_text(json.dumps({
        "status": "PASS",
        "file_count": len(paths),
        "files": {
            str(path.relative_to(root)): {
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
            for path in paths
        },
    }), encoding="utf-8")


def test_valid_lifecycle_and_aggregate_artifacts_pass(tmp_path: Path):
    root = tmp_path / "valid"
    _valid_artifact(root)
    report = AUDITOR.audit(root)
    assert report["status"] == "PASS"
    assert report["planned_lifecycle_cells"] == 5
    assert report["paired_contrast_rows_reconstructed"] == 8


def test_auditor_rejects_missing_lifecycle_cell(tmp_path: Path):
    root = tmp_path / "missing"
    _valid_artifact(root)
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    # Preserve row count and uniqueness while replacing one planned key with an invalid key.
    ledger.loc[ledger.index[-1], "validation_size"] = 11
    ledger.to_csv(root / "replicate_ledger.csv", index=False)
    _refresh_manifest(root)
    with pytest.raises(AssertionError, match="lifecycle key set mismatch"):
        AUDITOR.audit(root)


def test_auditor_rejects_summary_median_error_mismatch(tmp_path: Path):
    root = tmp_path / "bad-summary"
    _valid_artifact(root)
    summary = pd.read_csv(root / "replicate_posterior_summaries.csv")
    summary.loc[0, "median_error"] += 0.5
    summary.to_csv(root / "replicate_posterior_summaries.csv", index=False)
    _refresh_manifest(root)
    with pytest.raises(AssertionError, match="median error"):
        AUDITOR.audit(root)

def test_auditor_rejects_tampered_manifest_file(tmp_path: Path):
    root = tmp_path / "tampered"
    _valid_artifact(root)
    with (root / "aggregate_metrics.csv").open("a", encoding="utf-8") as stream:
        stream.write("tampered\\n")
    with pytest.raises(AssertionError, match="manifest size mismatch|byte-count mismatch|SHA-256 mismatch"):
        AUDITOR.audit(root)
