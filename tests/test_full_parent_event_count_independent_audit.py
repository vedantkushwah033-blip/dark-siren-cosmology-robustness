from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from run_full_parent_event_count_calibration import run

MODULE_PATH = ROOT / "scripts" / "audit_full_parent_event_count_calibration.py"
SPEC = importlib.util.spec_from_file_location("event_count_independent_audit", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def build_artifact(root: Path, replicates=2):
    raw, summary, paired = run(
        replicates=replicates, master_seed=2026111931, n_parent=300,
        h0_true=70.0, prior_min=20.0, prior_max=300.0, h0_points=41,
    )
    root.mkdir(parents=True, exist_ok=True)
    raw.to_csv(root / "replicates.csv", index=False, float_format="%.17g")
    summary.to_csv(root / "summary.csv", index=False, float_format="%.17g")
    paired.to_csv(root / "paired_contrast.csv", index=False, float_format="%.17g")
    design = {
        "replicates": replicates, "master_seed": 2026111931, "n_parent": 300,
        "event_counts": [4, 30], "event_design": "nested prefix of the same 30-event realization",
        "h0_true": 70.0, "uniform_h0_prior": [20.0, 300.0],
        "h0_grid_points": 41, "status": "COMPLETE",
    }
    (root / "design.json").write_text(json.dumps(design, sort_keys=True))
    manifest = {
        name: hashlib.sha256((root / name).read_bytes()).hexdigest()
        for name in ("replicates.csv", "summary.csv", "paired_contrast.csv", "design.json")
    }
    (root / "manifest.json").write_text(json.dumps(manifest, sort_keys=True))
    return raw


def refresh_manifest(root: Path):
    manifest = {
        name: hashlib.sha256((root / name).read_bytes()).hexdigest()
        for name in ("replicates.csv", "summary.csv", "paired_contrast.csv", "design.json")
    }
    (root / "manifest.json").write_text(json.dumps(manifest, sort_keys=True))


def test_independent_audit_reconstructs_both_event_count_arms(tmp_path):
    root = tmp_path / "artifact"
    build_artifact(root, replicates=2)
    report = MODULE.audit(root, expected_replicates=2, output=tmp_path / "audit.json")
    assert report["status"] == "PASS"
    assert report["raw_rows_reconstructed"] == 4
    assert report["parent_populations_regenerated"] == 2
    assert report["posterior_curves_reconstructed"] == 4


def test_audit_rejects_tampered_posterior_median_even_with_refreshed_hash(tmp_path):
    root = tmp_path / "artifact"
    build_artifact(root, replicates=2)
    raw = pd.read_csv(root / "replicates.csv")
    raw.loc[0, "h0_median"] += 0.25
    raw.to_csv(root / "replicates.csv", index=False, float_format="%.17g")
    refresh_manifest(root)
    with pytest.raises(AssertionError, match="h0_median"):
        MODULE.audit(root, expected_replicates=2)


def test_audit_rejects_missing_event_count_arm(tmp_path):
    root = tmp_path / "artifact"
    build_artifact(root, replicates=2)
    raw = pd.read_csv(root / "replicates.csv")
    raw = raw[raw.event_count == 4]
    raw.to_csv(root / "replicates.csv", index=False, float_format="%.17g")
    refresh_manifest(root)
    with pytest.raises(AssertionError, match="incomplete"):
        MODULE.audit(root, expected_replicates=2)


def test_audit_rejects_manifest_tampering(tmp_path):
    root = tmp_path / "artifact"
    build_artifact(root, replicates=2)
    with (root / "summary.csv").open("a") as handle:
        handle.write("\n")
    with pytest.raises(AssertionError, match="manifest mismatch"):
        MODULE.audit(root, expected_replicates=2)
