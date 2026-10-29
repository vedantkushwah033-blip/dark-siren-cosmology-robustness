from __future__ import annotations

import csv
import hashlib
import json
import shutil
from pathlib import Path

import pytest

from dark_siren_robustness.selection_conditioned_h0_dependent_noise import (
    # Importing the runner function is intentional: the test fixture creates an
    # authentic artifact schema, while the auditor must independently reconstruct it.
)
from scripts.run_h0_dependent_noise_sbc import run_experiment
from scripts.audit_h0_dependent_noise_sbc import audit


@pytest.fixture(scope="module")
def valid_artifact(tmp_path_factory):
    out = tmp_path_factory.mktemp("h0-sbc-valid")
    rows, summary = run_experiment(replicates=30, seed=991827)
    csv_path = out / "replicates.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary_path = out / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    files = {}
    for path in (csv_path, summary_path):
        files[path.name] = {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size,
        }
    manifest = {
        "schema_version": 1,
        "git_commit": summary["git_commit"],
        "replicate_count": len(rows),
        "files": files,
        "configuration": summary["configuration"],
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return out


def _copy_artifact(source: Path, destination: Path):
    shutil.copytree(source, destination)
    return destination


def _rehash(artifact: Path, *names: str):
    manifest_path = artifact / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for name in names:
        path = artifact / name
        manifest["files"][name] = {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size,
        }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def test_independent_auditor_accepts_replayed_and_reaggregated_artifact(valid_artifact, tmp_path):
    artifact = _copy_artifact(valid_artifact, tmp_path / "valid-copy")
    report = audit(artifact)
    assert report["status"] == "PASS"
    assert report["raw_records"] == 120
    assert report["summary_aggregates_reconstructed"] is True


def test_independent_auditor_rejects_rehashed_but_tampered_summary(valid_artifact, tmp_path):
    artifact = _copy_artifact(valid_artifact, tmp_path / "tampered-summary")
    summary_path = artifact / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["summaries"][0]["bias"] += 10.0
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _rehash(artifact, "summary.json")
    with pytest.raises(ValueError, match="summary aggregate bias mismatch"):
        audit(artifact)


def test_independent_auditor_rejects_rehashed_protocol_configuration_change(valid_artifact, tmp_path):
    artifact = _copy_artifact(valid_artifact, tmp_path / "tampered-config")
    summary_path = artifact / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["configuration"]["detection_distance_min_mpc"] = 301.0
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest_path = artifact / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["configuration"] = summary["configuration"]
    manifest["files"]["summary.json"] = {
        "sha256": hashlib.sha256(summary_path.read_bytes()).hexdigest(),
        "bytes": summary_path.stat().st_size,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="detection threshold"):
        audit(artifact)


def test_independent_auditor_rejects_rehashed_seed_mismatch(valid_artifact, tmp_path):
    artifact = _copy_artifact(valid_artifact, tmp_path / "tampered-seed")
    csv_path = artifact / "replicates.csv"
    rows = list(csv.DictReader(csv_path.open(newline="", encoding="utf-8")))
    rows[0]["seed"] = str(int(rows[0]["seed"]) + 1)
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    _rehash(artifact, "replicates.csv")
    with pytest.raises(ValueError, match="seed does not match"):
        audit(artifact)
