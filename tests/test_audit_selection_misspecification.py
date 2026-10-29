from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_selection_misspecification import audit_artifacts


def _run_smoke(output_dir: Path):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "run_selection_misspecification.py"),
         "--replicates", "2", "--n-parent", "1200", "--n-detected", "8",
         "--seed", "55129", "--h0-points", "51", "--z-points", "201",
         "--output-dir", str(output_dir)],
        cwd=ROOT, env=env, text=True, capture_output=True, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_runner_writes_manifest_and_independent_audit_passes(tmp_path):
    folder = tmp_path / "pilot"
    _run_smoke(folder)
    report = audit_artifacts(folder)
    assert report["status"] == "PASS"
    assert report["raw_records"] == 12
    assert report["replicates"] == 2
    assert report["posterior_grids_reconstructed"] == 12
    manifest = json.loads((folder / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["files"]["replicate_metrics.csv"]["sha256"] == report["replicate_metrics_sha256"]
    assert len(report["summary_sha256"]) == 64


def test_audit_rejects_tampered_posterior_even_if_manifest_is_rehashed(tmp_path):
    folder = tmp_path / "pilot"
    _run_smoke(folder)
    raw = folder / "replicate_metrics.csv"
    with raw.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        columns = list(rows[0])
    density = json.loads(rows[0]["posterior_density_grid_json"])
    density[len(density) // 2] *= 1.25
    rows[0]["posterior_density_grid_json"] = json.dumps(density, separators=(",", ":"))
    with raw.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    manifest_path = folder / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["files"]["replicate_metrics.csv"]["sha256"] = hashlib.sha256(raw.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + chr(10), encoding="utf-8")
    with pytest.raises(ValueError, match="reconstructed"):
        audit_artifacts(folder)


def test_audit_rejects_missing_arm_even_when_manifest_is_rehashed(tmp_path):
    folder = tmp_path / "pilot"
    _run_smoke(folder)
    raw = folder / "replicate_metrics.csv"
    with raw.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        columns = list(rows[0])
    rows = [row for row in rows if not (row["replicate"] == "1" and row["arm"] == "oracle")]
    with raw.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    # Rehash the modified CSV so the audit reaches structural validation
    # instead of passing merely because the manifest detects tampering.
    manifest_path = folder / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["files"]["replicate_metrics.csv"]["sha256"] = hashlib.sha256(raw.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + chr(10), encoding="utf-8")
    with pytest.raises(ValueError, match="raw record count mismatch"):
        audit_artifacts(folder)



def test_audit_rejects_tampered_primary_estimand_even_if_summary_is_rehashed(tmp_path):
    folder = tmp_path / "pilot"
    _run_smoke(folder)
    summary_path = folder / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["primary_estimand"]["mean"] += 1.0
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + chr(10), encoding="utf-8")
    manifest_path = folder / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["files"]["summary.json"]["sha256"] = hashlib.sha256(summary_path.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + chr(10), encoding="utf-8")
    with pytest.raises(ValueError, match="primary estimand mean mismatch"):
        audit_artifacts(folder)


def test_audit_rejects_impossible_weight_ess_even_if_manifest_is_rehashed(tmp_path):
    folder = tmp_path / "pilot"
    _run_smoke(folder)
    raw = folder / "replicate_metrics.csv"
    with raw.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        columns = list(rows[0])
    # Catalogue-only has at most n_catalogue_observed effective samples.
    row = next(row for row in rows if row["arm"] == "catalogue_only")
    row["weight_ess"] = str(int(row["n_catalogue_observed"]) + 10)
    with raw.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    manifest_path = folder / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["files"]["replicate_metrics.csv"]["sha256"] = hashlib.sha256(raw.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + chr(10), encoding="utf-8")
    with pytest.raises(ValueError, match="weight ESS"):
        audit_artifacts(folder)


def test_audit_rejects_invalid_assumed_inclusion_probability_even_if_manifest_is_rehashed(tmp_path):
    folder = tmp_path / "pilot"
    _run_smoke(folder)
    raw = folder / "replicate_metrics.csv"
    with raw.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        columns = list(rows[0])
    rows[0]["min_assumed_inclusion_probability"] = "0"
    with raw.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    manifest_path = folder / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["files"]["replicate_metrics.csv"]["sha256"] = hashlib.sha256(raw.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + chr(10), encoding="utf-8")
    with pytest.raises(ValueError, match="minimum assumed inclusion probability"):
        audit_artifacts(folder)



def test_audit_rejects_invalid_detection_efficiency_even_if_manifest_is_rehashed(tmp_path):
    folder = tmp_path / "pilot"
    _run_smoke(folder)
    raw = folder / "replicate_metrics.csv"
    with raw.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        columns = list(rows[0])
    rows[0]["min_detection_efficiency"] = "0"
    with raw.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    manifest_path = folder / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["files"]["replicate_metrics.csv"]["sha256"] = hashlib.sha256(raw.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + chr(10), encoding="utf-8")
    with pytest.raises(ValueError, match="minimum detection efficiency"):
        audit_artifacts(folder)
