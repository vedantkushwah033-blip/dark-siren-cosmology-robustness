from __future__ import annotations

import csv
import hashlib
import json
import shutil
from pathlib import Path

import pytest

from scripts.diagnose_oracle_bias import run_campaign, write_artifacts
from scripts.audit_oracle_bias import audit


@pytest.fixture(scope="module")
def valid_artifact(tmp_path_factory):
    folder = tmp_path_factory.mktemp("oracle-bias-artifact")
    summary, records = run_campaign(
        modes=["70", "60", "80", "prior"],
        reps=30,
        seed=991827,
        z_points=101,
        h0_points=31,
    )
    paths = write_artifacts(summary, records, folder / "oracle_bias_exact_density.json")
    return paths["summary"]


def _copy_artifact(summary_path: Path, destination: Path):
    destination.mkdir(parents=True, exist_ok=True)
    for path in (
        summary_path,
        summary_path.with_name(summary_path.stem + "_replicates.csv"),
        summary_path.with_name(summary_path.stem + "_manifest.json"),
    ):
        shutil.copy2(path, destination / path.name)
    return destination / summary_path.name


def _rehash(summary_path: Path, *names: str):
    manifest_path = summary_path.with_name(summary_path.stem + "_manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for name in names:
        path = summary_path.with_name(name)
        manifest["files"][name] = {
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "bytes": path.stat().st_size,
        }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def test_independent_oracle_auditor_replays_and_reaggregates_all_rows(valid_artifact):
    result = audit(valid_artifact)
    assert result["status"] == "PASS"
    assert result["raw_records"] == 120
    assert result["rows_replayed"] == 120
    assert result["posteriors_independently_reconstructed"] == 120
    assert result["summary_aggregates_reconstructed"] is True


def test_oracle_auditor_rejects_rehashed_but_tampered_summary(valid_artifact, tmp_path):
    summary_path = _copy_artifact(valid_artifact, tmp_path / "tampered-summary")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    summary["results"][0]["mean_bias"] += 1.0
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    _rehash(summary_path, summary_path.name)
    with pytest.raises(ValueError, match="summary aggregate mean_bias mismatch"):
        audit(summary_path)


def test_oracle_auditor_rejects_rehashed_but_tampered_data(valid_artifact, tmp_path):
    summary_path = _copy_artifact(valid_artifact, tmp_path / "tampered-data")
    raw_path = summary_path.with_name(summary_path.stem + "_replicates.csv")
    with raw_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    distances = json.loads(rows[0]["distances_json"])
    distances[0] += 5.0
    rows[0]["distances_json"] = json.dumps(distances, separators=(",", ":"))
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    _rehash(summary_path, raw_path.name)
    with pytest.raises(ValueError, match="data replay mismatch|posterior"):
        audit(summary_path)


def test_oracle_auditor_rejects_rehashed_nonfinite_metrics(valid_artifact, tmp_path):
    summary_path = _copy_artifact(valid_artifact, tmp_path / "nonfinite")
    raw_path = summary_path.with_name(summary_path.stem + "_replicates.csv")
    with raw_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    rows[0]["posterior_mean"] = "nan"
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    _rehash(summary_path, raw_path.name)
    with pytest.raises(ValueError, match="non-finite saved posterior metric"):
        audit(summary_path)


def test_oracle_auditor_rejects_rehashed_wrong_mode_seed(valid_artifact, tmp_path):
    summary_path = _copy_artifact(valid_artifact, tmp_path / "seed-mismatch")
    raw_path = summary_path.with_name(summary_path.stem + "_replicates.csv")
    with raw_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    rows[0]["mode_seed"] = str(int(rows[0]["mode_seed"]) + 1)
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    _rehash(summary_path, raw_path.name)
    with pytest.raises(ValueError, match="mode seed mismatch"):
        audit(summary_path)
