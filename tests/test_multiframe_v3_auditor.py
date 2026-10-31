from __future__ import annotations

import csv
import hashlib
import json

import pytest

from scripts.audit_multiframe_v3_pilot import audit_artifact
from scripts.run_multiframe_v3_pilot import run_multiframe_v3
from dark_siren_robustness.multiframe_design import MANIFEST_SCHEMA


def manifest():
    return {
        "schema": MANIFEST_SCHEMA, "root_seed": 20261031,
        "population_size": 80, "frames_per_scenario": 2,
        "outer_repetitions": 3, "bootstrap_replicates": 25,
        "nominal_level": 0.95, "sample_sizes": [10, 30],
        "methods": [{"method_id": "wald_srswor"},
                    {"method_id": "ordinary_percentile_bootstrap"},
                    {"method_id": "studentized_pseudopop_srswor"}],
        "scenarios": [{"scenario_id": 3, "name": "Normal", "family": "normal",
                       "params": {"loc": 0.0, "scale": 1.0}}],
    }


def refresh_checksum(directory, filename, key):
    digest = hashlib.sha256((directory / filename).read_bytes()).hexdigest()
    path = directory / "manifest.json"
    payload = json.loads(path.read_text())
    payload[key] = digest
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def make_pilot(tmp_path):
    out = tmp_path / "pilot"
    run_multiframe_v3(manifest(), out)
    return out


def test_v3_independent_auditor_passes_runner_artifact(tmp_path):
    out = make_pilot(tmp_path)
    report = audit_artifact(out)
    assert report["status"] == "PASS"
    assert report["frames_verified"] == 2
    assert report["replicate_rows_independently_reconstructed"] == 36


def test_auditor_detects_frame_tampering_even_after_checksum_refresh(tmp_path):
    out = make_pilot(tmp_path)
    path = out / "frames.csv"
    rows = list(csv.DictReader(path.open(newline="")))
    rows[0]["value"] = str(float(rows[0]["value"]) + 0.25)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    refresh_checksum(out, "frames.csv", "frame_csv_sha256")
    with pytest.raises(ValueError, match="finite-frame values"):
        audit_artifact(out)


def test_auditor_detects_sample_index_tampering_even_after_checksum_refresh(tmp_path):
    out = make_pilot(tmp_path)
    path = out / "replicates.csv"
    rows = list(csv.DictReader(path.open(newline="")))
    indices = json.loads(rows[0]["sample_indices"])
    indices[0] = (indices[0] + 1) % 80
    rows[0]["sample_indices"] = json.dumps(indices, separators=(",", ":"))
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    refresh_checksum(out, "replicates.csv", "replicate_csv_sha256")
    with pytest.raises(ValueError, match="sample indices"):
        audit_artifact(out)


def test_auditor_detects_endpoint_tampering_even_after_checksum_refresh(tmp_path):
    out = make_pilot(tmp_path)
    path = out / "replicates.csv"
    rows = list(csv.DictReader(path.open(newline="")))
    rows[0]["interval_lower"] = str(float(rows[0]["interval_lower"]) + 0.1)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    refresh_checksum(out, "replicates.csv", "replicate_csv_sha256")
    with pytest.raises(ValueError, match="endpoint"):
        audit_artifact(out)


def test_auditor_detects_missing_rows_even_after_checksum_refresh(tmp_path):
    out = make_pilot(tmp_path)
    path = out / "replicates.csv"
    rows = list(csv.DictReader(path.open(newline="")))
    rows.pop()
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    refresh_checksum(out, "replicates.csv", "replicate_csv_sha256")
    with pytest.raises(ValueError, match="incomplete"):
        audit_artifact(out)
