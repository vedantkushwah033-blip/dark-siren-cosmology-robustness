"""Regression tests for the independent multi-seed aggregate audit."""
import csv
import json

import pytest

from scripts.audit_multiseed_grid_sensitivity import audit_multiseed_sensitivity
from scripts.run_multiseed_grid_sensitivity import run_multiseed_sensitivity


def test_multiseed_audit_reconstructs_aggregate_and_rejects_tampering(tmp_path):
    run_multiseed_sensitivity(
        grid_points=[1201, 2401],
        seeds=[301, 302],
        n_replicates=100,
        output_dir=tmp_path,
    )
    report = json.loads((tmp_path / "multiseed_audit_report.json").read_text(encoding="utf-8"))
    assert report["status"] == "PASS"
    assert report["per_seed_pairs_complete"] is True
    assert report["seed_level_audits_reconstructed"] is True
    assert report["aggregate_statistics_reconstructed"] is True
    assert audit_multiseed_sensitivity(tmp_path) == report

    path = tmp_path / "multiseed_grid_sensitivity.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        fields = list(rows[0])
    rows[0]["mean_pit_mean_across_seeds"] = str(
        float(rows[0]["mean_pit_mean_across_seeds"]) + 0.01
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="mean_pit_mean_across_seeds"):
        audit_multiseed_sensitivity(tmp_path)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("n_seeds", "3", "incorrect aggregate seed count"),
        ("total_replicates_across_seeds", "999", "incorrect total replicate count"),
        ("seed_list", "301;999", "seed list does not match"),
    ],
)
def test_multiseed_audit_rejects_tampered_metadata(tmp_path, field, value, message):
    run_multiseed_sensitivity(
        grid_points=[1201, 2401],
        seeds=[301, 302],
        n_replicates=100,
        output_dir=tmp_path,
    )
    path = tmp_path / "multiseed_grid_sensitivity.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        fields = list(rows[0])
    rows[0][field] = value
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match=message):
        audit_multiseed_sensitivity(tmp_path)


def test_multiseed_audit_rejects_tampered_per_seed_input(tmp_path):
    run_multiseed_sensitivity(
        grid_points=[1201, 2401],
        seeds=[302, 301],
        n_replicates=100,
        output_dir=tmp_path,
    )
    path = tmp_path / "per_seed_paired_contrasts.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        fields = list(rows[0])
    rows[0]["pit_mean"] = str(float(rows[0]["pit_mean"]) + 0.01)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="aggregate input rows differ from audited comparison table"):
        audit_multiseed_sensitivity(tmp_path)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("seeds", [301, 999], "metadata seeds do not match"),
        ("grid_points", [1201, 9999], "metadata grid_points do not match"),
        ("replicates_per_seed_per_grid", 999, "metadata replicate count does not match"),
        ("independent_aggregate_audit", "FAIL", "does not confirm passing aggregate audit"),
    ],
)
def test_multiseed_audit_rejects_tampered_run_metadata(tmp_path, field, value, message):
    run_multiseed_sensitivity(
        grid_points=[1201, 2401],
        seeds=[301, 302],
        n_replicates=100,
        output_dir=tmp_path,
    )
    path = tmp_path / "multiseed_metadata.json"
    metadata = json.loads(path.read_text(encoding="utf-8"))
    metadata[field] = value
    path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        audit_multiseed_sensitivity(tmp_path)
