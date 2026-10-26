"""Regression tests for the independent grid-sensitivity comparison audit."""
import csv
import json

import pytest

from scripts.audit_calibration_grid_sensitivity import audit_grid_sensitivity
from scripts.run_calibration_grid_sensitivity import run_sensitivity


def test_audit_reconstructs_grid_comparison_and_rejects_tampering(tmp_path):
    run_sensitivity([1201, 2401], n_replicates=100, seed=24680, output_dir=tmp_path)
    result = audit_grid_sensitivity(tmp_path)
    assert result["status"] == "PASS"
    assert result["paired_raw_records_identical"] is True
    assert result["comparison_statistics_reconstructed"] is True

    path = tmp_path / "grid_sensitivity.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        fields = list(rows[0])
    rows[0]["rms_pit_delta_vs_finest"] = str(float(rows[0]["rms_pit_delta_vs_finest"]) + 0.01)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="rms_pit_delta_vs_finest"):
        audit_grid_sensitivity(tmp_path)


def test_audit_rejects_duplicate_grid_rows(tmp_path):
    run_sensitivity([1201, 2401], n_replicates=100, seed=24681, output_dir=tmp_path)
    path = tmp_path / "grid_sensitivity.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        fields = list(rows[0])
    rows[1]["grid_points"] = rows[0]["grid_points"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="unique"):
        audit_grid_sensitivity(tmp_path)
