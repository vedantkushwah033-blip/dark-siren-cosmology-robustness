from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "render_known_host_calibration_report.py"
SPEC = importlib.util.spec_from_file_location("render_known_host_calibration_report", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def fixture_summary():
    rows = []
    for truth in MODULE.TRUTHS:
        for noise in MODULE.NOISES:
            rows.append({
                "true_h0": str(truth),
                "fractional_distance_error": str(noise),
                "n": "100",
                "coverage_68": "0.68",
                "coverage_68_wilson_low": "0.58",
                "coverage_68_wilson_high": "0.77",
                "coverage_95": "0.95",
                "coverage_95_wilson_low": "0.88",
                "coverage_95_wilson_high": "0.98",
                "mean_median_error": "0.2",
                "rmse_median": "2.0",
                "mean_width_68": "5.0",
                "mean_width_95": "10.0",
                "median_near_prior_edge_fraction": "0.01",
            })
    metadata = {
        "cell_count": 36,
        "raw_row_count": 3600,
        "replicates_per_cell": 100,
        "prior": [50.0, 90.0],
        "original_commit": "abc123",
        "seed": 20261025,
    }
    return rows, metadata


def test_report_contains_all_cells_and_boundary_comparison():
    rows, metadata = fixture_summary()
    report = MODULE.render_report(rows, metadata)
    assert report.count("| 50 | 5% |") == 1
    assert report.count("| 90 | 20% |") == 1
    assert "Boundary (50, 90)" in report
    assert "Interior (55–85)" in report
    assert "not an observational gravitational-wave result" in report
    assert "Wilson" in report


def test_report_rejects_missing_or_duplicate_cells():
    rows, metadata = fixture_summary()
    with pytest.raises(ValueError, match="expected 36 design cells"):
        MODULE.validate_summary(rows[:-1], metadata)
    rows, metadata = fixture_summary()
    rows[-1] = dict(rows[0])
    with pytest.raises(ValueError, match="duplicate summary cell"):
        MODULE.validate_summary(rows, metadata)


def test_report_rejects_invalid_coverage_or_wilson_interval():
    rows, metadata = fixture_summary()
    rows[0]["coverage_68_wilson_low"] = "0.9"
    with pytest.raises(ValueError, match="invalid coverage/Wilson"):
        MODULE.validate_summary(rows, metadata)


def test_report_rejects_metadata_count_mismatch():
    rows, metadata = fixture_summary()
    metadata["raw_row_count"] = 1
    with pytest.raises(ValueError, match="metadata counts"):
        MODULE.validate_summary(rows, metadata)
