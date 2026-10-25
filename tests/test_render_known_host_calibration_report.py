from __future__ import annotations

import importlib.util
import hashlib
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
        "original_commit": "b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f",
        "original_repository": "vedantkushwah033-blip/bayesian-dark-siren-cosmology",
        "true_h0_values": list(MODULE.TRUTHS),
        "fractional_distance_errors": list(MODULE.NOISES),
        "posterior_grid_points": 1601,
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
    assert "Mean posterior-median error" in report
    assert "RMSE of posterior median" in report


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


def test_report_requires_matching_audit_pass_and_exact_input_hashes(tmp_path):
    metadata = {"original_commit": "abc123"}
    contents = {
        "per_run_results.csv": b"raw\n",
        "cell_summary.csv": b"summary\n",
        "cell_summary.json": b"{}\n",
    }
    hashes = {}
    for name, data in contents.items():
        (tmp_path / name).write_bytes(data)
        hashes[name] = hashlib.sha256(data).hexdigest()
    audit = {
        "audit": "PASS",
        "posterior_summaries_independently_reconstructed_from_raw_observations": True,
        "raw_records": 3600,
        "summary_cells": 36,
        "source_simulator_commit": "abc123",
        "input_sha256": hashes,
    }
    MODULE.validate_audit_report(tmp_path, audit, metadata, 3600)
    (tmp_path / "cell_summary.csv").write_text("tampered\n", encoding="utf-8")
    with pytest.raises(ValueError, match="changed since independent audit"):
        MODULE.validate_audit_report(tmp_path, audit, metadata, 3600)


def test_report_refuses_audit_without_independent_posterior_reconstruction(tmp_path):
    audit = {"audit": "PASS"}
    with pytest.raises(ValueError, match="did not independently reconstruct"):
        MODULE.validate_audit_report(tmp_path, audit, {"original_commit": "abc123"}, 3600)


def test_report_rejects_replicate_count_or_design_metadata_mismatch():
    rows, metadata = fixture_summary()
    rows[0]["n"] = "99"
    with pytest.raises(ValueError, match="does not match prespecified"):
        MODULE.validate_summary(rows, metadata)
    rows, metadata = fixture_summary()
    metadata["fractional_distance_errors"] = [0.05, 0.10, 0.15, 0.25]
    with pytest.raises(ValueError, match="noise grid"):
        MODULE.validate_summary(rows, metadata)
