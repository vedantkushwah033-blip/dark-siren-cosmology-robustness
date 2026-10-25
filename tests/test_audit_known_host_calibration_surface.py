import importlib.util
import math
from pathlib import Path

import pytest


PATH = Path(__file__).resolve().parents[1] / "scripts" / "audit_known_host_calibration_surface.py"
SPEC = importlib.util.spec_from_file_location("audit_known_host_calibration_surface", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def synthetic_row(run, median, standardized_noise):
    truth = 70.0
    frac = 0.1
    z = 0.03
    distance = 299792.458 * z / truth
    sigma = frac * distance
    return {
        "run": str(run),
        "true_h0": str(truth),
        "fractional_distance_error": str(frac),
        "host_redshift": str(z),
        "true_distance_mpc": str(distance),
        "observed_distance_mpc": str(distance + sigma * standardized_noise),
        "distance_sigma_mpc": str(sigma),
        "median": str(median),
        "mean": str(median),
        "posterior_sd": "2.0",
        "lower_68": str(median - 1.0),
        "upper_68": str(median + 1.0),
        "lower_95": str(median - 2.0),
        "upper_95": str(median + 2.0),
        "width_68": "2.0",
        "width_95": "4.0",
        "pull_median": str((median - truth) / 2.0),
    }


def test_auditor_wilson_interval_handles_zero_and_all_successes():
    lo0, hi0 = MODULE.wilson(0, 100)
    lo1, hi1 = MODULE.wilson(100, 100)
    assert lo0 == pytest.approx(0.0)
    assert 0 < hi0 < 0.05
    assert 0.95 < lo1 < 1
    assert hi1 == pytest.approx(1.0)


def test_auditor_independently_recomputes_coverage_and_errors():
    rows = [
        synthetic_row(0, 69.0, -1.0),
        synthetic_row(1, 70.0, 0.0),
        synthetic_row(2, 71.0, 1.0),
        synthetic_row(3, 70.0, 0.5),
    ]
    result = MODULE.recompute(rows)
    assert result["n"] == 4
    assert result["coverage_68"] == pytest.approx(1.0)
    assert result["coverage_95"] == pytest.approx(1.0)
    assert result["median_bias"] == pytest.approx(0.0)
    assert result["mean_median_error"] == pytest.approx(0.0)
    assert result["mean_width_68"] == pytest.approx(2.0)
    assert result["mean_width_95"] == pytest.approx(4.0)
    assert result["truth_on_prior_boundary"] is False


def test_auditor_rejects_disagreement_in_summary_numbers():
    with pytest.raises(ValueError, match="stored"):
        MODULE.close(0.5, 0.6, "coverage")


def test_auditor_rejects_credible_interval_that_excludes_its_median():
    rows = [synthetic_row(0, 70.0, 0.0)]
    rows[0]["lower_68"] = "71.0"
    rows[0]["upper_68"] = "72.0"
    with pytest.raises(ValueError, match="median"):
        MODULE.recompute(rows)


def test_auditor_rejects_nonpositive_distance_uncertainty():
    rows = [synthetic_row(0, 70.0, 0.0)]
    rows[0]["distance_sigma_mpc"] = "0.0"
    with pytest.raises(ValueError):
        MODULE.recompute(rows)


def test_auditor_rejects_nonfinite_raw_distance():
    rows = [synthetic_row(0, 70.0, 0.0)]
    rows[0]["observed_distance_mpc"] = "nan"
    with pytest.raises(ValueError, match="non-finite"):
        MODULE.recompute(rows)


def test_auditor_end_to_end_recomputes_and_writes_provenance(tmp_path, monkeypatch):
    import csv
    import json
    import sys

    truths = [50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0]
    noises = [0.05, 0.10, 0.15, 0.20]
    replicates = 100
    base_seed = 20261025
    raw_rows = []
    summary_rows = []
    for truth in truths:
        for frac in noises:
            cell = []
            for run in range(replicates):
                z = 0.03 + run * 0.00001
                distance = 299792.458 * z / truth
                sigma = frac * distance
                standardized_noise = ((run % 9) - 4) / 3.0
                median = truth
                lower68, upper68 = max(50.0, median - 1.0), min(90.0, median + 1.0)
                lower95, upper95 = max(50.0, median - 2.0), min(90.0, median + 2.0)
                row = {
                    "run": run, "seed": base_seed + run, "true_h0": truth,
                    "fractional_distance_error": frac, "host_redshift": z,
                    "true_distance_mpc": distance,
                    "observed_distance_mpc": distance + sigma * standardized_noise,
                    "distance_sigma_mpc": sigma, "median": median, "mean": median,
                    "posterior_sd": 2.0, "lower_68": lower68, "upper_68": upper68,
                    "lower_95": lower95, "upper_95": upper95,
                    "width_68": upper68 - lower68, "width_95": upper95 - lower95,
                    "pull_median": 0.0,
                }
                # Build a self-consistent fixture using the independently specified likelihood.
                row.update(MODULE.recompute_posterior_from_record({k: str(v) for k, v in row.items()}))
                row["width_68"] = row["upper_68"] - row["lower_68"]
                row["width_95"] = row["upper_95"] - row["lower_95"]
                row["pull_median"] = (row["median"] - truth) / row["posterior_sd"]
                raw_rows.append(row)
                cell.append({k: str(v) for k, v in row.items()})
            summary_rows.append(MODULE.recompute(cell))

    with (tmp_path / "per_run_results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(raw_rows[0]))
        writer.writeheader()
        writer.writerows(raw_rows)
    with (tmp_path / "cell_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    metadata = {
        "true_h0_values": truths, "fractional_distance_errors": noises,
        "seed": base_seed, "prior": [50.0, 90.0], "posterior_grid_points": 1601,
        "raw_row_count": len(raw_rows), "cell_count": len(summary_rows),
        "replicates_per_cell": replicates,
        "original_repository": "vedantkushwah033-blip/bayesian-dark-siren-cosmology",
        "original_commit": "b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f",
        "cells": summary_rows,
    }
    (tmp_path / "cell_summary.json").write_text(json.dumps(metadata), encoding="utf-8")
    monkeypatch.setattr(sys, "argv", [
        "audit_known_host_calibration_surface.py", str(tmp_path), "--replicates", str(replicates)
    ])
    monkeypatch.setenv("GITHUB_RUN_ID", "test-run-123")
    monkeypatch.setenv("GITHUB_SHA", "test-sha")
    monkeypatch.setenv("GITHUB_REF", "refs/heads/test")
    MODULE.main()

    report = json.loads((tmp_path / "audit_report.json").read_text(encoding="utf-8"))
    assert report["audit"] == "PASS"
    assert report["raw_records"] == 9 * 4 * replicates
    assert report["summary_cells"] == 36
    assert report["workflow_run_id"] == "test-run-123"
    assert report["workflow_commit_sha"] == "test-sha"
    assert len(report["input_sha256"]["per_run_results.csv"]) == 64
    assert len(report["auditor_script_sha256"]) == 64


def test_auditor_rejects_mixed_truth_within_cell():
    rows = [synthetic_row(0, 70.0, 0.0), synthetic_row(1, 70.0, 0.0)]
    rows[1]["true_h0"] = "71.0"
    with pytest.raises(ValueError, match="truth changes"):
        MODULE.recompute(rows)


def test_auditor_rejects_unprespecified_noise_level():
    rows = [synthetic_row(0, 70.0, 0.0)]
    rows[0]["fractional_distance_error"] = "0.12"
    with pytest.raises(ValueError, match="unexpected fractional"):
        MODULE.recompute(rows)
