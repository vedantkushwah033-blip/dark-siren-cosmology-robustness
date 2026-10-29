"""Regression tests for the independent distance-uncertainty grid-stability check."""
from __future__ import annotations

import csv

import pytest

from scripts.check_distance_uncertainty_grid_stability import check_grid_stability
from dark_siren_robustness.uncertainty_sensitivity import run_experiment


def _write_raw_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def test_grid_stability_reconstructs_runner_output(tmp_path):
    raw, _ = run_experiment(
        replicates=2,
        seed=731,
        grid_points=1601,
        truths=(70.0,),
        fractional_errors=(0.05,),
    )
    raw_path = tmp_path / "replicates.csv"
    _write_raw_csv(raw_path, raw)

    report = check_grid_stability(raw_path, fine_grid_points=3201)

    assert report["status"] == "PASS"
    assert report["raw_record_count"] == 6
    assert report["method"].startswith("independent Gaussian-likelihood")
    assert all(item["coverage_unchanged"] for item in report["coverage_68_by_arm"].values())
    assert all(
        item["maximum_absolute_change"] <= 1e-3
        for item in report["quantile_changes"].values()
    )


def test_grid_stability_rejects_material_coarse_quantile_error(tmp_path):
    raw, _ = run_experiment(
        replicates=2,
        seed=732,
        grid_points=1601,
        truths=(70.0,),
        fractional_errors=(0.05,),
    )
    raw[0]["median"] = str(float(raw[0]["median"]) + 0.01)
    raw_path = tmp_path / "replicates.csv"
    _write_raw_csv(raw_path, raw)

    report = check_grid_stability(raw_path, fine_grid_points=3201)

    assert report["status"] == "FAIL"
    assert report["grid_violation_count"] >= 1
    assert any(item["summary"] == "median" for item in report["grid_violations"])



def test_grid_stability_rejects_duplicate_raw_record(tmp_path):
    raw, _ = run_experiment(
        replicates=2,
        seed=733,
        grid_points=1601,
        truths=(70.0,),
        fractional_errors=(0.05,),
    )
    raw.append(dict(raw[0]))
    raw_path = tmp_path / "replicates.csv"
    _write_raw_csv(raw_path, raw)

    with pytest.raises(ValueError, match="duplicate raw record key"):
        check_grid_stability(raw_path, fine_grid_points=3201)

def test_grid_stability_rejects_missing_arm(tmp_path):
    raw, _ = run_experiment(
        replicates=2, seed=734, grid_points=1601,
        truths=(70.0,), fractional_errors=(0.05,),
    )
    raw = [row for row in raw if row["arm"] != "C_fixed_reference_sigma" or int(row["replicate"]) != 0]
    raw_path = tmp_path / "replicates.csv"
    _write_raw_csv(raw_path, raw)

    with pytest.raises(ValueError, match="incomplete or unexpected likelihood arms"):
        check_grid_stability(raw_path, fine_grid_points=3201)


def test_grid_stability_rejects_missing_replicate_from_one_cell(tmp_path):
    raw, _ = run_experiment(
        replicates=2, seed=735, grid_points=1601,
        truths=(55.0, 70.0), fractional_errors=(0.05,),
    )
    raw = [row for row in raw if not (
        float(row["true_h0"]) == 55.0 and int(row["replicate"]) == 0
    )]
    raw_path = tmp_path / "replicates.csv"
    _write_raw_csv(raw_path, raw)

    with pytest.raises(ValueError, match="replicate IDs must be contiguous|inconsistent replicate IDs"):
        check_grid_stability(raw_path, fine_grid_points=3201)


def test_grid_stability_rejects_even_fine_grid(tmp_path):
    with pytest.raises(ValueError, match="odd integer"):
        check_grid_stability(tmp_path / "not_read.csv", fine_grid_points=3200)


def test_grid_stability_rejects_missing_raw_columns(tmp_path):
    raw_path = tmp_path / "incomplete.csv"
    raw_path.write_text("true_h0,arm\n70,A_fixed_realized_sigma\n", encoding="utf-8")

    with pytest.raises(ValueError, match="raw CSV missing columns"):
        check_grid_stability(raw_path, fine_grid_points=3201)



def test_posterior_quantile_reconstruction_flags_boundary_mass():
    import numpy as np

    from scripts.check_distance_uncertainty_grid_stability import _posterior_quantiles

    c = 299_792.458
    z = 0.02
    truth_at_boundary_distance = c * z / 50.0
    row = {
        "host_redshift": str(z),
        "observed_distance_mpc": str(truth_at_boundary_distance),
        "fractional_error": "0.004",
        "generated_sigma_mpc": str(0.004 * truth_at_boundary_distance),
        "arm": "A_fixed_realized_sigma",
    }
    summary = _posterior_quantiles(row, np.linspace(50.0, 90.0, 3201))

    assert summary["boundary_mass_within_one_coarse_spacing"] > 0.05


def test_multi_resolution_sweep_reports_adjacent_grid_changes(tmp_path):
    from scripts.check_distance_uncertainty_grid_stability import compare_grid_resolutions

    raw, _ = run_experiment(
        replicates=2,
        seed=736,
        grid_points=1601,
        truths=(55.0, 70.0, 85.0),
        fractional_errors=(0.05,),
    )
    raw_path = tmp_path / "replicates.csv"
    _write_raw_csv(raw_path, raw)

    report = compare_grid_resolutions(raw_path, grid_points=(801, 1601, 3201, 6401))

    assert report["raw_record_count"] == 18
    assert report["grid_points"] == [801, 1601, 3201, 6401]
    assert len(report["comparisons"]) == 3
    assert report["comparisons"][-1]["coarse_grid_points"] == 3201
    assert report["comparisons"][-1]["fine_grid_points"] == 6401
    assert set(report["comparisons"][0]["summary_changes"]) == {
        "median", "lower_68", "upper_68", "lower_95", "upper_95"
    }
    assert "changed_68_coverage_classifications" in report["comparisons"][0]["coverage_classification_changes_by_arm"]["A_fixed_realized_sigma"]


def test_multi_resolution_sweep_rejects_unsorted_or_even_grids(tmp_path):
    from scripts.check_distance_uncertainty_grid_stability import compare_grid_resolutions

    with pytest.raises(ValueError, match="strictly increasing"):
        compare_grid_resolutions(tmp_path / "missing.csv", grid_points=(1601, 801))
    with pytest.raises(ValueError, match="odd integers"):
        compare_grid_resolutions(tmp_path / "missing.csv", grid_points=(801, 1600))
