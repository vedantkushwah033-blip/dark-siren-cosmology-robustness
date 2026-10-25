import importlib.util
import math
from pathlib import Path

import pytest


PATH = Path(__file__).resolve().parents[1] / "scripts" / "audit_known_host_calibration_surface.py"
SPEC = importlib.util.spec_from_file_location("audit_known_host_calibration_surface", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


RUNNER_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_known_host_calibration_surface.py"
RUNNER_SPEC = importlib.util.spec_from_file_location("run_known_host_calibration_surface", RUNNER_PATH)
RUNNER_MODULE = importlib.util.module_from_spec(RUNNER_SPEC)
RUNNER_SPEC.loader.exec_module(RUNNER_MODULE)


def test_runner_wilson_interval_has_exact_endpoint_bounds():
    lower_zero, upper_zero = RUNNER_MODULE.wilson_interval(0, 300)
    lower_all, upper_all = RUNNER_MODULE.wilson_interval(300, 300)
    assert lower_zero == 0.0
    assert 0.0 < upper_zero < 0.05
    assert 0.95 < lower_all < 1.0
    assert upper_all == 1.0


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
    assert lo0 == 0.0
    assert 0 < hi0 < 0.05
    assert 0.95 < lo1 < 1
    assert hi1 == 1.0


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
    for truth_index, truth in enumerate(truths):
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
    assert report["posterior_summaries_independently_reconstructed_from_raw_observations"] is True

    # Deliberately corrupt one posterior median and confirm the raw-data audit fails.
    raw_path = tmp_path / "per_run_results.csv"
    with raw_path.open(newline="", encoding="utf-8") as handle:
        corrupted_rows = list(csv.DictReader(handle))
        fieldnames = list(corrupted_rows[0])
    corrupted_rows[0]["median"] = str(float(corrupted_rows[0]["median"]) + 0.25)
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(corrupted_rows)
    with pytest.raises(ValueError, match="independent posterior median"):
        MODULE.main()


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


def test_independent_posterior_matches_separate_scipy_reference():
    import numpy as np
    from scipy.stats import norm

    row = synthetic_row(0, 70.0, 1.25)
    row["observed_distance_mpc"] = str(float(row["true_distance_mpc"]) * 1.13)
    reconstructed = MODULE.recompute_posterior_from_record(row)

    # Separate reference implementation using the simulator's normal-PDF
    # likelihood form rather than the auditor's log-likelihood calculation.
    grid = np.linspace(50.0, 90.0, 1601)
    predicted_distance = 299792.458 * float(row["host_redshift"]) / grid
    likelihood = norm.pdf(
        float(row["observed_distance_mpc"]),
        loc=predicted_distance,
        scale=float(row["distance_sigma_mpc"]),
    )
    likelihood /= np.trapezoid(likelihood, grid)
    mean = float(np.trapezoid(grid * likelihood, grid))
    sd = float(np.sqrt(np.trapezoid((grid - mean) ** 2 * likelihood, grid)))
    cdf = np.concatenate([[0.0], np.cumsum((likelihood[1:] + likelihood[:-1]) * np.diff(grid) / 2)])
    cdf /= cdf[-1]
    reference = {
        "mean": mean,
        "median": float(np.interp(0.5, cdf, grid)),
        "posterior_sd": sd,
        "lower_68": float(np.interp(0.16, cdf, grid)),
        "upper_68": float(np.interp(0.84, cdf, grid)),
        "lower_95": float(np.interp(0.025, cdf, grid)),
        "upper_95": float(np.interp(0.975, cdf, grid)),
    }
    for field, expected in reference.items():
        assert reconstructed[field] == pytest.approx(expected, rel=2e-6, abs=2e-6)


def test_auditor_rejects_host_redshift_outside_protocol_range():
    rows = [synthetic_row(0, 70.0, 0.0)]
    rows[0]["host_redshift"] = "0.2"
    with pytest.raises(ValueError, match="host redshift"):
        MODULE.recompute(rows)


def test_auditor_rejects_distance_sigma_inconsistent_with_noise_design():
    rows = [synthetic_row(0, 70.0, 0.0)]
    rows[0]["distance_sigma_mpc"] = str(float(rows[0]["distance_sigma_mpc"]) * 1.1)
    with pytest.raises(ValueError, match="distance sigma"):
        MODULE.recompute(rows)


@pytest.mark.parametrize("truth", [50.0, 90.0])
def test_independent_posterior_respects_prior_at_boundary_truths(truth):
    row = synthetic_row(0, truth, 0.0)
    row["true_h0"] = str(truth)
    row["true_distance_mpc"] = str(299792.458 * float(row["host_redshift"]) / truth)
    row["distance_sigma_mpc"] = str(0.2 * float(row["true_distance_mpc"]))
    row["observed_distance_mpc"] = row["true_distance_mpc"]
    result = MODULE.recompute_posterior_from_record(row)
    assert 50.0 <= result["lower_95"] <= result["lower_68"] <= result["median"]
    assert result["median"] <= result["upper_68"] <= result["upper_95"] <= 90.0
    assert result["posterior_sd"] > 0.0


@pytest.mark.parametrize("field,value", [
    ("observed_distance_mpc", "0"),
    ("observed_distance_mpc", "-1"),
    ("distance_sigma_mpc", "0"),
    ("distance_sigma_mpc", "-0.1"),
    ("host_redshift", "nan"),
])
def test_independent_posterior_rejects_invalid_raw_observation(field, value):
    row = synthetic_row(0, 70.0, 0.0)
    row[field] = value
    with pytest.raises(ValueError, match="invalid raw observation"):
        MODULE.recompute_posterior_from_record(row)


def test_independent_posterior_rejects_nonfinite_observed_distance():
    row = synthetic_row(0, 70.0, 0.0)
    row["observed_distance_mpc"] = "inf"
    with pytest.raises(ValueError, match="invalid raw observation"):
        MODULE.recompute_posterior_from_record(row)


def test_metadata_cell_count_validation_accepts_exact_json_types():
    expected = {
        "n": 100,
        "coverage_68_successes": 68,
        "coverage_95_successes": 95,
        "truth_on_prior_boundary": False,
    }
    metadata_cell = dict(expected)
    MODULE.validate_metadata_cell_counts(metadata_cell, expected, (70.0, 0.1))


@pytest.mark.parametrize(
    "field,value,message",
    [
        ("n", True, "must be an integer"),
        ("coverage_68_successes", 68.0, "must be an integer"),
        ("coverage_95_successes", "95", "must be an integer"),
        ("truth_on_prior_boundary", "false", "must be boolean"),
        ("truth_on_prior_boundary", 0, "must be boolean"),
    ],
)
def test_metadata_cell_validation_rejects_coerced_json_types(field, value, message):
    expected = {
        "n": 100,
        "coverage_68_successes": 68,
        "coverage_95_successes": 95,
        "truth_on_prior_boundary": False,
    }
    metadata_cell = dict(expected)
    metadata_cell[field] = value
    with pytest.raises(ValueError, match=message):
        MODULE.validate_metadata_cell_counts(metadata_cell, expected, (70.0, 0.1))


def test_metadata_cell_validation_rejects_wrong_integer_count():
    expected = {
        "n": 100,
        "coverage_68_successes": 68,
        "coverage_95_successes": 95,
        "truth_on_prior_boundary": False,
    }
    metadata_cell = dict(expected)
    metadata_cell["coverage_68_successes"] = 67
    with pytest.raises(ValueError, match="count coverage_68_successes mismatch"):
        MODULE.validate_metadata_cell_counts(metadata_cell, expected, (70.0, 0.1))


def test_base_seed_validation_accepts_nonnegative_integer():
    assert MODULE.validate_base_seed(20261025) == 20261025
    assert MODULE.validate_base_seed(0) == 0


@pytest.mark.parametrize("value", [-1, True, 1.0, "20261025", None])
def test_base_seed_validation_rejects_negative_or_coercible_nonintegers(value):
    with pytest.raises(ValueError, match="non-negative integer"):
        MODULE.validate_base_seed(value)


@pytest.mark.parametrize("value", ["0", "12", "0007"])
def test_csv_integer_parser_accepts_ascii_nonnegative_integers(value):
    assert MODULE.parse_nonnegative_integer(value, "test field") == int(value)


@pytest.mark.parametrize("value", ["-1", "+1", "1.0", "1e2", " 1", "1 ", "١", "", None, 1, True])
def test_csv_integer_parser_rejects_coercible_or_malformed_values(value):
    with pytest.raises(ValueError, match="must be a non-negative integer"):
        MODULE.parse_nonnegative_integer(value, "test field")


def valid_protocol_metadata():
    return {
        "seed": 20261025,
        "raw_row_count": 3600,
        "cell_count": 36,
        "replicates_per_cell": 100,
        "prior": [50.0, 90.0],
        "posterior_grid_points": 1601,
        "original_repository": "vedantkushwah033-blip/bayesian-dark-siren-cosmology",
        "original_commit": "b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f",
        "true_h0_values": [50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0],
        "fractional_distance_errors": [0.05, 0.10, 0.15, 0.20],
        "cells": [{} for _ in range(36)],
    }


def test_protocol_metadata_validator_accepts_valid_design():
    assert MODULE.validate_metadata_seed_and_design(valid_protocol_metadata(), 100) == 20261025


@pytest.mark.parametrize(
    "field,value,message",
    [
        ("seed", True, "base seed"),
        ("seed", "20261025", "base seed"),
        ("raw_row_count", True, "raw_row_count must be an integer"),
        ("cell_count", 36.0, "cell_count must be an integer"),
        ("replicates_per_cell", "100", "replicates_per_cell must be an integer"),
        ("prior", [0.0, 90.0], "unexpected prior"),
        ("posterior_grid_points", 1600, "unexpected prior"),
        ("original_commit", "un-pinned", "pinned protocol"),
        ("true_h0_values", [50, 55, 60, 65, 70, 75, 80, 85, 86], "truth grid"),
        ("fractional_distance_errors", [0.05, 0.10, 0.15, 0.21], "noise grid"),
        ("cells", [], "exactly 36 cell summaries"),
    ],
)
def test_protocol_metadata_validator_rejects_malformed_or_changed_design(field, value, message):
    metadata = valid_protocol_metadata()
    metadata[field] = value
    with pytest.raises(ValueError, match=message):
        MODULE.validate_metadata_seed_and_design(metadata, 100)


def test_protocol_metadata_validator_rejects_requested_replicate_mismatch():
    with pytest.raises(ValueError, match="row/cell counts"):
        MODULE.validate_metadata_seed_and_design(valid_protocol_metadata(), 300)


def test_protocol_metadata_validator_rejects_non_object_metadata():
    with pytest.raises(ValueError, match="metadata must be a JSON object"):
        MODULE.validate_metadata_seed_and_design([], 100)


def test_protocol_metadata_validator_rejects_non_object_cells():
    metadata = valid_protocol_metadata()
    metadata["cells"] = [{} for _ in range(35)] + ["not-an-object"]
    with pytest.raises(ValueError, match="each JSON metadata cell must be an object"):
        MODULE.validate_metadata_seed_and_design(metadata, 100)


def test_protocol_metadata_validator_rejects_malformed_cell_entry():
    metadata = valid_protocol_metadata()
    metadata["cells"][0] = "not-an-object"
    with pytest.raises(ValueError, match="each JSON metadata cell must be an object"):
        MODULE.validate_metadata_seed_and_design(metadata, 100)

def test_auditor_rejects_true_distance_inconsistent_with_low_z_truth():
    rows = [synthetic_row(0, 70.0, 0.0)]
    rows[0]["true_distance_mpc"] = str(float(rows[0]["true_distance_mpc"]) * 1.01)
    # Keep the uncertainty tied to the altered distance to isolate this check.
    rows[0]["distance_sigma_mpc"] = str(
        0.1 * float(rows[0]["true_distance_mpc"])
    )
    with pytest.raises(ValueError, match="low-z truth distance"):
        MODULE.recompute(rows)

def paired_rows(noise_by_fraction):
    truth = 70.0
    z = 0.03
    distance = 299792.458 * z / truth
    rows = []
    for frac, noise in zip((0.05, 0.10, 0.15, 0.20), noise_by_fraction):
        sigma = frac * distance
        rows.append({
            "host_redshift": str(z),
            "true_distance_mpc": str(distance),
            "observed_distance_mpc": str(distance + sigma * noise),
            "distance_sigma_mpc": str(sigma),
        })
    return rows


def test_paired_group_accepts_shared_host_and_standardized_noise():
    MODULE.validate_paired_group((70.0, 3), paired_rows((-1.25, -1.25, -1.25, -1.25)))


def test_paired_group_rejects_unpaired_standardized_noise():
    with pytest.raises(ValueError, match="standardized noise is not paired"):
        MODULE.validate_paired_group((70.0, 3), paired_rows((-1.25, -1.25, 0.5, -1.25)))


def test_paired_group_rejects_unpaired_host_population():
    rows = paired_rows((0.25, 0.25, 0.25, 0.25))
    rows[2]["host_redshift"] = str(float(rows[2]["host_redshift"]) + 0.001)
    with pytest.raises(ValueError, match="host/population not paired"):
        MODULE.validate_paired_group((70.0, 3), rows)


def test_paired_group_rejects_missing_uncertainty_level():
    with pytest.raises(ValueError, match="has 3 noise levels"):
        MODULE.validate_paired_group((70.0, 3), paired_rows((0.25, 0.25, 0.25, 0.25))[:3])

@pytest.mark.parametrize(
    "field,value,message",
    [
        ("host_redshift", "nan", "non-finite observation values"),
        ("true_distance_mpc", "inf", "non-finite observation values"),
        ("observed_distance_mpc", "nan", "non-finite observation values"),
        ("distance_sigma_mpc", "inf", "non-finite observation values"),
        ("distance_sigma_mpc", "0", "non-positive distance uncertainty"),
        ("distance_sigma_mpc", "-0.1", "non-positive distance uncertainty"),
    ],
)
def test_paired_group_rejects_nonfinite_or_invalid_observation_values(field, value, message):
    rows = paired_rows((0.25, 0.25, 0.25, 0.25))
    rows[1][field] = value
    with pytest.raises(ValueError, match=message):
        MODULE.validate_paired_group((70.0, 3), rows)


# The runner and the independent auditor intentionally implement their own
# summary paths. Keep regression coverage for the runner's cell statistics too.
RUNNER_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_known_host_calibration_surface.py"
RUNNER_SPEC = importlib.util.spec_from_file_location("run_known_host_calibration_surface", RUNNER_PATH)
RUNNER = importlib.util.module_from_spec(RUNNER_SPEC)
RUNNER_SPEC.loader.exec_module(RUNNER)


def runner_summary_row(median, mean, lower68, upper68, lower95, upper95, sd=2.0):
    return {
        "median": median, "mean": mean,
        "lower_68": lower68, "upper_68": upper68,
        "lower_95": lower95, "upper_95": upper95,
        "posterior_sd": sd,
    }


def test_runner_summarize_cell_reports_coverage_and_point_error():
    rows = [
        runner_summary_row(69.0, 69.2, 68.0, 70.0, 67.0, 71.0),
        runner_summary_row(70.0, 70.1, 69.0, 71.0, 68.0, 72.0),
        runner_summary_row(72.0, 71.8, 71.0, 73.0, 70.0, 74.0),
        runner_summary_row(71.0, 70.8, 70.0, 72.0, 69.0, 73.0),
    ]
    result = RUNNER.summarize_cell(rows, true_h0=70.0, fractional_error=0.1)
    assert result["n"] == 4
    assert result["coverage_68_successes"] == 3
    assert result["coverage_68"] == pytest.approx(0.75)
    assert result["coverage_95_successes"] == 4
    assert result["coverage_95"] == pytest.approx(1.0)
    assert result["mean_median_error"] == pytest.approx(0.5)
    assert result["rmse_median"] == pytest.approx(math.sqrt(1.5))
    assert result["mean_width_68"] == pytest.approx(2.0)
    assert result["mean_width_95"] == pytest.approx(4.0)


@pytest.mark.parametrize(
    "bad_row",
    [
        runner_summary_row(float("nan"), 70.0, 69.0, 71.0, 68.0, 72.0),
        runner_summary_row(70.0, float("inf"), 69.0, 71.0, 68.0, 72.0),
        runner_summary_row(70.0, 70.0, 72.0, 71.0, 68.0, 73.0),
        runner_summary_row(70.0, 70.0, 69.0, 71.0, 71.0, 72.0),
        runner_summary_row(70.0, 70.0, 49.0, 71.0, 48.0, 72.0),
        runner_summary_row(70.0, 70.0, 69.0, 71.0, 68.0, 72.0, sd=0.0),
    ],
)
def test_runner_summarize_cell_rejects_invalid_posterior_summaries(bad_row):
    with pytest.raises(ValueError):
        RUNNER.summarize_cell([bad_row], true_h0=70.0, fractional_error=0.1)


def test_runner_summarize_cell_rejects_empty_cell():
    with pytest.raises(ValueError, match="must not be empty"):
        RUNNER.summarize_cell([], true_h0=70.0, fractional_error=0.1)
