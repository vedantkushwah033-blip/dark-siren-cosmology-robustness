"""Regression tests for the selected-likelihood PIT calibration runner."""
import csv
import json

import numpy as np
import pytest

from scripts.audit_selected_likelihood_sbc import audit

from scripts.run_selected_likelihood_sbc import run_calibration, wilson_interval


def test_wilson_interval_has_valid_order_and_contains_observed_fraction():
    low, high = wilson_interval(50, 100)
    assert 0.0 <= low <= 0.5 <= high <= 1.0


def test_calibration_runner_records_finite_reproducible_pits_and_complete_bins():
    rows, summary = run_calibration(n_replicates=100, seed=8675309)

    assert len(rows) == 100
    assert [row["replicate"] for row in rows] == list(range(100))
    assert all(row["seed"] == 8675309 for row in rows)
    assert all(0.0 <= row["pit"] <= 1.0 for row in rows)
    assert all(row["rejection_attempts"] >= 1 for row in rows)
    assert np.isfinite(summary["pit_mean"])
    assert np.isfinite(summary["pit_mean_monte_carlo_se"])
    assert summary["pit_mean_monte_carlo_se"] > 0.0
    assert sum(item["count"] for item in summary["decile_bins_with_wilson_95_intervals"]) == 100
    assert all(
        0.0 <= item["wilson_95_low"] <= item["frequency"] <= item["wilson_95_high"] <= 1.0
        for item in summary["decile_bins_with_wilson_95_intervals"]
    )

    repeated_rows, repeated_summary = run_calibration(n_replicates=100, seed=8675309)
    assert [row["pit"] for row in repeated_rows] == [row["pit"] for row in rows]
    assert repeated_summary["pit_mean"] == summary["pit_mean"]
    assert repeated_summary["ks_statistic"] == summary["ks_statistic"]


def test_calibration_runner_rejects_too_small_replicate_count():
    with pytest.raises(ValueError, match="at least 100"):
        run_calibration(n_replicates=99, seed=1)


def test_calibration_runner_rejects_invalid_grid_resolution():
    with pytest.raises(ValueError, match="at least two points"):
        run_calibration(n_replicates=100, seed=1, grid_points=1)



def test_independent_audit_reconstructs_raw_pits_and_rejects_tampering(tmp_path):
    rows, summary = run_calibration(n_replicates=100, seed=246810)
    raw_path = tmp_path / "pit_replicates.csv"
    summary_path = tmp_path / "pit_summary.json"

    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    result = audit(raw_path, summary_path)
    assert result["status"] == "PASS"
    assert result["n_replicates_audited"] == 100

    # The audit must verify the seeded generation process, not only the PIT.
    altered_truth_rows = [dict(row) for row in rows]
    altered_truth_rows[0]["truth_h0"] = float(altered_truth_rows[0]["truth_h0"]) + 0.01
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(altered_truth_rows[0]))
        writer.writeheader()
        writer.writerows(altered_truth_rows)
    with pytest.raises(ValueError, match="generation replay mismatch"):
        audit(raw_path, summary_path)

    # Restore the original records before challenging the PIT itself.
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    tampered = [dict(row) for row in rows]
    tampered[0]["pit"] = min(1.0, float(tampered[0]["pit"]) + 0.02)
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(tampered[0]))
        writer.writeheader()
        writer.writerows(tampered)

    with pytest.raises(ValueError, match="PIT mismatch"):
        audit(raw_path, summary_path)

    # Restore untampered raw rows and independently challenge the summary.
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    altered_summary = dict(summary)
    altered_summary["pit_mean"] = float(summary["pit_mean"]) + 0.05
    summary_path.write_text(json.dumps(altered_summary), encoding="utf-8")
    with pytest.raises(ValueError, match="summary field pit_mean"):
        audit(raw_path, summary_path)

    # A duplicate replicate must also be rejected even if the summary count
    # is edited to match the longer file.
    duplicate_rows = [dict(row) for row in rows] + [dict(rows[0])]
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(duplicate_rows[0]))
        writer.writeheader()
        writer.writerows(duplicate_rows)
    altered_summary["n_replicates"] = len(duplicate_rows)
    summary_path.write_text(json.dumps(altered_summary), encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate replicate id"):
        audit(raw_path, summary_path)

    # The protocol itself is part of the provenance contract.
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    invalid_protocol = dict(summary)
    invalid_protocol["generation_protocol"] = dict(summary["generation_protocol"])
    invalid_protocol["generation_protocol"]["truth_draw"] = "redraw_after_rejection"
    summary_path.write_text(json.dumps(invalid_protocol), encoding="utf-8")
    with pytest.raises(ValueError, match="generation protocol metadata"):
        audit(raw_path, summary_path)



@pytest.mark.parametrize(
    ("field", "mutate"),
    [
        ("host_index", lambda value: (int(value) + 1) % 2),
        ("observed_distance_mpc", lambda value: float(value) + 0.01),
        ("rejection_attempts", lambda value: int(value) + 1),
    ],
)
def test_independent_audit_replays_host_distance_and_attempt_count(tmp_path, field, mutate):
    rows, summary = run_calibration(n_replicates=100, seed=1357911)
    rows[0][field] = mutate(rows[0][field])
    if field == "host_index":
        # Keep the host-index/redshift pair internally consistent so this
        # mutation reaches the seeded-generation replay check.
        rows[0]["host_redshift"] = summary["host_redshifts"][int(rows[0]["host_index"])]
    raw_path = tmp_path / "pit_replicates.csv"
    summary_path = tmp_path / "pit_summary.json"
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary_path.write_text(json.dumps(summary), encoding="utf-8")

    with pytest.raises(ValueError, match="generation replay mismatch"):
        audit(raw_path, summary_path)



@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda summary: summary["h0_prior"].update({"min": float("nan")}), "prior and posterior-grid support"),
        (lambda summary: summary["posterior_grid"].update({"integration": "simpson"}), "integration metadata"),
        (lambda summary: summary.update({"host_redshifts": [0.02, float("nan")]}), "host redshifts"),
        (lambda summary: summary.update({"distance_sigma_mpc": -1.0}), "noise or detection-interval"),
        (lambda summary: summary.update({"detection_interval_mpc_open": [260.0, 50.0]}), "noise or detection-interval"),
    ],
)
def test_independent_audit_rejects_invalid_numerical_metadata(tmp_path, mutate, message):
    rows, summary = run_calibration(n_replicates=100, seed=97531)
    raw_path = tmp_path / "pit_replicates.csv"
    summary_path = tmp_path / "pit_summary.json"
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    altered_summary = json.loads(json.dumps(summary))
    mutate(altered_summary)
    summary_path.write_text(json.dumps(altered_summary), encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        audit(raw_path, summary_path)


@pytest.mark.parametrize("grid_points", [1201, 2401, 4801])
def test_calibration_grid_resolution_is_recorded_and_auditable(tmp_path, grid_points):
    rows, summary = run_calibration(n_replicates=100, seed=86420, grid_points=grid_points)
    assert summary["posterior_grid"]["points"] == grid_points
    raw_path = tmp_path / f"pit_replicates_{grid_points}.csv"
    summary_path = tmp_path / f"pit_summary_{grid_points}.json"
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary_path.write_text(json.dumps(summary), encoding="utf-8")
    result = audit(raw_path, summary_path)
    assert result["status"] == "PASS"
    assert result["n_replicates_audited"] == 100


@pytest.mark.parametrize("seed", [-1, -100, 1.5, "seed", True])
def test_calibration_runner_rejects_invalid_seed(seed):
    with pytest.raises(ValueError, match="seed must be a non-negative integer"):
        run_calibration(n_replicates=100, seed=seed)


@pytest.mark.parametrize("n_replicates", [99, 100.5, True, "100"])
def test_calibration_runner_rejects_invalid_replicate_count(n_replicates):
    with pytest.raises(ValueError, match="n_replicates must be an integer of at least 100"):
        run_calibration(n_replicates=n_replicates, seed=1)


@pytest.mark.parametrize("grid_points", [1, 1201.5, True, "1201"])
def test_calibration_runner_rejects_invalid_grid_points(grid_points):
    with pytest.raises(ValueError, match="grid_points must be an integer of at least two"):
        run_calibration(n_replicates=100, seed=1, grid_points=grid_points)
