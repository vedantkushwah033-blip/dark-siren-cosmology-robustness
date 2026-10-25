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

    tampered = [dict(row) for row in rows]
    tampered[0]["pit"] = min(1.0, float(tampered[0]["pit"]) + 0.02)
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(tampered[0]))
        writer.writeheader()
        writer.writerows(tampered)

    with pytest.raises(ValueError, match="PIT mismatch"):
        audit(raw_path, summary_path)
