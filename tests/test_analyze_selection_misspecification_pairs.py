from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_selection_misspecification_pairs import analyze_records, ARMS


def _rows(n=3):
    true = 70.0
    error = {
        "analytic_truth": [1.0, -1.0, 2.0],
        "oracle": [0.8, -1.2, 2.2],
        "catalogue_only": [-10.0, -8.0, -9.0],
        "correct_weighting": [-1.0, 2.0, 1.0],
        "shifted_transition": [-2.0, 1.0, 0.0],
        "shallow_slope": [-4.0, 3.0, 0.0],
    }
    rows = []
    for rep in range(n):
        for arm in ARMS:
            rows.append({
                "replicate": str(rep), "seed": str(101 + rep), "arm": arm,
                "h0_true": str(true), "posterior_mean": str(true + error[arm][rep]),
                "ci68_width": "7.0", "covered_68": "true",
                "density_ise": "0.1", "n_parent": "1000",
                "n_catalogue_observed": "700", "n_detected": "8",
            })
    return rows


def test_primary_estimand_direction_and_analysis_determinism():
    first = analyze_records(_rows(), n_replicates=3, seed_start=101, h0_true=70.0, n_boot=1000)
    second = analyze_records(_rows(), n_replicates=3, seed_start=101, h0_true=70.0, n_boot=1000)
    primary = first["primary_estimand"]
    assert primary["mean"] == pytest.approx((10 - 1 + 8 - 2 + 9 - 1) / 3)
    assert primary["positive_favors"] == "correct_weighting"
    assert primary["bootstrap_95_interval"] == second["primary_estimand"]["bootstrap_95_interval"]
    assert len(first["secondary_paired_absolute_error_contrasts"]) == 7


def test_missing_arm_and_inconsistent_pairing_are_rejected():
    rows = _rows()
    with pytest.raises(ValueError, match="raw record count"):
        analyze_records(rows[:-1], n_replicates=3, seed_start=101, h0_true=70.0, n_boot=1000)
    rows = _rows()
    rows[-1]["seed"] = "999"
    with pytest.raises(ValueError, match="seed mismatch"):
        analyze_records(rows, n_replicates=3, seed_start=101, h0_true=70.0, n_boot=1000)
