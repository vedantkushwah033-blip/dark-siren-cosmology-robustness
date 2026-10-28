import csv
from pathlib import Path

import pytest

from scripts.audit_full_run_coverage_heterogeneity import audit, wilson

ARMS = [
    "A_fixed_realized_sigma",
    "B_fractional_sigma_H0",
    "C_fixed_reference_sigma",
]


def write_csv(path: Path, rows):
    columns = [
        "true_h0", "fractional_error", "replicate", "arm",
        "lower_68", "upper_68", "lower_95", "upper_95",
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def test_wilson_interval_contains_estimate():
    low, high = wilson(7, 10)
    assert low < 0.7 < high


def test_audit_separates_boundary_and_interior_truths(tmp_path):
    rows = []
    for arm in ARMS:
        for truth in (50, 55, 90):
            for replicate in range(2):
                # Exact boundary truths fall outside the central interval; 55 is covered.
                lower = 50.1 if truth in (50, 90) else 54.0
                upper = 89.9 if truth in (50, 90) else 56.0
                rows.append({
                    "true_h0": truth, "fractional_error": 0.1,
                    "replicate": replicate, "arm": arm,
                    "lower_68": lower, "upper_68": upper,
                    "lower_95": lower, "upper_95": upper,
                })
    path = tmp_path / "replicates.csv"
    write_csv(path, rows)
    report = audit(path)
    assert report["raw_record_count"] == 18
    for arm in ARMS:
        assert report["overall_by_arm"][arm]["coverage_68"]["estimate"] == pytest.approx(1 / 3)
        assert report["interior_by_arm"][arm]["coverage_68"]["estimate"] == 1.0
        assert report["boundary_by_arm"][arm]["coverage_68"]["estimate"] == 0.0


def test_audit_rejects_duplicate_records(tmp_path):
    rows = []
    for arm in ARMS:
        rows.append({
            "true_h0": 55, "fractional_error": 0.1, "replicate": 0, "arm": arm,
            "lower_68": 54, "upper_68": 56, "lower_95": 53, "upper_95": 57,
        })
    rows.append(dict(rows[0]))
    path = tmp_path / "replicates.csv"
    write_csv(path, rows)
    with pytest.raises(ValueError, match="duplicate"):
        audit(path)
