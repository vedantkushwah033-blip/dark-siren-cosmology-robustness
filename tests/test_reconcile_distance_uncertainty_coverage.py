from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "reconcile_distance_uncertainty_coverage",
    ROOT / "scripts" / "reconcile_distance_uncertainty_coverage.py",
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def _write_raw(path: Path) -> list[dict]:
    rows = []
    for truth in MODULE.EXPECTED_TRUTHS:
        for fraction in MODULE.EXPECTED_FRACTIONS:
            for replicate in range(3):
                for arm in MODULE.ARMS:
                    if truth in (50.0, 90.0):
                        lower68, upper68 = truth + 0.01, truth + 1.0
                    else:
                        lower68, upper68 = truth - 0.5, truth + 0.5
                    lower95, upper95 = truth - 1.0, truth + 1.0
                    rows.append({
                        "true_h0": truth,
                        "fractional_error": fraction,
                        "replicate": replicate,
                        "arm": arm,
                        "lower_68": lower68,
                        "upper_68": upper68,
                        "lower_95": lower95,
                        "upper_95": upper95,
                    })
    with path.open("w", encoding="utf-8", newline="") as stream:
        import csv
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def test_reconciliation_separates_boundary_and_interior_denominators(tmp_path):
    raw = tmp_path / "raw.csv"
    rows = _write_raw(raw)
    result = MODULE.reconcile(raw, bootstrap_resamples=200)
    assert result["status"] == "PASS"
    assert result["raw_record_count"] == 9 * 4 * 3 * 3
    assert result["interior_truths"] == [55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0]
    for arm in MODULE.ARMS:
        item = result["results"][arm]
        assert item["all_nine_truth_grid"]["68"]["n"] == 9 * 4 * 3
        assert item["all_nine_truth_grid"]["68"]["coverage"] == pytest.approx(7 / 9)
        assert item["interior_truth_grid"]["68"]["n"] == 7 * 4 * 3
        assert item["interior_truth_grid"]["68"]["coverage"] == pytest.approx(1.0)
        assert item["stratified_bootstrap"]["68"]["percentile_95_interval"] == [1.0, 1.0]


def test_reconciliation_is_deterministic_and_rejects_duplicate_keys(tmp_path):
    raw = tmp_path / "raw.csv"
    rows = _write_raw(raw)
    first = MODULE.reconcile(raw, bootstrap_seed=99, bootstrap_resamples=200)
    second = MODULE.reconcile(raw, bootstrap_seed=99, bootstrap_resamples=200)
    assert first == second

    rows.append(dict(rows[0]))
    with raw.open("w", encoding="utf-8", newline="") as stream:
        import csv
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="duplicate raw key"):
        MODULE.reconcile(raw, bootstrap_resamples=200)


def test_reconciliation_rejects_missing_design_cell(tmp_path):
    raw = tmp_path / "raw.csv"
    rows = _write_raw(raw)
    rows.pop()
    with raw.open("w", encoding="utf-8", newline="") as stream:
        import csv
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="raw design incomplete"):
        MODULE.reconcile(raw, bootstrap_resamples=200)
