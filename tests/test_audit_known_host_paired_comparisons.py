from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load_script(name: str, module_name: str):
    path = ROOT / "scripts" / name
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


ANALYSIS = load_script("analyze_known_host_paired_comparisons.py", "paired_analysis_for_audit_test")
AUDIT = load_script("audit_known_host_paired_comparisons.py", "paired_audit_under_test")


def make_rows(replicates=100):
    rows = []
    for truth_index, truth in enumerate(ANALYSIS.TRUTHS):
        for run in range(replicates):
            for noise in ANALYSIS.NOISES:
                median = truth + 0.2 * (noise / 0.05)
                half_width = noise * 100
                rows.append({
                    "run": str(run),
                    "seed": str(20261025 + truth_index * replicates + run),
                    "true_h0": str(truth),
                    "fractional_distance_error": str(noise),
                    "host_redshift": "0.03",
                    "true_distance_mpc": str(299792.458 * 0.03 / truth),
                    "median": str(median),
                    "lower_68": str(median - half_width),
                    "upper_68": str(median + half_width),
                    "lower_95": str(median - 2 * half_width),
                    "upper_95": str(median + 2 * half_width),
                })
    return rows


def table_as_csv_rows(rows):
    return [{key: str(value) for key, value in row.items()} for row in rows]


def test_independent_auditor_reconstructs_all_paired_contrasts():
    raw = make_rows()
    produced = ANALYSIS.analyze(raw, replicates=100, n_boot=1000, seed=20261025)
    result = AUDIT.validate_pair_table(
        raw, table_as_csv_rows(produced), replicates=100,
        bootstrap_resamples=1000, base_seed=20261025,
    )
    assert result == {
        "raw_records": 9 * 4 * 100,
        "paired_comparisons": 9 * 3 * 5,
        "paired_groups": 9 * 100,
    }


def test_independent_auditor_rejects_changed_mean_difference():
    raw = make_rows()
    produced = table_as_csv_rows(
        ANALYSIS.analyze(raw, replicates=100, n_boot=1000, seed=20261025)
    )
    produced[0]["mean_difference"] = str(float(produced[0]["mean_difference"]) + 0.1)
    with pytest.raises(ValueError, match="mean_difference mismatch"):
        AUDIT.validate_pair_table(
            raw, produced, replicates=100,
            bootstrap_resamples=1000, base_seed=20261025,
        )


def test_independent_auditor_rejects_missing_comparison():
    raw = make_rows()
    produced = table_as_csv_rows(
        ANALYSIS.analyze(raw, replicates=100, n_boot=1000, seed=20261025)
    )
    produced.pop()
    with pytest.raises(ValueError, match="expected 135 paired rows"):
        AUDIT.validate_pair_table(
            raw, produced, replicates=100,
            bootstrap_resamples=1000, base_seed=20261025,
        )
