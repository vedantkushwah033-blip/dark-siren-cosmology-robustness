from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "analyze_known_host_paired_comparisons.py"
SPEC = importlib.util.spec_from_file_location("paired_comparisons", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def make_rows(replicates=100):
    rows = []
    for truth in MODULE.TRUTHS:
        for run in range(replicates):
            for noise in MODULE.NOISES:
                # Same host quantities across noise levels; wider intervals at higher noise.
                half_width = noise * 100
                median = truth + (noise - 0.05)
                rows.append({
                    "run": str(run), "seed": str(20261025 + run),
                    "true_h0": str(truth), "fractional_distance_error": str(noise),
                    "host_redshift": "0.03", "true_distance_mpc": str(299792.458 * 0.03 / truth),
                    "median": str(median),
                    "lower_68": str(median - half_width), "upper_68": str(median + half_width),
                    "lower_95": str(median - 2 * half_width), "upper_95": str(median + 2 * half_width),
                })
    return rows


def test_analyze_emits_all_prespecified_paired_contrasts():
    results = MODULE.analyze(make_rows(), replicates=100, n_boot=1000, seed=20261025)
    assert len(results) == 9 * 3 * len(MODULE.METRICS)
    row = next(r for r in results if r["true_h0"] == 70.0
               and r["comparison_fractional_error"] == 0.10 and r["metric"] == "width_68")
    assert row["mean_difference"] == pytest.approx(10.0)
    assert row["n_paired_replicates"] == 100
    assert row["bootstrap_mean_difference_ci95_low"] == pytest.approx(10.0)
    assert row["bootstrap_mean_difference_ci95_high"] == pytest.approx(10.0)


def test_analyze_rejects_missing_pair_member():
    rows = make_rows()
    rows.pop()
    with pytest.raises(ValueError, match="complete paired"):
        MODULE.analyze(rows, replicates=100, n_boot=1000, seed=20261025)


def test_analyze_rejects_duplicate_pair_member():
    rows = make_rows()
    rows.append(dict(rows[0]))
    with pytest.raises(ValueError, match="duplicate"):
        MODULE.analyze(rows, replicates=100, n_boot=1000, seed=20261025)


def test_analyze_rejects_mismatched_host_pairing():
    rows = make_rows()
    target = next(r for r in rows if float(r["true_h0"]) == 70.0 and int(r["run"]) == 0
                  and float(r["fractional_distance_error"]) == 0.10)
    target["host_redshift"] = "0.031"
    with pytest.raises(ValueError, match="pairing mismatch"):
        MODULE.analyze(rows, replicates=100, n_boot=1000, seed=20261025)


def test_audit_gate_requires_pass_and_hashes(tmp_path):
    with pytest.raises(ValueError, match="audit=PASS"):
        MODULE.verify_audit(tmp_path, {"audit": "FAIL"})
