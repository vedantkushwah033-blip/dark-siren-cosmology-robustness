from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "analyze_known_host_calibration_pairs.py"
SPEC = importlib.util.spec_from_file_location("analyze_known_host_calibration_pairs", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC and SPEC.loader
SPEC.loader.exec_module(MODULE)


def test_paired_analysis_uses_within_replicate_differences():
    rows = []
    for run in range(100):
        for noise in MODULE.NOISES:
            shift = {0.05: 0.0, 0.10: 1.0, 0.15: 2.0, 0.20: 3.0}[noise]
            rows.append({
                "run": str(run), "seed": str(500 + run), "true_h0": "70.0",
                "fractional_distance_error": str(noise), "host_redshift": "0.03",
                "true_distance_mpc": "128.0", "observed_distance_mpc": str(128 + noise * (run % 5)),
                "distance_sigma_mpc": str(noise * 128.0),
                "median": str(70.0 + shift), "lower_68": str(68.0 + shift),
                "upper_68": str(72.0 + shift), "lower_95": str(66.0 + shift),
                "upper_95": str(74.0 + shift),
            })
    # Exercise analysis directly with a complete synthetic stratum.
    output = MODULE.analyze(rows, 100, 200, 42)
    subset = [r for r in output if r["true_h0"] == 70.0 and r["metric"] == "absolute_median_error"]
    by_noise = {r["comparison_fractional_error"]: r for r in subset}
    assert by_noise[0.10]["mean_paired_difference"] == pytest.approx(1.0)
    assert by_noise[0.15]["mean_paired_difference"] == pytest.approx(2.0)
    assert by_noise[0.20]["mean_paired_difference"] == pytest.approx(3.0)
    assert len(output) == 8 * 3 * 5


def test_percentile_interpolates_and_handles_endpoints():
    assert MODULE.percentile([0.0, 10.0], 0.25) == pytest.approx(2.5)
    assert MODULE.percentile([1.0, 2.0, 3.0], 0.0) == 1.0
    assert MODULE.percentile([1.0, 2.0, 3.0], 1.0) == 3.0


def test_pair_analysis_rejects_incomplete_blocks(tmp_path):
    # Minimal files to ensure audit gating is checked before attempting CSV parsing.
    (tmp_path / "audit_report.json").write_text(json.dumps({"audit": "FAIL"}))
    (tmp_path / "cell_summary.json").write_text("{}")
    (tmp_path / "per_run_results.csv").write_text("run\n0\n")
    with pytest.raises(ValueError, match="audit PASS"):
        MODULE.load_and_validate(tmp_path, None)
