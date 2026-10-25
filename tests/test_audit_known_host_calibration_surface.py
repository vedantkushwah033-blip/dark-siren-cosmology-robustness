import importlib.util
import math
from pathlib import Path

import pytest


PATH = Path(__file__).resolve().parents[1] / "scripts" / "audit_known_host_calibration_surface.py"
SPEC = importlib.util.spec_from_file_location("audit_known_host_calibration_surface", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


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
    assert lo0 == pytest.approx(0.0)
    assert 0 < hi0 < 0.05
    assert 0.95 < lo1 < 1
    assert hi1 == pytest.approx(1.0)


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
