import importlib.util
from pathlib import Path

import pytest


PATH = Path(__file__).resolve().parents[1] / "scripts" / "analyze_known_host_oracle_contrasts.py"
SPEC = importlib.util.spec_from_file_location("analyze_known_host_oracle_contrasts", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def row(run, scenario, median, seed=2026):
    return {
        "run": str(run),
        "scenario": scenario,
        "seed": str(seed),
        "host_redshift": "0.03",
        "observed_distance_mpc": "130.0",
        "distance_sigma_mpc": "19.5",
        "H0_true": "70.0",
        "median": str(median),
        "lower_68": "60.0",
        "upper_68": "80.0",
        "lower_95": "55.0",
        "upper_95": "85.0",
    }


def test_paired_oracle_contrast_validates_shared_event_inputs():
    full = [row(0, "random70_full_population_oracle", 72.0),
            row(1, "random70_full_population_oracle", 74.0, 2027)]
    known = [row(0, "known_true_host_oracle", 71.0),
             row(1, "known_true_host_oracle", 75.0, 2027)]
    result = MODULE.paired_delta(full, known)
    assert result["median_error"] == [-1.0, 1.0]
    assert result["coverage_68"] == [1.0, 1.0]


def test_paired_oracle_contrast_rejects_unmatched_event_inputs():
    full = [row(0, "random70_full_population_oracle", 72.0)]
    known = [row(0, "known_true_host_oracle", 71.0)]
    known[0]["observed_distance_mpc"] = "131.0"
    with pytest.raises(ValueError, match="shared event input observed_distance_mpc differs"):
        MODULE.paired_delta(full, known)


def test_paired_oracle_contrast_rejects_duplicate_run_indices():
    full = [row(0, "random70_full_population_oracle", 72.0),
            row(0, "random70_full_population_oracle", 73.0)]
    known = [row(0, "known_true_host_oracle", 71.0),
             row(0, "known_true_host_oracle", 72.0)]
    with pytest.raises(ValueError, match="run indices must be unique"):
        MODULE.paired_delta(full, known)
