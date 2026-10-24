"""Smoke tests for the independent factorial-artifact validator."""
import importlib.util
from pathlib import Path

import pytest

PATH = Path(__file__).resolve().parents[1] / "scripts" / "validate_prior_threshold_factorial.py"
SPEC = importlib.util.spec_from_file_location("validate_factorial_artifacts", PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_validator_imports_and_declares_expected_artifact_contract():
    assert module.RAW_KEY[0] == "true_H0"
    assert "posterior_median" in module.SUMMARY_KEYS
    assert "metric" in module.CONTRAST_KEY


def test_validator_rejects_empty_artifact(tmp_path):
    empty = tmp_path / "empty.csv"
    empty.write_text("true_H0,replicate\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Empty CSV"):
        module.read_csv(empty)


def test_validator_rejects_non_finite_numeric_artifact_value(tmp_path):
    bad = tmp_path / "bad.csv"
    bad.write_text("true_H0,rmse\n70,nan\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Non-finite numeric field"):
        module.read_csv(bad)


def test_independent_bootstrap_recomputation_matches_export_algorithm_on_fixture():
    experiment = module.factorial
    rows = []
    for truth in experiment.TRUTHS:
        for rep in range(2):
            for pmin, pmax in experiment.PRIORS:
                for threshold in experiment.ASSUMED_MAXIMA:
                    for method in experiment.METHODS:
                        shift = (pmin - 50.0) * 0.01 + (threshold - 120.0) * 0.001
                        if method == "corrected":
                            shift += 0.15
                        median = truth + (-1.0 if rep == 0 else 1.0) + shift
                        width = 4.0 + 0.01 * (pmax - pmin) + 0.001 * (threshold - 120.0)
                        rows.append({
                            "true_H0": truth, "replicate": rep, "prior_min": pmin,
                            "prior_max": pmax, "assumed_detection_max_mpc": threshold,
                            "method": method, "posterior_median": median,
                            "lower_68": median - width / 2, "upper_68": median + width / 2,
                            "lower_95": median - width, "upper_95": median + width,
                        })
    expected = experiment.paired_contrasts(rows, bootstrap_draws=101, seed=20261029)
    actual = module.independently_recompute_paired_bootstrap(
        rows, reps=2, bootstrap_draws=101, seed=20261029
    )
    assert len(actual) == len(expected)
    for left, right in zip(actual, expected):
        assert (left["true_H0"], left["method"], left["metric"],
                left["comparison_prior"], left["comparison_threshold_mpc"]) == (
                right["true_H0"], right["method"], right["metric"],
                right["comparison_prior"], right["comparison_threshold_mpc"])
        assert left["mean_difference_vs_baseline"] == pytest.approx(
            right["mean_difference_vs_baseline"], abs=1e-12
        )
        assert left["bootstrap_95_low"] == pytest.approx(right["bootstrap_95_low"], abs=1e-12)
        assert left["bootstrap_95_high"] == pytest.approx(right["bootstrap_95_high"], abs=1e-12)


def test_validator_rejects_observation_outside_declared_detection_interval():
    row = {"host_index": 0, "observed_distance_mpc": 120.0}
    with pytest.raises(ValueError, match="outside true detection interval"):
        module.validate_observation_fields(row, "fixture")


def test_validator_rejects_unknown_host_index():
    row = {"host_index": 9, "observed_distance_mpc": 42.0}
    with pytest.raises(ValueError, match="Invalid host_index"):
        module.validate_observation_fields(row, "fixture")
