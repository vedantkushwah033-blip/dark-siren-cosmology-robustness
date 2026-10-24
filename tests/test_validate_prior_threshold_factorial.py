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
