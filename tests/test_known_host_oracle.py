import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_known_host_oracle.py"
SPEC = importlib.util.spec_from_file_location("run_known_host_oracle", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_known_host_oracle_median_bias_uses_median_of_errors():
    rows = [{"median": 69.0}, {"median": 71.0}, {"median": 80.0}]
    assert MODULE.median_bias_from_rows(rows, true_h0=70.0) == 1.0


def test_known_host_oracle_median_bias_rejects_empty_rows():
    with pytest.raises(ValueError, match="must not be empty"):
        MODULE.median_bias_from_rows([], true_h0=70.0)
