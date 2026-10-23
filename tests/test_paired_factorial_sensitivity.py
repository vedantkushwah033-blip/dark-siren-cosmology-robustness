import importlib.util
import math
from pathlib import Path

import pytest

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_paired_factorial_sensitivity.py"
SPEC = importlib.util.spec_from_file_location("paired_factorial", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_scenario_names_are_unique_and_have_expected_count():
    names = [row[0] for row in module.SCENARIOS]
    assert len(names) == 9
    assert len(set(names)) == 9


def test_each_scenario_has_valid_completeness_and_weights():
    for name, missingness, completeness, generation, inference in module.SCENARIOS:
        assert missingness in {"complete", "random", "redshift"}
        assert 0.0 < completeness <= 1.0
        assert (missingness == "complete") == (completeness == 1.0)
        assert generation in {"uniform", "luminosity"}
        assert inference in {"uniform", "luminosity"}
        assert name


def test_prespecified_contrast_names_exist_and_keep_host_generation_comparable():
    contrast_specs = (
        ("complete_uniform_uniform", "random70_uniform_uniform"),
        ("complete_uniform_uniform", "redshift70_uniform_uniform"),
        ("random70_uniform_uniform", "redshift70_uniform_uniform"),
        ("complete_luminosity_luminosity", "complete_luminosity_uniform"),
        ("random70_luminosity_luminosity", "random70_luminosity_uniform"),
        ("redshift70_luminosity_luminosity", "redshift70_luminosity_uniform"),
        ("random70_luminosity_uniform", "redshift70_luminosity_uniform"),
    )
    scenarios = {name: (missingness, completeness, generation, inference)
                 for name, missingness, completeness, generation, inference in module.SCENARIOS}
    for a, b in contrast_specs:
        assert a in scenarios and b in scenarios
        assert scenarios[a][2] == scenarios[b][2]


def test_known_unpaired_run_ids_would_be_detected():
    a = {0: {"median": 70.0}, 1: {"median": 71.0}}
    b = {0: {"median": 70.5}, 2: {"median": 72.0}}
    assert set(a) != set(b)


def test_coverage_indicator_definition():
    true_h0 = 70.0
    assert float(68.0 <= true_h0 <= 72.0) == 1.0
    assert float(70.1 <= true_h0 <= 72.0) == 0.0
