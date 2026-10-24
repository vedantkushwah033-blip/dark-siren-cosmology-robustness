from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_posterior_grid_stability.py"
SPEC = spec_from_file_location("posterior_grid_stability", MODULE_PATH)
module = module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_reported_posterior_area_is_normalized_on_prior_support():
    cases = [
        module.evaluate_case("baseline", np.linspace(50.0, 90.0, 1601), 50.0, 90.0),
        module.evaluate_case("wider_grid_same_prior", np.linspace(30.0, 110.0, 6401), 50.0, 90.0),
        module.evaluate_case("wider_grid_wider_prior", np.linspace(30.0, 110.0, 6401), 30.0, 110.0),
    ]
    for case in cases:
        assert np.isclose(case["posterior_area"], 1.0, rtol=0.0, atol=1e-10), case


def test_wider_numerical_grid_does_not_change_reported_area_for_fixed_prior():
    baseline = module.evaluate_case("baseline", np.linspace(50.0, 90.0, 1601), 50.0, 90.0)
    wider = module.evaluate_case("wider_grid_same_prior", np.linspace(30.0, 110.0, 6401), 50.0, 90.0)
    assert abs(baseline["posterior_area"] - wider["posterior_area"]) < 1e-10
