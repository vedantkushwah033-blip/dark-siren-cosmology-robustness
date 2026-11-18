import numpy as np

from dark_siren_robustness.generative_selection import (
    calibrate_quadratic_logit_intercept,
    generate_quadratic_logit_inclusion_samples,
    quadratic_logit_inclusion_probability,
)
from scripts.run_quadratic_curvature_h0_smoke import simulate_detected_events


def test_quadratic_inclusion_calibration_preserves_parent_average():
    values = [
        calibrate_quadratic_logit_intercept(gamma=g)
        for g in (-2.0, -1.0, 0.0, 1.0, 2.0)
    ]
    targets = np.asarray([x["target_rate"] for x in values])
    achieved = np.asarray([x["achieved_rate"] for x in values])
    assert np.ptp(targets) < 1e-14
    assert np.max(np.abs(achieved - targets)) < 1e-11


def test_quadratic_generator_separates_observable_and_truth_columns():
    data = generate_quadratic_logit_inclusion_samples(
        population_seed=71, validation_seed=72, n_parent=500,
        n_validation_frame=300, n_validation=120, gamma=0.8,
    )
    assert set(data["validation_observed"].columns) == {
        "validation_id", "redshift", "retained"
    }
    assert "inclusion_probability" not in data["validation_observed"]
    assert "inclusion_probability" in data["validation_frame_truth"]
    assert set(data["catalogue_observed"].columns) == {"object_id", "redshift"}
    assert data["config"]["population_seed"] != data["config"]["validation_seed"]


def test_quadratic_inclusion_probability_matches_declared_logit():
    calibration = calibrate_quadratic_logit_intercept(gamma=1.0)
    z = np.asarray([0.03, 0.10, 0.20, 0.30, 0.39])
    p = quadratic_logit_inclusion_probability(
        z, gamma=1.0, beta=-3.9, alpha=calibration["alpha"]
    )
    x = (z - calibration["coordinate_center"]) / calibration["coordinate_scale"]
    expected = 1.0 / (1.0 + np.exp(-(
        calibration["alpha"] + calibration["beta"] * x
        + calibration["gamma"] * (x**2 - calibration["x2_center"])
    )))
    assert np.allclose(p, expected, rtol=1e-13, atol=1e-14)


def test_detected_event_generator_uses_parent_population_and_threshold():
    data = generate_quadratic_logit_inclusion_samples(
        population_seed=81, validation_seed=82, n_parent=500,
        n_validation_frame=300, n_validation=120, gamma=0.5,
    )
    observed, truth, attempts = simulate_detected_events(
        data["parent_truth"], seed=83, n_events=3,
    )
    assert len(observed) == len(truth) == 3
    assert np.all(observed["distance_mpc"].to_numpy() >= 250.0)
    assert np.all(truth["host_object_id"].to_numpy() < 500)
    assert attempts >= 3
