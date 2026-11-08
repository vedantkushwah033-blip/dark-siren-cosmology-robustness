from __future__ import annotations
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd
import pytest
from dark_siren_robustness.generative_selection import catalogue_inclusion_probability, generate_synthetic_population, latent_galaxy_intensity

def test_latent_intensity_is_positive_and_nonconstant() -> None:
    z = np.linspace(0.01, 0.40, 1000)
    intensity = latent_galaxy_intensity(z)
    assert np.all(np.isfinite(intensity)) and np.all(intensity > 0)
    assert float(np.ptp(intensity)) > 1.0

def test_catalogue_selection_decreases_with_redshift() -> None:
    p = catalogue_inclusion_probability(np.linspace(0.01, 0.40, 500))
    assert np.all((p > 0) & (p <= 1))
    assert np.all(np.diff(p) < 0)

def test_generator_is_reproducible_and_separates_truth() -> None:
    kwargs = {"seed": 9127, "n_parent": 1800, "n_detected": 35}
    first, second = generate_synthetic_population(**kwargs), generate_synthetic_population(**kwargs)
    for key in ("parent_truth", "catalogue_observed", "gw_events_observed", "gw_event_truth", "gw_generation_trace_truth"):
        pd.testing.assert_frame_equal(first[key], second[key])
    assert list(first["catalogue_observed"].columns) == ["object_id", "redshift"]
    assert list(first["gw_events_observed"].columns) == ["event_id", "distance_mpc", "detected"]
    assert "host_redshift" not in first["gw_events_observed"].columns
    assert "host_object_id" not in first["gw_events_observed"].columns
    assert len(first["gw_events_observed"]) == 35

def test_detected_events_obey_threshold_and_truth_ids_match_parent() -> None:
    result = generate_synthetic_population(seed=44, n_parent=1500, n_detected=40, detection_threshold_mpc=250.0)
    assert np.all(result["gw_events_observed"]["distance_mpc"] >= 250.0)
    truth, parent = result["gw_event_truth"], result["parent_truth"]
    assert np.all(truth["host_object_id"].isin(parent["object_id"]))
    expected = 299_792.458 * truth["host_redshift"] / result["config"]["h0_true_km_s_mpc"]
    np.testing.assert_allclose(truth["latent_mean_distance_mpc"], expected)
    assert len(result["catalogue_observed"]) < len(parent)
    config = result["config"]
    assert config["n_generation_attempts"] == (
        config["n_detected_candidates"] + config["n_rejected_candidates"]
    )
    assert config["n_detected_candidates"] == (
        config["n_detected_events"] + config["n_excess_detected_candidates_discarded"]
    )
    assert config["n_rejected_candidates"] > 0
    trace = result["gw_generation_trace_truth"]
    assert len(trace) == config["n_generation_attempts"]
    assert int(trace["detected"].sum()) == config["n_detected_candidates"]
    assert int((trace["detected"] == 0).sum()) == config["n_rejected_candidates"]
    assert int(trace["selected_as_event"].sum()) == config["n_detected_events"]
    assert int((trace["detected"] == 1).sum() - trace["selected_as_event"].sum()) == config["n_excess_detected_candidates_discarded"]

def test_invalid_configuration_is_rejected() -> None:
    with pytest.raises(ValueError):
        generate_synthetic_population(seed=-1)
    with pytest.raises(ValueError):
        generate_synthetic_population(n_parent=20)
    with pytest.raises(ValueError):
        catalogue_inclusion_probability(np.array([0.1, 0.2]), width=0.0)

def test_runner_writes_separate_observable_and_truth_artifacts(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    output_dir = tmp_path / "selection-output"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src") + os.pathsep + env.get("PYTHONPATH", "")
    completed = subprocess.run(
        [sys.executable, str(root / "scripts" / "run_synthetic_selection_experiment.py"),
         "--seed", "120", "--n-parent", "1200", "--n-detected", "25",
         "--output-dir", str(output_dir)],
        cwd=root, env=env, text=True, capture_output=True, check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    config = json.loads((output_dir / "config.json").read_text())
    assert config["n_parent"] == 1200
    assert config["n_detected_events"] == 25
    assert config["observable_files"] == ["observable/catalogue_observed.csv", "observable/gw_events_observed.csv"]
    assert config["truth_only_files"] == ["truth/parent_truth.csv", "truth/gw_event_truth.csv", "truth/gw_generation_trace_truth.csv"]
    observed_events = pd.read_csv(output_dir / "observable" / "gw_events_observed.csv")
    truth_events = pd.read_csv(output_dir / "truth" / "gw_event_truth.csv")
    trace = pd.read_csv(output_dir / "truth" / "gw_generation_trace_truth.csv")
    assert list(observed_events.columns) == ["event_id", "distance_mpc", "detected"]
    assert "host_redshift" in truth_events.columns
    assert len(observed_events) == len(truth_events) == 25
    assert int(trace["selected_as_event"].sum()) == 25
    assert not (output_dir / "observable" / "gw_event_truth.csv").exists()
    assert not (output_dir / "truth" / "gw_events_observed.csv").exists()


def test_zero_floor_is_exactly_logistic_and_positive_floor_is_misspecified():
    z = np.linspace(0.01, 0.40, 101)
    midpoint, width = 0.22, 0.05
    exact_logistic = catalogue_inclusion_probability(
        z, floor=0.0, midpoint=midpoint, width=width
    )
    expected = 1.0 / (1.0 + np.exp((z - midpoint) / width))
    np.testing.assert_allclose(exact_logistic, expected, rtol=1e-14, atol=1e-14)

    floor_misspecified = catalogue_inclusion_probability(
        z, floor=0.10, midpoint=midpoint, width=width
    )
    assert np.all(floor_misspecified > exact_logistic)
    assert np.all(floor_misspecified <= 1.0)
    assert not np.allclose(
        np.log(floor_misspecified / (1.0 - floor_misspecified)),
        np.log(exact_logistic / (1.0 - exact_logistic)),
    )
