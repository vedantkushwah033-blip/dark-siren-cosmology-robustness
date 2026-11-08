from __future__ import annotations

import numpy as np
import pandas as pd

from dark_siren_robustness.generative_selection import (
    generate_independent_validation_sample,
    generate_synthetic_population,
)
from dark_siren_robustness.inclusion_posterior import fit_inclusion_posterior_quadrature
from dark_siren_robustness.h0_posterior import normalize_h0_posterior
from dark_siren_robustness.shared_inclusion_likelihood import (
    shared_inclusion_selected_log_likelihood,
)


def test_end_to_end_independent_validation_to_selected_h0_posterior():
    """Exercise fixed analysis catalogue plus a genuinely independent validation frame."""
    generated = generate_synthetic_population(
        seed=20261108, n_parent=2400, n_detected=8, h0_true=70.0,
        sigma_mpc=140.0, detection_threshold_mpc=250.0,
    )
    validation_generated = generate_independent_validation_sample(
        seed=20261109, n_validation_frame=4000, n_validation=300,
    )
    catalogue = generated["catalogue_observed"]
    events_frame = generated["gw_events_observed"]
    config = generated["config"]
    validation = validation_generated["validation_observed"]

    # Distinct generation seeds and independent parent-frame provenance.
    assert generated["config"]["seed"] != validation_generated["config"]["seed"]
    assert validation_generated["config"]["independent_frame"] is True
    assert list(validation.columns) == ["validation_id", "redshift", "retained"]
    assert len(validation) == 300
    assert validation["retained"].nunique() == 2
    assert "inclusion_probability" not in validation.columns
    assert "inclusion_probability" not in catalogue.columns
    assert "host_object_id" not in events_frame.columns
    assert "host_redshift" not in events_frame.columns
    assert np.all(events_frame["distance_mpc"].to_numpy() >= config["detection_threshold_mpc"])

    inclusion_posterior = fit_inclusion_posterior_quadrature(
        validation["redshift"].to_numpy(), validation["retained"].to_numpy(),
        degree=1, order=7,
    )
    assert np.isclose(inclusion_posterior.weights.sum(), 1.0)

    z_catalogue = catalogue["redshift"].to_numpy()
    z_min, z_max = config["redshift_domain"]
    merger_base_weight = 1.0 + 1.5 * (z_catalogue - z_min) / (z_max - z_min)
    events = [
        {"gw_distance": float(row.distance_mpc), "gw_sigma": float(config["sigma_mpc"])}
        for row in events_frame.itertuples(index=False)
    ]
    h0_grid = np.linspace(50.0, 90.0, 81)
    selected = shared_inclusion_selected_log_likelihood(
        h0_grid, z_catalogue, events, inclusion_posterior,
        detection_threshold_mpc=float(config["detection_threshold_mpc"]),
        population_base_weights=merger_base_weight,
    )
    posterior = normalize_h0_posterior(
        h0_grid, selected.marginal_log_likelihood, prior_min=50.0, prior_max=90.0,
    )
    assert selected.node_log_selected_likelihood.shape == (
        inclusion_posterior.nodes.shape[0], h0_grid.size
    )
    assert np.all(np.isfinite(selected.marginal_log_likelihood))
    assert np.all(np.isfinite(selected.node_log_detection_efficiency))
    assert np.allclose(selected.node_catalogue_weights.sum(axis=0), 1.0)
    assert np.isclose(np.trapezoid(posterior.density, posterior.grid), 1.0)
    assert 50.0 <= posterior.quantile(0.16) < posterior.quantile(0.84) <= 90.0


def test_independent_validation_frame_is_reproducible_and_hides_truth():
    a = generate_independent_validation_sample(
        seed=9921, n_validation_frame=1200, n_validation=150,
    )
    b = generate_independent_validation_sample(
        seed=9921, n_validation_frame=1200, n_validation=150,
    )
    pd.testing.assert_frame_equal(a["validation_observed"], b["validation_observed"])
    assert "inclusion_probability" not in a["validation_observed"].columns
    assert "inclusion_probability" in a["validation_frame_truth"].columns
    assert len(a["validation_frame_truth"]) == 1200
    assert a["validation_frame_truth"]["selected_for_validation"].sum() == 150


def test_independent_validation_h0_runner_writes_auditable_artifacts(tmp_path):
    import hashlib
    import json
    import os
    import subprocess
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    output = tmp_path / "independent-validation-smoke"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src") + os.pathsep + env.get("PYTHONPATH", "")
    completed = subprocess.run(
        [
            sys.executable,
            str(root / "scripts" / "run_independent_validation_inclusion_h0_smoke.py"),
            "--seed", "8811",
            "--validation-seed", "9922",
            "--n-parent", "1200",
            "--n-detected", "4",
            "--n-validation-frame", "1000",
            "--n-validation", "150",
            "--orders", "5,7,9,11",
            "--h0-points", "41",
            "--catalogue-floor", "0.0",
            "--output-dir", str(output),
        ],
        cwd=root, env=env, text=True, capture_output=True, check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
    config = json.loads((output / "config.json").read_text())
    convergence = json.loads((output / "quadrature_convergence.json").read_text())
    prior_sensitivity = json.loads(
        (output / "prior_scale_sensitivity.json").read_text()
    )
    prior_density = pd.read_csv(output / "prior_scale_sensitivity_density.csv")
    assert prior_sensitivity["slope_multipliers"] == [0.5, 1.0, 2.0]
    assert len(prior_sensitivity["records"]) == 3
    assert list(prior_density.columns) == [
        "H0",
        "density_slope_multiplier_0.5",
        "density_slope_multiplier_1",
        "density_slope_multiplier_2",
    ]
    assert config["analysis_generator_config"]["seed"] != config["validation_generator_config"]["seed"]
    assert config["truth_firewall"].startswith("Only validation")
    assert config["analysis_generator_config"]["catalogue_selection"]["floor"] == 0.0
    assert config["validation_generator_config"]["catalogue_selection"]["floor"] == 0.0
    assert convergence["status"] in {"CONVERGED", "NON_CONVERGED"}
    assert len(convergence["posterior_by_order"]) >= 2

    validation = pd.read_csv(output / "observable" / "validation_observed.csv")
    catalogue = pd.read_csv(output / "observable" / "analysis_catalogue.csv")
    events = pd.read_csv(output / "observable" / "gw_events_observed.csv")
    assert list(validation.columns) == ["validation_id", "redshift", "retained"]
    assert list(catalogue.columns) == ["object_id", "redshift"]
    assert list(events.columns) == ["event_id", "distance_mpc", "detected"]
    assert "inclusion_probability" not in validation.columns
    assert "host_object_id" not in events.columns
    density = pd.read_csv(output / "h0_posterior_density_by_order.csv")
    assert "H0" in density.columns and len(density) == 41

    manifest = json.loads((output / "manifest.json").read_text())
    for relative_path, expected_digest in manifest.items():
        digest = hashlib.sha256((output / relative_path).read_bytes()).hexdigest()
        assert digest == expected_digest

    audit = subprocess.run(
        [
            sys.executable,
            str(root / "scripts" / "audit_independent_validation_inclusion_h0_smoke.py"),
            str(output),
        ],
        cwd=root, env=env, text=True, capture_output=True, check=False,
    )
    assert audit.returncode == 0, audit.stdout + audit.stderr
    audit_report = json.loads(audit.stdout)
    assert audit_report["status"] == "PASS"
    assert audit_report["generation_attempts_reconstructed"] > 0
    assert audit_report["selected_events_reconstructed"] == 4
    assert len(audit_report["independent_selected_likelihood_reconstruction_max_abs_density_differences"]) >= 2
    assert max(audit_report["independent_selected_likelihood_reconstruction_max_abs_density_differences"].values()) < 2e-8
    assert len(audit_report["prior_scale_sensitivity_reconstruction_max_abs_density_differences"]) == 3
    assert max(audit_report["prior_scale_sensitivity_reconstruction_max_abs_density_differences"].values()) < 2e-8
