import csv
import json
from pathlib import Path
import platform
import sys
from importlib.metadata import version

import numpy as np

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_distance_uncertainty_sensitivity import audit_artifacts
from dark_siren_robustness.uncertainty_sensitivity import ARMS, paired_ab_contrasts, run_experiment


def _write_artifacts(folder, replicates=2, seed=123, grid_points=101):
    raw, summaries = run_experiment(replicates, seed, grid_points)
    folder.mkdir(parents=True)
    raw_path = folder / "replicates.csv"
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(raw[0]))
        writer.writeheader()
        writer.writerows(raw)
    bootstrap_resamples = 100
    metadata = {
        "status": "synthetic sensitivity diagnostic; not an observational GW analysis",
        "code_revision": "test-fixture",
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": version("scipy"),
        "failed_run_count": 0,
        "likelihood_equations": {
            "A_fixed_realized_sigma": "logL=-0.5*((d_obs-c*z/H0)/sigma_realized)^2-log(sigma_realized)-0.5*log(2*pi)",
            "B_fractional_sigma_H0": "sigma(H0)=f*c*z/H0; logL=-0.5*((d_obs-c*z/H0)/sigma(H0))^2-log(sigma(H0))-0.5*log(2*pi)",
            "C_fixed_reference_sigma": "sigma_ref=f*c*z/70; logL=-0.5*((d_obs-c*z/H0)/sigma_ref)^2-log(sigma_ref)-0.5*log(2*pi)",
        },
        "seed": seed,
        "replicates_per_cell": replicates,
        "grid_points": grid_points,
        "bootstrap_resamples": bootstrap_resamples,
        "prior": [50.0, 90.0],
        "truths": [50, 55, 60, 65, 70, 75, 80, 85, 90],
        "fractional_errors": [0.05, 0.10, 0.15, 0.20],
        "arms": list(ARMS),
        "raw_record_count": len(raw),
        "summary_cell_count": len(summaries),
        "boundary_truths_must_be_interpreted_separately": True,
        "summary": summaries,
        "paired_ab_contrasts": paired_ab_contrasts(raw, seed=seed, bootstrap_resamples=bootstrap_resamples),
        "paired_ab_contrast_count": 9 * 4 * 5,
    }
    (folder / "summary.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def test_independent_audit_reconstructs_raw_rows_and_summaries(tmp_path):
    folder = tmp_path / "artifacts"
    _write_artifacts(folder)
    report = audit_artifacts(folder)
    assert report["status"] == "PASS"
    assert report["raw_record_count"] == 9 * 4 * 2 * 3
    assert report["summary_cell_count"] == 9 * 4 * 3
    assert report["a_b_pairs_verified"] == 9 * 4 * 2
    assert report["posterior_rows_reconstructed"] == report["raw_record_count"]
    assert report["paired_ab_contrast_count"] == 9 * 4 * 5
    assert report["paired_ab_bootstrap_intervals_reconstructed"] == 9 * 4 * 5


def test_independent_audit_rejects_tampered_observation(tmp_path):
    folder = tmp_path / "artifacts"
    _write_artifacts(folder)
    raw_path = folder / "replicates.csv"
    with raw_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        columns = list(rows[0])
    rows[0]["observed_distance_mpc"] = str(float(rows[0]["observed_distance_mpc"]) + 1.0)
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="observed distance"):
        audit_artifacts(folder)


def test_independent_audit_rejects_tampered_summary(tmp_path):
    folder = tmp_path / "artifacts"
    _write_artifacts(folder)
    path = folder / "summary.json"
    metadata = json.loads(path.read_text(encoding="utf-8"))
    metadata["summary"][0]["coverage_68"] = 1.0
    path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="coverage_68 mismatch"):
        audit_artifacts(folder)


def test_independent_audit_rejects_tampered_paired_contrast(tmp_path):
    folder = tmp_path / "artifacts"
    _write_artifacts(folder)
    path = folder / "summary.json"
    metadata = json.loads(path.read_text(encoding="utf-8"))
    metadata["paired_ab_contrasts"][0]["mean_difference"] += 1.0
    path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="mean difference mismatch"):
        audit_artifacts(folder)


def test_independent_audit_rejects_changed_likelihood_equation_metadata(tmp_path):
    folder = tmp_path / "artifacts"
    _write_artifacts(folder)
    path = folder / "summary.json"
    metadata = json.loads(path.read_text(encoding="utf-8"))
    metadata["likelihood_equations"]["B_fractional_sigma_H0"] = "omit normalization"
    path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="likelihood equation metadata"):
        audit_artifacts(folder)
