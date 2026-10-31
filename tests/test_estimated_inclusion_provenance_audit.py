from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

MODULE_PATH = Path(__file__).parents[1] / "scripts" / "audit_estimated_inclusion_provenance.py"
SPEC = importlib.util.spec_from_file_location("estimated_inclusion_provenance_audit", MODULE_PATH)
AUDITOR = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(AUDITOR)


def _write_fixture(tmp_path: Path) -> tuple[Path, Path]:
    from dark_siren_robustness.generative_selection import generate_synthetic_population

    master_seed = 20261030
    replicate = 0
    sim_seed, catalogue_seed, audit_seed = AUDITOR._component_seeds(master_seed, replicate)
    n_parent = 500
    n_detected = 30
    z_min, z_max = 0.01, 0.40
    h0_true, sigma_mpc, threshold = 70.0, 140.0, 250.0
    generated = generate_synthetic_population(
        seed=sim_seed,
        n_parent=n_parent,
        n_detected=n_detected,
        h0_true=h0_true,
        sigma_mpc=sigma_mpc,
        detection_threshold_mpc=threshold,
        z_min=z_min,
        z_max=z_max,
    )
    parent_z = generated["parent_truth"]["redshift"].to_numpy(dtype=float)
    x = (parent_z - (z_min + z_max) / 2) / ((z_max - z_min) / 2)
    from scipy.special import expit

    true_p = expit(0.25 - 1.2 * x + 0.3 * (x**2 - 1 / 3))
    retained = np.random.default_rng(catalogue_seed).random(n_parent) < true_p
    order = np.random.default_rng(audit_seed).permutation(n_parent)
    m = 100
    inference = {
        "validation_redshift": parent_z[order[:m]].tolist(),
        "validation_included": retained[order[:m]].astype(int).tolist(),
        "catalogue_redshift": parent_z[retained].tolist(),
        "observed_distances_mpc": generated["gw_events_observed"]["distance_mpc"].tolist(),
    }
    scoring = {
        "gw_event_host_object_id": generated["gw_event_truth"]["host_object_id"].tolist(),
        "gw_event_host_redshift": generated["gw_event_truth"]["host_redshift"].tolist(),
        "gw_event_latent_mean_distance_mpc": generated["gw_event_truth"]["latent_mean_distance_mpc"].tolist(),
        "parent_redshift": parent_z.tolist(),
        "true_inclusion_probability": true_p.tolist(),
    }
    rows = [{
        "replicate": replicate,
        "replicate_seed": sim_seed,
        "catalogue_seed": catalogue_seed,
        "audit_seed": audit_seed,
        "arm": f"estimated_{m}",
        "validation_size": m,
        "status": "success",
        "inference_audit_payload_json": json.dumps(inference),
        "scoring_audit_payload_json": json.dumps(scoring),
    }]
    replicate_path = tmp_path / "replicates.csv"
    pd.DataFrame(rows).to_csv(replicate_path, index=False)
    summary_path = tmp_path / "summary.json"
    summary_path.write_text(json.dumps({
        "seed": master_seed,
        "n_parent": n_parent,
        "n_detected": n_detected,
        "h0_true": h0_true,
        "sigma_mpc": sigma_mpc,
        "detection_threshold_mpc": threshold,
        "redshift_grid": [z_min, z_max, 401],
        "true_inclusion_model": {
            "intercept": 0.25, "slope": -1.2, "quadratic": 0.3,
        },
    }), encoding="utf-8")
    return replicate_path, summary_path

def test_seed_level_provenance_audit_reconstructs_parent_catalogue_and_validation(tmp_path):
    replicates, summary = _write_fixture(tmp_path)
    report = AUDITOR.audit_provenance(replicates, summary)
    assert report["status"] == "PASS"
    assert report["replicates_reconstructed"] == 1
    assert report["estimated_arm_rows_reconstructed"] == 1


def test_seed_level_provenance_audit_rejects_tampered_validation_label(tmp_path):
    replicates, summary = _write_fixture(tmp_path)
    frame = pd.read_csv(replicates)
    payload = json.loads(frame.loc[0, "inference_audit_payload_json"])
    payload["validation_included"][0] = 1 - payload["validation_included"][0]
    frame.loc[0, "inference_audit_payload_json"] = json.dumps(payload)
    frame.to_csv(replicates, index=False)
    with pytest.raises(AssertionError, match="validation labels mismatch"):
        AUDITOR.audit_provenance(replicates, summary)


def test_seed_level_provenance_audit_rejects_wrong_component_seed(tmp_path):
    replicates, summary = _write_fixture(tmp_path)
    frame = pd.read_csv(replicates)
    frame.loc[0, "audit_seed"] = int(frame.loc[0, "audit_seed"]) + 1
    frame.to_csv(replicates, index=False)
    with pytest.raises(AssertionError, match="seed schedule"):
        AUDITOR.audit_provenance(replicates, summary)



def test_seed_level_provenance_audit_rejects_tampered_event_distance(tmp_path):
    replicates, summary = _write_fixture(tmp_path)
    frame = pd.read_csv(replicates)
    payload = json.loads(frame.loc[0, "inference_audit_payload_json"])
    payload["observed_distances_mpc"][0] += 0.25
    frame.loc[0, "inference_audit_payload_json"] = json.dumps(payload)
    frame.to_csv(replicates, index=False)
    with pytest.raises(AssertionError, match="detected GW distance observations"):
        AUDITOR.audit_provenance(replicates, summary)


def test_seed_level_provenance_audit_rejects_tampered_gw_host_identity(tmp_path):
    replicates, summary = _write_fixture(tmp_path)
    frame = pd.read_csv(replicates)
    payload = json.loads(frame.loc[0, "scoring_audit_payload_json"])
    payload["gw_event_host_object_id"][0] = (payload["gw_event_host_object_id"][0] + 1) % 500
    frame.loc[0, "scoring_audit_payload_json"] = json.dumps(payload)
    frame.to_csv(replicates, index=False)
    with pytest.raises(AssertionError, match="detected GW host IDs"):
        AUDITOR.audit_provenance(replicates, summary)


def test_seed_level_provenance_audit_reconstructs_failed_estimated_arm(tmp_path):
    replicates, summary = _write_fixture(tmp_path)
    frame = pd.read_csv(replicates)
    frame.loc[0, "status"] = "failed"
    frame.loc[0, "error"] = "synthetic injected failure after data generation"
    frame.to_csv(replicates, index=False)
    report = AUDITOR.audit_provenance(replicates, summary)
    assert report["status"] == "PASS"
    assert report["failed_estimated_arm_rows_reconstructed"] == 1
    assert report["estimated_arm_rows_reconstructed"] == 1
