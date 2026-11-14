"""End-to-end tests for the lightweight bandwidth pilot and independent audit."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_parent_density_bandwidth_pilot import run_bandwidth_pilot, write_artifacts
from audit_parent_density_bandwidth_pilot import audit_bandwidth_pilot
from select_bandwidth_train_holdout import select_bandwidth
from audit_bandwidth_train_holdout import audit_selection


def test_parent_density_bandwidth_pilot_audit_and_selection(tmp_path):
    source = tmp_path / "source"
    attempts, samples, truth, config = run_bandwidth_pilot(
        n_replicates=8, n_parent=250, seed=7713, grid_points=201,
        bandwidth_multipliers=(0.05, 0.10, 0.20),
    )
    write_artifacts(attempts, samples, truth, config, source)
    report = audit_bandwidth_pilot(source)
    assert report["status"] == "PASS"
    assert report["bandwidth_attempts_reconstructed"] == 24

    selected = select_bandwidth(
        pd.read_csv(source / "attempt_metrics.csv"),
        candidate_multipliers=config["bandwidth_multipliers"],
        n_replicates=8, split_seed=22, training_fraction=0.5,
        source_manifest_sha256=__import__("hashlib").sha256(
            (source / "manifest.json").read_bytes()
        ).hexdigest(),
    )
    selected["source_attempts_sha256"] = __import__("hashlib").sha256(
        (source / "attempt_metrics.csv").read_bytes()
    ).hexdigest()
    selected["source_config_sha256"] = __import__("hashlib").sha256(
        (source / "config.json").read_bytes()
    ).hexdigest()
    selection_path = tmp_path / "selection.json"
    import json
    selection_path.write_text(json.dumps(selected), encoding="utf-8")
    selection_audit = audit_selection(source, selection_path)
    assert selection_audit["status"] == "PASS"
    assert selection_audit["training_replicates"] == 4
    assert selection_audit["holdout_replicates"] == 4
