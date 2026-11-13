"""Tests for leakage-free KDE bandwidth training/holdout selection."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from select_bandwidth_train_holdout import select_bandwidth
from audit_bandwidth_train_holdout import audit_selection


def _fixture_attempts():
    rows = []
    for replicate in range(12):
        for multiplier in (0.05, 0.10, 0.20):
            # 0.10 is best on training replicates. Holdout is deliberately
            # worse for it, proving selection is not re-optimized on holdout.
            if replicate < 8:
                iae = {0.05: 0.12, 0.10: 0.08, 0.20: 0.16}[multiplier] + replicate * 1e-4
            else:
                iae = {0.05: 0.10, 0.10: 0.30, 0.20: 0.14}[multiplier] + replicate * 1e-4
            rows.append({
                "replicate": replicate, "arm": "parent_sample_kde", "degree": 0,
                "gamma": 0.0, "target_retention": 0.25, "n_validation": 300,
                "bandwidth_multiplier": multiplier, "status": "SUCCESS", "iae": iae,
            })
    return pd.DataFrame(rows)


def test_selection_uses_training_only_and_reports_holdout_without_reselection():
    result = select_bandwidth(
        _fixture_attempts(), candidate_multipliers=[0.05, 0.10, 0.20],
        n_replicates=12, split_seed=17, training_fraction=0.67,
    )
    assert result["selected_multiplier"] == 0.10
    assert set(result["training_replicate_ids"]).isdisjoint(result["holdout_replicate_ids"])
    assert len(result["training_replicate_ids"]) + len(result["holdout_replicate_ids"]) == 12
    selected_holdout = result["selected_holdout_score"]["holdout_mean_iae"]
    assert selected_holdout > next(
        row["holdout_mean_iae"] for row in result["holdout_scores"]
        if row["bandwidth_multiplier"] == 0.05
    )


def test_independent_auditor_reconstructs_selection_and_hashes(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    attempts = _fixture_attempts()
    attempts.to_csv(source / "attempt_metrics.csv", index=False)
    config = {
        "n_replicates": 12,
        "bandwidth_multipliers": [0.05, 0.10, 0.20],
    }
    (source / "config.json").write_text(json.dumps(config), encoding="utf-8")
    import hashlib
    files = {}
    for name in ("attempt_metrics.csv", "config.json"):
        path = source / name
        files[name] = {
            "path": name, "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
    (source / "manifest.json").write_text(json.dumps({"files": files}), encoding="utf-8")
    manifest_hash = hashlib.sha256((source / "manifest.json").read_bytes()).hexdigest()
    result = select_bandwidth(
        attempts, candidate_multipliers=config["bandwidth_multipliers"],
        n_replicates=12, split_seed=17, training_fraction=0.67,
        source_manifest_sha256=manifest_hash,
    )
    result["source_attempts_sha256"] = hashlib.sha256((source / "attempt_metrics.csv").read_bytes()).hexdigest()
    result["source_config_sha256"] = hashlib.sha256((source / "config.json").read_bytes()).hexdigest()
    selection_path = tmp_path / "selection.json"
    selection_path.write_text(json.dumps(result), encoding="utf-8")
    report = audit_selection(source, selection_path)
    assert report["status"] == "PASS"
    assert report["selected_multiplier"] == 0.10
