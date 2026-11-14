"""Independently audit a training/holdout KDE bandwidth selection artifact."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit_selection(source_dir: Path, selection_path: Path) -> dict:
    source_dir = source_dir.resolve()
    selection = json.loads(selection_path.read_text(encoding="utf-8"))
    manifest_path = source_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for record in manifest["files"].values():
        candidate = source_dir / record["path"]
        if not candidate.is_file() or candidate.stat().st_size != int(record["bytes"]):
            raise AssertionError("source manifest file set or size mismatch")
        if _hash(candidate) != record["sha256"]:
            raise AssertionError("source manifest checksum mismatch")

    config = json.loads((source_dir / "config.json").read_text(encoding="utf-8"))
    attempts = pd.read_csv(source_dir / "attempt_metrics.csv")
    candidates = sorted(set(float(x) for x in config["bandwidth_multipliers"]))
    n_replicates = int(config["n_replicates"])
    split_seed = int(selection["split_seed"])
    fraction = float(selection["training_fraction"])
    n_train = min(max(int(round(n_replicates * fraction)), 2), n_replicates - 2)
    rng = np.random.default_rng(split_seed)
    permutation = rng.permutation(np.arange(n_replicates, dtype=int))
    expected_train = sorted(int(x) for x in permutation[:n_train])
    expected_holdout = sorted(int(x) for x in permutation[n_train:])
    if selection["training_replicate_ids"] != expected_train:
        raise AssertionError("training replicate partition does not match frozen seed split")
    if selection["holdout_replicate_ids"] != expected_holdout:
        raise AssertionError("holdout replicate partition does not match frozen seed split")
    if set(expected_train) & set(expected_holdout):
        raise AssertionError("training and holdout replicate IDs overlap")

    gamma = float(selection["gamma"])
    retention = float(selection["target_retention"])
    validation_size = int(selection["validation_size"])
    data = attempts.loc[
        (attempts.arm == "parent_sample_kde")
        & (attempts.degree == 0)
        & np.isclose(attempts.gamma.astype(float), gamma)
        & np.isclose(attempts.target_retention.astype(float), retention)
        & (attempts.n_validation.astype(int) == validation_size)
        & attempts.bandwidth_multiplier.astype(float).isin(candidates)
    ].copy()
    if data.duplicated(["replicate", "bandwidth_multiplier"]).any():
        raise AssertionError("duplicate source rows")
    if set(data.replicate.astype(int)) != set(range(n_replicates)):
        raise AssertionError("source replicate coverage mismatch")
    if len(data) != n_replicates * len(candidates):
        raise AssertionError("source candidate grid is incomplete")
    data["iae"] = pd.to_numeric(data.iae, errors="coerce")
    if not np.isfinite(data.iae).all() or not data.status.eq("SUCCESS").all():
        raise AssertionError("source has failed or non-finite density-recovery attempts")

    expected_training = []
    expected_holdout_scores = []
    for multiplier in candidates:
        train_values = data.loc[
            data.replicate.isin(expected_train)
            & np.isclose(data.bandwidth_multiplier.astype(float), multiplier), "iae"
        ].to_numpy(float)
        hold_values = data.loc[
            data.replicate.isin(expected_holdout)
            & np.isclose(data.bandwidth_multiplier.astype(float), multiplier), "iae"
        ].to_numpy(float)
        expected_training.append({
            "bandwidth_multiplier": multiplier,
            "training_mean_iae": float(np.mean(train_values)),
            "training_sd_iae": float(np.std(train_values, ddof=1)),
            "training_replicates": len(train_values),
        })
        expected_holdout_scores.append({
            "bandwidth_multiplier": multiplier,
            "holdout_mean_iae": float(np.mean(hold_values)),
            "holdout_sd_iae": float(np.std(hold_values, ddof=1)) if len(hold_values) > 1 else None,
            "holdout_mcse_iae": float(np.std(hold_values, ddof=1) / np.sqrt(len(hold_values)))
                if len(hold_values) > 1 else None,
            "holdout_replicates": len(hold_values),
        })
    expected_selected = min(
        expected_training, key=lambda row: (row["training_mean_iae"], row["bandwidth_multiplier"])
    )["bandwidth_multiplier"]
    if selection["candidate_multipliers"] != candidates:
        raise AssertionError("candidate multiplier set mismatch")
    if selection["selected_multiplier"] != expected_selected:
        raise AssertionError("selected multiplier was not minimized on training data")
    if selection.get("schema_version") != 1:
        raise AssertionError("unsupported selection artifact schema")
    if selection.get("status") != "PASS_SELECTION_ONLY_NOT_CONFIRMATORY":
        raise AssertionError("selection artifact status is not the declared non-confirmatory status")

    def _check_score_rows(observed_rows, expected_rows, label):
        if not isinstance(observed_rows, list) or len(observed_rows) != len(expected_rows):
            raise AssertionError(f"{label} score row count mismatch")
        for index, (observed, expected) in enumerate(zip(observed_rows, expected_rows)):
            if set(observed) != set(expected):
                raise AssertionError(f"{label} score fields mismatch at row {index}")
            for key, expected_value in expected.items():
                observed_value = observed[key]
                if expected_value is None:
                    if observed_value is not None:
                        raise AssertionError(f"{label} {key} should be null at row {index}")
                elif isinstance(expected_value, (int, float)):
                    if not np.isclose(
                        float(observed_value), float(expected_value),
                        rtol=1e-12, atol=1e-14,
                    ):
                        raise AssertionError(f"{label} {key} mismatch at row {index}")
                elif observed_value != expected_value:
                    raise AssertionError(f"{label} {key} mismatch at row {index}")

    _check_score_rows(selection["training_scores"], expected_training, "training")
    _check_score_rows(selection["holdout_scores"], expected_holdout_scores, "holdout")
    if selection["source_manifest_sha256"] != _hash(manifest_path):
        raise AssertionError("source manifest hash mismatch")
    if selection["source_attempts_sha256"] != _hash(source_dir / "attempt_metrics.csv"):
        raise AssertionError("source attempts hash mismatch")
    if selection["source_config_sha256"] != _hash(source_dir / "config.json"):
        raise AssertionError("source config hash mismatch")
    selected_holdout = next(
        row for row in expected_holdout_scores
        if row["bandwidth_multiplier"] == expected_selected
    )
    observed_selected_holdout = selection["selected_holdout_score"]
    if set(observed_selected_holdout) != set(selected_holdout):
        raise AssertionError("selected holdout summary fields mismatch")
    for key, expected_value in selected_holdout.items():
        observed_value = observed_selected_holdout[key]
        if expected_value is None:
            if observed_value is not None:
                raise AssertionError(f"selected holdout {key} should be null")
        elif isinstance(expected_value, (int, float)):
            if not np.isclose(
                float(observed_value), float(expected_value), rtol=1e-12, atol=1e-14,
            ):
                raise AssertionError(f"selected holdout {key} mismatch")
        elif observed_value != expected_value:
            raise AssertionError(f"selected holdout {key} mismatch")
    return {
        "status": "PASS",
        "selected_multiplier": expected_selected,
        "training_replicates": len(expected_train),
        "holdout_replicates": len(expected_holdout),
        "candidate_count": len(candidates),
        "source_attempts_sha256": _hash(source_dir / "attempt_metrics.csv"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("selection_json", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit_selection(args.source_dir, args.selection_json)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
