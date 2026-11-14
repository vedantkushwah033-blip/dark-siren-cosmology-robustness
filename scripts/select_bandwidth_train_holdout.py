"""Select a KDE Scott multiplier on independent training replicates.

Selection uses analytic parent-density integrated absolute error (IAE) only.
It never reads H0 likelihoods or posterior outcomes. A frozen holdout split
estimates out-of-sample density-recovery performance.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _verify_source_manifest(source: Path) -> tuple[dict, str]:
    manifest_path = source / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for record in manifest["files"].values():
        path = source / record["path"]
        if not path.is_file() or path.stat().st_size != int(record["bytes"]):
            raise ValueError(f"source artifact missing or size mismatch: {record['path']}")
        if _sha256(path) != record["sha256"]:
            raise ValueError(f"source artifact checksum mismatch: {record['path']}")
    return manifest, _sha256(manifest_path)


def select_bandwidth(
    attempts: pd.DataFrame,
    *,
    candidate_multipliers: list[float],
    n_replicates: int,
    split_seed: int = 20261113,
    training_fraction: float = 0.70,
    gamma: float = 0.0,
    target_retention: float = 0.25,
    validation_size: int = 300,
    source_manifest_sha256: str = "",
) -> dict:
    """Choose a multiplier on training replicate IDs and score it on holdout IDs."""
    candidates = sorted(set(float(x) for x in candidate_multipliers))
    if len(candidates) < 2 or any(not np.isfinite(x) or x <= 0 for x in candidates):
        raise ValueError("at least two distinct positive finite bandwidth candidates are required")
    if n_replicates < 4:
        raise ValueError("at least four independent replicates are required for train/holdout")
    if not 0.5 <= training_fraction <= 0.9:
        raise ValueError("training_fraction must lie in [0.5, 0.9]")
    selected = attempts.loc[
        (attempts["arm"] == "parent_sample_kde")
        & (attempts["degree"] == 0)
        & np.isclose(attempts["gamma"].astype(float), gamma)
        & np.isclose(attempts["target_retention"].astype(float), target_retention)
        & (attempts["n_validation"].astype(int) == int(validation_size))
        & attempts["bandwidth_multiplier"].astype(float).isin(candidates)
    ].copy()
    expected_ids = set(range(n_replicates))
    observed_ids = set(selected["replicate"].astype(int).unique())
    if observed_ids != expected_ids:
        raise ValueError("selected parent-KDE rows do not cover the declared replicate IDs")
    if selected.duplicated(["replicate", "bandwidth_multiplier"]).any():
        raise ValueError("duplicate replicate/bandwidth rows in source attempts")
    if len(selected) != n_replicates * len(candidates):
        raise ValueError("incomplete replicate-by-bandwidth grid")
    selected["replicate"] = selected["replicate"].astype(int)
    selected["bandwidth_multiplier"] = selected["bandwidth_multiplier"].astype(float)
    selected["iae"] = pd.to_numeric(selected["iae"], errors="coerce")
    if not np.isfinite(selected["iae"]).all() or not selected["status"].eq("SUCCESS").all():
        raise ValueError("training/holdout selection requires every parent-KDE attempt to succeed")

    rng = np.random.default_rng(int(split_seed))
    shuffled = rng.permutation(np.array(sorted(expected_ids), dtype=int))
    n_train = int(round(n_replicates * training_fraction))
    n_train = min(max(n_train, 2), n_replicates - 2)
    training_ids = sorted(int(x) for x in shuffled[:n_train])
    holdout_ids = sorted(int(x) for x in shuffled[n_train:])
    training = selected.loc[selected.replicate.isin(training_ids)]
    holdout = selected.loc[selected.replicate.isin(holdout_ids)]

    training_scores = []
    for multiplier in candidates:
        values = training.loc[
            np.isclose(training.bandwidth_multiplier, multiplier), "iae"
        ].to_numpy(float)
        if len(values) != len(training_ids):
            raise ValueError("incomplete training observations for a bandwidth candidate")
        training_scores.append({
            "bandwidth_multiplier": multiplier,
            "training_mean_iae": float(np.mean(values)),
            "training_sd_iae": float(np.std(values, ddof=1)),
            "training_replicates": len(values),
        })
    selected_multiplier = min(
        training_scores, key=lambda row: (row["training_mean_iae"], row["bandwidth_multiplier"])
    )["bandwidth_multiplier"]

    holdout_scores = []
    for multiplier in candidates:
        values = holdout.loc[
            np.isclose(holdout.bandwidth_multiplier, multiplier), "iae"
        ].to_numpy(float)
        if len(values) != len(holdout_ids):
            raise ValueError("incomplete holdout observations for a bandwidth candidate")
        holdout_scores.append({
            "bandwidth_multiplier": multiplier,
            "holdout_mean_iae": float(np.mean(values)),
            "holdout_sd_iae": float(np.std(values, ddof=1)) if len(values) > 1 else None,
            "holdout_mcse_iae": float(np.std(values, ddof=1) / np.sqrt(len(values)))
                if len(values) > 1 else None,
            "holdout_replicates": len(values),
        })
    return {
        "schema_version": 1,
        "status": ("BOUNDARY_LIMITED_EXPAND_GRID" if selected_multiplier in (candidates[0], candidates[-1]) else "PASS_SELECTION_ONLY_NOT_CONFIRMATORY"),
        "selection_eligible": bool(selected_multiplier not in (candidates[0], candidates[-1])),
        "boundary_limited": bool(selected_multiplier in (candidates[0], candidates[-1])),
        "estimand": "analytic merger-weighted parent-host redshift density",
        "selection_metric": "integrated absolute error against analytic parent density",
        "selection_population": "parent_sample_kde only; no H0 likelihood or posterior metric",
        "gamma": float(gamma),
        "target_retention": float(target_retention),
        "validation_size": int(validation_size),
        "candidate_multipliers": candidates,
        "split_seed": int(split_seed),
        "training_fraction": float(training_fraction),
        "training_replicate_ids": training_ids,
        "holdout_replicate_ids": holdout_ids,
        "training_scores": training_scores,
        "selected_multiplier": float(selected_multiplier),
        "holdout_scores": holdout_scores,
        "selected_holdout_score": next(
            row for row in holdout_scores
            if row["bandwidth_multiplier"] == selected_multiplier
        ),
        "source_manifest_sha256": source_manifest_sha256,
        "limitations": [
            "Synthetic parent-density recovery only; not an H0 calibration or observational validation.",
            "The selected multiplier is conditional on this parent DGP, sample size, grid, and candidate set.",
            "Holdout results are descriptive and must be reported without re-selecting on holdout.",
            "Arm-specific weighted Scott widths may differ because effective sample size and weighted variance differ.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--split-seed", type=int, default=20261113)
    parser.add_argument("--training-fraction", type=float, default=0.70)
    parser.add_argument("--gamma", type=float, default=0.0)
    parser.add_argument("--target-retention", type=float, default=0.25)
    parser.add_argument("--validation-size", type=int, default=300)
    args = parser.parse_args()

    source = args.source_dir.resolve()
    manifest, manifest_hash = _verify_source_manifest(source)
    config = json.loads((source / "config.json").read_text(encoding="utf-8"))
    attempts = pd.read_csv(source / "attempt_metrics.csv")
    result = select_bandwidth(
        attempts,
        candidate_multipliers=[float(x) for x in config["bandwidth_multipliers"]],
        n_replicates=int(config["n_replicates"]),
        split_seed=args.split_seed,
        training_fraction=args.training_fraction,
        gamma=args.gamma,
        target_retention=args.target_retention,
        validation_size=args.validation_size,
        source_manifest_sha256=manifest_hash,
    )
    result["source_attempts_sha256"] = _sha256(source / "attempt_metrics.csv")
    result["source_config_sha256"] = _sha256(source / "config.json")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
