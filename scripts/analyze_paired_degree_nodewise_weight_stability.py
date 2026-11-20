"""Reconstruct posterior-mass-weighted node-wise IPW stability diagnostics.

Input is an extracted paired-degree pilot artifact directory containing
degree1/ and degree2/ child-run trees. The script verifies its reconstructed
worst-node diagnostics against the stored replicate diagnostic ledgers.
"""
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    """Left-continuous weighted empirical quantile, with normalized weights."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if values.ndim != 1 or weights.shape != values.shape or values.size == 0:
        raise ValueError("values and weights must be matching non-empty 1D arrays")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(weights)):
        raise ValueError("values and weights must be finite")
    if np.any(weights < 0) or not np.isfinite(weights.sum()) or weights.sum() <= 0:
        raise ValueError("weights must be non-negative with positive finite total")
    if not np.isfinite(q) or not 0 <= q <= 1:
        raise ValueError("q must lie in [0, 1]")
    normalized = weights / weights.sum()
    order = np.argsort(values, kind="mergesort")
    cumulative = np.cumsum(normalized[order])
    index = min(int(np.searchsorted(cumulative, q, side="left")), values.size - 1)
    return float(values[order[index]])


def compute_nodewise_diagnostics(
    catalogue_redshift: np.ndarray,
    nodes: np.ndarray,
    posterior_weights: np.ndarray,
    *,
    degree: int,
    coordinate_center: float,
    coordinate_scale: float,
    quadratic_center: float = 0.0,
) -> dict[str, float]:
    """Compute node-wise IPW ESS/tail mass and posterior-weighted summaries."""
    z = np.asarray(catalogue_redshift, dtype=float)
    nodes = np.asarray(nodes, dtype=float)
    posterior_weights = np.asarray(posterior_weights, dtype=float)
    if z.ndim != 1 or z.size < 1 or not np.all(np.isfinite(z)):
        raise ValueError("catalogue_redshift must be a finite non-empty vector")
    if degree not in (1, 2) or nodes.ndim != 2 or nodes.shape[1] != degree + 1:
        raise ValueError("nodes must have shape (n_nodes, degree + 1), degree 1 or 2")
    if not np.all(np.isfinite(nodes)):
        raise ValueError("quadrature nodes must be finite")
    if posterior_weights.shape != (nodes.shape[0],):
        raise ValueError("posterior_weights must match the quadrature-node count")
    if not np.all(np.isfinite(posterior_weights)) or np.any(posterior_weights < 0):
        raise ValueError("posterior_weights must be finite and non-negative")
    if posterior_weights.sum() <= 0:
        raise ValueError("posterior_weights must have positive total mass")
    if not np.isfinite(coordinate_center) or not np.isfinite(coordinate_scale) or coordinate_scale <= 0:
        raise ValueError("coordinate center/scale must be finite with positive scale")
    if not np.isfinite(quadratic_center):
        raise ValueError("quadratic_center must be finite")

    x = (z - coordinate_center) / coordinate_scale
    columns = [np.ones_like(x), x]
    if degree == 2:
        columns.append(x * x - quadratic_center)
    design = np.column_stack(columns)
    probability = expit(design @ nodes.T)
    if not np.all(np.isfinite(probability)) or np.any(probability <= 0):
        raise ArithmeticError("node-specific inclusion probabilities are numerically unusable")

    # This is the shared host-rate/base-population weighting convention used
    # by run_quadratic_five_arm_smoke.py in this synthetic benchmark.
    base_weight = 1.0 + 1.5 * (z - 0.01) / 0.39
    if np.any(base_weight <= 0) or not np.all(np.isfinite(base_weight)):
        raise ArithmeticError("base-population weights must be finite and positive")
    raw_weights = base_weight[:, None] / probability
    normalized = raw_weights / raw_weights.sum(axis=0, keepdims=True)
    node_ess = 1.0 / np.sum(normalized * normalized, axis=0)
    n_top = max(1, int(math.ceil(0.05 * z.size)))
    node_top5 = np.sort(normalized, axis=0)[-n_top:, :].sum(axis=0)
    qweights = posterior_weights / posterior_weights.sum()

    return {
        "minimum_node_kish_ess": float(np.min(node_ess)),
        "posterior_weighted_q05_node_kish_ess": weighted_quantile(node_ess, qweights, 0.05),
        "posterior_weighted_q50_node_kish_ess": weighted_quantile(node_ess, qweights, 0.50),
        "posterior_weighted_q90_node_kish_ess": weighted_quantile(node_ess, qweights, 0.90),
        "maximum_node_top_5_percent_mass": float(np.max(node_top5)),
        "posterior_weighted_q50_node_top_5_percent_mass": weighted_quantile(node_top5, qweights, 0.50),
        "posterior_weighted_q90_node_top_5_percent_mass": weighted_quantile(node_top5, qweights, 0.90),
        "quadrature_ess": float(1.0 / np.sum(qweights * qweights)),
    }


def analyze_artifact(artifact_dir: str | Path, *, tolerance: float = 1e-9) -> pd.DataFrame:
    """Analyze both degree trees and verify worst-node metrics against ledgers."""
    root = Path(artifact_dir)
    rows: list[dict[str, float | int]] = []
    for degree in (1, 2):
        degree_root = root / f"degree{degree}"
        ledger_path = degree_root / "replicate_weight_diagnostics.csv"
        if not ledger_path.is_file():
            raise FileNotFoundError(f"missing stored diagnostic ledger: {ledger_path}")
        ledger = pd.read_csv(ledger_path)
        ledger = ledger[ledger["arm"] == "E_uncertainty_integrated_ipw"]
        stored = {
            (int(r.replicate), float(r.gamma), int(r.validation_size)): r
            for r in ledger.itertuples(index=False)
        }
        seen: set[tuple[int, float, int]] = set()
        for child in sorted((degree_root / "replicates").glob("rep*")):
            match = re.fullmatch(r"rep(\d+)_g(m1|p0|p1)", child.name)
            if match is None:
                continue
            replicate = int(match.group(1))
            gamma = {"m1": -1.0, "p0": 0.0, "p1": 1.0}[match.group(2)]
            config = json.loads((child / "config.json").read_text())
            calibration = config["quadratic_calibration"]
            catalogue = pd.read_csv(child / "observable" / "catalogue_observed.csv")
            catalogue_z = catalogue["redshift"].to_numpy(float)
            with np.load(child / "results" / "log_likelihoods.npz") as archive:
                for n_validation in config["validation_sizes"]:
                    key = (replicate, gamma, int(n_validation))
                    if key not in stored:
                        raise ValueError(f"no stored E-arm diagnostic for {key} in degree {degree}")
                    diagnostics = compute_nodewise_diagnostics(
                        catalogue_z,
                        archive[f"posterior_nodes_n{n_validation}"],
                        archive[f"posterior_weights_n{n_validation}"],
                        degree=degree,
                        coordinate_center=float(calibration["coordinate_center"]),
                        coordinate_scale=float(calibration["coordinate_scale"]),
                        quadratic_center=float(calibration["x2_center"]) if degree == 2 else 0.0,
                    )
                    ledger_row = stored[key]
                    if abs(diagnostics["minimum_node_kish_ess"] - ledger_row.minimum_node_kish_ess) > tolerance:
                        raise ArithmeticError(f"minimum-node ESS mismatch for degree {degree}, cell {key}")
                    if abs(diagnostics["maximum_node_top_5_percent_mass"] - ledger_row.maximum_node_top_5_percent_mass) > tolerance:
                        raise ArithmeticError(f"maximum top-5% mass mismatch for degree {degree}, cell {key}")
                    if abs(diagnostics["quadrature_ess"] - ledger_row.quadrature_ess) > tolerance:
                        raise ArithmeticError(f"quadrature ESS mismatch for degree {degree}, cell {key}")
                    rows.append({
                        "degree": degree,
                        "replicate": replicate,
                        "gamma": gamma,
                        "validation_size": int(n_validation),
                        **diagnostics,
                    })
                    seen.add(key)
        if seen != set(stored):
            missing = sorted(set(stored) - seen)
            extra = sorted(seen - set(stored))
            raise ValueError(f"degree {degree} artifact/ledger cell mismatch; missing={missing[:5]}, extra={extra[:5]}")
    result = pd.DataFrame(rows).sort_values(["degree", "gamma", "validation_size", "replicate"])
    return result.reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-dir", required=True, type=Path,
                        help="Path to the extracted paired-degree pilot artifact")
    parser.add_argument("--output", type=Path, default=Path("nodewise_weight_stability.csv"))
    args = parser.parse_args()
    result = analyze_artifact(args.artifact_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(json.dumps({
        "status": "PASS",
        "rows": int(len(result)),
        "degrees": sorted(int(x) for x in result.degree.unique()),
        "output": str(args.output),
        "interpretation": "Posterior-mass-weighted node diagnostics; no calibration or estimator-superiority claim.",
    }, indent=2))


if __name__ == "__main__":
    main()
