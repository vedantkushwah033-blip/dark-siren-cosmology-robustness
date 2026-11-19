"""Pool two independently audited paired five-arm 30-event pilot batches.

Checks matching scientific configuration, complete lifecycle/audit denominators,
and disjoint deterministic seed pairs before constructing replicate-level
pooled metrics and paired contrasts. Descriptive synthetic pilot only.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm, t

ARMS = [
    "A_catalogue_only",
    "B_parent_oracle",
    "C_known_probability_ipw",
    "D_map_plugin_ipw",
    "E_uncertainty_integrated_ipw",
]
CONTRASTS = [
    ("E_uncertainty_integrated_ipw", "D_map_plugin_ipw", "E_minus_D"),
    ("E_uncertainty_integrated_ipw", "C_known_probability_ipw", "E_minus_C"),
    ("D_map_plugin_ipw", "C_known_probability_ipw", "D_minus_C"),
    ("B_parent_oracle", "C_known_probability_ipw", "B_minus_C"),
]
CONFIG_KEYS = [
    "experiment", "curvature_beta", "detection_threshold_mpc", "failure_bounds",
    "gammas", "h0_points", "h0_prior_support", "h0_true", "n_detected",
    "n_parent", "quadrature_order", "scott_multiplier", "sigma_mpc",
    "validation_frame", "validation_sizes", "replicates",
]


def wilson_interval(successes: int, n: int) -> tuple[float, float]:
    """Wilson 95% interval for a binomial proportion."""
    if n <= 0:
        return float("nan"), float("nan")
    z = float(norm.ppf(0.975))
    p = successes / n
    denominator = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denominator
    half = z * np.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / denominator
    return float(center - half), float(center + half)


def paired_t_interval(values: np.ndarray) -> tuple[float, float, float, float, int]:
    """Mean, SE, and two-sided 95% Student-t interval for paired differences."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    n = int(values.size)
    if n == 0:
        return (float("nan"),) * 4 + (0,)
    mean = float(values.mean())
    if n == 1:
        return mean, float("nan"), float("nan"), float("nan"), n
    standard_error = float(values.std(ddof=1) / np.sqrt(n))
    critical = float(t.ppf(0.975, n - 1))
    return (
        mean,
        standard_error,
        mean - critical * standard_error,
        mean + critical * standard_error,
        n,
    )


def load_and_validate(artifact_dir: Path):
    root = Path(artifact_dir)
    config = json.loads((root / "config.json").read_text(encoding="utf-8"))
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    summaries = pd.read_csv(root / "replicate_posterior_summaries.csv")
    audits = pd.read_csv(root / "replicate_audits.csv")
    independent = json.loads(
        (root / "independent_orchestration_audit.json").read_text(encoding="utf-8")
    )

    replicates = int(config["replicates"])
    expected_cells = replicates * len(config["gammas"]) * len(config["validation_sizes"]) * len(ARMS)
    expected_children = replicates * len(config["gammas"])

    assert len(ledger) == expected_cells, "planned lifecycle denominator mismatch"
    assert ledger["cell_id"].is_unique, "duplicate lifecycle cell ID"
    assert set(ledger["status"]) == {"SUCCESS"}, "batch has failed/non-evaluable cells"
    assert set(ledger["audit_status"]) == {"PASS"}, "batch has child audit failures"
    assert len(audits) == expected_children, "child audit denominator mismatch"
    assert set(audits["status"]) == {"PASS"}, "child audit table contains non-PASS status"
    assert len(summaries) == expected_cells, "posterior summary denominator mismatch"
    assert set(summaries["status"]) == {"SUCCESS"}, "posterior summaries contain non-success rows"
    assert set(summaries["audit_status"]) == {"PASS"}, "posterior summaries lack audit PASS"
    assert independent["status"] == "PASS", "independent orchestration audit did not pass"
    assert int(independent["planned_cells"]) == expected_cells
    assert int(independent["successful_child_runs"]) == expected_children

    seed_pairs = {
        (int(row["analysis_seed"]), int(row["validation_seed"]))
        for row in config["seed_schedule"]
    }
    assert len(seed_pairs) == replicates, "duplicate seed pair within batch"
    summary_keys = ["replicate", "gamma", "validation_size", "arm"]
    assert not summaries.duplicated(summary_keys).any(), "duplicate posterior summary key"
    assert len(summaries) == expected_cells
    return config, ledger, summaries, seed_pairs, independent


def build_pooled_tables(first: pd.DataFrame, second: pd.DataFrame, truth: float):
    pooled = pd.concat(
        [first.assign(batch=0), second.assign(batch=1)],
        ignore_index=True,
    )
    pooled["replicate_global"] = (
        pooled["replicate"] + pooled["batch"] * int(first["replicate"].max() + 1)
    )
    pooled["covered68"] = (
        (pooled["q16"] <= truth) & (truth <= pooled["q84"])
    )
    pooled["width68"] = pooled["q84"] - pooled["q16"]

    cell_rows = []
    for (gamma, validation_size, arm), group in pooled.groupby(
        ["gamma", "validation_size", "arm"], sort=True
    ):
        errors = group["median_error"].to_numpy(dtype=float)
        successes = int(group["covered68"].sum())
        n = int(len(group))
        low, high = wilson_interval(successes, n)
        cell_rows.append({
            "gamma": float(gamma),
            "validation_size": int(validation_size),
            "arm": str(arm),
            "planned_n": n,
            "successful_n": n,
            "bias": float(errors.mean()),
            "rmse": float(np.sqrt(np.mean(errors ** 2))),
            "coverage68": successes / n,
            "coverage_wilson95_low": low,
            "coverage_wilson95_high": high,
            "mean_width68": float(group["width68"].mean()),
            "mean_q025": float(group["q025"].mean()),
            "fraction_q025_below_50p5": float((group["q025"] < 50.5).mean()),
        })
    cells = pd.DataFrame(cell_rows)

    paired_rows = []
    for (gamma, validation_size), group in pooled.groupby(
        ["gamma", "validation_size"], sort=True
    ):
        median = group.pivot(
            index=["batch", "replicate"], columns="arm", values="median_error"
        )
        width = group.pivot(
            index=["batch", "replicate"], columns="arm", values="width68"
        )
        coverage = group.pivot(
            index=["batch", "replicate"], columns="arm", values="covered68"
        )
        for left, right, label in CONTRASTS:
            metrics = {
                "median_error": median[left] - median[right],
                "absolute_error": median[left].abs() - median[right].abs(),
                "squared_error": median[left] ** 2 - median[right] ** 2,
                "coverage68": coverage[left].astype(float) - coverage[right].astype(float),
                "width68": width[left] - width[right],
            }
            for metric, values in metrics.items():
                mean, standard_error, low, high, n = paired_t_interval(
                    values.to_numpy(dtype=float)
                )
                paired_rows.append({
                    "gamma": float(gamma),
                    "validation_size": int(validation_size),
                    "contrast": label,
                    "left_arm": left,
                    "right_arm": right,
                    "metric": metric,
                    "paired_n": n,
                    "mean_difference": mean,
                    "paired_se": standard_error,
                    "t95_low": low,
                    "t95_high": high,
                })
    return cells, pd.DataFrame(paired_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_a", type=Path)
    parser.add_argument("artifact_b", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    cfg_a, ledger_a, summaries_a, seeds_a, audit_a = load_and_validate(args.artifact_a)
    cfg_b, ledger_b, summaries_b, seeds_b, audit_b = load_and_validate(args.artifact_b)

    for key in CONFIG_KEYS:
        if cfg_a[key] != cfg_b[key]:
            raise AssertionError(f"scientific configuration mismatch: {key}")
    if seeds_a & seeds_b:
        raise AssertionError("master-seed streams overlap")

    cells, paired = build_pooled_tables(summaries_a, summaries_b, float(cfg_a["h0_true"]))
    expected_cells = len(cfg_a["gammas"]) * len(cfg_a["validation_sizes"]) * len(ARMS)
    expected_paired = (
        len(cfg_a["gammas"]) * len(cfg_a["validation_sizes"]) * len(CONTRASTS) * 5
    )
    assert len(cells) == expected_cells
    assert len(paired) == expected_paired
    assert set(paired["paired_n"]) == {2 * int(cfg_a["replicates"])}

    args.out.mkdir(parents=True, exist_ok=True)
    cells.to_csv(args.out / "pooled_cells.csv", index=False)
    paired.to_csv(args.out / "pooled_paired_contrasts.csv", index=False)
    report = {
        "status": "PASS",
        "master_seeds": [int(cfg_a["master_seed"]), int(cfg_b["master_seed"])],
        "replicates_per_batch": int(cfg_a["replicates"]),
        "pooled_replicates_per_gamma_validation_cell": 2 * int(cfg_a["replicates"]),
        "unique_seed_pairs": len(seeds_a | seeds_b),
        "seed_overlap": len(seeds_a & seeds_b),
        "planned_cells_per_batch": [len(ledger_a), len(ledger_b)],
        "successful_cells_per_batch": [
            int((ledger_a["status"] == "SUCCESS").sum()),
            int((ledger_b["status"] == "SUCCESS").sum()),
        ],
        "independently_audited_children_per_batch": [
            int(audit_a["successful_child_runs"]),
            int(audit_b["successful_child_runs"]),
        ],
        "configuration_keys_checked": CONFIG_KEYS,
        "pooled_cell_rows": len(cells),
        "paired_contrast_metric_rows": len(paired),
        "contrasts": [label for _, _, label in CONTRASTS],
        "coverage_intervals": "Wilson 95%",
        "paired_intervals": "two-sided Student-t 95%; exploratory and not multiplicity-adjusted",
        "interpretation": (
            "Pooled paired synthetic engineering-pilot readout only; H0 prior [50,90] "
            "and 281-point grid retained; no observational inference or universal estimator ranking."
        ),
    }
    (args.out / "pooled_audit.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
