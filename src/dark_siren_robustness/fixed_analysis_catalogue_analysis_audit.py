"""Independent audit of fixed-catalogue confirmatory contrast calculations.

This module deliberately does not import the production analysis function. It
reconstructs paired contrasts, intervals, p-values and Holm adjustment from the
raw attempt ledger, then compares them with the reported JSON.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

GAMMAS = (-1.2, 0.0, 1.2)
VALIDATION_SIZES = (100, 300, 1000)
BANDWIDTH = 0.15
ARMS = {
    "catalogue_only": 0,
    "known_inclusion_ipw": 0,
    "parent_sample_kde": 0,
    "estimated_degree_1_ipw": 1,
    "estimated_degree_2_ipw": 2,
}


def _holm(pvalues: np.ndarray) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    result = np.full(values.shape, np.nan)
    valid = np.flatnonzero(np.isfinite(values))
    order = valid[np.argsort(values[valid])]
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, min(1.0, (len(order) - rank) * values[idx]))
        result[idx] = running
    return result


def _close(actual, expected, label: str, *, atol: float = 1e-11) -> None:
    if expected is None:
        if actual is not None:
            raise ValueError(f"{label}: expected null, got {actual!r}")
        return
    if actual is None or not np.isclose(float(actual), float(expected), rtol=1e-10, atol=atol):
        raise ValueError(f"{label}: reported {actual!r}, reconstructed {expected!r}")


def _list_close(actual, expected, label: str) -> None:
    if expected is None:
        if actual is not None:
            raise ValueError(f"{label}: expected null, got {actual!r}")
        return
    if actual is None or len(actual) != len(expected):
        raise ValueError(f"{label}: list length mismatch")
    for index, (left, right) in enumerate(zip(actual, expected)):
        _close(left, right, f"{label}[{index}]")


def _paired_delta(frame: pd.DataFrame, gamma: float, n_validation: int, degree: int) -> pd.Series:
    estimate_arm = f"estimated_degree_{degree}_ipw"
    subset = frame[
        (frame["gamma"] == gamma)
        & (frame["n_validation"] == n_validation)
        & np.isclose(frame["bandwidth_multiplier"], BANDWIDTH)
        & frame["arm"].isin(("catalogue_only", estimate_arm))
    ].copy()
    subset = subset[subset["status"] == "SUCCESS"]
    wide = subset.pivot(index="replicate", columns="arm", values="iae")
    if "catalogue_only" not in wide or estimate_arm not in wide:
        return pd.Series(dtype=float)
    return (wide[estimate_arm] - wide["catalogue_only"]).dropna()


def _expected_failure_counts(frame: pd.DataFrame, gamma: float, arm: str) -> dict:
    cell = frame[
        (frame["gamma"] == gamma)
        & (frame["n_validation"] == 300)
        & np.isclose(frame["bandwidth_multiplier"], BANDWIDTH)
        & (frame["arm"] == arm)
    ]
    return {
        "n_planned": int(len(cell)),
        "n_success": int((cell["status"] == "SUCCESS").sum()),
        "n_failed": int((cell["status"] == "FAILED").sum()),
        "n_non_evaluable": int((cell["status"] == "NON_EVALUABLE").sum()),
        "n_not_run": int((cell["status"] == "NOT_RUN").sum()),
    }


def audit_confirmatory_analysis(
    attempts: pd.DataFrame,
    reported: dict,
    *,
    planned_replicates: int,
    bootstrap_resamples: int = 10000,
    bootstrap_seed: int = 20261110,
) -> dict:
    required = {
        "replicate", "gamma", "target_retention", "n_validation",
        "bandwidth_multiplier", "degree", "arm", "status", "iae",
    }
    missing = required.difference(attempts.columns)
    if missing:
        raise ValueError(f"Raw attempt ledger missing columns: {sorted(missing)}")
    if planned_replicates < 2 or bootstrap_resamples < 1000:
        raise ValueError("planned_replicates>=2 and bootstrap_resamples>=1000 required")

    frame = attempts.copy()
    if frame[list(required)].isna().any().any():
        # IAE may be null only for non-success lifecycle states.
        allowed_null = frame["iae"].isna() & frame["status"].isin(
            ("FAILED", "NON_EVALUABLE", "NOT_RUN")
        )
        other_nulls = frame.drop(columns=["iae"]).isna().any(axis=1)
        if other_nulls.any() or (frame["iae"].isna() & ~allowed_null).any():
            raise ValueError("Unexpected null values in raw attempt ledger")

    if set(frame["status"].unique()) - {"SUCCESS", "FAILED", "NON_EVALUABLE", "NOT_RUN"}:
        raise ValueError("Unknown lifecycle status in raw attempt ledger")
    success = frame["status"] == "SUCCESS"
    if not np.isfinite(frame.loc[success, "iae"].to_numpy(dtype=float)).all():
        raise ValueError("SUCCESS attempts must have finite IAE")
    if frame.loc[~success, "iae"].notna().any():
        raise ValueError("Non-success attempts must not carry a numeric IAE")

    expected_arms = set(ARMS)
    if set(frame["arm"].unique()) != expected_arms:
        raise ValueError("Unexpected estimator-arm set")
    if set(frame["gamma"].unique()) != set(GAMMAS):
        raise ValueError("Unexpected curvature grid")
    if set(frame["n_validation"].unique()) != set(VALIDATION_SIZES):
        raise ValueError("Unexpected validation-size grid")
    if not np.allclose(frame["bandwidth_multiplier"].unique(), [BANDWIDTH]):
        raise ValueError("Unexpected bandwidth grid")
    if set(frame["target_retention"].unique()) != {0.25}:
        raise ValueError("Unexpected target-retention grid")
    if set(frame["replicate"].unique()) != set(range(planned_replicates)):
        raise ValueError("Replicate IDs do not match the planned design")

    expected_rows = planned_replicates * len(GAMMAS) * len(VALIDATION_SIZES) * len(ARMS)
    if len(frame) != expected_rows:
        raise ValueError(f"Expected {expected_rows} attempt rows, found {len(frame)}")
    keys = ["replicate", "gamma", "target_retention", "n_validation",
            "bandwidth_multiplier", "degree", "arm"]
    if frame.duplicated(keys).any():
        raise ValueError("Duplicate attempt keys")
    for arm, degree in ARMS.items():
        if not (frame.loc[frame["arm"] == arm, "degree"] == degree).all():
            raise ValueError(f"Degree metadata inconsistent for {arm}")

    # Reconstruct the primary symmetric-curvature interaction.
    negative = _paired_delta(frame, -1.2, 300, 1)
    positive = _paired_delta(frame, 1.2, 300, 1)
    paired = pd.concat([negative.rename("negative"), positive.rename("positive")], axis=1).dropna()
    interaction = (paired["positive"] - paired["negative"]).to_numpy(dtype=float)
    n = len(interaction)
    primary = reported["primary_interaction"]
    expected_mean = float(interaction.mean()) if n else None
    expected_se = float(interaction.std(ddof=1) / np.sqrt(n)) if n > 1 else None
    _close(primary["mean_interaction"], expected_mean, "primary.mean_interaction")
    _close(primary["monte_carlo_se"], expected_se, "primary.monte_carlo_se")
    if primary["n_replicates_planned"] != planned_replicates:
        raise ValueError("primary.n_replicates_planned mismatch")
    if primary["n_complete_pairs"] != n or primary["n_incomplete_replicates"] != planned_replicates - n:
        raise ValueError("primary paired-replicate counts mismatch")

    expected_counts = {}
    for gamma in (-1.2, 1.2):
        for arm in ("catalogue_only", "estimated_degree_1_ipw"):
            expected_counts[f"gamma={gamma}:{arm}"] = _expected_failure_counts(frame, gamma, arm)
    if primary["failure_aware_cell_counts"] != expected_counts:
        raise ValueError("Primary failure-aware lifecycle counts mismatch")

    expected_boot = expected_t = None
    if n >= 2:
        rng = np.random.default_rng(bootstrap_seed)
        indices = rng.integers(0, n, size=(bootstrap_resamples, n))
        bootstrap = interaction[indices].mean(axis=1)
        expected_boot = [float(np.quantile(bootstrap, 0.025)), float(np.quantile(bootstrap, 0.975))]
        se = float(interaction.std(ddof=1) / np.sqrt(n))
        critical = float(student_t.ppf(0.975, n - 1))
        mean = float(interaction.mean())
        expected_t = [mean - critical * se, mean + critical * se]
    _list_close(primary["bootstrap_95_interval"], expected_boot, "primary.bootstrap_95_interval")
    _list_close(primary["paired_t_95_interval"], expected_t, "primary.paired_t_95_interval")

    # Reconstruct all 18 secondary contrasts and the Holm family.
    secondary_expected = []
    for gamma in GAMMAS:
        for n_validation in VALIDATION_SIZES:
            for degree in (1, 2):
                differences = _paired_delta(frame, gamma, n_validation, degree).to_numpy(dtype=float)
                count = len(differences)
                mean = float(differences.mean()) if count else None
                se = float(differences.std(ddof=1) / np.sqrt(count)) if count > 1 else None
                if count > 1:
                    critical = float(student_t.ppf(0.975, count - 1))
                    lower, upper = mean - critical * se, mean + critical * se
                    pvalue = float(2 * student_t.sf(abs(mean / se), count - 1)) if se > 0 else (
                        0.0 if mean else 1.0
                    )
                else:
                    lower = upper = pvalue = None
                secondary_expected.append({
                    "gamma": float(gamma),
                    "n_validation": int(n_validation),
                    "degree": int(degree),
                    "n_pairs": int(count),
                    "n_replicates_planned": int(planned_replicates),
                    "mean_iae_difference": mean,
                    "monte_carlo_se": se,
                    "paired_t_95_interval": [lower, upper] if lower is not None else None,
                    "p_value": pvalue,
                })
    adjusted = _holm(np.array([
        item["p_value"] if item["p_value"] is not None else np.nan
        for item in secondary_expected
    ]))
    for item, value in zip(secondary_expected, adjusted):
        item["holm_adjusted_p_value"] = float(value) if np.isfinite(value) else None

    reported_secondary = reported.get("secondary_contrasts", [])
    if len(reported_secondary) != len(secondary_expected):
        raise ValueError("Secondary contrast count mismatch")
    for index, (actual, expected) in enumerate(zip(reported_secondary, secondary_expected)):
        for key in ("gamma", "n_validation", "degree", "n_pairs", "n_replicates_planned"):
            if actual.get(key) != expected[key]:
                raise ValueError(f"secondary[{index}].{key} mismatch")
        for key in ("mean_iae_difference", "monte_carlo_se", "p_value", "holm_adjusted_p_value"):
            _close(actual.get(key), expected[key], f"secondary[{index}].{key}")
        _list_close(actual.get("paired_t_95_interval"), expected["paired_t_95_interval"],
                    f"secondary[{index}].paired_t_95_interval")

    return {
        "status": "PASS",
        "audit_scope": "raw attempt design, lifecycle ledger, primary interaction, bootstrap/t intervals, 18 secondary contrasts and Holm adjustment",
        "attempt_count": int(len(frame)),
        "primary_complete_pairs": int(n),
        "secondary_contrasts_reconstructed": int(len(secondary_expected)),
        "successful_attempts": int(success.sum()),
        "failed_attempts": int((frame["status"] == "FAILED").sum()),
        "non_evaluable_attempts": int((frame["status"] == "NON_EVALUABLE").sum()),
        "not_run_attempts": int((frame["status"] == "NOT_RUN").sum()),
        "limitations": [
            "This checks arithmetic and ledger consistency, not the scientific validity of the estimand.",
            "It does not replace the independent raw-artifact auditor.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--attempts-csv", type=Path, required=True)
    parser.add_argument("--analysis-json", type=Path, required=True)
    parser.add_argument("--planned-replicates", type=int, default=100)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    attempts = pd.read_csv(args.attempts_csv)
    reported = json.loads(args.analysis_json.read_text(encoding="utf-8"))
    result = audit_confirmatory_analysis(attempts, reported,
                                         planned_replicates=args.planned_replicates)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
