"""Independent seed-level provenance audit for the estimated-inclusion pilot.

Reconstructs the parent redshift frame, catalogue-retention indicators, and
nested validation samples from the recorded component seeds. This is separate
from the posterior/artifact auditor and intentionally does not import the
production simulator or study implementation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit


def _parent_redshifts(
    seed: int, n_parent: int, z_min: float, z_max: float
) -> np.ndarray:
    """Independent inverse-CDF reconstruction of the declared parent frame."""
    rng = np.random.default_rng(seed)
    grid = np.linspace(z_min, z_max, 20001)
    intensity = (
        1.0
        + 2.2 * np.exp(-0.5 * ((grid - 0.105) / 0.025) ** 2)
        + 1.6 * np.exp(-0.5 * ((grid - 0.275) / 0.040) ** 2)
    )
    increments = 0.5 * (intensity[:-1] + intensity[1:]) * np.diff(grid)
    cdf = np.concatenate(([0.0], np.cumsum(increments)))
    cdf /= cdf[-1]
    return np.interp(rng.random(n_parent), cdf, grid)


def _component_seeds(master_seed: int, replicate: int) -> tuple[int, int, int]:
    sequence = np.random.SeedSequence([int(master_seed), int(replicate)])
    children = sequence.spawn(3)
    return tuple(
        int(child.generate_state(1, dtype=np.uint32)[0]) for child in children
    )


def _same_array(actual: object, expected: np.ndarray, label: str) -> None:
    values = np.asarray(actual, dtype=float)
    if values.shape != expected.shape or not np.all(np.isfinite(values)):
        raise AssertionError(f"{label}: shape or finite-value mismatch")
    if not np.allclose(values, expected, rtol=0.0, atol=2e-14):
        difference = float(np.max(np.abs(values - expected))) if values.size else 0.0
        raise AssertionError(f"{label}: reconstruction mismatch (max abs diff={difference:.3g})")


def audit_provenance(
    replicates_path: Path, summary_path: Path
) -> dict[str, object]:
    frame = pd.read_csv(replicates_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    required_summary = {
        "seed", "n_parent", "true_inclusion_model", "redshift_grid"
    }
    if not required_summary.issubset(summary):
        raise AssertionError("summary lacks required generation metadata")
    if frame.empty or not {"replicate", "replicate_seed", "catalogue_seed", "audit_seed", "arm", "status"}.issubset(frame):
        raise AssertionError("replicate table is empty or lacks provenance columns")

    n_parent = int(summary["n_parent"])
    z_min, z_max = map(float, summary["redshift_grid"][:2])
    model = summary["true_inclusion_model"]
    intercept = float(model["intercept"])
    slope = float(model["slope"])
    quadratic = float(model["quadratic"])
    midpoint = (z_min + z_max) / 2.0
    half_range = (z_max - z_min) / 2.0

    checked_replicates = 0
    checked_estimated_rows = 0
    failures_without_payload: list[int] = []
    for replicate, group in frame.groupby("replicate", sort=True):
        replicate = int(replicate)
        expected_sim, expected_catalogue, expected_audit = _component_seeds(
            int(summary["seed"]), replicate
        )
        for column, expected in (
            ("replicate_seed", expected_sim),
            ("catalogue_seed", expected_catalogue),
            ("audit_seed", expected_audit),
        ):
            observed = pd.to_numeric(group[column], errors="raise").to_numpy(dtype=np.int64)
            if np.any(observed != expected):
                raise AssertionError(
                    f"replicate {replicate}: {column} does not match the declared seed schedule"
                )

        estimated = group[group["arm"].astype(str).str.startswith("estimated_")]
        successful = estimated[estimated["status"] == "success"]
        if successful.empty:
            failures_without_payload.append(replicate)
            continue

        parent_z = _parent_redshifts(expected_sim, n_parent, z_min, z_max)
        x = (parent_z - midpoint) / half_range
        true_probability = expit(
            intercept + slope * x + quadratic * (x**2 - 1.0 / 3.0)
        )
        retained = np.random.default_rng(expected_catalogue).random(n_parent) < true_probability
        audit_order = np.random.default_rng(expected_audit).permutation(n_parent)

        for _, row in successful.iterrows():
            try:
                inference = json.loads(row["inference_audit_payload_json"])
                scoring = json.loads(row["scoring_audit_payload_json"])
                validation_z = np.asarray(inference["validation_redshift"], dtype=float)
                validation_y = np.asarray(inference["validation_included"], dtype=int)
                stored_parent = np.asarray(scoring["parent_redshift"], dtype=float)
                stored_probability = np.asarray(scoring["true_inclusion_probability"], dtype=float)
                catalogue_z = np.asarray(inference["catalogue_redshift"], dtype=float)
            except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
                raise AssertionError(
                    f"replicate {replicate}, arm {row['arm']}: incomplete provenance payload"
                ) from exc

            _same_array(stored_parent, parent_z, f"replicate {replicate} parent frame")
            _same_array(stored_probability, true_probability, f"replicate {replicate} true inclusion probabilities")
            _same_array(catalogue_z, parent_z[retained], f"replicate {replicate} retained catalogue")
            m = int(row["validation_size"])
            if validation_z.shape != (m,) or validation_y.shape != (m,):
                raise AssertionError(
                    f"replicate {replicate}, arm {row['arm']}: validation size mismatch"
                )
            _same_array(validation_z, parent_z[audit_order[:m]], f"replicate {replicate} validation redshifts")
            if not np.array_equal(validation_y, retained[audit_order[:m]].astype(int)):
                raise AssertionError(
                    f"replicate {replicate}, arm {row['arm']}: validation labels mismatch"
                )
            checked_estimated_rows += 1
        checked_replicates += 1

    if failures_without_payload:
        raise AssertionError(
            "Cannot establish seed-level provenance for replicates with no successful "
            "estimated arm because failure rows do not persist the raw validation frame: "
            + ", ".join(map(str, failures_without_payload))
        )
    if checked_replicates != frame["replicate"].nunique():
        raise AssertionError("not every replicate received a provenance reconstruction")

    return {
        "status": "PASS",
        "replicates_reconstructed": checked_replicates,
        "estimated_arm_rows_reconstructed": checked_estimated_rows,
        "checks": [
            "component seed schedule",
            "parent redshift frame",
            "true quadratic-logistic inclusion probabilities",
            "catalogue-retention indicators",
            "nested validation redshifts and retention labels",
        ],
        "scope": (
            "Seed-level reconstruction of parent/catalogue/validation inputs only. "
            "GW event generation is not reconstructed by this audit."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Independently reconstruct estimated-inclusion pilot provenance."
    )
    parser.add_argument("--replicates", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = audit_provenance(args.replicates, args.summary)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
