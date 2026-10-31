"""Independent seed-level provenance audit for the estimated-inclusion pilot.

Reconstructs parent/catalogue/validation data and detected GW event truth from
recorded component seeds. This is separate from the posterior/artifact auditor
and intentionally does not import the production simulator or study code.
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



def _detected_event_distances(
    seed: int,
    n_parent: int,
    parent_redshift: np.ndarray,
    n_detected: int,
    h0_true: float,
    sigma_mpc: float,
    detection_threshold_mpc: float,
    z_min: float,
    z_max: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Replay the production generator's event RNG stream without importing it.

    The production generator consumes one parent-redshift draw per parent
    galaxy, then one internal catalogue-thinning draw per galaxy, before
    sampling hosts and noisy distances. The internal catalogue is not used by
    the estimated-inclusion study, but its RNG consumption is part of the
    event-generation provenance and must be replayed.
    """
    rng = np.random.default_rng(int(seed))
    rng.random(int(n_parent))  # parent inverse-CDF uniforms
    rng.random(int(n_parent))  # generator's own catalogue-thinning uniforms
    host_weight = 1.0 + 1.5 * (parent_redshift - z_min) / (z_max - z_min)
    host_probability = host_weight / host_weight.sum()
    accepted_ids: list[np.ndarray] = []
    accepted_distances: list[np.ndarray] = []
    n_accepted = 0
    attempts = 0
    max_attempts = max(100_000, int(n_detected) * 2_000)
    while n_accepted < n_detected and attempts < max_attempts:
        batch = min(max(1000, 2 * (n_detected - n_accepted)), max_attempts - attempts)
        host_ids = rng.choice(n_parent, size=batch, replace=True, p=host_probability)
        host_z = parent_redshift[host_ids]
        measured = rng.normal(299_792.458 * host_z / h0_true, sigma_mpc)
        detected = measured >= detection_threshold_mpc
        if np.any(detected):
            take = min(n_detected - n_accepted, int(np.sum(detected)))
            accepted_ids.append(host_ids[detected][:take])
            accepted_distances.append(measured[detected][:take])
            n_accepted += take
        attempts += batch
    if n_accepted != n_detected:
        raise AssertionError(
            f"event replay generated {n_accepted}/{n_detected} detections"
        )
    host_ids = np.concatenate(accepted_ids).astype(np.int64)
    distances = np.concatenate(accepted_distances)
    host_redshift = parent_redshift[host_ids]
    latent_mean = 299_792.458 * host_redshift / h0_true
    return host_ids, host_redshift, latent_mean, distances


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
        "seed", "n_parent", "n_detected", "n_replicates_requested",
        "h0_true", "sigma_mpc", "detection_threshold_mpc",
        "true_inclusion_model", "redshift_grid"
    }
    if not required_summary.issubset(summary):
        raise AssertionError("summary lacks required generation metadata")
    required_columns = {
        "replicate", "replicate_seed", "catalogue_seed", "audit_seed",
        "arm", "status",
    }
    if frame.empty or not required_columns.issubset(frame):
        raise AssertionError("replicate table is empty or lacks provenance columns")

    requested_replicates = summary["n_replicates_requested"]
    if (
        isinstance(requested_replicates, bool)
        or not isinstance(requested_replicates, int)
        or requested_replicates < 1
    ):
        raise AssertionError("summary.n_replicates_requested must be a positive integer")
    replicate_values = pd.to_numeric(frame["replicate"], errors="raise").to_numpy(dtype=float)
    if not np.all(np.isfinite(replicate_values)) or not np.all(
        replicate_values == np.floor(replicate_values)
    ):
        raise AssertionError("replicate IDs must be finite integers")
    if frame.duplicated(["replicate", "arm"]).any():
        raise AssertionError("replicate table contains duplicate replicate/arm rows")
    observed_replicates = set(replicate_values.astype(int).tolist())
    expected_replicates = set(range(requested_replicates))
    if observed_replicates != expected_replicates:
        missing = sorted(expected_replicates - observed_replicates)
        unexpected = sorted(observed_replicates - expected_replicates)
        raise AssertionError(
            "replicate IDs do not match requested replicate range "
            f"(missing={missing}, unexpected={unexpected})"
        )

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
    checked_event_rows = 0
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
        if estimated.empty:
            raise AssertionError(f"replicate {replicate}: no estimated-inclusion arms were recorded")

        parent_z = _parent_redshifts(expected_sim, n_parent, z_min, z_max)
        x = (parent_z - midpoint) / half_range
        true_probability = expit(
            intercept + slope * x + quadratic * (x**2 - 1.0 / 3.0)
        )
        expected_host_ids, expected_host_z, expected_latent_mean, expected_distances = _detected_event_distances(
            expected_sim, n_parent, parent_z, int(summary["n_detected"]),
            float(summary["h0_true"]), float(summary["sigma_mpc"]),
            float(summary["detection_threshold_mpc"]), z_min, z_max,
        )
        retained = np.random.default_rng(expected_catalogue).random(n_parent) < true_probability
        audit_order = np.random.default_rng(expected_audit).permutation(n_parent)

        for _, row in estimated.iterrows():
            try:
                inference = json.loads(row["inference_audit_payload_json"])
                scoring = json.loads(row["scoring_audit_payload_json"])
                validation_z = np.asarray(inference["validation_redshift"], dtype=float)
                validation_y = np.asarray(inference["validation_included"], dtype=int)
                stored_parent = np.asarray(scoring["parent_redshift"], dtype=float)
                stored_probability = np.asarray(scoring["true_inclusion_probability"], dtype=float)
                catalogue_z = np.asarray(inference["catalogue_redshift"], dtype=float)
                observed_distances = np.asarray(
                    inference["observed_distances_mpc"], dtype=float
                )
                stored_host_ids = np.asarray(scoring["gw_event_host_object_id"], dtype=np.int64)
                stored_host_z = np.asarray(scoring["gw_event_host_redshift"], dtype=float)
                stored_latent_mean = np.asarray(
                    scoring["gw_event_latent_mean_distance_mpc"], dtype=float
                )
            except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
                raise AssertionError(
                    f"replicate {replicate}, arm {row['arm']}: incomplete provenance payload"
                ) from exc

            _same_array(stored_parent, parent_z, f"replicate {replicate} parent frame")
            _same_array(stored_probability, true_probability, f"replicate {replicate} true inclusion probabilities")
            _same_array(catalogue_z, parent_z[retained], f"replicate {replicate} retained catalogue")
            if not np.array_equal(stored_host_ids, expected_host_ids):
                raise AssertionError(
                    f"replicate {replicate}: detected GW host IDs do not match event RNG replay"
                )
            _same_array(stored_host_z, expected_host_z, f"replicate {replicate} detected GW host redshifts")
            _same_array(stored_latent_mean, expected_latent_mean, f"replicate {replicate} latent mean distances")
            _same_array(
                observed_distances, expected_distances,
                f"replicate {replicate} detected GW distance observations",
            )
            checked_event_rows += 1
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

    if checked_replicates != frame["replicate"].nunique():
        raise AssertionError("not every replicate received a provenance reconstruction")

    return {
        "status": "PASS",
        "replicates_reconstructed": checked_replicates,
        "estimated_arm_rows_reconstructed": checked_estimated_rows,
        "failed_estimated_arm_rows_reconstructed": int((frame.loc[frame["arm"].astype(str).str.startswith("estimated_"), "status"] != "success").sum()),
        "event_observation_rows_reconstructed": checked_event_rows,
        "checks": [
            "component seed schedule",
            "parent redshift frame",
            "true quadratic-logistic inclusion probabilities",
            "catalogue-retention indicators",
            "nested validation redshifts and retention labels",
            "detected GW host IDs, host redshifts, latent distances, observed distances, and event RNG stream",
        ],
        "scope": (
            "Seed-level reconstruction of parent redshifts, catalogue retention, "
            "validation samples, detected GW host identities, host redshifts, "
            "latent mean distances, and observed detected distances for all "
            "estimated-inclusion rows, including failed fits. This validates deterministic synthetic "
            "provenance, not astrophysical realism or real detector selection."
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
