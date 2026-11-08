"""Independently reconstruct an archived inclusion-posterior smoke artifact.

The auditor reads observable inputs to reconstruct the validation posterior,
node-specific host weights, selected multi-event likelihood and H0 posterior.
Truth-only files are consulted only for provenance checks after reconstruction.
This remains a smoke-artifact audit, not a calibration study.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit, log_ndtr, logsumexp
from scipy.stats import norm

from dark_siren_robustness.dataset_likelihood import C_KM_S
from dark_siren_robustness.inclusion_posterior import (
    fit_inclusion_posterior_quadrature,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def audit(output_dir: Path) -> dict[str, object]:
    manifest_path = output_dir / "sha256_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    hash_checks = {}
    for relative, expected in manifest.items():
        path = output_dir / relative
        _require(path.is_file(), f"manifested file missing: {relative}")
        actual = _sha256(path)
        hash_checks[relative] = actual == expected
        _require(actual == expected, f"checksum mismatch: {relative}")

    summary = json.loads(
        (output_dir / "analysis" / "summary.json").read_text(encoding="utf-8")
    )
    config = summary["analysis_config"]
    validation_config = summary["validation_config"]
    _require(
        summary["analysis_seed"] != summary["validation_seed"],
        "analysis and validation seeds are not distinct",
    )
    _require(
        validation_config["seed"] == summary["validation_seed"],
        "validation seed does not match validation configuration",
    )
    _require(
        config["seed"] == summary["analysis_seed"],
        "analysis seed does not match analysis configuration",
    )

    # Reconstruction inputs: observable files only.
    catalogue = pd.read_csv(output_dir / "observable" / "catalogue_observed.csv")
    events_table = pd.read_csv(output_dir / "observable" / "gw_events_observed.csv")
    validation = pd.read_csv(output_dir / "observable" / "validation_observed.csv")
    _require(len(catalogue) >= 2, "observable catalogue has fewer than two hosts")
    _require(len(events_table) == config["n_detected_events"], "event count mismatch")
    _require(len(validation) == validation_config["n_validation"], "validation count mismatch")
    _require(
        set(validation.columns) == {"validation_id", "redshift", "retained"},
        "unexpected validation observable schema",
    )
    _require(
        set(events_table.columns) == {"event_id", "distance_mpc", "detected"},
        "unexpected GW observable schema",
    )
    _require(
        np.all(events_table["distance_mpc"].to_numpy() >= config["detection_threshold_mpc"]),
        "an observed event violates the one-sided selection threshold",
    )

    fit = fit_inclusion_posterior_quadrature(
        validation["redshift"].to_numpy(dtype=float),
        validation["retained"].to_numpy(dtype=int),
        degree=1,
        order=9,
    )
    z = catalogue["redshift"].to_numpy(dtype=float)
    z_min, z_max = config["redshift_domain"]
    base = 1.0 + 1.5 * (z - z_min) / (z_max - z_min)
    x = (z - fit.center) / fit.scale
    design = np.column_stack([x ** power for power in range(fit.degree + 1)])
    log_probability = -np.logaddexp(0.0, -(design @ fit.nodes.T))
    log_raw = np.log(base)[:, None] - log_probability
    log_host_weights = log_raw - logsumexp(log_raw, axis=0, keepdims=True)

    h0 = np.linspace(50.0, 90.0, 81)
    means = C_KM_S * z[None, :] / h0[:, None]
    sigma = float(config["sigma_mpc"])
    threshold = float(config["detection_threshold_mpc"])
    node_log_joint = np.zeros((fit.nodes.shape[0], h0.size), dtype=float)
    for event in events_table.itertuples(index=False):
        event_terms = norm.logpdf(
            float(event.distance_mpc), loc=means, scale=sigma
        )
        for node in range(fit.nodes.shape[0]):
            node_log_joint[node] += logsumexp(
                event_terms + log_host_weights[:, node][None, :], axis=1
            )

    log_detection_by_host = log_ndtr((means - threshold) / sigma)
    node_log_beta = np.empty_like(node_log_joint)
    for node in range(fit.nodes.shape[0]):
        node_log_beta[node] = logsumexp(
            log_detection_by_host + log_host_weights[:, node][None, :],
            axis=1,
        )
    node_log_selected = (
        node_log_joint - len(events_table) * node_log_beta
    )
    marginal_log_likelihood = logsumexp(
        node_log_selected + np.log(fit.weights)[:, None], axis=0
    )

    # Independently normalize the uniform-prior H0 posterior on [50, 90].
    shifted = marginal_log_likelihood - np.max(marginal_log_likelihood)
    unnormalized = np.exp(shifted)
    area = float(np.trapezoid(unnormalized, h0))
    _require(np.isfinite(area) and area > 0, "independent H0 normalization failed")
    density = unnormalized / area
    cdf = np.concatenate((
        [0.0], np.cumulative_trapezoid(density, h0)
    ))
    cdf /= cdf[-1]

    archived = pd.read_csv(output_dir / "analysis" / "h0_posterior.csv")
    archived_grid = archived["H0_km_s_Mpc"].to_numpy(dtype=float)
    archived_density = archived["posterior_density"].to_numpy(dtype=float)
    archived_cdf = archived["posterior_cdf"].to_numpy(dtype=float)
    _require(
        np.array_equal(archived_grid, h0),
        "archived H0 grid differs from independent audit grid",
    )
    max_density_difference = float(np.max(np.abs(archived_density - density)))
    max_cdf_difference = float(np.max(np.abs(archived_cdf - cdf)))
    _require(
        max_density_difference < 2e-8,
        f"posterior density reconstruction mismatch: {max_density_difference}",
    )
    _require(
        max_cdf_difference < 2e-8,
        f"posterior CDF reconstruction mismatch: {max_cdf_difference}",
    )

    # Truth-only files are opened only after the observable-only inference
    # reconstruction above. They are used for provenance/scoring checks only.
    parent_truth = pd.read_csv(output_dir / "truth" / "analysis_parent_truth.csv")
    event_truth = pd.read_csv(output_dir / "truth" / "gw_event_truth.csv")
    validation_truth = pd.read_csv(
        output_dir / "truth" / "validation_frame_truth.csv"
    )
    _require(
        len(parent_truth) == config["n_parent"],
        "analysis parent truth size mismatch",
    )
    _require(
        event_truth["host_object_id"].isin(parent_truth["object_id"]).all(),
        "event truth contains unknown parent host IDs",
    )
    selected_ids = validation["validation_id"].to_numpy(dtype=int)
    truth_selected = validation_truth.set_index("frame_id").loc[selected_ids]
    _require(
        np.allclose(
            validation["redshift"].to_numpy(),
            truth_selected["redshift"].to_numpy(),
            rtol=0.0,
            atol=1e-14,
        ),
        "validation observable redshifts do not match selected truth frame",
    )
    _require(
        np.array_equal(
            validation["retained"].to_numpy(dtype=int),
            truth_selected["retained"].to_numpy(dtype=int),
        ),
        "validation observable labels do not match selected truth frame",
    )
    _require(
        int(validation_truth["selected_for_validation"].sum()) == len(validation),
        "validation truth selection indicator count mismatch",
    )

    return {
        "status": "PASS",
        "audit_type": "independent_observable_input_likelihood_and_H0_reconstruction",
        "checksums_verified": len(hash_checks),
        "observable_only_reconstruction": True,
        "truth_used_only_after_reconstruction": True,
        "quadrature_order": 9,
        "quadrature_nodes": int(fit.nodes.shape[0]),
        "max_abs_posterior_density_difference": max_density_difference,
        "max_abs_posterior_cdf_difference": max_cdf_difference,
        "validation_rows_crosschecked_against_truth": int(len(validation)),
        "event_host_ids_crosschecked": int(len(event_truth)),
        "limitations": [
            "one synthetic replicate only",
            "posterior quadrature nodes are regenerated with the production fitter",
            "likelihood and beta equations are independently reconstructed here",
            "no SBC, calibration, or frequentist coverage claim",
            "no real survey or detector validation",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    result = audit(args.output_dir)
    destination = args.output_dir / "analysis" / "independent_audit.json"
    destination.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
