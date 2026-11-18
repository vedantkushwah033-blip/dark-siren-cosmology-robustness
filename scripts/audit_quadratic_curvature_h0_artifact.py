"""Independent reconstruction audit for a quadratic-curvature H0 smoke artifact.

This auditor intentionally does not import project likelihood or posterior code.
It rebuilds catalogue weights, event numerators, one-sided selection efficiencies,
shared-node marginalization, and normalized H0 curves from persisted raw arrays.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import log_expit, log_ndtr, logsumexp
from scipy.stats import norm

C_KM_S = 299_792.458
ATOL = 2e-10


def _assert_close(name: str, actual: np.ndarray, expected: np.ndarray, atol: float = ATOL):
    if actual.shape != expected.shape:
        raise AssertionError(f"{name}: shape mismatch {actual.shape} != {expected.shape}")
    if not np.all(np.isfinite(actual)) or not np.all(np.isfinite(expected)):
        raise AssertionError(f"{name}: non-finite values")
    error = float(np.max(np.abs(actual - expected))) if actual.size else 0.0
    if error > atol:
        raise AssertionError(f"{name}: max absolute difference {error:.6g} > {atol}")
    return error


def _weights(z, theta, degree, base, center, scale):
    x = (z - center) / scale
    design = np.column_stack([x ** k for k in range(degree + 1)])
    log_pi = log_expit(design @ theta)
    log_raw = np.log(base) - log_pi
    log_norm = logsumexp(log_raw)
    log_w = log_raw - log_norm
    weights = np.exp(log_w)
    if np.any(weights <= 0) or not np.isclose(weights.sum(), 1.0, atol=1e-12):
        raise AssertionError("reconstructed host weights lost support or normalization")
    return weights, log_pi


def _node_likelihood(z, events, weights, h0, threshold):
    log_w = np.log(weights)
    means = C_KM_S * z[None, :] / h0[:, None]
    sigma = float(events.sigma_mpc.iloc[0])
    if not np.allclose(events.sigma_mpc.to_numpy(float), sigma, rtol=1e-12, atol=0):
        raise AssertionError("event sigma is not common, contrary to toy-model contract")
    event_log_numerator = np.zeros((weights.shape[1], h0.size))
    for row in events.itertuples(index=False):
        terms = norm.logpdf(float(row.distance_mpc), loc=means, scale=sigma)
        for j in range(weights.shape[1]):
            event_log_numerator[j] += logsumexp(terms + log_w[:, j][None, :], axis=1)
    log_detection = log_ndtr((means - threshold) / sigma)
    log_beta = np.column_stack([
        logsumexp(log_detection + log_w[:, j][None, :], axis=1)
        for j in range(weights.shape[1])
    ]).T
    selected = event_log_numerator - len(events) * log_beta
    return event_log_numerator, log_beta, selected


def _normalize(loglike, h0):
    shifted = np.exp(loglike - np.max(loglike))
    area = float(np.trapezoid(shifted, h0))
    if not np.isfinite(area) or area <= 0:
        raise AssertionError("invalid posterior normalization integral")
    return shifted / area


def audit(root: Path):
    root = root.resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    expected_paths = set(manifest)
    actual_paths = {
        str(p.relative_to(root)) for p in root.rglob("*")
        if p.is_file() and p.name != "manifest.json"
    }
    if expected_paths != actual_paths:
        raise AssertionError(
            f"manifest path mismatch: missing={sorted(expected_paths-actual_paths)}, "
            f"unlisted={sorted(actual_paths-expected_paths)}"
        )
    for rel, meta in manifest.items():
        raw = (root / rel).read_bytes()
        digest = hashlib.sha256(raw).hexdigest()
        if len(raw) != int(meta["bytes"]) or digest != meta["sha256"]:
            raise AssertionError(f"manifest integrity failure: {rel}")

    summary = json.loads((root / "summary.json").read_text())
    cfg = summary["data_generating_config"]
    catalogue = pd.read_csv(root / "observable/catalogue.csv")
    validation = pd.read_csv(root / "observable/validation.csv")
    events = pd.read_csv(root / "observable/events.csv")
    parent = pd.read_csv(root / "truth/parent_truth.csv")
    validation_frame = pd.read_csv(root / "truth/validation_frame_truth.csv")
    event_hosts = pd.read_csv(root / "truth/event_hosts.csv")
    curves = pd.read_csv(root / "h0_posterior_curves.csv")
    arrays = np.load(root / "node_level_likelihoods.npz")

    if summary["status"] != "synthetic_smoke_not_measurement":
        raise AssertionError("unexpected artifact status")
    if len(catalogue) != summary["n_catalogue"] or len(validation) != summary["n_validation"]:
        raise AssertionError("observable row counts disagree with summary")
    if len(events) != summary["n_detected_events"] or len(event_hosts) != len(events):
        raise AssertionError("event row counts disagree with summary")
    if int(cfg["n_parent"]) != len(parent) or int(cfg["n_validation_frame"]) != len(validation_frame):
        raise AssertionError("truth frame row counts disagree with DGP configuration")

    # Reconstruct the declared quadratic inclusion law from truth-only DGP parameters.
    zmin, zmax = float(parent.redshift.min()), float(parent.redshift.max())
    # The generator's requested support is recorded by the actual parent range;
    # the DGP x coordinate is fixed by the declared coordinate center/scale.
    x_parent = (parent.redshift.to_numpy(float) - float(cfg["coordinate_center"])) / float(cfg["coordinate_scale"])
    eta_parent = (
        float(cfg["alpha"]) + float(cfg["beta"]) * x_parent
        + float(cfg["gamma"]) * (x_parent * x_parent - float(cfg["x2_center"]))
    )
    _assert_close(
        "parent inclusion probabilities",
        parent.inclusion_probability.to_numpy(float),
        1.0 / (1.0 + np.exp(-eta_parent)),
        atol=2e-12,
    )
    retained_parent = parent.loc[parent.retained == 1, ["object_id", "redshift"]].reset_index(drop=True)
    _assert_close("catalogue object IDs", catalogue.object_id.to_numpy(float), retained_parent.object_id.to_numpy(float), atol=0)
    _assert_close("catalogue redshifts", catalogue.redshift.to_numpy(float), retained_parent.redshift.to_numpy(float), atol=2e-15)

    selected_validation = validation_frame.loc[
        validation_frame.selected_for_validation == 1,
        ["frame_id", "redshift", "retained"],
    ].reset_index(drop=True)
    _assert_close("validation IDs", validation.validation_id.to_numpy(float), selected_validation.frame_id.to_numpy(float), atol=0)
    _assert_close("validation redshifts", validation.redshift.to_numpy(float), selected_validation.redshift.to_numpy(float), atol=2e-15)
    _assert_close("validation labels", validation.retained.to_numpy(float), selected_validation.retained.to_numpy(float), atol=0)

    if np.any(events.distance_mpc.to_numpy(float) < float(summary.get("threshold_mpc", 250.0))):
        # The threshold is stored in run metadata below for all new artifacts.
        threshold = float(summary["threshold_mpc"])
        if np.any(events.distance_mpc.to_numpy(float) < threshold):
            raise AssertionError("an observed event fails the declared detection threshold")
    threshold = float(summary["threshold_mpc"])
    if not np.all(events.distance_mpc.to_numpy(float) >= threshold):
        raise AssertionError("detected event violates one-sided threshold")
    host_ids = event_hosts.host_object_id.to_numpy(int)
    if np.any(host_ids < 0) or np.any(host_ids >= len(parent)):
        raise AssertionError("event host ID outside parent population")
    _assert_close(
        "event host redshifts",
        event_hosts.host_redshift.to_numpy(float),
        parent.redshift.to_numpy(float)[host_ids],
        atol=2e-15,
    )
    _assert_close("event IDs", events.event_id.to_numpy(float), event_hosts.event_id.to_numpy(float), atol=0)

    # Replay the event-generation RNG independently from the full parent population.
    rng = np.random.default_rng(int(summary["event_seed"]))
    parent_z = parent.redshift.to_numpy(float)
    merger = 1.0 + 1.5 * (parent_z - float(summary["z_min"])) / (
        float(summary["z_max"]) - float(summary["z_min"])
    )
    host_probability = merger / merger.sum()
    replayed = []
    replay_attempts = 0
    target_events = len(events)
    max_attempts = max(100_000, target_events * 10_000)
    while len(replayed) < target_events and replay_attempts < max_attempts:
        batch = max(256, 2 * (target_events - len(replayed)))
        indices = rng.choice(len(parent_z), size=batch, replace=True, p=host_probability)
        distances = rng.normal(
            C_KM_S * parent_z[indices] / float(summary["truth_h0"]),
            float(events.sigma_mpc.iloc[0]),
        )
        for index, distance in zip(indices, distances):
            replay_attempts += 1
            if distance >= threshold:
                replayed.append((int(index), float(distance)))
                if len(replayed) == target_events:
                    break
            if replay_attempts >= max_attempts:
                break
    if replay_attempts != int(summary["event_generation_attempts"]):
        raise AssertionError("event-generation attempt count does not replay")
    _assert_close(
        "replayed event host IDs",
        event_hosts.host_object_id.to_numpy(float),
        np.asarray([item[0] for item in replayed], dtype=float),
        atol=0,
    )
    _assert_close(
        "replayed observed distances",
        events.distance_mpc.to_numpy(float),
        np.asarray([item[1] for item in replayed], dtype=float),
        atol=2e-12,
    )

    z = catalogue.redshift.to_numpy(float)
    base = 1.0 + 1.5 * (z - threshold * 0 + float(summary["z_min"])) / (
        float(summary["z_max"]) - float(summary["z_min"])
    )
    h0 = curves.H0.to_numpy(float)
    if np.any(np.diff(h0) <= 0):
        raise AssertionError("H0 grid is not strictly increasing")
    errors = {}
    reconstructed_curves = {}

    for degree in (1, 2):
        for order in summary["orders"]:
            prefix = f"degree{degree}_order{order}"
            nodes = arrays[f"{prefix}_nodes"]
            node_weights = arrays[f"{prefix}_weights"]
            saved_host_weights = arrays[f"{prefix}_host_weights"]
            saved_beta = arrays[f"{prefix}_node_log_beta"]
            saved_selected = arrays[f"{prefix}_node_log_selected"]
            if nodes.shape[1] != degree + 1 or node_weights.shape != (nodes.shape[0],):
                raise AssertionError(f"{prefix}: node dimensions inconsistent")
            if not np.isclose(node_weights.sum(), 1.0, atol=1e-12) or np.any(node_weights < 0):
                raise AssertionError(f"{prefix}: invalid posterior node weights")
            host_weights = np.column_stack([
                _weights(z, theta, degree, base, float(cfg["coordinate_center"]),
                         float(cfg["coordinate_scale"]))[0]
                for theta in nodes
            ])
            errors[f"{prefix}_host_weights"] = _assert_close(
                f"{prefix} host weights", saved_host_weights, host_weights, atol=2e-12
            )
            numerator, beta, selected = _node_likelihood(
                z, events, host_weights, h0, threshold
            )
            errors[f"{prefix}_beta"] = _assert_close(
                f"{prefix} selection normalization", saved_beta, beta, atol=2e-10
            )
            errors[f"{prefix}_selected"] = _assert_close(
                f"{prefix} selected node likelihood", saved_selected, selected, atol=2e-10
            )
            marginal = logsumexp(selected + np.log(node_weights)[:, None], axis=0)
            density = _normalize(marginal, h0)
            curve_name = f"density_d{degree}_o{order}"
            errors[f"{prefix}_posterior_curve"] = _assert_close(
                f"{prefix} normalized H0 curve", curves[curve_name].to_numpy(float),
                density, atol=2e-10
            )
            reconstructed_curves[curve_name] = density

    # Independently reconstruct catalogue-only and oracle arms.
    constant_theta = np.zeros((1, 3))
    alpha = float(cfg["alpha"])
    beta_coef = float(cfg["beta"])
    gamma = float(cfg["gamma"])
    c_g = float(cfg["x2_center"])
    oracle_theta = np.asarray([[alpha - gamma * c_g, beta_coef, gamma]])
    for label, theta in (("catalogue_only", constant_theta), ("oracle", oracle_theta)):
        weights, _ = _weights(z, theta[0], 2, base, float(cfg["coordinate_center"]),
                              float(cfg["coordinate_scale"]))
        numerator, log_beta, selected = _node_likelihood(
            z, events, weights[:, None], h0, threshold
        )
        density = _normalize(selected[0], h0)
        curve_name = "density_catalogue_only" if label == "catalogue_only" else "density_oracle_true_inclusion"
        errors[f"{label}_posterior_curve"] = _assert_close(
            f"{label} normalized H0 curve", curves[curve_name].to_numpy(float),
            density, atol=2e-10
        )
        reconstructed_curves[curve_name] = density

    report = {
        "status": "PASS",
        "files_manifest_verified": len(manifest),
        "parent_rows": len(parent),
        "catalogue_rows": len(catalogue),
        "validation_rows": len(validation),
        "detected_events": len(events),
        "coefficient_nodes_reconstructed": int(sum(arrays[f"degree{d}_order{o}_nodes"].shape[0] for d in (1, 2) for o in summary["orders"])),
        "max_abs_reconstruction_error": max(errors.values()) if errors else 0.0,
        "checks": list(errors),
        "limitations": [
            "auditor reconstructs likelihood and H0 normalization from persisted coefficient nodes; it does not independently refit the validation posterior",
            "single synthetic realization; no SBC, coverage, or estimator superiority claim",
        ],
    }
    (root / "independent_audit.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    audit(parser.parse_args().artifact_dir)


if __name__ == "__main__":
    main()
