"""Independently reconstruct the independent-validation H0 smoke artifacts."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def audit(artifact_dir: Path) -> dict[str, object]:
    root = artifact_dir.resolve()
    config = json.loads((root / "config.json").read_text(encoding="utf-8"))
    convergence = json.loads((root / "quadrature_convergence.json").read_text(encoding="utf-8"))
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))

    checksum_failures = []
    for relative, expected in manifest.items():
        path = root / relative
        if not path.is_file():
            checksum_failures.append(relative)
            continue
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            checksum_failures.append(relative)
    if checksum_failures:
        raise AssertionError(f"artifact checksum failures: {checksum_failures}")

    analysis_cfg = config["analysis_generator_config"]
    validation_cfg = config["validation_generator_config"]
    if analysis_cfg["seed"] == validation_cfg["seed"]:
        raise AssertionError("analysis and validation seeds must be distinct")
    if not validation_cfg.get("independent_frame", False):
        raise AssertionError("validation frame is not declared independent")

    analysis_selection = analysis_cfg["catalogue_selection"]
    validation_selection = validation_cfg["catalogue_selection"]
    for parameter in ("floor", "midpoint", "width"):
        if not np.isclose(
            float(analysis_selection[parameter]),
            float(validation_selection[parameter]),
            rtol=0.0,
            atol=0.0,
        ):
            raise AssertionError(
                f"analysis and validation catalogue-selection {parameter} differ"
            )

    catalogue = pd.read_csv(root / "observable" / "analysis_catalogue.csv")
    events = pd.read_csv(root / "observable" / "gw_events_observed.csv")
    validation = pd.read_csv(root / "observable" / "validation_observed.csv")
    parent = pd.read_csv(root / "truth" / "analysis_parent_truth.csv")
    event_truth = pd.read_csv(root / "truth" / "gw_event_truth.csv")
    trace = pd.read_csv(root / "truth" / "gw_generation_trace_truth.csv")
    validation_truth = pd.read_csv(root / "truth" / "validation_frame_truth.csv")

    floor = float(analysis_selection["floor"])
    midpoint = float(analysis_selection["midpoint"])
    width = float(analysis_selection["width"])
    if not (0.0 <= floor <= 1.0 and midpoint > 0.0 and width > 0.0):
        raise AssertionError("catalogue-selection parameters are outside their domain")
    expected_parent_probability = floor + (1.0 - floor) / (
        1.0 + np.exp(
            (parent["redshift"].to_numpy(dtype=float) - midpoint) / width
        )
    )
    expected_validation_probability = floor + (1.0 - floor) / (
        1.0 + np.exp(
            (validation_truth["redshift"].to_numpy(dtype=float) - midpoint) / width
        )
    )
    if not np.allclose(
        parent["inclusion_probability"].to_numpy(dtype=float),
        expected_parent_probability,
        rtol=1e-13,
        atol=1e-14,
    ):
        raise AssertionError("analysis parent inclusion probabilities do not match configured DGP")
    if not np.allclose(
        validation_truth["inclusion_probability"].to_numpy(dtype=float),
        expected_validation_probability,
        rtol=1e-13,
        atol=1e-14,
    ):
        raise AssertionError("validation inclusion probabilities do not match configured DGP")

    if list(catalogue.columns) != ["object_id", "redshift"]:
        raise AssertionError("analysis catalogue exposes unexpected columns")
    if list(events.columns) != ["event_id", "distance_mpc", "detected"]:
        raise AssertionError("GW observable table exposes unexpected columns")
    if list(validation.columns) != ["validation_id", "redshift", "retained"]:
        raise AssertionError("validation observable table exposes unexpected columns")
    if "inclusion_probability" in validation.columns or "host_object_id" in events.columns:
        raise AssertionError("truth fields leaked into observable inputs")

    expected_catalogue = parent.loc[parent["retained"] == 1, ["object_id", "redshift"]]
    expected_catalogue = expected_catalogue.sort_values("object_id").reset_index(drop=True)
    observed_catalogue = catalogue.sort_values("object_id").reset_index(drop=True)
    pd.testing.assert_frame_equal(
        observed_catalogue, expected_catalogue, check_dtype=False, atol=1e-14, rtol=1e-14
    )

    frame_ids = validation_truth.loc[
        validation_truth["selected_for_validation"] == 1, "frame_id"
    ].to_numpy(dtype=int)
    if not np.array_equal(np.sort(frame_ids), np.sort(validation["validation_id"].to_numpy(dtype=int))):
        raise AssertionError("validation SRS IDs do not match truth-frame selection flags")
    selected_validation_truth = validation_truth.set_index("frame_id").loc[
        validation["validation_id"].to_numpy(dtype=int)
    ].reset_index()
    if not np.allclose(selected_validation_truth["redshift"], validation["redshift"], rtol=0, atol=1e-14):
        raise AssertionError("validation redshifts do not reconstruct from raw frame")
    if not np.array_equal(selected_validation_truth["retained"].to_numpy(), validation["retained"].to_numpy()):
        raise AssertionError("validation retention labels do not reconstruct from raw frame")

    threshold = float(analysis_cfg["detection_threshold_mpc"])
    if not np.all(np.isfinite(trace["distance_mpc"])):
        raise AssertionError("generation trace has non-finite distances")
    expected_detected = (trace["distance_mpc"].to_numpy() >= threshold).astype(int)
    if not np.array_equal(expected_detected, trace["detected"].to_numpy(dtype=int)):
        raise AssertionError("generation trace detection flags do not match threshold rule")
    if np.any(trace["selected_as_event"].to_numpy(dtype=int) > expected_detected):
        raise AssertionError("a rejected candidate was marked as a selected event")
    if trace["attempt_id"].duplicated().any() or not np.array_equal(
        trace["attempt_id"].to_numpy(dtype=int), np.arange(len(trace), dtype=int)
    ):
        raise AssertionError("attempt IDs are not a contiguous unique sequence")
    if not set(trace["host_object_id"]).issubset(set(parent["object_id"])):
        raise AssertionError("generation trace references unknown parent objects")

    n_attempts = len(trace)
    n_detected = int(expected_detected.sum())
    n_rejected = int((expected_detected == 0).sum())
    n_selected = int(trace["selected_as_event"].sum())
    n_excess = n_detected - n_selected
    reconstructed = {
        "n_generation_attempts": n_attempts,
        "n_detected_candidates": n_detected,
        "n_rejected_candidates": n_rejected,
        "n_excess_detected_candidates_discarded": n_excess,
        "n_detected_events": n_selected,
    }
    for key, value in reconstructed.items():
        if int(analysis_cfg[key]) != value:
            raise AssertionError(f"generator config mismatch for {key}: {analysis_cfg[key]} != {value}")
    if n_attempts != n_detected + n_rejected:
        raise AssertionError("attempt count does not equal detected plus rejected candidates")

    selected_trace = trace.loc[trace["selected_as_event"] == 1].reset_index(drop=True)
    if len(selected_trace) != len(events) or len(event_truth) != len(events):
        raise AssertionError("selected candidate, observed event, and event truth counts differ")
    if not np.array_equal(selected_trace["host_object_id"], event_truth["host_object_id"]):
        raise AssertionError("selected generation hosts do not match event truth table")
    if not np.allclose(selected_trace["distance_mpc"], events["distance_mpc"], rtol=0, atol=1e-12):
        raise AssertionError("selected candidate distances do not match observed GW events")
    if not np.array_equal(events["event_id"].to_numpy(), event_truth["event_id"].to_numpy()):
        raise AssertionError("event IDs do not align between observed and truth tables")
    if not np.all(events["distance_mpc"].to_numpy() >= threshold):
        raise AssertionError("an observed event violates the GW threshold")

    density = pd.read_csv(root / "h0_posterior_density_by_order.csv")
    orders = convergence["orders_evaluated"]
    expected_columns = ["H0"] + [f"density_order_{order}" for order in orders]
    if list(density.columns) != expected_columns:
        raise AssertionError("posterior density columns do not match evaluated quadrature orders")
    h0 = density["H0"].to_numpy(dtype=float)
    if not np.all(np.isfinite(h0)) or np.any(np.diff(h0) <= 0):
        raise AssertionError("H0 grid is invalid")
    for order in orders:
        values = density[f"density_order_{order}"].to_numpy(dtype=float)
        if not np.all(np.isfinite(values)) or np.any(values < 0):
            raise AssertionError(f"invalid H0 posterior density at quadrature order {order}")
        if not np.isclose(np.trapezoid(values, h0), 1.0, rtol=1e-8, atol=1e-8):
            raise AssertionError(f"H0 posterior at order {order} is not normalized")

    # Reconstruct each shared-node selected likelihood from observable inputs.
    # The inclusion quadrature is regenerated with the production fitter; the
    # event-mixture and one-sided beta calculations below are independently
    # implemented here rather than calling the production likelihood function.
    from scipy.special import log_ndtr, logsumexp
    from scipy.stats import norm
    from dark_siren_robustness.dataset_likelihood import C_KM_S
    from dark_siren_robustness.inclusion_posterior import (
        fit_inclusion_posterior_quadrature,
    )

    coefficient_prior = config["inference"].get("coefficient_prior")
    if not isinstance(coefficient_prior, dict) or coefficient_prior.get("family") != "independent_normal":
        raise AssertionError("coefficient prior family is missing or unsupported")
    prior_mean = np.asarray(coefficient_prior.get("mean"), dtype=float)
    prior_sd = np.asarray(coefficient_prior.get("sd"), dtype=float)
    expected_dimension = int(config["inference"]["degree"]) + 1
    if prior_mean.shape != (expected_dimension,) or prior_sd.shape != (expected_dimension,):
        raise AssertionError("coefficient prior vectors do not match fitted model dimension")
    if not np.all(np.isfinite(prior_mean)) or not np.all(np.isfinite(prior_sd)) or np.any(prior_sd <= 0):
        raise AssertionError("coefficient prior parameters are invalid")

    z_catalogue = catalogue["redshift"].to_numpy(dtype=float)
    z_min, z_max = analysis_cfg["redshift_domain"]
    base_weights = 1.0 + 1.5 * (z_catalogue - z_min) / (z_max - z_min)
    h0_grid = density["H0"].to_numpy(dtype=float)
    means = C_KM_S * z_catalogue[None, :] / h0_grid[:, None]
    sigma = float(analysis_cfg["sigma_mpc"])
    node_reconstruction_differences = {}

    for order in orders:
        fit = fit_inclusion_posterior_quadrature(
            validation["redshift"].to_numpy(dtype=float),
            validation["retained"].to_numpy(dtype=int),
            degree=int(config["inference"]["degree"]),
            prior_mean=prior_mean,
            prior_sd=prior_sd,
            order=int(order),
        )
        x = (z_catalogue - fit.center) / fit.scale
        design = np.column_stack(
            [x ** power for power in range(fit.degree + 1)]
        )
        log_inclusion = -np.logaddexp(
            0.0, -(design @ fit.nodes.T)
        )
        log_raw = np.log(base_weights)[:, None] - log_inclusion
        log_host_weights = log_raw - logsumexp(
            log_raw, axis=0, keepdims=True
        )

        node_log_joint = np.zeros(
            (fit.nodes.shape[0], h0_grid.size), dtype=float
        )
        for event in events.itertuples(index=False):
            event_terms = norm.logpdf(
                float(event.distance_mpc), loc=means, scale=sigma
            )
            for node in range(fit.nodes.shape[0]):
                node_log_joint[node] += logsumexp(
                    event_terms + log_host_weights[:, node][None, :],
                    axis=1,
                )

        log_detect_by_host = log_ndtr(
            (means - threshold) / sigma
        )
        node_log_beta = np.empty_like(node_log_joint)
        for node in range(fit.nodes.shape[0]):
            node_log_beta[node] = logsumexp(
                log_detect_by_host + log_host_weights[:, node][None, :],
                axis=1,
            )
        node_log_selected = (
            node_log_joint - len(events) * node_log_beta
        )

        # Verify persisted posterior nodes/weights and nodewise likelihood arrays
        # against this independent reconstruction from the raw validation sample.
        quadrature_table = pd.read_csv(
            root / f"inclusion_quadrature_order_{order}.csv"
        )
        theta_columns = [
            f"theta_{index}" for index in range(fit.nodes.shape[1])
        ]
        archived_nodes = quadrature_table[theta_columns].to_numpy(dtype=float)
        archived_weights = quadrature_table["posterior_weight"].to_numpy(dtype=float)
        if not np.allclose(archived_nodes, fit.nodes, rtol=1e-11, atol=1e-12):
            raise AssertionError(f"archived quadrature nodes differ at order {order}")
        if not np.allclose(archived_weights, fit.weights, rtol=1e-11, atol=1e-14):
            raise AssertionError(f"archived quadrature weights differ at order {order}")

        archived_likelihood = pd.read_csv(
            root / f"selected_likelihood_order_{order}.csv"
        )
        expected_rows = fit.nodes.shape[0] * h0_grid.size
        if len(archived_likelihood) != expected_rows:
            raise AssertionError(f"archived likelihood row count differs at order {order}")
        expected_node_ids = np.repeat(
            np.arange(fit.nodes.shape[0]), h0_grid.size
        )
        expected_h0 = np.tile(h0_grid, fit.nodes.shape[0])
        if not np.array_equal(
            archived_likelihood["node_index"].to_numpy(dtype=int),
            expected_node_ids,
        ) or not np.allclose(
            archived_likelihood["H0"].to_numpy(dtype=float),
            expected_h0, rtol=0, atol=1e-12,
        ):
            raise AssertionError(f"archived likelihood index does not match at order {order}")
        archived_node_selected = archived_likelihood[
            "node_log_selected_likelihood"
        ].to_numpy(dtype=float).reshape(fit.nodes.shape[0], h0_grid.size)
        archived_node_beta = archived_likelihood[
            "node_log_detection_efficiency"
        ].to_numpy(dtype=float).reshape(fit.nodes.shape[0], h0_grid.size)
        if not np.allclose(
            archived_node_selected, node_log_selected,
            rtol=1e-10, atol=1e-10,
        ):
            raise AssertionError(f"archived nodewise selected likelihood differs at order {order}")
        if not np.allclose(
            archived_node_beta, node_log_beta,
            rtol=1e-10, atol=1e-10,
        ):
            raise AssertionError(f"archived nodewise detection efficiency differs at order {order}")

        marginal_log_likelihood = logsumexp(
            node_log_selected + np.log(fit.weights)[:, None],
            axis=0,
        )
        shifted = marginal_log_likelihood - np.max(marginal_log_likelihood)
        reconstructed_unnormalized = np.exp(shifted)
        reconstructed_area = float(
            np.trapezoid(reconstructed_unnormalized, h0_grid)
        )
        if not np.isfinite(reconstructed_area) or reconstructed_area <= 0:
            raise AssertionError(
                f"independent H0 normalization failed at order {order}"
            )
        reconstructed_density = (
            reconstructed_unnormalized / reconstructed_area
        )
        archived_marginal = archived_likelihood[
            "marginal_log_likelihood"
        ].to_numpy(dtype=float).reshape(fit.nodes.shape[0], h0_grid.size)[0]
        if not np.allclose(
            archived_marginal, marginal_log_likelihood,
            rtol=1e-10, atol=1e-10,
        ):
            raise AssertionError(
                f"archived marginal log likelihood differs at order {order}"
            )
        archived_density = density[f"density_order_{order}"].to_numpy(
            dtype=float
        )
        difference = float(
            np.max(np.abs(reconstructed_density - archived_density))
        )
        node_reconstruction_differences[str(order)] = difference
        if difference >= 2e-8:
            raise AssertionError(
                f"independent selected-likelihood reconstruction mismatch "
                f"at order {order}: {difference}"
            )

    # Independently reconstruct the slope-prior sensitivity outputs.
    prior_summary = json.loads(
        (root / "prior_scale_sensitivity.json").read_text(encoding="utf-8")
    )
    prior_density = pd.read_csv(
        root / "prior_scale_sensitivity_density.csv"
    )
    prior_cfg = config["inference"]["coefficient_prior"]
    prior_mean = np.asarray(prior_cfg["mean"], dtype=float)
    prior_sd = np.asarray(prior_cfg["sd"], dtype=float)
    multipliers = [float(value) for value in prior_summary["slope_multipliers"]]
    expected_density_columns = ["H0"] + [
        f"density_slope_multiplier_{value:g}" for value in multipliers
    ]
    if list(prior_density.columns) != expected_density_columns:
        raise AssertionError("prior sensitivity density columns do not match summary")
    prior_h0 = prior_density["H0"].to_numpy(dtype=float)
    if not np.allclose(prior_h0, h0_grid, rtol=0, atol=1e-12):
        raise AssertionError("prior sensitivity H0 grid differs from inference grid")
    if prior_summary["quadrature_order"] not in orders:
        raise AssertionError("prior sensitivity order was not in the convergence run")
    prior_sensitivity_differences = {}
    prior_records = {
        float(record["slope_multiplier"]): record
        for record in prior_summary["records"]
    }
    if set(prior_records) != set(multipliers):
        raise AssertionError("prior sensitivity records do not match multipliers")

    for multiplier in multipliers:
        scaled_sd = prior_sd.copy()
        scaled_sd[1:] *= multiplier
        fit = fit_inclusion_posterior_quadrature(
            validation["redshift"].to_numpy(dtype=float),
            validation["retained"].to_numpy(dtype=int),
            degree=int(config["inference"]["degree"]),
            prior_mean=prior_mean,
            prior_sd=scaled_sd,
            order=int(prior_summary["quadrature_order"]),
        )
        x = (z_catalogue - fit.center) / fit.scale
        design = np.column_stack(
            [x ** power for power in range(fit.degree + 1)]
        )
        log_inclusion = -np.logaddexp(0.0, -(design @ fit.nodes.T))
        log_raw = np.log(base_weights)[:, None] - log_inclusion
        log_host_weights = log_raw - logsumexp(
            log_raw, axis=0, keepdims=True
        )
        node_log_joint = np.zeros(
            (fit.nodes.shape[0], h0_grid.size), dtype=float
        )
        for event in events.itertuples(index=False):
            event_terms = norm.logpdf(
                float(event.distance_mpc), loc=means, scale=sigma
            )
            for node in range(fit.nodes.shape[0]):
                node_log_joint[node] += logsumexp(
                    event_terms + log_host_weights[:, node][None, :],
                    axis=1,
                )
        log_detect_by_host = log_ndtr((means - threshold) / sigma)
        node_log_beta = np.empty_like(node_log_joint)
        for node in range(fit.nodes.shape[0]):
            node_log_beta[node] = logsumexp(
                log_detect_by_host + log_host_weights[:, node][None, :],
                axis=1,
            )
        node_log_selected = (
            node_log_joint - len(events) * node_log_beta
        )
        marginal = logsumexp(
            node_log_selected + np.log(fit.weights)[:, None],
            axis=0,
        )
        shifted = marginal - np.max(marginal)
        unnormalized = np.exp(shifted)
        area = float(np.trapezoid(unnormalized, h0_grid))
        if not np.isfinite(area) or area <= 0:
            raise AssertionError(
                f"prior sensitivity posterior normalization failed for {multiplier}"
            )
        reconstructed_density = unnormalized / area
        archived_density = prior_density[
            f"density_slope_multiplier_{multiplier:g}"
        ].to_numpy(dtype=float)
        difference = float(np.max(np.abs(
            reconstructed_density - archived_density
        )))
        prior_sensitivity_differences[str(multiplier)] = difference
        if difference >= 2e-8:
            raise AssertionError(
                f"prior sensitivity reconstruction mismatch for {multiplier}: {difference}"
            )
        if not np.allclose(
            prior_records[multiplier]["prior_sd"], scaled_sd,
            rtol=0, atol=1e-14,
        ):
            raise AssertionError(
                f"recorded coefficient prior scales mismatch for {multiplier}"
            )

    return {
        "status": "PASS",
        "artifact_dir": str(root),
        "manifest_files_verified": len(manifest),
        "analysis_parent_count": int(len(parent)),
        "analysis_catalogue_count": int(len(catalogue)),
        "validation_frame_count": int(len(validation_truth)),
        "validation_srs_count": int(len(validation)),
        "generation_attempts_reconstructed": n_attempts,
        "detected_candidates_reconstructed": n_detected,
        "rejected_candidates_reconstructed": n_rejected,
        "selected_events_reconstructed": n_selected,
        "quadrature_orders": orders,
        "h0_posterior_normalization_checks": len(orders),
        "independent_selected_likelihood_reconstruction_max_abs_density_differences": node_reconstruction_differences,
        "likelihood_reconstruction": "shared-node event mixtures and one-sided detection efficiency independently recomputed from observable inputs; quadrature regenerated with production fitter",
        "quadrature_status": convergence["status"],
        "prior_scale_sensitivity_reconstruction_max_abs_density_differences": prior_sensitivity_differences,
        "scientific_claim_limit": "Artifact and arithmetic integrity only; not calibration, SBC, coverage, or real-data validation.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    report = audit(args.artifact_dir)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
