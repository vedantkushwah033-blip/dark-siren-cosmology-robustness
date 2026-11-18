"""Curvature-matched synthetic H0 propagation smoke.

The validation fitter sees only independently sampled validation redshifts and
retention labels. Parent inclusion probabilities and latent event hosts are
written to truth-only outputs and never enter the likelihood construction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from dataclasses import replace

from dark_siren_robustness.generative_selection import (
    generate_quadratic_logit_inclusion_samples,
)
from dark_siren_robustness.inclusion_posterior import fit_inclusion_posterior_quadrature
from dark_siren_robustness.shared_inclusion_likelihood import (
    shared_inclusion_selected_log_likelihood,
)
from dark_siren_robustness.h0_posterior import normalize_h0_posterior

C_KM_S = 299_792.458


def simulate_detected_events(parent_truth: pd.DataFrame, *, seed: int,
                             n_events: int, h0_true: float = 70.0,
                             sigma_mpc: float = 140.0,
                             threshold_mpc: float = 250.0,
                             z_min: float = 0.01, z_max: float = 0.40):
    """Draw hosts from the full parent population and retain threshold-passing distances."""
    rng = np.random.default_rng(seed)
    z = parent_truth["redshift"].to_numpy(float)
    merger_weight = 1.0 + 1.5 * (z - z_min) / (z_max - z_min)
    host_probability = merger_weight / merger_weight.sum()
    accepted = []
    attempts = 0
    max_attempts = max(100_000, n_events * 10_000)
    while len(accepted) < n_events and attempts < max_attempts:
        batch = max(256, 2 * (n_events - len(accepted)))
        ids = rng.choice(len(z), size=batch, replace=True, p=host_probability)
        distances = rng.normal(C_KM_S * z[ids] / h0_true, sigma_mpc)
        for idx, distance in zip(ids, distances):
            attempts += 1
            if distance >= threshold_mpc:
                accepted.append((int(idx), float(distance)))
                if len(accepted) == n_events:
                    break
            if attempts >= max_attempts:
                break
    if len(accepted) != n_events:
        raise RuntimeError(f"generated only {len(accepted)}/{n_events} detected events")
    truth = pd.DataFrame({
        "event_id": np.arange(n_events, dtype=np.int64),
        "host_object_id": [x[0] for x in accepted],
        "host_redshift": [z[x[0]] for x in accepted],
    })
    observed = pd.DataFrame({
        "event_id": np.arange(n_events, dtype=np.int64),
        "distance_mpc": [x[1] for x in accepted],
        "sigma_mpc": np.full(n_events, sigma_mpc),
    })
    return observed, truth, attempts


def _posterior_curve(z_validation, retained, z_catalogue, events, base_weights,
                     h0_grid, *, degree, order, center, scale, threshold):
    fit = fit_inclusion_posterior_quadrature(
        z_validation, retained, degree=degree, order=order,
        coordinate_center=center, coordinate_scale=scale,
    )
    likelihood = shared_inclusion_selected_log_likelihood(
        h0_grid, z_catalogue, events, fit,
        detection_threshold_mpc=threshold,
        population_base_weights=base_weights,
    )
    posterior = normalize_h0_posterior(
        h0_grid, likelihood.marginal_log_likelihood,
        prior_min=float(h0_grid[0]), prior_max=float(h0_grid[-1]),
    )
    return fit, posterior, likelihood


def run(args):
    data = generate_quadratic_logit_inclusion_samples(
        population_seed=args.seed, validation_seed=args.validation_seed,
        n_parent=args.n_parent, n_validation_frame=args.n_validation_frame,
        n_validation=args.n_validation, gamma=args.gamma, beta=args.beta,
        z_min=args.z_min, z_max=args.z_max,
    )
    z_min, z_max = args.z_min, args.z_max
    events_observed, events_truth, attempts = simulate_detected_events(
        data["parent_truth"], seed=args.event_seed, n_events=args.n_events,
        h0_true=args.h0_true, sigma_mpc=args.sigma_mpc,
        threshold_mpc=args.threshold, z_min=z_min, z_max=z_max,
    )
    catalogue = data["catalogue_observed"]
    validation = data["validation_observed"]
    z_catalogue = catalogue["redshift"].to_numpy(float)
    z_validation = validation["redshift"].to_numpy(float)
    retained = validation["retained"].to_numpy(int)
    base = 1.0 + 1.5 * (z_catalogue - z_min) / (z_max - z_min)
    events = [
        {"gw_distance": float(row.distance_mpc), "gw_sigma": float(row.sigma_mpc)}
        for row in events_observed.itertuples(index=False)
    ]
    h0_grid = np.linspace(args.h0_min, args.h0_max, args.h0_points)
    center, scale = 0.205, 0.195
    results = {}
    curves = {"H0": h0_grid}
    fits = {}
    node_payloads = {}
    for degree in (1, 2):
        previous = None
        max_density_difference = None
        max_quantile_difference = None
        for order in args.orders:
            fit, posterior, likelihood = _posterior_curve(
                z_validation, retained, z_catalogue, events, base, h0_grid,
                degree=degree, order=order, center=center, scale=scale,
                threshold=args.threshold,
            )
            curves[f"density_d{degree}_o{order}"] = np.interp(
                h0_grid, posterior.grid, posterior.density
            )
            if previous is not None:
                max_density_difference = float(np.trapezoid(
                    np.abs(posterior.density - previous.density), h0_grid
                ))
                max_quantile_difference = float(max(
                    abs(posterior.quantile(q) - previous.quantile(q))
                    for q in (0.05, 0.5, 0.95)
                ))
            previous = posterior
            fits[(degree, order)] = fit
            prefix = f"degree{degree}_order{order}"
            node_payloads[f"{prefix}_nodes"] = fit.nodes
            node_payloads[f"{prefix}_weights"] = fit.weights
            node_payloads[f"{prefix}_node_log_selected"] = likelihood.node_log_selected_likelihood
            node_payloads[f"{prefix}_node_log_beta"] = likelihood.node_log_detection_efficiency
            node_payloads[f"{prefix}_host_weights"] = likelihood.node_catalogue_weights
        results[f"degree_{degree}"] = {
            "last_order": args.orders[-1],
            "posterior_mean_coefficients": fits[(degree, args.orders[-1])].posterior_mean.tolist(),
            "quadrature_ess": float(fits[(degree, args.orders[-1])].quadrature_ess),
            "maximum_node_mass": float(fits[(degree, args.orders[-1])].maximum_normalized_weight),
            "minimum_log_inclusion_probability": float(np.min(
                likelihood.node_minimum_log_inclusion_probability
            )),
            "h0_q05": float(previous.quantile(0.05)),
            "h0_median": float(previous.quantile(0.5)),
            "h0_q95": float(previous.quantile(0.95)),
            "successive_order_density_l1": max_density_difference,
            "successive_order_max_quantile_shift": max_quantile_difference,
        }

    # Catalogue-only reference: a constant inclusion function cancels on normalization.
    fitted = fits[(2, args.orders[-1])]
    constant = replace(fitted, nodes=np.zeros((1, 3)), weights=np.ones(1))
    no_correction = shared_inclusion_selected_log_likelihood(
        h0_grid, z_catalogue, events, constant,
        detection_threshold_mpc=args.threshold, population_base_weights=base,
    )
    truth_theta = np.asarray([
        data["config"]["alpha"] - args.gamma * data["config"]["x2_center"],
        args.beta, args.gamma,
    ])
    oracle = replace(fitted, nodes=truth_theta[None, :], weights=np.ones(1))
    oracle_likelihood = shared_inclusion_selected_log_likelihood(
        h0_grid, z_catalogue, events, oracle,
        detection_threshold_mpc=args.threshold, population_base_weights=base,
    )
    node_payloads["catalogue_only_node_log_selected"] = no_correction.node_log_selected_likelihood
    node_payloads["catalogue_only_node_log_beta"] = no_correction.node_log_detection_efficiency
    node_payloads["catalogue_only_host_weights"] = no_correction.node_catalogue_weights
    node_payloads["oracle_node_log_selected"] = oracle_likelihood.node_log_selected_likelihood
    node_payloads["oracle_node_log_beta"] = oracle_likelihood.node_log_detection_efficiency
    node_payloads["oracle_host_weights"] = oracle_likelihood.node_catalogue_weights
    for label, loglike in (
        ("catalogue_only", no_correction.marginal_log_likelihood),
        ("oracle_true_inclusion", oracle_likelihood.marginal_log_likelihood),
    ):
        post = normalize_h0_posterior(
            h0_grid, loglike, prior_min=args.h0_min, prior_max=args.h0_max
        )
        curves[f"density_{label}"] = np.interp(h0_grid, post.grid, post.density)
        results[label] = {
            "h0_q05": float(post.quantile(0.05)),
            "h0_median": float(post.quantile(0.5)),
            "h0_q95": float(post.quantile(0.95)),
        }

    output = Path(args.output_dir)
    (output / "observable").mkdir(parents=True, exist_ok=True)
    (output / "truth").mkdir(parents=True, exist_ok=True)
    catalogue.to_csv(output / "observable/catalogue.csv", index=False, float_format="%.17g")
    validation.to_csv(output / "observable/validation.csv", index=False, float_format="%.17g")
    events_observed.to_csv(output / "observable/events.csv", index=False, float_format="%.17g")
    data["parent_truth"].to_csv(output / "truth/parent_truth.csv", index=False, float_format="%.17g")
    data["validation_frame_truth"].to_csv(output / "truth/validation_frame_truth.csv", index=False, float_format="%.17g")
    events_truth.to_csv(output / "truth/event_hosts.csv", index=False, float_format="%.17g")
    pd.DataFrame(curves).to_csv(output / "h0_posterior_curves.csv", index=False, float_format="%.17g")
    np.savez_compressed(output / "node_level_likelihoods.npz", **node_payloads)
    summary = {
        "status": "synthetic_smoke_not_measurement",
        "data_generating_config": data["config"],
        "truth_h0": args.h0_true,
        "gamma": args.gamma,
        "n_catalogue": int(len(catalogue)),
        "n_validation": int(len(validation)),
        "n_detected_events": int(len(events_observed)),
        "event_generation_attempts": int(attempts),
        "results": results,
        "limitations": [
            "single synthetic realization; not coverage or SBC",
            "degree-one and degree-two posterior comparisons are model-conditional",
            "oracle true-inclusion arm is diagnostic only",
            "common event noise and one-sided threshold are deliberate toy-model restrictions",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    manifest = {}
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            manifest[str(path.relative_to(output))] = {
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=2026111901)
    parser.add_argument("--validation-seed", type=int, default=2026111902)
    parser.add_argument("--event-seed", type=int, default=2026111903)
    parser.add_argument("--n-parent", type=int, default=1200)
    parser.add_argument("--n-validation-frame", type=int, default=1000)
    parser.add_argument("--n-validation", type=int, default=300)
    parser.add_argument("--n-events", type=int, default=4)
    parser.add_argument("--gamma", type=float, default=1.0)
    parser.add_argument("--beta", type=float, default=-3.9)
    parser.add_argument("--h0-true", type=float, default=70.0)
    parser.add_argument("--sigma-mpc", type=float, default=140.0)
    parser.add_argument("--threshold", type=float, default=250.0)
    parser.add_argument("--z-min", type=float, default=0.01)
    parser.add_argument("--z-max", type=float, default=0.40)
    parser.add_argument("--h0-min", type=float, default=55.0)
    parser.add_argument("--h0-max", type=float, default=85.0)
    parser.add_argument("--h0-points", type=int, default=61)
    parser.add_argument("--orders", type=lambda s: tuple(int(x) for x in s.split(",")), default=(5, 7))
    parser.add_argument("--output-dir", default="artifacts/quadratic_curvature_h0_smoke")
    args = parser.parse_args()
    if args.seed == args.validation_seed or args.seed == args.event_seed or args.validation_seed == args.event_seed:
        parser.error("population, validation, and event seeds must be distinct")
    if len(args.orders) < 2 or any(x < 3 for x in args.orders) or any(b <= a for a, b in zip(args.orders, args.orders[1:])):
        parser.error("--orders must contain at least two increasing integers >= 3")
    summary = run(args)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
