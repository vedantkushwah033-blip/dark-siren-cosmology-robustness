"""Paired degree-one/degree-two five-arm inclusion quadrature sensitivity pilot.

Uses identical synthetic inputs across inclusion-model degrees and quadrature
orders, independently audits every child artifact, and measures convergence of
the integrated-IPW H0 posterior. This is a numerical engineering pilot, not
confirmatory calibration or an observational analysis.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _posterior(grid: np.ndarray, loglike: np.ndarray) -> tuple[np.ndarray, dict[str, float]]:
    shifted = np.asarray(loglike, float) - float(np.max(loglike))
    density = np.exp(shifted)
    area = float(trapezoid(density, grid))
    if not np.isfinite(area) or area <= 0:
        raise AssertionError("posterior normalization failed")
    density /= area
    cdf = np.concatenate(([0.0], cumulative_trapezoid(density, grid)))
    cdf /= cdf[-1]
    quantiles = {
        f"q{int(round(p * 1000)):03d}": float(np.interp(p, cdf, grid))
        for p in (0.025, 0.16, 0.50, 0.84, 0.975)
    }
    return density, quantiles


def _assert_expected_config(config: dict, expected: dict, label: str) -> None:
    """Fail closed if a paired child changes any frozen simulation setting."""
    for key, value in expected.items():
        if config.get(key) != value:
            raise AssertionError(
                f"{label} config mismatch for {key}: {config.get(key)!r} != {value!r}"
            )


def _compare_inputs(reference: Path, current: Path) -> int:
    pairs = [
        ("observable/catalogue_observed.csv", "catalogue"),
        ("observable/gw_events_observed.csv", "events"),
        ("truth/parent_truth.csv", "parent"),
        ("truth/gw_event_truth.csv", "event truth"),
        ("truth/validation_frame_truth.csv", "validation frame"),
    ]
    checked = 0
    for rel, label in pairs:
        left = pd.read_csv(reference / rel)
        right = pd.read_csv(current / rel)
        try:
            pd.testing.assert_frame_equal(left, right, check_dtype=False, check_exact=True)
        except AssertionError as exc:
            raise AssertionError(f"degree/order pairing mismatch in {label}: {exc}") from exc
        checked += 1
    for file in reference.glob("observable/validation_n*.csv"):
        rel = file.relative_to(reference)
        left = pd.read_csv(file)
        right = pd.read_csv(current / rel)
        pd.testing.assert_frame_equal(left, right, check_dtype=False, check_exact=True)
        checked += 1
    return checked


def run(args: argparse.Namespace) -> dict:
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError("output directory must be empty")
    if args.orders != [5, 7, 9]:
        raise ValueError("the preregistered pilot requires quadrature orders exactly 5,7,9")
    if (not args.validation_sizes or args.validation_sizes[0] < 5
        or sorted(set(args.validation_sizes)) != args.validation_sizes
        or args.validation_sizes[-1] > args.validation_frame):
        raise ValueError("validation sizes must be unique, increasing, >=5, and within frame")
    if args.n_detected < 1 or args.n_parent < 100:
        raise ValueError("invalid parent/event dimensions")
    if (not np.isfinite(args.gamma) or not np.isfinite(args.curvature_beta)
        or not np.isfinite(args.h0_true) or not np.isfinite(args.h0_prior_min)
        or not np.isfinite(args.h0_prior_max) or args.h0_prior_min >= args.h0_prior_max
        or not args.h0_prior_min <= args.h0_true <= args.h0_prior_max):
        raise ValueError("invalid curvature or H0 prior/truth configuration")

    runner = Path(__file__).with_name("run_quadratic_five_arm_smoke.py")
    auditor = Path(__file__).with_name("audit_quadratic_five_arm_smoke.py")
    records = []
    diagnostic_parts = []
    reference = None
    input_checks = 0
    posterior_store: dict[tuple[int, int, str], tuple[np.ndarray, dict[str, float]]] = {}
    likelihood_store: dict[tuple[int, int, str], np.ndarray] = {}
    h0_grid = None

    for degree in (1, 2):
        for order in args.orders:
            label = f"degree{degree}_order{order}"
            child = out / label
            command = [
                sys.executable, str(runner),
                "--seed", str(args.seed),
                "--validation-seed", str(args.validation_seed),
                "--n-parent", str(args.n_parent),
                "--n-detected", str(args.n_detected),
                "--validation-frame", str(args.validation_frame),
                "--validation-sizes", ",".join(map(str, args.validation_sizes)),
                "--curvature-gamma", str(args.gamma),
                "--curvature-beta", str(args.curvature_beta),
                "--inclusion-degree", str(degree),
                "--h0-true", str(args.h0_true),
                "--h0-prior-min", str(args.h0_prior_min),
                "--h0-prior-max", str(args.h0_prior_max),
                "--h0-points", str(args.h0_points),
                "--scott-multiplier", str(args.scott_multiplier),
                "--quadrature-order", str(order),
                "--sigma-mpc", str(args.sigma_mpc),
                "--detection-threshold-mpc", str(args.detection_threshold_mpc),
                "--output-dir", str(child),
            ]
            proc = subprocess.run(command, capture_output=True, text=True)
            (out / f"{label}.stdout.txt").write_text(proc.stdout)
            (out / f"{label}.stderr.txt").write_text(proc.stderr)
            if proc.returncode:
                raise RuntimeError(f"{label} runner failed: {proc.stderr[-3000:]}")
            audit_path = child / "audit_report.json"
            audited = subprocess.run(
                [sys.executable, str(auditor), str(child), "--output", str(audit_path)],
                capture_output=True, text=True,
            )
            (out / f"{label}.audit.stdout.txt").write_text(audited.stdout)
            (out / f"{label}.audit.stderr.txt").write_text(audited.stderr)
            if audited.returncode:
                raise RuntimeError(f"{label} independent audit failed: {audited.stderr[-3000:]}")
            audit_report = json.loads(audit_path.read_text())
            if audit_report.get("status") != "PASS":
                raise AssertionError(f"{label} independent audit did not PASS")

            config = json.loads((child / "config.json").read_text())
            expected_config = {
                "seed": args.seed,
                "validation_seed": args.validation_seed,
                "n_parent": args.n_parent,
                "n_detected": args.n_detected,
                "validation_frame": args.validation_frame,
                "validation_sizes": args.validation_sizes,
                "catalogue_floor": 0.10,
                "catalogue_model": "quadratic_logit",
                "inclusion_model_degree": degree,
                "curvature_gamma": args.gamma,
                "curvature_beta": args.curvature_beta,
                "h0_true": args.h0_true,
                "h0_prior_support": [args.h0_prior_min, args.h0_prior_max],
                "h0_grid_points": args.h0_points,
                "scott_multiplier": args.scott_multiplier,
                "quadrature_order": order,
                "sigma_mpc": args.sigma_mpc,
                "detection_threshold_mpc": args.detection_threshold_mpc,
                "truth_used_for_non_oracle_inference": False,
                "truth_used_for_arm_c_only": True,
            }
            _assert_expected_config(config, expected_config, label)
            if reference is None:
                reference = child
            else:
                input_checks += _compare_inputs(reference, child)
                reference_config = json.loads((reference / "config.json").read_text())
                for key in (
                    "n_catalogue_observed", "quadratic_calibration",
                    "redshift_grid_points", "bandwidth_selection_status",
                ):
                    if config.get(key) != reference_config.get(key):
                        raise AssertionError(f"{label} paired configuration mismatch for {key}")

            child_diagnostics = pd.read_csv(child / "results" / "weight_diagnostics.csv")
            child_diagnostics.insert(0, "quadrature_order", order)
            child_diagnostics.insert(0, "degree", degree)
            diagnostic_parts.append(child_diagnostics)

            with np.load(child / "results" / "log_likelihoods.npz") as archive:
                grid = np.asarray(archive["H0_grid"], float)
                if h0_grid is None:
                    h0_grid = grid
                elif not np.array_equal(h0_grid, grid):
                    raise AssertionError("H0 grid differs across paired runs")
                for nval in args.validation_sizes:
                    for arm in (
                        "A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
                        "D_map_plugin_ipw", "E_uncertainty_integrated_ipw",
                    ):
                        curve = np.asarray(archive[f"{arm}_n{nval}"], float)
                        likelihood_store[(degree, order, f"{arm}_n{nval}")] = curve
                        density, quantiles = _posterior(grid, curve)
                        posterior_store[(degree, order, f"{arm}_n{nval}")] = (density, quantiles)
                        records.append({
                            "degree": degree, "quadrature_order": order,
                            "validation_size": nval, "arm": arm,
                            **quantiles,
                            "quadrature_nodes": int(archive[f"posterior_nodes_n{nval}"].shape[0]),
                            "posterior_quadrature_ess": float(
                                1.0 / np.sum(archive[f"posterior_weights_n{nval}"] ** 2)
                            ),
                            "child_audit_status": "PASS",
                        })

    invariant_rows = []
    order_invariant_arms = (
        "A_catalogue_only", "B_parent_oracle",
        "C_known_probability_ipw", "D_map_plugin_ipw",
    )
    for degree in (1, 2):
        for nval in args.validation_sizes:
            for arm in order_invariant_arms:
                key = f"{arm}_n{nval}"
                reference_curve = likelihood_store[(degree, 5, key)]
                for order in (7, 9):
                    candidate = likelihood_store[(degree, order, key)]
                    difference = float(np.max(np.abs(candidate-reference_curve)))
                    if not np.allclose(candidate, reference_curve, rtol=1e-11, atol=1e-11):
                        raise AssertionError(
                            f"{key} should be quadrature-order invariant for degree {degree}; "
                            f"max absolute difference={difference}"
                        )
                    invariant_rows.append({
                        "invariant": "fixed_degree_across_quadrature_order",
                        "degree": degree, "quadrature_order": order,
                        "validation_size": nval, "arm": arm,
                        "max_abs_log_likelihood_difference": difference,
                    })
    for order in args.orders:
        for nval in args.validation_sizes:
            for arm in (
                "A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
            ):
                key = f"{arm}_n{nval}"
                left = likelihood_store[(1, order, key)]
                right = likelihood_store[(2, order, key)]
                difference = float(np.max(np.abs(left-right)))
                if not np.allclose(left, right, rtol=1e-11, atol=1e-11):
                    raise AssertionError(
                        f"{key} must be degree-invariant; max absolute difference={difference}"
                    )
                invariant_rows.append({
                    "invariant": "fixed_arm_across_model_degree",
                    "degree": 0, "quadrature_order": order,
                    "validation_size": nval, "arm": arm,
                    "max_abs_log_likelihood_difference": difference,
                })
    pd.DataFrame(invariant_rows).to_csv(out / "likelihood_invariance_checks.csv", index=False)

    posterior_rows = pd.DataFrame(records)
    posterior_rows.to_csv(out / "posterior_summaries.csv", index=False)
    diagnostics = pd.concat(diagnostic_parts, ignore_index=True)
    diagnostics.to_csv(out / "weight_diagnostics.csv", index=False)
    convergence = []
    for degree in (1, 2):
        for nval in args.validation_sizes:
            key = f"E_uncertainty_integrated_ipw_n{nval}"
            for low, high in ((5, 7), (7, 9), (5, 9)):
                p_low, q_low = posterior_store[(degree, low, key)]
                p_high, q_high = posterior_store[(degree, high, key)]
                convergence.append({
                    "degree": degree, "validation_size": nval,
                    "lower_order": low, "higher_order": high,
                    "posterior_density_l1": float(trapezoid(np.abs(p_low-p_high), h0_grid)),
                    "median_abs_shift": abs(q_high["q500"]-q_low["q500"]),
                    "q025_abs_shift": abs(q_high["q025"]-q_low["q025"]),
                    "q16_abs_shift": abs(q_high["q160"]-q_low["q160"]),
                    "q84_abs_shift": abs(q_high["q840"]-q_low["q840"]),
                    "q975_abs_shift": abs(q_high["q975"]-q_low["q975"]),
                })
    convergence_df = pd.DataFrame(convergence)
    convergence_df.to_csv(out / "integrated_ipw_quadrature_convergence.csv", index=False)

    degree_contrasts = []
    for order in args.orders:
        for nval in args.validation_sizes:
            for arm in (
                "A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
                "D_map_plugin_ipw", "E_uncertainty_integrated_ipw",
            ):
                one = posterior_store[(1, order, f"{arm}_n{nval}")][1]
                two = posterior_store[(2, order, f"{arm}_n{nval}")][1]
                degree_contrasts.append({
                    "quadrature_order": order, "validation_size": nval, "arm": arm,
                    "degree1_median": one["q500"], "degree2_median": two["q500"],
                    "degree2_minus_degree1_median": two["q500"]-one["q500"],
                    "degree1_q025": one["q025"], "degree2_q025": two["q025"],
                    "degree1_q975": one["q975"], "degree2_q975": two["q975"],
                })
    pd.DataFrame(degree_contrasts).to_csv(out / "degree_model_posterior_contrasts.csv", index=False)

    config = {
        "experiment": "paired_degree1_degree2_quadratic_logit_inclusion_quadrature_sensitivity",
        "status": "engineering_numerical_pilot_not_calibration_not_observational",
        "seed": args.seed, "validation_seed": args.validation_seed,
        "n_parent": args.n_parent, "n_detected": args.n_detected,
        "validation_frame": args.validation_frame, "validation_sizes": args.validation_sizes,
        "gamma": args.gamma, "curvature_beta": args.curvature_beta,
        "degrees": [1, 2], "quadrature_orders": args.orders,
        "h0_true": args.h0_true, "h0_prior_support": [args.h0_prior_min, args.h0_prior_max],
        "h0_points": args.h0_points, "scott_multiplier": args.scott_multiplier,
        "input_pairing_checks": input_checks,
        "independent_child_audits": len(args.orders)*2,
        "claim_limit": "Numerical quadrature and paired-model sensitivity only; not SBC, calibration, or real-data inference.",
    }
    (out / "config.json").write_text(json.dumps(config, indent=2, sort_keys=True)+"\n")
    files = sorted(p for p in out.rglob("*") if p.is_file() and p.name != "manifest.json")
    manifest = {"status":"PASS","file_count":len(files),"files":{
        str(p.relative_to(out)):{"bytes":p.stat().st_size,"sha256":_sha(p)} for p in files
    }}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n")
    return {
        "status":"PASS", "child_runs":len(args.orders)*2,
        "child_audits_passed":len(args.orders)*2,
        "input_pairing_checks":input_checks,
        "posterior_rows":len(posterior_rows), "weight_diagnostic_rows":len(diagnostics),
        "likelihood_invariance_checks":len(invariant_rows),
        "convergence_rows":len(convergence_df),
        "maximum_9_vs_7_density_l1":float(convergence_df.loc[
            (convergence_df.lower_order==7)&(convergence_df.higher_order==9),
            "posterior_density_l1"].max()),
        "output_dir":str(out),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--seed", type=int, default=20261120)
    p.add_argument("--validation-seed", type=int, default=20261121)
    p.add_argument("--n-parent", type=int, default=500)
    p.add_argument("--n-detected", type=int, default=4)
    p.add_argument("--validation-frame", type=int, default=300)
    p.add_argument("--validation-sizes", type=lambda s:[int(x) for x in s.split(",")], default=[100,300])
    p.add_argument("--gamma", type=float, default=1.0)
    p.add_argument("--curvature-beta", type=float, default=-3.9)
    p.add_argument("--h0-true", type=float, default=70.0)
    p.add_argument("--h0-prior-min", type=float, default=30.0)
    p.add_argument("--h0-prior-max", type=float, default=110.0)
    p.add_argument("--h0-points", type=int, default=41)
    p.add_argument("--scott-multiplier", type=float, default=0.7)
    p.add_argument("--sigma-mpc", type=float, default=140.0)
    p.add_argument("--detection-threshold-mpc", type=float, default=250.0)
    p.add_argument("--orders", type=lambda s:[int(x) for x in s.split(",")], default=[5,7,9])
    p.add_argument("--output-dir", type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(run(a), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
