"""Generate one immutable synthetic replicate/gamma input bundle without inference.

This is the input-preparation boundary for independent arm × validation-size
workers. It deliberately performs no inclusion-model fit, KDE, H0 likelihood,
or posterior calculation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from dark_siren_robustness.generative_selection import (
    generate_quadratic_logit_inclusion_samples,
    generate_synthetic_population,
)
from dark_siren_robustness.shared_input_bundle import write_shared_input_manifest


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run(args: argparse.Namespace) -> dict:
    output = args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("output directory must be empty; shared input bundles are immutable")
    output.mkdir(parents=True, exist_ok=True)
    observable = output / "observable"
    truth = output / "truth"
    observable.mkdir()
    truth.mkdir()
    sizes = list(args.validation_sizes)
    if not sizes or sizes != sorted(set(sizes)) or sizes[0] < 5:
        raise ValueError("validation sizes must be unique, increasing, and >=5")
    if sizes[-1] > args.validation_frame:
        raise ValueError("largest validation size exceeds validation frame")
    if args.inclusion_degree not in (1, 2):
        raise ValueError("inclusion degree must be 1 or 2")
    if not np.isfinite(args.h0_true) or args.h0_true <= 0:
        raise ValueError("H0 truth must be finite and positive")

    population = generate_synthetic_population(
        seed=args.seed, event_seed=args.event_seed,
        n_parent=args.n_parent, n_detected=args.n_detected,
        h0_true=args.h0_true, sigma_mpc=args.sigma_mpc,
        detection_threshold_mpc=args.detection_threshold_mpc,
        z_min=0.01, z_max=0.40, catalogue_floor=0.10,
        catalogue_midpoint=0.22, catalogue_width=0.05,
    )
    quadratic = generate_quadratic_logit_inclusion_samples(
        population_seed=args.seed, validation_seed=args.validation_seed,
        validation_sample_seed=args.validation_sample_seed,
        n_parent=args.n_parent, n_validation_frame=args.validation_frame,
        n_validation=sizes[-1], gamma=args.curvature_gamma,
        beta=args.curvature_beta, z_min=0.01, z_max=0.40,
    )
    parent = quadratic["parent_truth"]
    if not np.array_equal(
        population["parent_truth"]["redshift"].to_numpy(float),
        parent["redshift"].to_numpy(float),
    ):
        raise RuntimeError("GW and catalogue generators do not share the same parent population")
    catalogue = quadratic["catalogue_observed"]
    validation_observed = quadratic["validation_observed"]
    validation_truth = quadratic["validation_frame_truth"]

    catalogue.to_csv(observable / "catalogue_observed.csv", index=False)
    population["gw_events_observed"].to_csv(observable / "gw_events_observed.csv", index=False)
    parent.to_csv(truth / "parent_truth.csv", index=False)
    population["gw_event_truth"].to_csv(truth / "gw_event_truth.csv", index=False)
    population["gw_generation_trace_truth"].to_csv(
        truth / "gw_generation_trace_truth.csv", index=False
    )
    validation_truth.to_csv(truth / "validation_frame_truth.csv", index=False)
    for size in sizes:
        validation_observed.iloc[:size].copy().to_csv(
            observable / f"validation_n{size}.csv", index=False
        )

    config = {
        "experiment": "quadratic_logit_shared_input_bundle",
        "status": "shared_inputs_only_no_inference",
        "seed": int(args.seed),
        "event_seed": args.event_seed,
        "validation_seed": int(args.validation_seed),
        "validation_sample_seed": args.validation_sample_seed,
        "n_parent": int(args.n_parent),
        "n_catalogue_observed": int(len(catalogue)),
        "n_detected": int(args.n_detected),
        "validation_frame": int(args.validation_frame),
        "validation_sizes": sizes,
        "catalogue_model": "quadratic_logit",
        "inclusion_model_degree": int(args.inclusion_degree),
        "curvature_gamma": float(args.curvature_gamma),
        "curvature_beta": float(args.curvature_beta),
        "quadratic_calibration": quadratic["config"],
        "h0_true": float(args.h0_true),
        "h0_prior_support": [float(args.h0_prior_min), float(args.h0_prior_max)],
        "h0_grid_points": int(args.h0_points),
        "redshift_grid_points": 501,
        "scott_multiplier": float(args.scott_multiplier),
        "quadrature_order": int(args.quadrature_order),
        "sigma_mpc": float(args.sigma_mpc),
        "detection_threshold_mpc": float(args.detection_threshold_mpc),
        "truth_used_for_non_oracle_inference": False,
        "claim_limit": "Input-generation artifact only; no arm inference or calibration.",
    }
    config_path = output / "config.json"
    config_path.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    provenance = {
        "seed": int(args.seed), "event_seed": args.event_seed,
        "validation_seed": int(args.validation_seed),
        "validation_sample_seed": args.validation_sample_seed,
        "n_parent": int(args.n_parent), "n_detected": int(args.n_detected),
        "h0_true": float(args.h0_true), "sigma_mpc": float(args.sigma_mpc),
        "detection_threshold_mpc": float(args.detection_threshold_mpc),
        "n_validation_frame": int(args.validation_frame),
        "curvature_gamma": float(args.curvature_gamma),
        "curvature_beta": float(args.curvature_beta),
        "inclusion_model_degree": int(args.inclusion_degree),
        "quadratic_basis_center": float(quadratic["config"]["x2_center"]),
    }
    shared = write_shared_input_manifest(
        output, validation_sizes=sizes, provenance=provenance,
    )
    files = sorted(
        path for path in output.rglob("*")
        if path.is_file() and path.name not in {"manifest.json"}
    )
    manifest = {
        "status": "PASS",
        "file_count": len(files),
        "files": {
            str(path.relative_to(output)): {
                "bytes": path.stat().st_size, "sha256": _sha256(path),
            } for path in files
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return {
        "status": "PASS", "output_dir": str(output),
        "shared_input_files": len(shared["files"]),
        "validation_sizes": sizes,
        "inference_performed": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--event-seed", type=int, default=None)
    parser.add_argument("--validation-seed", type=int, required=True)
    parser.add_argument("--validation-sample-seed", type=int, default=None)
    parser.add_argument("--n-parent", type=int, default=1200)
    parser.add_argument("--n-detected", type=int, default=30)
    parser.add_argument("--validation-frame", type=int, default=500)
    parser.add_argument("--validation-sizes", type=lambda x: [int(v) for v in x.split(",")], default=[100, 300, 500])
    parser.add_argument("--curvature-gamma", type=float, required=True)
    parser.add_argument("--curvature-beta", type=float, default=-3.9)
    parser.add_argument("--inclusion-degree", type=int, choices=(1, 2), default=1)
    parser.add_argument("--h0-true", type=float, required=True)
    parser.add_argument("--h0-prior-min", type=float, default=50.0)
    parser.add_argument("--h0-prior-max", type=float, default=90.0)
    parser.add_argument("--h0-points", type=int, default=1121)
    parser.add_argument("--scott-multiplier", type=float, default=0.7)
    parser.add_argument("--quadrature-order", type=int, default=5)
    parser.add_argument("--sigma-mpc", type=float, default=140.0)
    parser.add_argument("--detection-threshold-mpc", type=float, default=250.0)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.h0_prior_min >= args.h0_prior_max or not args.h0_prior_min <= args.h0_true <= args.h0_prior_max:
        parser.error("H0 truth must lie within finite increasing prior support")
    print(json.dumps(run(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
