"""Prepare immutable synthetic inputs without fitting or evaluating any inference arm.

This command creates only the shared observable/truth bundle consumed by
independent arm-cell workers. It deliberately imports no inference, KDE,
posterior-fitting, or likelihood modules.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile

import pandas as pd

from dark_siren_robustness.generative_selection import (
    generate_quadratic_logit_inclusion_samples,
    generate_synthetic_population,
)
from dark_siren_robustness.shared_input_bundle import write_shared_input_manifest


def _parse_sizes(raw: str) -> list[int]:
    try:
        values = [int(item.strip()) for item in raw.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("validation sizes must be comma-separated integers") from exc
    if not values or values != sorted(set(values)) or any(value < 5 for value in values):
        raise argparse.ArgumentTypeError("validation sizes must be unique, increasing integers >= 5")
    return values


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def prepare(args: argparse.Namespace) -> dict[str, object]:
    destination = args.output_dir.resolve()
    if destination.exists():
        raise FileExistsError(f"refusing to overwrite existing bundle: {destination}")
    sizes = args.validation_sizes
    if args.n_parent < 20 or args.n_detected < 1 or args.validation_frame < max(sizes):
        raise ValueError("population/event/frame sizes are inconsistent with the requested design")
    if args.inclusion_degree not in (1, 2):
        raise ValueError("inclusion degree must be 1 or 2")

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{destination.name}.", dir=destination.parent))
    try:
        observable = temporary / "observable"
        truth = temporary / "truth"
        observable.mkdir()
        truth.mkdir()

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
        if not (population["parent_truth"]["redshift"].to_numpy(float)
                == quadratic["parent_truth"]["redshift"].to_numpy(float)).all():
            raise RuntimeError("GW and catalogue generators do not share identical parent redshifts")

        catalogue = quadratic["catalogue_observed"]
        validation_full = quadratic["validation_observed"]
        parent_truth = quadratic["parent_truth"]
        events = population["gw_events_observed"]

        catalogue.to_csv(observable / "catalogue_observed.csv", index=False)
        events.to_csv(observable / "gw_events_observed.csv", index=False)
        parent_truth.to_csv(truth / "parent_truth.csv", index=False)
        population["gw_event_truth"].to_csv(truth / "gw_event_truth.csv", index=False)
        population["gw_generation_trace_truth"].to_csv(
            truth / "gw_generation_trace_truth.csv", index=False
        )
        quadratic["validation_frame_truth"].to_csv(
            truth / "validation_frame_truth.csv", index=False
        )
        for n in sizes:
            validation_full.iloc[:n].copy().to_csv(
                observable / f"validation_n{n}.csv", index=False
            )

        basis_center = float(quadratic["config"]["x2_center"])
        config = {
            "experiment": "preparation_only_shared_input_bundle",
            "status": "FROZEN_INPUTS_ONLY_NO_INFERENCE_EXECUTED",
            "seed": int(args.seed),
            "event_seed": args.event_seed,
            "validation_seed": int(args.validation_seed),
            "validation_sample_seed": args.validation_sample_seed,
            "n_parent": int(args.n_parent),
            "n_catalogue_observed": int(len(catalogue)),
            "n_detected": int(args.n_detected),
            "validation_frame": int(args.validation_frame),
            "validation_sizes": sizes,
            "catalogue_floor": 0.10,
            "catalogue_model": "quadratic_logit",
            "inclusion_model_degree": int(args.inclusion_degree),
            "fitted_quadratic_basis_center": basis_center if args.inclusion_degree == 2 else 0.0,
            "fitted_coefficient_basis": "[1, x, x^2 - c_g]" if args.inclusion_degree == 2 else "[1, x]",
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
            "truth_used_for_arm_c_only": True,
            "claim_limit": "Synthetic input preparation only; no inference, calibration, or observational claim.",
        }
        (temporary / "config.json").write_text(
            json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        frozen = write_shared_input_manifest(
            temporary, validation_sizes=sizes,
            provenance={
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
                "quadratic_basis_center": basis_center,
            },
        )
        files = sorted(p for p in temporary.rglob("*") if p.is_file())
        manifest = {
            "schema_version": 1, "status": "PASS",
            "preparation_only": True, "inference_executed": False,
            "shared_input_manifest_sha256": _sha256(temporary / "shared_input_manifest.json"),
            "file_count": len(files),
            "files": {
                str(p.relative_to(temporary)): {
                    "bytes": p.stat().st_size, "sha256": _sha256(p)
                } for p in files
            },
        }
        (temporary / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.rename(temporary, destination)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise

    return {
        "status": "PASS", "output_dir": str(destination),
        "preparation_only": True, "inference_executed": False,
        "shared_input_files": len(frozen["files"]),
        "validation_sizes": sizes, "n_parent": int(args.n_parent),
        "n_catalogue_observed": int(len(catalogue)), "n_detected": int(args.n_detected),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261122)
    parser.add_argument("--event-seed", type=int, default=None)
    parser.add_argument("--validation-seed", type=int, default=20261123)
    parser.add_argument("--validation-sample-seed", type=int, default=None)
    parser.add_argument("--n-parent", type=int, default=1200)
    parser.add_argument("--n-detected", type=int, default=4)
    parser.add_argument("--validation-frame", type=int, default=1500)
    parser.add_argument("--validation-sizes", type=_parse_sizes, default=_parse_sizes("100,300,1000"))
    parser.add_argument("--curvature-gamma", type=float, default=0.0)
    parser.add_argument("--curvature-beta", type=float, default=-3.9)
    parser.add_argument("--inclusion-degree", type=int, choices=(1, 2), default=1)
    parser.add_argument("--h0-true", type=float, default=70.0)
    parser.add_argument("--h0-prior-min", type=float, default=50.0)
    parser.add_argument("--h0-prior-max", type=float, default=90.0)
    parser.add_argument("--h0-points", type=int, default=1121)
    parser.add_argument("--scott-multiplier", type=float, default=0.7)
    parser.add_argument("--quadrature-order", type=int, default=5)
    parser.add_argument("--sigma-mpc", type=float, default=140.0)
    parser.add_argument("--detection-threshold-mpc", type=float, default=250.0)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
