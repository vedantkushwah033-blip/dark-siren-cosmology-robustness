"""Execute and persist exactly one (replicate, gamma, validation_size, arm) cell.

The input directory is an immutable output bundle from the quadratic five-arm
simulation generator. This worker reads only the inputs required by the chosen
arm, verifies their hashes against the generator manifest, and writes one
self-contained numerical result. It does not update the lifecycle ledger.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dark_siren_robustness.arm_cell_execution import ARM_NAMES, evaluate_arm_cell

INPUTS = {
    "config": "config.json",
    "manifest": "manifest.json",
    "catalogue": "observable/catalogue_observed.csv",
    "events": "observable/gw_events_observed.csv",
    "parent": "truth/parent_truth.csv",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_verified(root: Path, relative: str, manifest: dict) -> Path:
    path = root / relative
    if not path.is_file():
        raise FileNotFoundError(f"required shared input is missing: {relative}")
    entry = manifest.get("files", {}).get(relative)
    if not isinstance(entry, dict):
        raise ValueError(f"shared input is not listed in the immutable manifest: {relative}")
    if path.stat().st_size != int(entry.get("bytes", -1)):
        raise ValueError(f"shared input byte count differs from manifest: {relative}")
    if _sha256(path) != entry.get("sha256"):
        raise ValueError(f"shared input SHA-256 differs from manifest: {relative}")
    return path


def _write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


def run(args: argparse.Namespace) -> dict:
    root = args.input_dir.resolve()
    out = args.output_dir.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError("cell output directory must be absent or empty")
    out.mkdir(parents=True, exist_ok=True)

    manifest_path = root / INPUTS["manifest"]
    config_path = root / INPUTS["config"]
    if not manifest_path.is_file() or not config_path.is_file():
        raise FileNotFoundError("shared input bundle requires config.json and manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS":
        raise ValueError("shared input bundle manifest is not PASS")
    if config.get("catalogue_model") != "quadratic_logit":
        raise ValueError("input bundle is not the expected quadratic-logit design")
    if args.arm not in ARM_NAMES:
        raise ValueError(f"unknown arm: {args.arm}")
    if args.validation_size not in config["validation_sizes"]:
        raise ValueError("requested validation size is not present in the frozen input bundle")

    required = [INPUTS["config"], INPUTS["manifest"], INPUTS["events"]]
    if args.arm != "B_parent_oracle":
        required.append(INPUTS["catalogue"])
    if args.arm in {"B_parent_oracle", "C_known_probability_ipw"}:
        required.append(INPUTS["parent"])
    validation_relative = f"observable/validation_n{args.validation_size}.csv"
    if args.arm in {"D_map_plugin_ipw", "E_uncertainty_integrated_ipw"}:
        required.append(validation_relative)
    verified = {relative: _load_verified(root, relative, manifest) for relative in required}

    events_frame = pd.read_csv(verified[INPUTS["events"]])
    events = [
        {"gw_distance": float(row.distance_mpc), "gw_sigma": float(config["sigma_mpc"])}
        for row in events_frame.itertuples(index=False)
    ]
    parent = catalogue = validation_z = validation_retained = true_probability = None
    parent_base = catalogue_base = None
    if INPUTS["parent"] in verified:
        parent_frame = pd.read_csv(verified[INPUTS["parent"]])
        parent = parent_frame["redshift"].to_numpy(float)
        parent_base = 1.0 + 1.5 * (parent - 0.01) / 0.39
    if INPUTS["catalogue"] in verified:
        catalogue_frame = pd.read_csv(verified[INPUTS["catalogue"]])
        catalogue = catalogue_frame["redshift"].to_numpy(float)
        catalogue_base = 1.0 + 1.5 * (catalogue - 0.01) / 0.39
    if args.arm == "C_known_probability_ipw":
        probability_by_id = parent_frame.set_index("object_id")["inclusion_probability"]
        true_probability = catalogue_frame["object_id"].map(probability_by_id).to_numpy(float)
        if not np.all(np.isfinite(true_probability)):
            raise ValueError("catalogue IDs could not be joined to parent inclusion probabilities")
    if args.arm in {"D_map_plugin_ipw", "E_uncertainty_integrated_ipw"}:
        validation = pd.read_csv(verified[validation_relative])
        validation_z = validation["redshift"].to_numpy(float)
        validation_retained = validation["retained"].to_numpy(np.int8)

    prior_min, prior_max = map(float, config["h0_prior_support"])
    h0_grid = np.linspace(prior_min, prior_max, int(config["h0_grid_points"]))
    redshift_grid = np.linspace(0.01, 0.40, int(config["redshift_grid_points"]))
    result = evaluate_arm_cell(
        args.arm,
        parent_redshift=parent,
        catalogue_redshift=catalogue,
        validation_redshift=validation_z,
        validation_retained=validation_retained,
        true_catalogue_inclusion_probability=true_probability,
        events=events,
        H0_grid=h0_grid,
        redshift_grid=redshift_grid,
        merger_weight=np.ones_like(redshift_grid),
        h0_true=float(config["h0_true"]),
        h0_prior_min=prior_min,
        h0_prior_max=prior_max,
        detection_threshold_mpc=float(config["detection_threshold_mpc"]),
        scott_multiplier=float(config["scott_multiplier"]),
        quadrature_order=int(config["quadrature_order"]),
        inclusion_degree=int(config["inclusion_model_degree"]),
        coordinate_center=0.205,
        coordinate_scale=0.195,
        quadratic_center=float(config["quadratic_calibration"]["x2_center"])
            if int(config["inclusion_model_degree"]) == 2 else 0.0,
        parent_base_weights=parent_base,
        catalogue_base_weights=catalogue_base,
    )
    summary = {
        "replicate": args.replicate,
        "gamma": float(config["curvature_gamma"]),
        "validation_size": args.validation_size,
        "arm": args.arm,
        "status": result.status,
        "h0_true": float(config["h0_true"]),
        "median_error": float(result.posterior_summary["posterior_median"] - config["h0_true"]),
        **{key: float(value) for key, value in result.posterior_summary.items()},
        "weight_diagnostics": {key: float(value) for key, value in result.weight_diagnostics.items()},
        "realized_bandwidths": [float(value) for value in result.realized_bandwidths],
        "input_sha256": {relative: _sha256(path) for relative, path in verified.items()},
    }
    np.savez_compressed(
        out / "arm_cell_arrays.npz",
        H0_grid=h0_grid,
        redshift_grid=redshift_grid,
        raw_weights=result.raw_weights,
        component_weights=result.component_weights,
        mixture_weights=result.mixture_weights,
        component_density=result.component_density,
        component_log_likelihood=result.component_log_likelihood,
        marginal_log_likelihood=result.marginal_log_likelihood,
        component_log_detection_efficiency=result.component_log_detection_efficiency,
        realized_bandwidths=result.realized_bandwidths,
        inclusion_mode=(np.asarray([]) if result.inclusion_mode is None else result.inclusion_mode),
        inclusion_nodes=(np.asarray([[]]) if result.inclusion_nodes is None else result.inclusion_nodes),
        inclusion_quadrature_weights=(np.asarray([]) if result.inclusion_quadrature_weights is None else result.inclusion_quadrature_weights),
    )
    _write_json(out / "summary.json", summary)
    provenance = {
        "schema_version": 1,
        "replicate": args.replicate,
        "gamma": float(config["curvature_gamma"]),
        "validation_size": args.validation_size,
        "arm": args.arm,
        "source_config_sha256": _sha256(config_path),
        "source_manifest_sha256": _sha256(manifest_path),
        "verified_inputs": summary["input_sha256"],
        "claim_limit": "One independently executed synthetic arm cell; not calibration or observational inference.",
    }
    _write_json(out / "provenance.json", provenance)
    files = sorted(path for path in out.iterdir() if path.is_file())
    output_manifest = {
        "status": "PASS",
        "file_count": len(files),
        "files": {
            path.name: {"bytes": path.stat().st_size, "sha256": _sha256(path)}
            for path in files
        },
    }
    _write_json(out / "manifest.json", output_manifest)
    return {"status": "PASS", "arm": args.arm, "validation_size": args.validation_size,
            "output_dir": str(out), "files": output_manifest["file_count"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--validation-size", type=int, required=True)
    parser.add_argument("--arm", choices=ARM_NAMES, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
