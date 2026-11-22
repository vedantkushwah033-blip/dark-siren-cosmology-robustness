"""Evaluate exactly one arm × validation-size cell from a frozen shared-input bundle.

Execution primitive only: no lifecycle transition and no independent audit.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile

import numpy as np
import pandas as pd

from dark_siren_robustness.arm_cell_execution import ARM_NAMES, evaluate_arm_cell
from dark_siren_robustness.shared_input_bundle import load_shared_input_manifest


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_value(value):
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def run_cell(bundle_dir: Path, output_dir: Path, *, replicate: int, gamma: float,
             validation_size: int, arm: str) -> dict:
    bundle = Path(bundle_dir).resolve()
    destination = Path(output_dir).resolve()
    if destination.exists():
        raise FileExistsError(f"cell output already exists: {destination}")
    if arm not in ARM_NAMES:
        raise ValueError(f"unknown arm: {arm}")

    shared = load_shared_input_manifest(bundle)
    sizes = shared["validation_sizes"]
    if validation_size not in sizes:
        raise ValueError("validation_size is not part of the frozen shared-input plan")
    provenance = shared["provenance"]
    if float(provenance.get("curvature_gamma")) != float(gamma):
        raise ValueError("requested gamma differs from frozen shared-input provenance")

    config = json.loads((bundle / "config.json").read_text(encoding="utf-8"))
    if int(config["seed"]) != int(provenance["seed"]):
        raise ValueError("config/shared-input seed mismatch")
    if list(config["validation_sizes"]) != list(sizes):
        raise ValueError("config/shared-input validation-size mismatch")

    catalogue = pd.read_csv(bundle / "observable/catalogue_observed.csv")
    parent = pd.read_csv(bundle / "truth/parent_truth.csv")
    events_frame = pd.read_csv(bundle / "observable/gw_events_observed.csv")
    validation = pd.read_csv(bundle / f"observable/validation_n{validation_size}.csv")

    catalogue_z = catalogue["redshift"].to_numpy(float)
    parent_z = parent["redshift"].to_numpy(float)
    catalogue_base = 1.0 + 1.5 * (catalogue_z - 0.01) / 0.39
    parent_base = 1.0 + 1.5 * (parent_z - 0.01) / 0.39
    true_probability = catalogue.merge(
        parent[["object_id", "inclusion_probability"]],
        on="object_id", validate="one_to_one",
    )["inclusion_probability"].to_numpy(float)
    events = [
        {"gw_distance": float(row.distance_mpc), "gw_sigma": float(config["sigma_mpc"])}
        for row in events_frame.itertuples(index=False)
    ]
    prior_min, prior_max = map(float, config["h0_prior_support"])
    h0_grid = np.linspace(prior_min, prior_max, int(config["h0_grid_points"]))
    redshift_grid = np.linspace(0.01, 0.40, int(config["redshift_grid_points"]))
    merger_weight = np.ones_like(redshift_grid)
    degree = int(config["inclusion_model_degree"])

    result = evaluate_arm_cell(
        arm,
        parent_redshift=parent_z if arm == "B_parent_oracle" else None,
        catalogue_redshift=catalogue_z if arm != "B_parent_oracle" else None,
        validation_redshift=(validation["redshift"].to_numpy(float)
                             if arm in {"D_map_plugin_ipw", "E_uncertainty_integrated_ipw"} else None),
        validation_retained=(validation["retained"].to_numpy(np.int8)
                             if arm in {"D_map_plugin_ipw", "E_uncertainty_integrated_ipw"} else None),
        true_catalogue_inclusion_probability=(true_probability
                                              if arm == "C_known_probability_ipw" else None),
        events=events, H0_grid=h0_grid, redshift_grid=redshift_grid,
        merger_weight=merger_weight, h0_true=float(config["h0_true"]),
        h0_prior_min=prior_min, h0_prior_max=prior_max,
        detection_threshold_mpc=float(config["detection_threshold_mpc"]),
        scott_multiplier=float(config["scott_multiplier"]),
        quadrature_order=int(config["quadrature_order"]),
        inclusion_degree=degree, coordinate_center=0.205, coordinate_scale=0.195,
        quadratic_center=float(config.get("fitted_quadratic_basis_center", 0.0)),
        parent_base_weights=parent_base if arm == "B_parent_oracle" else None,
        catalogue_base_weights=catalogue_base if arm != "B_parent_oracle" else None,
    )

    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{destination.name}.", dir=destination.parent))
    try:
        np.savez_compressed(
            temporary / "cell_arrays.npz",
            H0_grid=h0_grid, redshift_grid=redshift_grid,
            raw_weights=result.raw_weights,
            component_weights=result.component_weights,
            mixture_weights=result.mixture_weights,
            component_density=result.component_density,
            component_log_likelihood=result.component_log_likelihood,
            marginal_log_likelihood=result.marginal_log_likelihood,
            component_log_detection_efficiency=result.component_log_detection_efficiency,
            realized_bandwidths=result.realized_bandwidths,
            inclusion_mode=(np.asarray([]) if result.inclusion_mode is None else result.inclusion_mode),
            inclusion_nodes=(np.asarray([]) if result.inclusion_nodes is None else result.inclusion_nodes),
            inclusion_quadrature_weights=(np.asarray([]) if result.inclusion_quadrature_weights is None else result.inclusion_quadrature_weights),
        )
        metadata = {
            "schema_version": 1,
            "status": "COMPUTED_NOT_AUDITED",
            "cell_key": {"replicate": int(replicate), "gamma": float(gamma),
                         "validation_size": int(validation_size), "arm": arm},
            "bundle_manifest_sha256": _sha256(bundle / "shared_input_manifest.json"),
            "config_sha256": _sha256(bundle / "config.json"),
            "posterior_summary": _json_value(dict(result.posterior_summary)),
            "weight_diagnostics": _json_value(dict(result.weight_diagnostics)),
            "realized_bandwidths": _json_value(result.realized_bandwidths),
            "inclusion_mode": _json_value(result.inclusion_mode),
            "inclusion_nodes": _json_value(result.inclusion_nodes),
            "inclusion_quadrature_weights": _json_value(result.inclusion_quadrature_weights),
            "limitations": [
                "Numerical cell execution only; no independent audit has run.",
                "This command does not mutate the lifecycle checkpoint.",
                "No confirmatory calibration or estimator-superiority claim is implied.",
            ],
        }
        (temporary / "cell_result.json").write_text(
            json.dumps(metadata, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        files = {}
        for path in sorted(temporary.iterdir()):
            if path.is_file():
                files[path.name] = {"bytes": path.stat().st_size, "sha256": _sha256(path)}
        manifest = {"schema_version": 1, "status": "COMPUTED_NOT_AUDITED", "files": files}
        (temporary / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        if destination.exists():
            raise FileExistsError(f"cell output already exists: {destination}")
        os.rename(temporary, destination)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    return {"status": "COMPUTED_NOT_AUDITED", "output_dir": str(destination),
            "cell_key": metadata["cell_key"], "posterior_summary": metadata["posterior_summary"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--gamma", type=float, required=True)
    parser.add_argument("--validation-size", type=int, required=True)
    parser.add_argument("--arm", choices=ARM_NAMES, required=True)
    args = parser.parse_args()
    print(json.dumps(run_cell(
        args.bundle_dir, args.output_dir, replicate=args.replicate,
        gamma=args.gamma, validation_size=args.validation_size, arm=args.arm,
    ), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
