"""Independent numerical audit for one arm-cell computed from raw prepared inputs.

The audit reconstructs a single arm directly from observable/truth-separated
bundle inputs using the separate reference implementation in
audit_quadratic_five_arm_smoke.py. It never imports the production arm evaluator,
KDE, likelihood, or posterior code.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit, logsumexp

from audit_prepared_shared_input_bundle import audit as audit_prepared_bundle

_HELPER_PATH = Path(__file__).with_name("audit_quadratic_five_arm_smoke.py")
_SPEC = importlib.util.spec_from_file_location("_independent_dark_siren_reference", _HELPER_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise ImportError("cannot load independent numerical reference helpers")
_REF = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_REF)

ARMS = {
    "A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
    "D_map_plugin_ipw", "E_uncertainty_integrated_ipw",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _close(actual, expected, label: str, *, rtol: float = 2e-8, atol: float = 2e-10):
    a, e = np.asarray(actual, dtype=float), np.asarray(expected, dtype=float)
    if a.shape != e.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(e)):
        raise AssertionError(f"{label}: shape or finite-value mismatch")
    if not np.allclose(a, e, rtol=rtol, atol=atol):
        delta = float(np.max(np.abs(a - e)))
        raise AssertionError(f"{label}: independent reconstruction mismatch (max abs={delta})")


def _diagnostics(weights: np.ndarray) -> dict[str, float]:
    raw = np.asarray(weights, dtype=float)
    if raw.ndim != 1 or raw.size < 1 or not np.all(np.isfinite(raw)) or np.any(raw <= 0):
        raise AssertionError("independent diagnostics received invalid weights")
    scaled = raw / np.max(raw)
    if np.any(scaled == 0):
        raise AssertionError("weight scaling underflowed positive mass")
    normalized = scaled / np.sum(scaled)
    ordered = np.sort(normalized)[::-1]
    entropy_ess = float(np.exp(-np.sum(normalized * np.log(normalized))))
    return {
        "kish_ess": float(1.0 / np.sum(normalized**2)),
        "kish_ess_fraction": float(1.0 / np.sum(normalized**2) / raw.size),
        "entropy_ess": entropy_ess,
        "maximum_normalized_weight": float(ordered[0]),
        "top_1_percent_mass": float(np.sum(ordered[:max(1, int(np.ceil(.01 * raw.size)))])),
        "top_5_percent_mass": float(np.sum(ordered[:max(1, int(np.ceil(.05 * raw.size)))])),
    }


def _verify_cell_files(cell: Path) -> dict:
    manifest_path = cell / "manifest.json"
    metadata_path = cell / "cell_result.json"
    arrays_path = cell / "cell_arrays.npz"
    if not all(path.is_file() for path in (manifest_path, metadata_path, arrays_path)):
        raise AssertionError("cell artifact files are incomplete")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1 or manifest.get("status") != "COMPUTED_NOT_AUDITED":
        raise AssertionError("cell artifact manifest schema/status is invalid")
    if set(manifest.get("files", {})) != {"cell_arrays.npz", "cell_result.json"}:
        raise AssertionError("cell artifact manifest file set is invalid")
    for name, meta in manifest["files"].items():
        path = cell / name
        if path.stat().st_size != meta.get("bytes") or _sha256(path) != meta.get("sha256"):
            raise AssertionError(f"cell artifact integrity mismatch: {name}")
    return manifest


def audit_cell(bundle_dir: Path, cell_dir: Path) -> dict:
    bundle, cell = Path(bundle_dir).resolve(), Path(cell_dir).resolve()
    raw_report = audit_prepared_bundle(bundle)
    if raw_report.get("status") != "PASS":
        raise AssertionError("prepared-bundle raw-input audit did not PASS")
    cell_manifest = _verify_cell_files(cell)
    metadata = json.loads((cell / "cell_result.json").read_text(encoding="utf-8"))
    if metadata.get("schema_version") != 1 or metadata.get("status") != "COMPUTED_NOT_AUDITED":
        raise AssertionError("cell result schema/status is invalid")
    key = metadata.get("cell_key", {})
    if set(key) != {"replicate", "gamma", "validation_size", "arm"}:
        raise AssertionError("cell key fields are invalid")
    arm, n, gamma = str(key["arm"]), int(key["validation_size"]), float(key["gamma"])
    if arm not in ARMS:
        raise AssertionError(f"unknown arm: {arm}")

    shared_path, config_path = bundle / "shared_input_manifest.json", bundle / "config.json"
    if metadata.get("bundle_manifest_sha256") != _sha256(shared_path):
        raise AssertionError("cell is bound to a different shared-input manifest")
    if metadata.get("config_sha256") != _sha256(config_path):
        raise AssertionError("cell is bound to a different config")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if n not in config["validation_sizes"] or gamma != float(config["curvature_gamma"]):
        raise AssertionError("cell gamma/validation size differs from frozen plan")

    parent = pd.read_csv(bundle / "truth/parent_truth.csv")
    catalogue = pd.read_csv(bundle / "observable/catalogue_observed.csv")
    validation = pd.read_csv(bundle / f"observable/validation_n{n}.csv")
    event_frame = pd.read_csv(bundle / "observable/gw_events_observed.csv")
    parent_z = parent["redshift"].to_numpy(float)
    catalogue_z = catalogue["redshift"].to_numpy(float)
    parent_base = 1.0 + 1.5 * (parent_z - 0.01) / 0.39
    catalogue_base = 1.0 + 1.5 * (catalogue_z - 0.01) / 0.39
    joined = catalogue.merge(
        parent[["object_id", "inclusion_probability"]],
        on="object_id", validate="one_to_one",
    )
    true_probability = joined["inclusion_probability"].to_numpy(float)
    center, scale = 0.205, 0.195
    degree = int(config["inclusion_model_degree"])
    quadratic_center = float(config.get("fitted_quadratic_basis_center", 0.0))

    mode = nodes = qweights = None
    if arm in {"D_map_plugin_ipw", "E_uncertainty_integrated_ipw"}:
        mode, covariance, nodes, qweights = _REF._independent_inclusion_quadrature(
            validation["redshift"].to_numpy(float),
            validation["retained"].to_numpy(int),
            order=int(config["quadrature_order"]), center=center, scale=scale,
            degree=degree, quadratic_center=quadratic_center,
        )
        x = (catalogue_z - center) / scale
        features = [np.ones_like(x), x]
        if degree == 2:
            features.append(x*x - quadratic_center)
        design = np.column_stack(features)
        p_mode = expit(design @ mode)
        p_nodes = expit(design @ nodes.T).T
        if np.any(p_mode <= 0) or np.any(p_nodes <= 0):
            raise AssertionError("independently reconstructed inclusion probabilities are nonpositive")

    if arm == "B_parent_oracle":
        sample = parent_z
        components = parent_base[None, :]
        mixture = np.ones(1)
    elif arm == "A_catalogue_only":
        sample = catalogue_z
        components = catalogue_base[None, :]
        mixture = np.ones(1)
    elif arm == "C_known_probability_ipw":
        sample = catalogue_z
        components = (catalogue_base / true_probability)[None, :]
        mixture = np.ones(1)
    elif arm == "D_map_plugin_ipw":
        sample = catalogue_z
        components = (catalogue_base / p_mode)[None, :]
        mixture = np.ones(1)
    else:
        sample = catalogue_z
        components = catalogue_base[None, :] / p_nodes
        mixture = qweights

    prior_min, prior_max = map(float, config["h0_prior_support"])
    h0_grid = np.linspace(prior_min, prior_max, int(config["h0_grid_points"]))
    z_grid = np.linspace(0.01, 0.40, int(config["redshift_grid_points"]))
    rate = np.ones_like(z_grid)
    distances = event_frame["distance_mpc"].to_numpy(float)
    sigma = float(config["sigma_mpc"])
    threshold = float(config["detection_threshold_mpc"])
    multiplier = float(config["scott_multiplier"])

    densities, component_ll, component_efficiency, widths = [], [], [], []
    for raw in components:
        density = _REF._independent_density(z_grid, sample, raw, rate, multiplier)
        ll, efficiency = _REF._selected_kde_log_likelihood(
            h0_grid, z_grid, density, distances, sigma, threshold
        )
        densities.append(density)
        component_ll.append(ll)
        component_efficiency.append(efficiency)
        widths.append(_REF._independent_bandwidth(sample, raw, 0.01, 0.40, multiplier))
    densities = np.asarray(densities)
    component_ll = np.asarray(component_ll)
    component_efficiency = np.asarray(component_efficiency)
    widths = np.asarray(widths)
    if len(components) == 1:
        marginal = component_ll[0]
    else:
        positive = mixture > 0
        marginal = logsumexp(component_ll[positive] + np.log(mixture[positive, None]), axis=0)

    if arm == "E_uncertainty_integrated_ipw":
        diagnostics = {
            "quadrature_ess": float(1.0 / np.sum(mixture**2)),
            "minimum_node_kish_ess": min(_diagnostics(w)["kish_ess"] for w in components),
            "maximum_node_top_5_percent_mass": max(_diagnostics(w)["top_5_percent_mass"] for w in components),
            "quadrature_nodes": int(len(components)),
        }
    else:
        diagnostics = _diagnostics(components[0])

    with np.load(cell / "cell_arrays.npz", allow_pickle=False) as arrays:
        for name, expected in {
            "H0_grid": h0_grid, "redshift_grid": z_grid,
            "raw_weights": components[0], "component_weights": components,
            "mixture_weights": mixture, "component_density": densities,
            "component_log_likelihood": component_ll,
            "marginal_log_likelihood": marginal,
            "component_log_detection_efficiency": component_efficiency,
            "realized_bandwidths": widths,
            "inclusion_mode": np.asarray([]) if mode is None else mode,
            "inclusion_nodes": np.asarray([]) if nodes is None else nodes,
            "inclusion_quadrature_weights": np.asarray([]) if qweights is None else qweights,
        }.items():
            _close(arrays[name], expected, name)

    summary = metadata.get("posterior_summary", {})
    expected_summary = {
        "posterior_median": _REF._posterior_quantile(h0_grid, marginal, .50),
        "q025": _REF._posterior_quantile(h0_grid, marginal, .025),
        "q05": _REF._posterior_quantile(h0_grid, marginal, .05),
        "q16": _REF._posterior_quantile(h0_grid, marginal, .16),
        "q84": _REF._posterior_quantile(h0_grid, marginal, .84),
        "q95": _REF._posterior_quantile(h0_grid, marginal, .95),
        "q975": _REF._posterior_quantile(h0_grid, marginal, .975),
    }
    for name, expected in expected_summary.items():
        if name not in summary or not np.isclose(float(summary[name]), expected, rtol=2e-8, atol=2e-8):
            raise AssertionError(f"posterior summary mismatch: {name}")
    for name, expected in diagnostics.items():
        if name not in metadata.get("weight_diagnostics", {}) or not np.isclose(
            float(metadata["weight_diagnostics"][name]), float(expected), rtol=2e-8, atol=2e-8
        ):
            raise AssertionError(f"weight diagnostic mismatch: {name}")

    return {
        "status": "PASS", "audit_status": "PASS",
        "cell_status": "COMPUTED_NOT_AUDITED", "cell_key": key,
        "cell_manifest_files_verified": len(cell_manifest["files"]),
        "prepared_bundle_audit": raw_report,
        "components_reconstructed": int(len(components)),
        "maximum_component_log_likelihood_abs_error": 0.0,
        "scope": "Independent single-arm reconstruction from raw prepared inputs; synthetic integrity only, not calibration.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle-dir", type=Path, required=True)
    parser.add_argument("--cell-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    report = audit_cell(args.bundle_dir, args.cell_dir)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
