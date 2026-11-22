"""Independently cross-audit one arm-cell artifact against a frozen bundle.

The audit first invokes the standalone full-bundle auditor, which independently
reconstructs the five-arm weights, KDEs, selected likelihoods, posterior
summaries, and bandwidths without importing the production inference path.
It then verifies that the cell artifact agrees with the independently audited
arm/validation-size slice. It never imports the arm-cell evaluator.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

import importlib.util

_audit_path = Path(__file__).with_name("audit_quadratic_five_arm_smoke.py")
_audit_spec = importlib.util.spec_from_file_location("_independent_five_arm_auditor", _audit_path)
if _audit_spec is None or _audit_spec.loader is None:
    raise ImportError("could not load the independent five-arm artifact auditor")
_audit_module = importlib.util.module_from_spec(_audit_spec)
_audit_spec.loader.exec_module(_audit_module)
audit_bundle = _audit_module.audit


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _close(actual, expected, label: str, *, rtol: float = 2e-8, atol: float = 2e-10) -> None:
    a = np.asarray(actual, dtype=float)
    e = np.asarray(expected, dtype=float)
    if a.shape != e.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(e)):
        raise AssertionError(f"{label} has invalid shape or non-finite values")
    if not np.allclose(a, e, rtol=rtol, atol=atol):
        raise AssertionError(f"{label} differs from independent bundle reconstruction")


def _verify_cell_manifest(cell_dir: Path) -> dict:
    manifest_path = cell_dir / "manifest.json"
    if not manifest_path.is_file():
        raise AssertionError("cell artifact manifest is missing")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1 or manifest.get("status") != "COMPUTED_NOT_AUDITED":
        raise AssertionError("cell artifact manifest schema/status is invalid")
    expected = {"cell_arrays.npz", "cell_result.json"}
    if set(manifest.get("files", {})) != expected:
        raise AssertionError("cell artifact manifest file set is invalid")
    for relative, metadata in manifest["files"].items():
        path = cell_dir / relative
        if not path.is_file():
            raise AssertionError(f"cell artifact is missing: {relative}")
        if path.stat().st_size != int(metadata["bytes"]):
            raise AssertionError(f"cell artifact byte count mismatch: {relative}")
        if _sha256(path) != metadata["sha256"]:
            raise AssertionError(f"cell artifact SHA-256 mismatch: {relative}")
    return manifest


def audit_cell(bundle_dir: Path, cell_dir: Path) -> dict:
    bundle = Path(bundle_dir).resolve()
    cell = Path(cell_dir).resolve()
    cell_manifest = _verify_cell_manifest(cell)
    metadata = json.loads((cell / "cell_result.json").read_text(encoding="utf-8"))
    if metadata.get("schema_version") != 1 or metadata.get("status") != "COMPUTED_NOT_AUDITED":
        raise AssertionError("cell result schema/status is invalid")
    key = metadata.get("cell_key", {})
    required = {"replicate", "gamma", "validation_size", "arm"}
    if set(key) != required:
        raise AssertionError("cell key must contain replicate, gamma, validation_size, and arm")
    arm = key["arm"]
    validation_size = int(key["validation_size"])
    if arm not in {
        "A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
        "D_map_plugin_ipw", "E_uncertainty_integrated_ipw",
    }:
        raise AssertionError("unknown arm in cell key")

    shared_manifest = bundle / "shared_input_manifest.json"
    config_path = bundle / "config.json"
    if metadata.get("bundle_manifest_sha256") != _sha256(shared_manifest):
        raise AssertionError("cell is bound to a different shared-input manifest")
    if metadata.get("config_sha256") != _sha256(config_path):
        raise AssertionError("cell is bound to a different bundle config")

    independent_report = audit_bundle(bundle)
    if independent_report.get("status") != "PASS":
        raise AssertionError("independent full-bundle audit did not PASS")

    config = json.loads(config_path.read_text(encoding="utf-8"))
    if validation_size not in config["validation_sizes"]:
        raise AssertionError("cell validation size is not in the audited bundle plan")
    if not np.isfinite(float(key["gamma"])) or float(key["gamma"]) != float(config["curvature_gamma"]):
        raise AssertionError("cell gamma differs from audited bundle config")

    with np.load(cell / "cell_arrays.npz", allow_pickle=False) as cell_arrays, \
         np.load(bundle / "results/log_likelihoods.npz", allow_pickle=False) as legacy:
        suffix = f"_n{validation_size}"
        expected_ll_key = f"{arm}{suffix}"
        if expected_ll_key not in legacy:
            raise AssertionError(f"independent bundle archive lacks {expected_ll_key}")
        _close(cell_arrays["H0_grid"], legacy["H0_grid"], "H0 grid")
        _close(cell_arrays["redshift_grid"], legacy["redshift_grid"], "redshift grid")
        _close(cell_arrays["marginal_log_likelihood"], legacy[expected_ll_key],
               "marginal selected log likelihood")

        if arm == "E_uncertainty_integrated_ipw":
            expected_density = legacy[f"density_E_nodes{suffix}"]
            expected_component_ll = legacy[f"E_node_log_likelihood{suffix}"]
            expected_efficiency = legacy[f"E_node_log_efficiency{suffix}"]
            _close(cell_arrays["component_log_detection_efficiency"], expected_efficiency,
                   "E component detection efficiencies")
            _close(cell_arrays["inclusion_nodes"], legacy[f"posterior_nodes{suffix}"],
                   "E inclusion posterior nodes")
            _close(cell_arrays["inclusion_quadrature_weights"], legacy[f"posterior_weights{suffix}"],
                   "E inclusion posterior quadrature weights")
        else:
            expected_density = legacy[f"density_{arm}{suffix}"][None, :]
            expected_component_ll = legacy[expected_ll_key][None, :]
            if arm in {"D_map_plugin_ipw"}:
                _close(cell_arrays["inclusion_mode"], legacy[f"posterior_mode{suffix}"],
                       "D inclusion posterior mode")
                _close(cell_arrays["inclusion_nodes"], legacy[f"posterior_nodes{suffix}"],
                       "D inclusion posterior nodes")
                _close(cell_arrays["inclusion_quadrature_weights"], legacy[f"posterior_weights{suffix}"],
                       "D inclusion posterior quadrature weights")
        _close(cell_arrays["component_density"], expected_density, "component KDE density")
        _close(cell_arrays["component_log_likelihood"], expected_component_ll,
               "component selected log likelihood")

    summary = pd.read_csv(bundle / "results/posterior_summaries.csv")
    summary = summary.loc[
        (summary["validation_size"] == validation_size) & (summary["arm"] == arm)
    ]
    if len(summary) != 1 or summary.iloc[0]["status"] != "SUCCESS":
        raise AssertionError("independent bundle posterior summary row is missing or invalid")
    for field, value in metadata["posterior_summary"].items():
        if field not in summary or not np.isclose(
            float(value), float(summary.iloc[0][field]), rtol=2e-8, atol=2e-10
        ):
            raise AssertionError(f"cell posterior summary mismatch: {field}")

    diagnostics = pd.read_csv(bundle / "results/weight_diagnostics.csv")
    diagnostics = diagnostics.loc[
        (diagnostics["validation_size"] == validation_size) & (diagnostics["arm"] == arm)
    ]
    if len(diagnostics) != 1:
        raise AssertionError("independent bundle weight-diagnostics row is missing")
    for field, value in metadata["weight_diagnostics"].items():
        if field not in diagnostics or not np.isclose(
            float(value), float(diagnostics.iloc[0][field]), rtol=2e-8, atol=2e-10
        ):
            raise AssertionError(f"cell weight diagnostic mismatch: {field}")

    bandwidths = pd.read_csv(bundle / "results/realized_bandwidths.csv")
    bandwidths = bandwidths.loc[
        (bandwidths["validation_size"] == validation_size) & (bandwidths["arm"] == arm)
    ].sort_values("node")
    _close(
        metadata["realized_bandwidths"],
        bandwidths["realized_bandwidth_redshift"].to_numpy(float),
        "realized KDE bandwidths",
    )

    return {
        "status": "PASS",
        "cell_status": "COMPUTED_NOT_AUDITED",
        "audit_status": "PASS",
        "cell_key": key,
        "cell_manifest_files_verified": len(cell_manifest["files"]),
        "independent_bundle_audit": independent_report,
        "scope": (
            "Cell output cross-checked against a separate full-bundle reconstruction. "
            "This is not yet the production coordinator or a four-state lifecycle transition."
        ),
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
