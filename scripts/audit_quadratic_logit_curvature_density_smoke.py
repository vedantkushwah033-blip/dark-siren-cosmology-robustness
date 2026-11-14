"""Independent raw-input auditor for the quadratic-logit density smoke artifact."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dark_siren_robustness.generative_selection import (
    generate_quadratic_logit_inclusion_samples,
    latent_galaxy_intensity,
)
from scipy.integrate import trapezoid
from scipy.special import expit

# Reuse the separately implemented numerical reference routines used by the
# floor-family raw-input auditor; these do not call the production fitter/KDE.
_reference_path = Path(__file__).with_name("audit_absolute_bandwidth_density_smoke.py")
_reference_spec = importlib.util.spec_from_file_location("independent_floor_density_reference", _reference_path)
_reference = importlib.util.module_from_spec(_reference_spec)
assert _reference_spec.loader is not None
_reference_spec.loader.exec_module(_reference)
_independent_logistic_quadrature = _reference._independent_logistic_quadrature
_independent_weight_diagnostics = _reference._independent_weight_diagnostics
_independent_scott_width = _reference._independent_scott_width
_independent_kde = _reference._independent_kde

Z_MIN, Z_MAX = 0.01, 0.40
BOUNDARY_LOWER, BOUNDARY_UPPER = 0.049, 0.361


def _hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _assert_frame_close(actual: pd.DataFrame, expected: pd.DataFrame, label: str) -> None:
    if list(actual.columns) != list(expected.columns) or len(actual) != len(expected):
        raise AssertionError(f"{label} schema or row-count mismatch")
    for column in actual.columns:
        left = actual[column].to_numpy()
        right = expected[column].to_numpy()
        if np.issubdtype(left.dtype, np.number) and np.issubdtype(right.dtype, np.number):
            if not np.allclose(left.astype(float), right.astype(float), rtol=1e-13, atol=1e-14, equal_nan=True):
                raise AssertionError(f"{label} numerical mismatch: {column}")
        elif left.tolist() != right.tolist():
            raise AssertionError(f"{label} value mismatch: {column}")


def _metric(density: np.ndarray, truth: np.ndarray, grid: np.ndarray) -> dict[str, float]:
    delta = np.asarray(density, dtype=float) - truth
    left = grid <= BOUNDARY_LOWER
    right = grid >= BOUNDARY_UPPER
    return {
        "iae": float(trapezoid(np.abs(delta), grid)),
        "ise": float(trapezoid(delta**2, grid)),
        "boundary_iae": float(
            trapezoid(np.abs(delta[left]), grid[left])
            + trapezoid(np.abs(delta[right]), grid[right])
        ),
    }


def audit(root: Path) -> dict[str, object]:
    root = Path(root).resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    checked = 0
    for relative, metadata in manifest["files"].items():
        path = root / relative
        if (
            not path.is_file()
            or path.stat().st_size != int(metadata["bytes"])
            or _hash(path) != metadata["sha256"]
        ):
            raise AssertionError(f"manifest verification failed: {relative}")
        checked += 1
    if checked != int(manifest["file_count"]):
        raise AssertionError("manifest file count mismatch")

    cfg = json.loads((root / "config.json").read_text())
    if cfg.get("experiment") != "quadratic_logit_curvature_five_arm_density_smoke":
        raise AssertionError("unexpected experiment identifier")
    if cfg.get("status") != "ENGINEERING_SMOKE_NOT_SELECTION_NOT_CONFIRMATORY":
        raise AssertionError("unexpected experiment claim status")
    if cfg.get("analysis_model") != "degree_one_logistic":
        raise AssertionError("analysis model mismatch")
    if cfg.get("arm_C") != "true quadratic inclusion probabilities; oracle diagnostic only":
        raise AssertionError("Arm C oracle labeling missing")

    regenerated = generate_quadratic_logit_inclusion_samples(
        population_seed=int(cfg["population_seed"]),
        validation_seed=int(cfg["validation_seed"]),
        n_parent=int(cfg["n_parent"]),
        n_validation_frame=int(cfg["validation_frame"]),
        n_validation=max(cfg["validation_sizes"]),
        gamma=float(cfg["gamma"]),
        beta=float(cfg["beta"]),
        z_min=Z_MIN,
        z_max=Z_MAX,
    )
    parent = pd.read_csv(root / "truth" / "parent_truth.csv")
    catalogue = pd.read_csv(root / "observable" / "catalogue_observed.csv")
    validation_max = pd.read_csv(root / "observable" / "validation_max.csv")
    validation_truth = pd.read_csv(root / "truth" / "validation_frame_truth.csv")
    _assert_frame_close(parent, regenerated["parent_truth"], "regenerated parent truth")
    _assert_frame_close(catalogue, regenerated["catalogue_observed"], "regenerated observable catalogue")
    _assert_frame_close(validation_max, regenerated["validation_observed"], "regenerated observable validation")
    _assert_frame_close(validation_truth, regenerated["validation_frame_truth"], "regenerated validation truth")

    pz = parent["redshift"].to_numpy(float)
    cz = catalogue["redshift"].to_numpy(float)
    p_catalogue = parent.loc[parent["retained"].astype(bool), "inclusion_probability"].to_numpy(float)
    if len(p_catalogue) != len(catalogue):
        raise AssertionError("truth probability/catalogue row alignment mismatch")
    parent_base = 1.0 + 1.5 * (pz - Z_MIN) / (Z_MAX - Z_MIN)
    catalogue_base = 1.0 + 1.5 * (cz - Z_MIN) / (Z_MAX - Z_MIN)
    grid = np.linspace(Z_MIN, Z_MAX, int(cfg["grid_points"]))
    merger = 1.0 + 1.5 * (grid - Z_MIN) / (Z_MAX - Z_MIN)
    truth_raw = latent_galaxy_intensity(grid) * merger
    truth_density = truth_raw / float(trapezoid(truth_raw, grid))
    widths = [float(x) for x in cfg["absolute_bandwidth_candidates"]]

    expected_metrics: list[dict[str, object]] = []
    expected_diagnostics: list[dict[str, object]] = []
    expected_scott: list[dict[str, object]] = []
    for n in cfg["validation_sizes"]:
        validation = pd.read_csv(root / "observable" / f"validation_n{int(n)}.csv")
        expected_prefix = validation_max.iloc[:int(n)].reset_index(drop=True)
        _assert_frame_close(validation, expected_prefix, f"nested validation prefix n={n}")
        mode, posterior_nodes, posterior_weights = _independent_logistic_quadrature(
            validation["redshift"].to_numpy(float),
            validation["retained"].to_numpy(int),
            cfg,
        )
        archive = np.load(root / "results" / "density_reconstruction_inputs.npz")
        stored_mode = archive[f"posterior_mode_n{int(n)}"]
        stored_nodes = archive[f"posterior_nodes_n{int(n)}"]
        stored_q = archive[f"posterior_weights_n{int(n)}"]
        if not np.allclose(stored_mode, mode, rtol=2e-7, atol=2e-7):
            raise AssertionError("stored inclusion posterior mode differs from independent refit")
        if not np.allclose(stored_nodes, posterior_nodes, rtol=2e-7, atol=2e-7):
            raise AssertionError("stored inclusion posterior nodes differ from independent refit")
        if not np.allclose(stored_q, posterior_weights, rtol=2e-7, atol=2e-7):
            raise AssertionError("stored inclusion posterior weights differ from independent refit")
        xcat = (cz - float(cfg["coordinate_center"])) / float(cfg["coordinate_scale"])
        design_cat = np.column_stack((np.ones(cz.size), xcat))
        p_map = expit(design_cat @ mode)
        p_nodes = expit(design_cat @ posterior_nodes.T).T
        node_weights = catalogue_base[None, :] / p_nodes
        raw_inputs = {
            "A_catalogue_only": (cz, catalogue_base),
            "B_parent_oracle": (pz, parent_base),
            "C_known_probability_ipw": (cz, catalogue_base / p_catalogue),
            "D_map_plugin_ipw": (cz, catalogue_base / p_map),
        }
        node_diag = [
            _independent_weight_diagnostics(w)
            for w in node_weights
        ]
        for name, (sample_z, weights) in raw_inputs.items():
            expected_diagnostics.append({
                "validation_size": int(n), "arm": name, "node": -1,
                **_independent_weight_diagnostics(weights),
            })
            expected_scott.append({
                "validation_size": int(n), "arm": name, "node": -1,
                "realized_bandwidth_redshift": _independent_scott_width(sample_z, weights, float(cfg["scott_multiplier"])),
            })
        for j, weights in enumerate(node_weights):
            expected_diagnostics.append({
                "validation_size": int(n), "arm": "E_uncertainty_integrated_ipw",
                "node": j, **node_diag[j],
            })
            expected_scott.append({
                "validation_size": int(n), "arm": "E_uncertainty_integrated_ipw",
                "node": j,
                "realized_bandwidth_redshift": _independent_scott_width(cz, weights, float(cfg["scott_multiplier"])),
            })
        expected_diagnostics.append({
            "validation_size": int(n), "arm": "E_uncertainty_integrated_ipw",
            "node": -1,
            "kish_ess": min(x["kish_ess"] for x in node_diag),
            "kish_ess_fraction": min(x["kish_ess_fraction"] for x in node_diag),
            "maximum_normalized_weight": max(x["maximum_normalized_weight"] for x in node_diag),
            "top_1_percent_mass": max(x["top_1_percent_mass"] for x in node_diag),
            "top_5_percent_mass": max(x["top_5_percent_mass"] for x in node_diag),
        })

        for width in widths:
            for name, (sample_z, weights) in raw_inputs.items():
                density = _independent_kde(grid, sample_z, weights, merger, width)
                expected_metrics.append({
                    "validation_size": int(n), "absolute_bandwidth_sd": width,
                    "arm": name, "node": -1, **_metric(density, truth_density, grid),
                    "status": "SUCCESS",
                    "kish_ess": _independent_weight_diagnostics(weights)["kish_ess"],
                })
            node_densities = np.asarray([
                _independent_kde(grid, cz, weights, merger, width)
                for weights in node_weights
            ])
            q = posterior_weights
            for j, weights in enumerate(node_weights):
                expected_metrics.append({
                    "validation_size": int(n), "absolute_bandwidth_sd": width,
                    "arm": "E_uncertainty_integrated_ipw", "node": j,
                    **_metric(node_densities[j], truth_density, grid),
                    "status": "SUCCESS", "kish_ess": node_diag[j]["kish_ess"],
                })
            mean_density = np.tensordot(q, node_densities, axes=(0, 0))
            expected_metrics.append({
                "validation_size": int(n), "absolute_bandwidth_sd": width,
                "arm": "E_posterior_mean_density", "node": -1,
                **_metric(mean_density, truth_density, grid),
                "status": "SUCCESS", "kish_ess": float("nan"),
            })
            node_metrics = [_metric(d, truth_density, grid) for d in node_densities]
            expected_metrics.append({
                "validation_size": int(n), "absolute_bandwidth_sd": width,
                "arm": "E_posterior_weighted_node_metric", "node": -1,
                **{
                    key: float(np.dot(q, [m[key] for m in node_metrics]))
                    for key in ("iae", "ise", "boundary_iae")
                },
                "status": "SUCCESS", "kish_ess": float("nan"),
            })

    actual_metrics = pd.read_csv(root / "results" / "density_metrics.csv")
    expected_metrics_df = pd.DataFrame(expected_metrics)
    metric_keys = ["validation_size", "absolute_bandwidth_sd", "arm", "node"]
    metric_values = ["iae", "ise", "boundary_iae", "kish_ess"]
    _compare_frames(actual_metrics, expected_metrics_df, metric_keys, metric_values, "density metric")

    actual_diag = pd.read_csv(root / "results" / "weight_diagnostics.csv")
    expected_diag_df = pd.DataFrame(expected_diagnostics)
    diag_keys = ["validation_size", "arm", "node"]
    diag_values = [
        "kish_ess", "kish_ess_fraction", "maximum_normalized_weight",
        "top_1_percent_mass", "top_5_percent_mass",
    ]
    _compare_frames(actual_diag, expected_diag_df, diag_keys, diag_values, "weight diagnostic")

    actual_scott = pd.read_csv(root / "results" / "scott_widths.csv")
    expected_scott_df = pd.DataFrame(expected_scott)
    _compare_frames(
        actual_scott, expected_scott_df,
        ["validation_size", "arm", "node"],
        ["realized_bandwidth_redshift"], "Scott width",
    )
    return {
        "status": "PASS",
        "manifest_files_verified": checked,
        "parent_rows_reconstructed": len(parent),
        "catalogue_rows_reconstructed": len(catalogue),
        "validation_rows_reconstructed": len(validation_max),
        "density_metric_rows_reconstructed": len(expected_metrics_df),
        "weight_diagnostic_rows_reconstructed": len(expected_diag_df),
        "scott_width_rows_reconstructed": len(expected_scott_df),
        "gamma": float(cfg["gamma"]),
        "scope": "Seed-regenerated raw-input audit plus independent reconstruction of five-arm density metrics, weight diagnostics, and Scott widths; no bandwidth selection or H0 inference.",
    }


def _compare_frames(
    actual: pd.DataFrame,
    expected: pd.DataFrame,
    keys: list[str],
    values: list[str],
    label: str,
) -> None:
    if len(actual) != len(expected):
        raise AssertionError(f"{label} row count mismatch")
    if actual.duplicated(keys).any() or expected.duplicated(keys).any():
        raise AssertionError(f"{label} duplicate key")
    actual = actual.sort_values(keys).reset_index(drop=True)
    expected = expected.sort_values(keys).reset_index(drop=True)
    for key in keys:
        if actual[key].tolist() != expected[key].tolist():
            raise AssertionError(f"{label} key mismatch: {key}")
    for column in values:
        if not np.allclose(
            actual[column].to_numpy(float),
            expected[column].to_numpy(float),
            rtol=2e-10, atol=2e-10, equal_nan=True,
        ):
            raise AssertionError(f"{label} numerical mismatch: {column}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.input_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
