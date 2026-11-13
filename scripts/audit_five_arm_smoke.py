"""Independent raw-input and likelihood audit for five-arm smoke artifacts.

This auditor does not import the five-arm runner, weight builder, likelihood
implementation, or H0 posterior normalizer. It reconstructs weights and the
selected likelihood directly from saved observable/truth inputs and persisted
validation-posterior nodes. It audits likelihood propagation, not the posterior
quadrature construction itself.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
import pandas as pd
from scipy.special import expit, log_ndtr, logsumexp
from scipy.stats import gaussian_kde, norm

C_KM_S = 299_792.458
ARM_KEYS = {
    "A_catalogue_only": "A_catalogue_only",
    "B_parent_oracle": "B_parent_oracle",
    "C_known_probability_ipw": "C_known_probability_ipw",
    "D_map_plugin_ipw": "D_map_plugin_ipw",
    "E_uncertainty_integrated_ipw": "E_uncertainty_integrated_ipw",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _assert_close(actual, expected, label: str, *, atol: float = 2e-10) -> None:
    if not np.allclose(actual, expected, rtol=2e-10, atol=atol, equal_nan=False):
        difference = float(np.max(np.abs(np.asarray(actual) - np.asarray(expected))))
        raise AssertionError(f"{label} differs; maximum absolute difference={difference}")


def _independent_density(
    z_grid: np.ndarray, sample: np.ndarray, probability: np.ndarray | None,
    rate: np.ndarray, bandwidth_factor: float,
) -> np.ndarray:
    """Independent reflected weighted KDE reconstruction; no production imports."""
    if probability is None:
        raw = np.ones(sample.size, dtype=float)
    else:
        raw = 1.0 / np.asarray(probability, dtype=float)
    scaled = raw / np.max(raw)
    reflected = np.concatenate((
        sample, 2.0 * z_grid[0] - sample, 2.0 * z_grid[-1] - sample,
    ))
    reflected_weights = np.tile(scaled / 3.0, 3)
    estimator = gaussian_kde(
        reflected, weights=reflected_weights, bw_method=float(bandwidth_factor)
    )
    unnormalized = estimator(z_grid) * rate
    area = float(trapezoid(unnormalized, z_grid))
    if not np.isfinite(area) or area <= 0:
        raise AssertionError("independent KDE has invalid normalization")
    density = unnormalized / area
    if not np.isclose(trapezoid(density, z_grid), 1.0, rtol=1e-10, atol=1e-10):
        raise AssertionError("independent KDE failed unit-integral check")
    return density


def _selected_kde_log_likelihood(
    h0_grid: np.ndarray, z_grid: np.ndarray, density: np.ndarray,
    distances: np.ndarray, sigma: float, threshold: float,
) -> tuple[np.ndarray, np.ndarray]:
    means = C_KM_S * z_grid[None, :] / h0_grid[:, None]
    event_log_mix, event_log_beta = [], []
    for distance in distances:
        numerator = density[None, :] * norm.pdf(distance, loc=means, scale=sigma)
        efficiency = density[None, :] * np.exp(log_ndtr((means - threshold) / sigma))
        event_log_mix.append(np.log(trapezoid(numerator, z_grid, axis=1)))
        event_log_beta.append(np.log(trapezoid(efficiency, z_grid, axis=1)))
    event_log_mix = np.asarray(event_log_mix)
    event_log_beta = np.asarray(event_log_beta)
    return np.sum(event_log_mix - event_log_beta, axis=0), event_log_beta

def _posterior_quantile(grid: np.ndarray, log_likelihood: np.ndarray, q: float) -> float:
    shifted = np.exp(log_likelihood - np.max(log_likelihood))
    area = float(trapezoid(shifted, grid))
    density = shifted / area
    increments = 0.5 * (density[:-1] + density[1:]) * np.diff(grid)
    cdf = np.concatenate(([0.0], np.cumsum(increments)))
    cdf /= cdf[-1]
    return float(np.interp(q, cdf, grid))


def audit(root: Path) -> dict:
    root = root.resolve()
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    checked_files = 0
    for relative, metadata in manifest["files"].items():
        path = root / relative
        if not path.is_file():
            raise AssertionError(f"manifest file is missing: {relative}")
        if path.stat().st_size != int(metadata["bytes"]):
            raise AssertionError(f"file size mismatch: {relative}")
        if _sha256(path) != metadata["sha256"]:
            raise AssertionError(f"SHA-256 mismatch: {relative}")
        checked_files += 1
    if checked_files != int(manifest["file_count"]):
        raise AssertionError("manifest file_count does not match entries")

    config = json.loads((root / "config.json").read_text(encoding="utf-8"))
    if config["status"] != "engineering_smoke_not_confirmatory_not_observational":
        raise AssertionError("artifact claim-status label is missing or unexpected")
    floor = float(config["catalogue_floor"])
    midpoint, width = 0.22, 0.05
    center, scale = 0.205, 0.195
    threshold = float(config["detection_threshold_mpc"])
    sigma = float(config["sigma_mpc"])
    h0_true = float(config["h0_true"])

    parent = pd.read_csv(root / "truth" / "parent_truth.csv")
    catalogue = pd.read_csv(root / "observable" / "catalogue_observed.csv")
    events = pd.read_csv(root / "observable" / "gw_events_observed.csv")
    validation_truth = pd.read_csv(root / "truth" / "validation_frame_truth.csv")
    event_truth = pd.read_csv(root / "truth" / "gw_event_truth.csv")
    generation_trace = pd.read_csv(root / "truth" / "gw_generation_trace_truth.csv")
    selected_trace = generation_trace.loc[
        generation_trace["selected_as_event"] == 1
    ].sort_values("attempt_id")
    if len(events) != int(config["n_detected"]) or len(event_truth) != len(events):
        raise AssertionError("detected-event observable/truth row counts differ")
    if np.any(events["distance_mpc"].to_numpy(float) < threshold):
        raise AssertionError("an observed event fails the declared GW detection threshold")
    if len(catalogue) != int(config["n_catalogue_observed"]):
        raise AssertionError("catalogue row count differs from config")
    expected_catalogue_ids = set(parent.loc[parent["retained"] == 1, "object_id"].astype(int))
    if set(catalogue["object_id"].astype(int)) != expected_catalogue_ids:
        raise AssertionError("observable catalogue does not match retained parent IDs")
    if not np.array_equal(events["event_id"].to_numpy(int), event_truth["event_id"].to_numpy(int)):
        raise AssertionError("GW event observable/truth IDs differ")
    if len(selected_trace) != len(events):
        raise AssertionError("selected generation-trace rows do not match detected-event count")
    _assert_close(
        events["distance_mpc"].to_numpy(float),
        selected_trace["distance_mpc"].to_numpy(float),
        "observed GW distances versus selected generation trace",
    )
    if not np.array_equal(
        event_truth["host_object_id"].to_numpy(int),
        selected_trace["host_object_id"].to_numpy(int),
    ):
        raise AssertionError("GW event host IDs disagree with selected generation trace")
    event_join = event_truth.merge(
        parent[["object_id", "redshift"]],
        left_on="host_object_id", right_on="object_id",
        suffixes=("_event", "_parent"), validate="many_to_one",
    )
    _assert_close(
        event_join["host_redshift"].to_numpy(float),
        event_join["redshift"].to_numpy(float),
        "GW host-redshift truth join",
    )
    _assert_close(
        event_join["latent_mean_distance_mpc"].to_numpy(float),
        C_KM_S * event_join["redshift"].to_numpy(float) / h0_true,
        "GW latent distance truth",
    )
    joined = catalogue.merge(
        parent[["object_id", "redshift", "inclusion_probability"]],
        on="object_id", suffixes=("_observed", "_truth"), validate="one_to_one",
    )
    _assert_close(
        joined["redshift_observed"].to_numpy(float),
        joined["redshift_truth"].to_numpy(float),
        "catalogue redshift join",
    )
    catalogue_z = joined["redshift_observed"].to_numpy(float)
    true_probability = floor + (1.0 - floor) / (
        1.0 + np.exp((catalogue_z - midpoint) / width)
    )
    _assert_close(
        joined["inclusion_probability"].to_numpy(float),
        true_probability,
        "catalogue truth inclusion probabilities",
    )
    parent_z = parent["redshift"].to_numpy(float)
    parent_base = 1.0 + 1.5 * (parent_z - 0.01) / 0.39
    catalogue_base = 1.0 + 1.5 * (catalogue_z - 0.01) / 0.39
    distances = events["distance_mpc"].to_numpy(float)
    grid = np.linspace(50.0, 90.0, int(config["h0_grid_points"]))
    z_grid = np.linspace(0.01, 0.40, int(config["redshift_grid_points"]))
    rate_grid = 1.0 + 1.5 * (z_grid - 0.01) / 0.39
    bandwidth_factor = float(config["bandwidth_factor"])
    archive = np.load(root / "results" / "log_likelihoods.npz")
    _assert_close(archive["H0_grid"], grid, "H0 grid")
    _assert_close(archive["redshift_grid"], z_grid, "redshift grid")
    _assert_close(archive["merger_weight_grid"], rate_grid, "merger weight grid")
    stored_summary = pd.read_csv(root / "results" / "posterior_summaries.csv")

    validation_sizes = [int(value) for value in config["validation_sizes"]]
    validation_ids = {}
    posterior_rows = 0
    max_likelihood_error = 0.0
    for n in validation_sizes:
        validation = pd.read_csv(root / "observable" / f"validation_n{n}.csv")
        validation_ids[n] = validation["validation_id"].to_numpy(int)
        if len(validation) != n:
            raise AssertionError(f"validation size mismatch at n={n}")
        validation_join = validation.merge(
            validation_truth[["frame_id", "redshift", "retained"]],
            left_on="validation_id", right_on="frame_id",
            suffixes=("_observable", "_truth"), validate="one_to_one",
        )
        _assert_close(
            validation_join["redshift_observable"].to_numpy(float),
            validation_join["redshift_truth"].to_numpy(float),
            f"validation redshift join n={n}",
        )
        _assert_close(
            validation_join["retained_observable"].to_numpy(int),
            validation_join["retained_truth"].to_numpy(int),
            f"validation retention join n={n}",
            atol=0.0,
        )
        y = validation["retained"].to_numpy(int)
        z_validation = validation["redshift"].to_numpy(float)
        x_validation = (z_validation - center) / scale
        design_validation = np.column_stack([np.ones(n), x_validation])
        nodes = archive[f"posterior_nodes_n{n}"]
        mode = archive[f"posterior_mode_n{n}"]
        quadrature = archive[f"posterior_weights_n{n}"]
        if nodes.shape[1] != 2 or mode.shape != (2,):
            raise AssertionError(f"unexpected inclusion model dimension at n={n}")
        if not np.isclose(np.sum(quadrature), 1.0) or np.any(quadrature < 0):
            raise AssertionError(f"invalid coefficient quadrature weights at n={n}")
        if not np.all(np.isfinite(design_validation @ mode)):
            raise AssertionError(f"non-finite validation MAP predictor at n={n}")

        x_catalogue = (catalogue_z - center) / scale
        design_catalogue = np.column_stack([np.ones(len(catalogue_z)), x_catalogue])
        p_map = expit(design_catalogue @ mode)
        p_nodes = expit(design_catalogue @ nodes.T).T
        if np.any(p_map <= 0) or np.any(p_nodes <= 0):
            raise AssertionError(f"non-positive estimated inclusion probability at n={n}")

        arm_inputs = {
            "A_catalogue_only": (catalogue_z, None),
            "B_parent_oracle": (parent_z, None),
            "C_known_probability_ipw": (catalogue_z, true_probability),
            "D_map_plugin_ipw": (catalogue_z, p_map),
        }
        reconstructed = {}
        for arm, (z_host, probability) in arm_inputs.items():
            p_host = _independent_density(
                z_grid, z_host, probability, rate_grid, bandwidth_factor
            )
            _assert_close(
                archive[f"density_{arm}_n{n}"], p_host,
                f"{arm} independently reconstructed KDE n={n}",
            )
            reconstructed[arm], _ = _selected_kde_log_likelihood(
                grid, z_grid, p_host, distances, sigma, threshold
            )
        node_likelihoods = []
        node_efficiencies = []
        node_densities = []
        for probability in p_nodes:
            p_host = _independent_density(
                z_grid, catalogue_z, probability, rate_grid, bandwidth_factor
            )
            node_densities.append(p_host)
            ll, efficiency = _selected_kde_log_likelihood(
                grid, z_grid, p_host, distances, sigma, threshold,
            )
            node_likelihoods.append(ll)
            node_efficiencies.append(efficiency)
        _assert_close(
            archive[f"density_E_nodes_n{n}"], np.asarray(node_densities),
            f"E node KDEs n={n}",
        )
        node_likelihoods = np.asarray(node_likelihoods)
        node_efficiencies = np.asarray(node_efficiencies)
        positive = quadrature > 0
        e_log_likelihood = logsumexp(
            node_likelihoods[positive] + np.log(quadrature[positive, None]),
            axis=0,
        )
        reconstructed["E_uncertainty_integrated_ipw"] = e_log_likelihood

        for arm, expected in reconstructed.items():
            key = f"{arm}_n{n}"
            _assert_close(archive[key], expected, f"{arm} log likelihood n={n}")
            max_likelihood_error = max(
                max_likelihood_error,
                float(np.max(np.abs(archive[key] - expected))),
            )
            row = stored_summary.loc[
                (stored_summary["validation_size"] == n)
                & (stored_summary["arm"] == arm)
            ]
            if len(row) != 1 or row.iloc[0]["status"] != "SUCCESS":
                raise AssertionError(f"posterior summary row missing/invalid: {arm}, n={n}")
            expected_quantiles = {
                "posterior_median": _posterior_quantile(grid, expected, 0.50),
                "q025": _posterior_quantile(grid, expected, 0.025),
                "q05": _posterior_quantile(grid, expected, 0.05),
                "q16": _posterior_quantile(grid, expected, 0.16),
                "q84": _posterior_quantile(grid, expected, 0.84),
                "q95": _posterior_quantile(grid, expected, 0.95),
                "q975": _posterior_quantile(grid, expected, 0.975),
            }
            for field, value in expected_quantiles.items():
                if not np.isclose(float(row.iloc[0][field]), value, rtol=2e-9, atol=2e-9):
                    raise AssertionError(f"{arm} n={n}: posterior summary {field} mismatch")
            if not np.isclose(float(row.iloc[0]["median_error"]), expected_quantiles["posterior_median"] - h0_true):
                raise AssertionError(f"{arm} n={n}: median error mismatch")
            posterior_rows += 1

        _assert_close(
            archive[f"E_node_log_likelihood_n{n}"],
            node_likelihoods,
            f"E node likelihoods n={n}",
        )
        _assert_close(
            archive[f"E_node_log_efficiency_n{n}"],
            node_efficiencies,
            f"E node detection efficiencies n={n}",
        )

    for smaller, larger in zip(validation_sizes, validation_sizes[1:]):
        if not np.array_equal(
            validation_ids[smaller], validation_ids[larger][:len(validation_ids[smaller])]
        ):
            raise AssertionError("validation samples are not nested prefixes")

    return {
        "status": "PASS",
        "manifest_files_verified": checked_files,
        "posterior_summary_rows_reconstructed": posterior_rows,
        "maximum_log_likelihood_absolute_difference": max_likelihood_error,
        "validation_sizes_verified": validation_sizes,
        "scope": (
            "Independently reconstructs input joins, inclusion probabilities, "
            "A-D/E weights, selected likelihoods, posterior summaries and file hashes. "
            "Uses stored coefficient posterior nodes; does not independently "
            "validate the logistic posterior/quadrature construction."
        ),
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
