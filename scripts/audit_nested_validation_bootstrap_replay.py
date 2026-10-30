"""Independent end-to-end replay audit for nested validation-bootstrap artifacts.

Regenerates each synthetic scenario from its persisted seeds using a separately
written simulation reconstruction, rebuilds the finite-parent target and outer
SRSWOR estimates with the independent posterior implementation, and delegates
inner-attempt replay to the separately implemented posterior-bootstrap auditor.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit

from dark_siren_robustness.posterior_validation_bootstrap_audit import (
    _independent_posterior,
    audit_posterior_validation_bootstrap,
)
from dark_siren_robustness.posterior_validation_bootstrap import (
    bootstrap_validation_h0_posterior,
)

C_KM_S = 299_792.458


def _require(condition: bool, message: str) -> None:
    if not bool(condition):
        raise AssertionError(message)


def _sha_float64(values) -> str:
    arr = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
    digest = hashlib.sha256()
    digest.update(np.asarray(arr.shape, dtype=np.int64).tobytes())
    digest.update(arr.tobytes())
    return digest.hexdigest()


def _reconstruct_generator(seed: int, n_parent: int, n_detected: int,
                           h0_true: float, sigma: float, threshold: float,
                           z_min: float, z_max: float) -> tuple[np.ndarray, np.ndarray]:
    """Independent reproduction of the generator's parent-z and detected-distance draws."""
    rng = np.random.default_rng(int(seed))
    grid = np.linspace(z_min, z_max, 20001)
    intensity = (
        1.0
        + 2.2 * np.exp(-0.5 * ((grid - 0.105) / 0.025) ** 2)
        + 1.6 * np.exp(-0.5 * ((grid - 0.275) / 0.040) ** 2)
    )
    increments = 0.5 * (intensity[:-1] + intensity[1:]) * np.diff(grid)
    cdf = np.r_[0.0, np.cumsum(increments)]
    cdf /= cdf[-1]
    parent_z = np.interp(rng.random(n_parent), cdf, grid)

    # The upstream generator draws catalogue retention before GW hosts. It is
    # intentionally consumed here even though this benchmark uses its own
    # independently seeded validation-retention mechanism.
    p_catalogue = 0.10 + 0.90 / (1.0 + np.exp((parent_z - 0.22) / 0.05))
    _ = rng.random(n_parent) < p_catalogue

    host_weight = 1.0 + 1.5 * (parent_z - z_min) / (z_max - z_min)
    host_probability = host_weight / host_weight.sum()
    accepted_distances = []
    n_accepted = 0
    attempts = 0
    max_attempts = max(100_000, n_detected * 2_000)
    while n_accepted < n_detected and attempts < max_attempts:
        batch = min(max(1000, 2 * (n_detected - n_accepted)), max_attempts - attempts)
        host_ids = rng.choice(n_parent, size=batch, replace=True, p=host_probability)
        measured = rng.normal(C_KM_S * parent_z[host_ids] / h0_true, sigma)
        detected = measured >= threshold
        if np.any(detected):
            accepted = measured[detected]
            take = min(n_detected - n_accepted, len(accepted))
            accepted_distances.append(accepted[:take])
            n_accepted += take
        attempts += batch
    _require(n_accepted == n_detected, "independent generator failed to reproduce detected event count")
    return parent_z, np.concatenate(accepted_distances)


def _seed(root_seed: int, *keys: int) -> int:
    return int(np.random.SeedSequence([int(root_seed), *map(int, keys)]).generate_state(
        1, dtype=np.uint32
    )[0])


def _reconstruct_result(row, method: str, n_parent: int, b: int) -> dict:
    mask = np.asarray(row.bootstrap_success_mask, dtype=bool)
    means = np.asarray(row.bootstrap_attempt_posterior_mean, dtype=float)
    successes = int(mask.sum())
    success_means = means[mask]
    raw_se = float(np.std(success_means, ddof=1)) if successes >= 2 else float("nan")
    interval = (
        np.quantile(success_means, [0.025, 0.975]).tolist()
        if successes >= 2 else [float("nan"), float("nan")]
    )
    fpc = float(np.sqrt(1.0 - int(row.validation_size) / n_parent))
    failures = row.bootstrap_failure_records
    return {
        "schema": "validation_h0_posterior_resampling_v2",
        "method": method,
        "design_aware_method_status": (
            "exploratory_pseudopopulation_resampling_not_studentized_inference"
            if method == "design_aware_srswor" else "ordinary_paired_record_bootstrap"
        ),
        "seed": int(row.bootstrap_seed),
        "n_attempted": int(row.n_bootstrap_attempted),
        "n_successful": int(row.n_bootstrap_successful),
        "n_failed": int(row.n_bootstrap_attempted) - int(row.n_bootstrap_successful),
        "success_mask": mask,
        "attempt_posterior_mean": means,
        "reference_posterior_mean_by_attempt": np.asarray(
            row.bootstrap_reference_posterior_mean_by_attempt, dtype=float
        ),
        "posterior_median_draws_by_attempt": np.asarray(
            row.bootstrap_posterior_median_draws_by_attempt, dtype=float
        ),
        "ci68_low_draws_by_attempt": np.asarray(
            row.bootstrap_ci68_low_draws_by_attempt, dtype=float
        ),
        "ci68_high_draws_by_attempt": np.asarray(
            row.bootstrap_ci68_high_draws_by_attempt, dtype=float
        ),
        "inverse_weight_ess_by_attempt": np.asarray(
            row.bootstrap_inverse_weight_ess_by_attempt, dtype=float
        ),
        "inverse_weight_maximum_by_attempt": np.asarray(
            row.bootstrap_inverse_weight_maximum_by_attempt, dtype=float
        ),
        "coefficient_boundary_limited_by_attempt": np.asarray(
            row.bootstrap_coefficient_boundary_limited_by_attempt, dtype=bool
        ),
        "posterior_density_draws": np.asarray(row.bootstrap_posterior_density_draws, dtype=float),
        "sample_indices": row.bootstrap_sample_indices,
        "pseudo_population_indices": row.bootstrap_pseudo_population_indices,
        "selected_pseudo_population_positions": row.bootstrap_selected_pseudo_population_positions,
        "failure_records": failures,
        "input_fingerprints": row.bootstrap_input_fingerprints,
        "fixed_full_validation_reference_mean": float(
            row.bootstrap_fixed_full_validation_reference_mean
        ),
        "raw_bootstrap_standard_error": raw_se,
        "fpc_factor": fpc,
        "fpc_adjusted_standard_error_diagnostic_only": raw_se * fpc,
        "posterior_mean_percentile_interval_95": interval,
    }


def audit(directory: Path) -> dict:
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    records = pd.read_csv(directory / "calibration_records.csv")
    attempts = pd.read_json(directory / "attempt_audit.jsonl", lines=True)
    # Status lives in the calibration table; join it explicitly rather than
    # assuming the attempt-level trace duplicates the row-level outcome.
    key = ["scenario", "validation_size", "outer_draw", "method"]
    attempts = attempts.merge(
        records[key + ["bootstrap_status"]],
        on=key, how="left", validate="one_to_one",
    )
    _require(attempts["bootstrap_status"].notna().all(),
             "attempt trace is missing calibration status rows")
    params = manifest["model_parameters"]
    n_parent = int(manifest["parent_size"])
    n_detected = int(manifest["detected_event_count"])
    root_seed = int(manifest["seed"])
    sigma = float(params["distance_noise_sigma_mpc"])
    threshold = float(params["detection_threshold_mpc"])
    h0_true = float(params["h0_true_km_s_mpc"])
    z_min = float(params["redshift_grid"]["start"])
    z_max = float(params["redshift_grid"]["stop"])
    z_grid = np.linspace(z_min, z_max, int(params["redshift_grid"]["points"]))
    h0_grid = np.linspace(
        float(params["h0_grid"]["start"]), float(params["h0_grid"]["stop"]),
        int(params["h0_grid"]["points"]),
    )
    rate = 1.0 + 1.5 * (z_grid-z_min)/(z_max-z_min)
    midpoint, half_range = (z_min+z_max)/2, (z_max-z_min)/2
    scenarios = int(manifest["scenarios"])
    draws = int(manifest["outer_draws"])
    sizes = tuple(sorted(int(v) for v in manifest["validation_sizes"]))
    b = int(manifest["bootstrap_replicates_per_method_per_outer_draw"])
    audited_outer = 0
    audited_bootstrap = 0
    audited_wrapper_failures = 0

    for scenario in range(scenarios):
        scenario_rows = attempts.loc[attempts.scenario.eq(scenario)]
        first = scenario_rows.iloc[0]
        root = np.random.SeedSequence([root_seed, scenario])
        sim_ss, thin_ss, outer_ss = root.spawn(3)
        sim_seed = int(sim_ss.generate_state(1, dtype=np.uint32)[0])
        thin_seed = int(thin_ss.generate_state(1, dtype=np.uint32)[0])
        outer_seed = int(outer_ss.generate_state(1, dtype=np.uint32)[0])
        _require(int(first.simulation_seed) == sim_seed, f"simulation seed mismatch in scenario {scenario}")
        _require(int(first.thinning_seed) == thin_seed, f"thinning seed mismatch in scenario {scenario}")
        _require(int(first.outer_seed) == outer_seed, f"outer seed mismatch in scenario {scenario}")

        parent_z, distances = _reconstruct_generator(
            sim_seed, n_parent, n_detected, h0_true, sigma, threshold, z_min, z_max
        )
        x = (parent_z-midpoint)/half_range
        p_true = expit(0.25-1.2*x)
        retained = np.random.default_rng(thin_seed).random(n_parent) < p_true
        observed_z = parent_z[retained]
        _require(_sha_float64(parent_z) == first.parent_frame_sha256,
                 f"parent-frame fingerprint mismatch in scenario {scenario}")
        _require(_sha_float64(observed_z) == first.catalogue_sha256,
                 f"catalogue fingerprint mismatch in scenario {scenario}")
        _require(_sha_float64(distances) == first.gw_distance_sha256,
                 f"GW-distance fingerprint mismatch in scenario {scenario}")
        _require(observed_z.size >= 5 and np.unique(observed_z).size >= 3
                 and np.unique(retained).size == 2,
                 f"scenario {scenario} violates finite-parent target preconditions")

        _, full_summary, _ = _independent_posterior(
            z_grid, h0_grid, observed_z, parent_z, retained.astype(int),
            distances, sigma, threshold, rate, None,
        )
        target = float(full_summary["posterior_mean"])
        scenario_records = records.loc[records.scenario.eq(scenario)]
        _require(np.allclose(scenario_records.finite_parent_target.to_numpy(float), target,
                             rtol=5e-6, atol=5e-7),
                 f"finite-parent target does not independently reconstruct in scenario {scenario}")

        outer_rng = np.random.default_rng(outer_seed)
        for draw in range(draws):
            permutation = outer_rng.permutation(n_parent)
            for size in sizes:
                expected_indices = permutation[:size]
                draw_rows = attempts.loc[
                    attempts.scenario.eq(scenario)
                    & attempts.outer_draw.eq(draw)
                    & attempts.validation_size.eq(size)
                ]
                _require(len(draw_rows) == 3, "outer draw must have exactly three method rows")
                for row in draw_rows.itertuples(index=False):
                    _require(np.array_equal(np.asarray(row.validation_indices, dtype=int), expected_indices),
                             f"outer SRSWOR indices fail independent replay: scenario={scenario}, draw={draw}, n={size}")
                vz, vy = parent_z[expected_indices], retained[expected_indices].astype(int)
                record_rows = records.loc[
                    records.scenario.eq(scenario)
                    & records.outer_draw.eq(draw)
                    & records.validation_size.eq(size)
                ]
                status = str(record_rows.outer_status.iloc[0])
                if np.unique(vy).size != 2:
                    _require(status == "failed", "one-class outer sample was not labelled failed")
                    audited_outer += 1
                    continue
                try:
                    _, outer_summary, _ = _independent_posterior(
                        z_grid, h0_grid, observed_z, vz, vy, distances,
                        sigma, threshold, rate, None,
                    )
                    rebuilt_outer_mean = float(outer_summary["posterior_mean"])
                except (ValueError, RuntimeError, ArithmeticError, FloatingPointError):
                    _require(status == "failed", "independently failed outer fit was labelled successful")
                    audited_outer += 1
                    continue
                _require(status == "success", "independently successful outer fit was labelled failed")
                stored_outer = float(record_rows.outer_posterior_mean.iloc[0])
                _require(np.isclose(stored_outer, rebuilt_outer_mean, rtol=5e-6, atol=5e-7),
                         "outer posterior mean mismatch")
                audited_outer += 1

                for method in ("paired_bootstrap", "design_aware_srswor"):
                    row = draw_rows.loc[draw_rows.method.eq(method)].iloc[0]
                    expected_seed = _seed(root_seed, scenario, draw, size,
                                          1 if method == "paired_bootstrap" else 2)
                    _require(int(row.bootstrap_seed) == expected_seed,
                             f"bootstrap seed mismatch for {method}, scenario={scenario}, draw={draw}")
                    if str(row.bootstrap_status) == "failed":
                        # Re-run the full wrapper with the same input arrays and
                        # deterministic seed, then compare the exact failure.
                        observed_class = str(row.bootstrap_wrapper_failure_class)
                        observed_message = str(row.bootstrap_wrapper_failure_message)
                        try:
                            replayed = bootstrap_validation_h0_posterior(
                                z_grid, h0_grid, observed_z, vz, vy, distances,
                                sigma_mpc=sigma,
                                detection_threshold_mpc=threshold,
                                parent_population_size=n_parent,
                                method=method, n_bootstrap=b,
                                seed=expected_seed, merger_weight=rate,
                            )
                            if int(replayed["n_successful"]) < 2:
                                raise RuntimeError(
                                    "fewer than two successful bootstrap replicates "
                                    f"({int(replayed['n_successful'])}/{b})"
                                )
                        except (ValueError, RuntimeError, ArithmeticError, FloatingPointError) as exc:
                            _require(type(exc).__name__ == observed_class,
                                     f"wrapper failure class mismatch for {method}: "
                                     f"{type(exc).__name__} != {observed_class}")
                            _require(str(exc) == observed_message,
                                     f"wrapper failure message mismatch for {method}")
                            audited_wrapper_failures += 1
                        else:
                            raise AssertionError(
                                f"recorded wrapper failure did not reproduce for {method}, "
                                f"scenario={scenario}, draw={draw}"
                            )
                        continue
                    if str(row.bootstrap_status) != "success":
                        continue
                    payload = _reconstruct_result(row, method, n_parent, b)
                    audit_posterior_validation_bootstrap(
                        payload, z_grid, h0_grid, observed_z, vz, vy, distances,
                        sigma_mpc=sigma, detection_threshold_mpc=threshold,
                        parent_population_size=n_parent, merger_weight=rate,
                        bandwidth_factor=None,
                    )
                    audited_bootstrap += 1

    return {
        "status": "INDEPENDENT_SCENARIO_AND_POSTERIOR_REPLAY_PASS",
        "scenarios_reconstructed": scenarios,
        "outer_fits_reconstructed": audited_outer,
        "successful_bootstrap_methods_replayed": audited_bootstrap,
        "wrapper_level_bootstrap_failures_replayed": audited_wrapper_failures,
        "scope": "synthetic parent/GW generation, catalogue thinning, finite-parent target, outer SRSWOR, successful inner posterior resampling",
        "limitations": "confirmatory-scale run not established; unexpected non-determinism remains possible across dependency/platform versions",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.artifact_directory), indent=2))


if __name__ == "__main__":
    main()
