"""Benchmark scalar versus batched event-kernel evaluation at production-like dimensions.

The benchmark isolates likelihood-kernel evaluation on already normalized host
densities. It reports timing and numerical agreement; it does not assert a
machine-independent speed threshold or include KDE density-construction time.
"""
from __future__ import annotations

import argparse
import json
import time

import numpy as np
from scipy.integrate import trapezoid
from scipy.special import log_ndtr
from scipy.stats import norm

from dark_siren_robustness.dataset_likelihood import C_KM_S
from dark_siren_robustness.five_arm_kde_likelihood import _batched_density_likelihood


def _scalar_reference(h0, z, densities, distances, sigmas, threshold):
    all_ll = np.empty((len(densities), len(h0)))
    all_eff = np.empty((len(densities), len(distances), len(h0)))
    means = C_KM_S * z[None, :] / h0[:, None]
    for density_index, density in enumerate(densities):
        for event_index, (distance, sigma) in enumerate(zip(distances, sigmas)):
            numerator = density[None, :] * norm.pdf(
                distance, loc=means, scale=sigma
            )
            efficiency = density[None, :] * np.exp(
                log_ndtr((means - threshold) / sigma)
            )
            log_num = np.log(trapezoid(numerator, z, axis=1))
            log_eff = np.log(trapezoid(efficiency, z, axis=1))
            all_ll[density_index] += 0.0 if event_index else 0.0
            all_eff[density_index, event_index] = log_eff
            if event_index == 0:
                all_ll[density_index] = log_num - log_eff
            else:
                all_ll[density_index] += log_num - log_eff
    return all_ll, all_eff


def benchmark(*, h0_points=321, z_points=801, n_events=30,
              n_densities=23, repeats=3, seed=20261114):
    if min(h0_points, z_points, n_events, n_densities, repeats) < 1:
        raise ValueError("all benchmark dimensions and repeats must be positive")
    rng = np.random.default_rng(seed)
    z = np.linspace(0.01, 0.40, z_points)
    h0 = np.linspace(55.0, 85.0, h0_points)
    densities = []
    for i in range(n_densities):
        center1 = 0.08 + 0.20 * rng.random()
        center2 = 0.12 + 0.24 * rng.random()
        width1 = 0.018 + 0.025 * rng.random()
        width2 = 0.020 + 0.035 * rng.random()
        raw = (
            0.25
            + 0.7 * np.exp(-0.5 * ((z - center1) / width1) ** 2)
            + 0.45 * np.exp(-0.5 * ((z - center2) / width2) ** 2)
        )
        densities.append(raw / float(trapezoid(raw, z)))
    densities = np.asarray(densities)
    distances = np.linspace(270.0, 1500.0, n_events)
    sigmas = np.linspace(120.0, 220.0, n_events)
    threshold = 250.0

    # Warm up both paths before collecting timings.
    reference_ll, reference_eff = _scalar_reference(
        h0, z, densities, distances, sigmas, threshold
    )
    batched_ll, batched_eff = _batched_density_likelihood(
        h0, z, densities, distances, sigmas, threshold
    )
    max_ll_difference = float(np.max(np.abs(reference_ll - batched_ll)))
    max_eff_difference = float(np.max(np.abs(reference_eff - batched_eff)))
    if max(max_ll_difference, max_eff_difference) > 2e-10:
        raise AssertionError("batched likelihood differs from scalar reference")

    scalar_times, batched_times = [], []
    for _ in range(repeats):
        start = time.perf_counter()
        _scalar_reference(h0, z, densities, distances, sigmas, threshold)
        scalar_times.append(time.perf_counter() - start)
        start = time.perf_counter()
        _batched_density_likelihood(
            h0, z, densities, distances, sigmas, threshold
        )
        batched_times.append(time.perf_counter() - start)
    scalar_median = float(np.median(scalar_times))
    batched_median = float(np.median(batched_times))
    return {
        "status": "PASS_NUMERICAL_EQUIVALENCE",
        "dimensions": {
            "h0_points": h0_points, "redshift_points": z_points,
            "events": n_events, "density_nodes": n_densities,
            "repeats": repeats,
        },
        "scalar_seconds": scalar_times,
        "batched_seconds": batched_times,
        "scalar_median_seconds": scalar_median,
        "batched_median_seconds": batched_median,
        "observed_speedup": scalar_median / batched_median,
        "max_log_likelihood_absolute_difference": max_ll_difference,
        "max_log_efficiency_absolute_difference": max_eff_difference,
        "limitations": [
            "Timing is hardware, BLAS, and load dependent; no universal speed claim.",
            "This isolates likelihood kernels on precomputed densities and excludes KDE construction.",
            "Synthetic benchmark only; it does not establish scientific calibration or estimator performance.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h0-points", type=int, default=321)
    parser.add_argument("--z-points", type=int, default=801)
    parser.add_argument("--events", type=int, default=30)
    parser.add_argument("--density-nodes", type=int, default=23)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20261114)
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()
    report = benchmark(
        h0_points=args.h0_points, z_points=args.z_points,
        n_events=args.events, n_densities=args.density_nodes,
        repeats=args.repeats, seed=args.seed,
    )
    with open(args.output, "w", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
