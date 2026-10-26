"""Independently audit selected-likelihood PIT raw records and summary.

This audit reconstructs the Gaussian host mixture and selection denominator
without importing the production posterior helper. It is a synthetic-model
consistency check, not evidence of astrophysical validity.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.special import ndtr
from scipy.stats import kstest, norm


def wilson_interval(successes: int, trials: int, z: float = 1.959963984540054):
    if trials <= 0 or successes < 0 or successes > trials:
        raise ValueError("invalid binomial counts")
    p = successes / trials
    denom = 1.0 + z * z / trials
    center = (p + z * z / (2.0 * trials)) / denom
    half = z * np.sqrt(p * (1.0 - p) / trials + z * z / (4.0 * trials * trials)) / denom
    return float(center - half), float(center + half)


def audit(raw_path: str | Path, summary_path: str | Path) -> dict:
    raw_path, summary_path = Path(raw_path), Path(summary_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    with raw_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    n = len(rows)
    declared_n = summary.get("n_replicates")
    if isinstance(declared_n, bool) or not isinstance(declared_n, int) or declared_n < 1:
        raise ValueError("summary replicate count must be a positive integer")
    if n != declared_n:
        raise ValueError("raw row count does not match summary replicate count")

    seed_value = summary.get("seed")
    if isinstance(seed_value, bool) or not isinstance(seed_value, int) or seed_value < 0:
        raise ValueError("summary seed must be a non-negative integer")

    prior = summary["h0_prior"]
    grid_meta = summary["posterior_grid"]
    grid_size_value = grid_meta.get("points")
    if isinstance(grid_size_value, bool) or not isinstance(grid_size_value, int) or grid_size_value < 2:
        raise ValueError("posterior grid point count must be an integer of at least two")
    lo, hi = float(prior["min"]), float(prior["max"])
    if (
        prior["distribution"] != "uniform"
        or not np.isfinite([lo, hi]).all()
        or not lo < hi
        or lo != float(grid_meta["min"])
        or hi != float(grid_meta["max"])
    ):
        raise ValueError("prior and posterior-grid support metadata disagree")
    if grid_meta.get("integration") != "trapezoidal":
        raise ValueError("unsupported or missing posterior-grid integration metadata")
    protocol = summary.get("generation_protocol", {})
    expected_protocol = {
        "truth_draw": "once_per_replicate_from_declared_prior",
        "selection_conditioning": "rejection_sample_host_and_distance_conditional_on_fixed_truth",
        "acceptance_rule": "detection_distance_min < observed_distance < detection_distance_max",
        "rng": "numpy.random.default_rng (PCG64 default bit generator)",
    }
    if protocol != expected_protocol:
        raise ValueError("generation protocol metadata missing or inconsistent")
    grid_size = grid_size_value
    grid = np.linspace(lo, hi, grid_size)
    redshifts = np.asarray(summary["host_redshifts"], dtype=float)
    weights = np.asarray(summary["host_weights"], dtype=float)
    if redshifts.ndim != 1 or weights.shape != redshifts.shape or len(redshifts) == 0:
        raise ValueError("invalid host-mixture metadata")
    if not np.isfinite(redshifts).all() or np.any(redshifts <= 0):
        raise ValueError("host redshifts must be finite and positive")
    if (
        np.any(weights < 0)
        or not np.isfinite(weights).all()
        or not np.isclose(weights.sum(), 1.0)
    ):
        raise ValueError("host weights must be finite, non-negative, and sum to one")
    sigma = float(summary["distance_sigma_mpc"])
    dmin, dmax = map(float, summary["detection_interval_mpc_open"])
    if not (
        np.isfinite(sigma)
        and sigma > 0
        and np.isfinite([dmin, dmax]).all()
        and 0 <= dmin < dmax
    ):
        raise ValueError("invalid noise or detection-interval metadata")

    # Replay the declared RNG protocol from the recorded seed. This verifies
    # the actual simulated truth, host, observed distance, and rejection count,
    # not merely that those values look plausible.
    rng = np.random.default_rng(seed_value)
    seen = set()
    pits = []
    truths = []
    for row in rows:
        replicate = int(row["replicate"])
        seed = int(row["seed"])
        truth = float(row["truth_h0"])
        host = int(row["host_index"])
        host_z = float(row["host_redshift"])
        distance = float(row["observed_distance_mpc"])
        attempts = int(row["rejection_attempts"])
        pit_recorded = float(row["pit"])
        if replicate in seen:
            raise ValueError(f"duplicate replicate id {replicate}")
        seen.add(replicate)
        if seed != seed_value:
            raise ValueError("raw seed does not match summary")
        if not lo <= truth <= hi or not 0 <= host < len(redshifts):
            raise ValueError("truth or host index outside declared support")
        if not np.isclose(host_z, redshifts[host], rtol=0.0, atol=1e-14):
            raise ValueError("host redshift does not match host index")
        if not dmin < distance < dmax or attempts < 1:
            raise ValueError("raw observation violates detection or attempt constraints")

        # Exact deterministic replay catches altered raw data and confirms the
        # fixed-truth rejection-sampling protocol used by the generator.
        expected_truth = float(rng.uniform(lo, hi))
        expected_attempts = 0
        while True:
            expected_attempts += 1
            expected_host = int(rng.choice(len(redshifts), p=weights))
            expected_mean = 299792.458 * redshifts[expected_host] / expected_truth
            expected_distance = float(rng.normal(expected_mean, sigma))
            if dmin < expected_distance < dmax:
                break
            if expected_attempts > 1_000_000:
                raise ValueError("deterministic replay exceeded rejection-sampling safety limit")
        if (
            truth != expected_truth
            or host != expected_host
            or distance != expected_distance
            or attempts != expected_attempts
        ):
            raise ValueError(f"generation replay mismatch at replicate {replicate}")
        truths.append(truth)

        if not 0.0 <= pit_recorded <= 1.0:
            raise ValueError("recorded PIT outside [0, 1]")

        # Independent direct Gaussian-mixture reconstruction of p(d | H0)
        # and beta(H0) = P(detection | H0).
        means = 299792.458 * redshifts[None, :] / grid[:, None]
        likelihood = np.sum(
            weights[None, :] * norm.pdf(distance, loc=means, scale=sigma),
            axis=1,
        )
        host_efficiency = ndtr((dmax - means) / sigma) - ndtr((dmin - means) / sigma)
        beta = host_efficiency @ weights
        if not np.all(np.isfinite(beta)) or np.any(beta <= 0):
            raise ValueError("invalid selection efficiency during independent reconstruction")
        posterior_unnormalized = likelihood / beta
        normalizer = trapezoid(posterior_unnormalized, grid)
        if not np.isfinite(normalizer) or normalizer <= 0:
            raise ValueError("invalid posterior normalization during independent reconstruction")
        posterior = posterior_unnormalized / normalizer
        increments = 0.5 * (posterior[1:] + posterior[:-1]) * np.diff(grid)
        cdf = np.concatenate(([0.0], np.cumsum(increments)))
        cdf /= cdf[-1]
        pit_reconstructed = float(np.interp(truth, grid, cdf))
        if not np.isclose(pit_recorded, pit_reconstructed, rtol=2e-8, atol=2e-10):
            raise ValueError(
                f"PIT mismatch at replicate {replicate}: recorded={pit_recorded:.12g}, "
                f"independent={pit_reconstructed:.12g}"
            )
        # Reconstructing the PIT above verifies each recorded value; use the
        # recorded values for summary reconstruction to avoid bin-edge changes
        # caused solely by harmless floating-point differences across libraries.
        pits.append(pit_recorded)

    if seen != set(range(n)):
        raise ValueError("replicate ids must be exactly 0..n-1")
    pits = np.asarray(pits, dtype=float)
    counts, edges = np.histogram(pits, bins=np.linspace(0.0, 1.0, 11))
    expected_mean = float(np.mean(pits))
    expected_se = float(np.std(pits, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
    expected_variance = float(np.var(pits, ddof=1)) if n > 1 else 0.0
    ks = kstest(pits, "uniform")
    truth_prior_ks = kstest((np.asarray(truths) - lo) / (hi - lo), "uniform")
    expected_truth_mean = float(np.mean(truths))
    checks = [
        ("truth_prior_ks_statistic", float(truth_prior_ks.statistic)),
        ("truth_prior_ks_p_value_descriptive_only", float(truth_prior_ks.pvalue)),
        ("truth_prior_mean", expected_truth_mean),
        ("truth_prior_mean_expected", float((lo + hi) / 2.0)),
        ("pit_mean", expected_mean),
        ("pit_mean_monte_carlo_se", expected_se),
        ("pit_variance", expected_variance),
        ("ks_statistic", float(ks.statistic)),
        ("ks_p_value_descriptive_only", float(ks.pvalue)),
        ("max_abs_decile_frequency_deviation", float(np.max(np.abs(counts / n - 0.1)))),
    ]
    for key, expected in checks:
        if not np.isclose(float(summary[key]), expected, rtol=2e-8, atol=2e-10):
            raise ValueError(f"summary field {key} disagrees with independently reconstructed raw PITs")

    reported_bins = summary["decile_bins_with_wilson_95_intervals"]
    if len(reported_bins) != 10:
        raise ValueError("summary must contain ten PIT decile bins")
    for i, item in enumerate(reported_bins):
        low, high = wilson_interval(int(counts[i]), n)
        expected = {
            "left": float(edges[i]),
            "right": float(edges[i + 1]),
            "count": int(counts[i]),
            "frequency": float(counts[i] / n),
            "wilson_95_low": low,
            "wilson_95_high": high,
        }
        for key, value in expected.items():
            if not np.isclose(float(item[key]), float(value), rtol=2e-8, atol=2e-10):
                raise ValueError(f"decile {i} field {key} disagrees with raw PITs")

    return {
        "status": "PASS",
        "n_replicates_audited": n,
        "independent_pit_reconstruction": "passed",
        "summary_reconstruction": "passed",
        "interpretation": "Synthetic consistency only; not astrophysical validation.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", default="artifacts/selected_likelihood_sbc/pit_replicates.csv")
    parser.add_argument("--summary", default="artifacts/selected_likelihood_sbc/pit_summary.json")
    parser.add_argument("--output", default="artifacts/selected_likelihood_sbc/audit_report.json")
    args = parser.parse_args()
    result = audit(args.raw, args.summary)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
