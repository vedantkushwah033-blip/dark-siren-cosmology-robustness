"""Audit Bayesian fixed-truth coverage against an inverted Gaussian distance interval.

This is an independent artifact-level diagnostic for Arm A. It does not import
the simulation runner or production posterior code.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import NormalDist

C_KM_S = 299_792.458
REQUIRED = {
    "true_h0", "fractional_error", "replicate", "arm", "host_redshift",
    "observed_distance_mpc", "generated_sigma_mpc",
    "lower_68", "upper_68", "lower_95", "upper_95",
}
ARM_A = "A_fixed_realized_sigma"


def wilson(successes: int, n: int, z: float = 1.959963984540054) -> list[float]:
    if n <= 0 or not 0 <= successes <= n:
        raise ValueError("Wilson interval requires n > 0 and 0 <= successes <= n")
    p = successes / n
    den = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / den
    half = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / den
    return [max(0.0, center - half), min(1.0, center + half)]


def _summarize(rows: list[dict], level: int) -> dict:
    bayes_key = f"bayes_covered_{level}"
    reference_key = f"reference_covered_{level}"
    n = len(rows)
    out = {"n": n}
    for label, key in (("bayesian_equal_tail", bayes_key), ("inverted_gaussian_distance", reference_key)):
        successes = sum(bool(row[key]) for row in rows)
        out[label] = {
            "successes": successes,
            "n": n,
            "coverage": successes / n,
            "wilson_95": wilson(successes, n),
        }
    return out


def audit(path: Path) -> dict:
    raw_bytes = path.read_bytes()
    sha256 = hashlib.sha256(raw_bytes).hexdigest()
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        missing = REQUIRED - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing required columns: {sorted(missing)}")
        parsed = []
        seen = set()
        for line, row in enumerate(reader, start=2):
            try:
                truth = float(row["true_h0"])
                frac = float(row["fractional_error"])
                replicate = int(row["replicate"])
                arm = row["arm"]
                z_host = float(row["host_redshift"])
                observed = float(row["observed_distance_mpc"])
                sigma = float(row["generated_sigma_mpc"])
                bounds = {
                    68: (float(row["lower_68"]), float(row["upper_68"])),
                    95: (float(row["lower_95"]), float(row["upper_95"])),
                }
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid numeric field at CSV line {line}") from exc
            numeric = [truth, frac, z_host, observed, sigma, *(x for pair in bounds.values() for x in pair)]
            if not all(math.isfinite(x) for x in numeric):
                raise ValueError(f"non-finite value at CSV line {line}")
            if truth <= 0 or frac <= 0 or z_host <= 0 or observed <= 0 or sigma <= 0:
                raise ValueError(f"non-positive physical/numerical input at CSV line {line}")
            if replicate < 0:
                raise ValueError(f"negative replicate index at CSV line {line}")
            for level, (lo, hi) in bounds.items():
                if lo > hi:
                    raise ValueError(f"reversed {level}% posterior interval at CSV line {line}")
                if level == 95 and (lo > bounds[68][0] or hi < bounds[68][1]):
                    raise ValueError(f"95% interval does not contain 68% interval at CSV line {line}")
            key = (truth, frac, replicate, arm)
            if key in seen:
                raise ValueError(f"duplicate record key at CSV line {line}: {key}")
            seen.add(key)
            if arm != ARM_A:
                continue

            # Invert a central Gaussian interval for d ~ Normal(c z/H0, sigma).
            coverage = {}
            for level, tail_probability in ((68, 0.16), (95, 0.025)):
                q = NormalDist().inv_cdf(1.0 - tail_probability)
                d_low = observed - q * sigma
                d_high = observed + q * sigma
                if d_high <= 0:
                    raise ValueError(f"non-positive upper distance endpoint at CSV line {line}")
                h_low = C_KM_S * z_host / d_high
                h_high = math.inf if d_low <= 0 else C_KM_S * z_host / d_low
                lo, hi = bounds[level]
                coverage[f"bayes_covered_{level}"] = lo <= truth <= hi
                coverage[f"reference_covered_{level}"] = h_low <= truth <= h_high
            parsed.append({
                "true_h0": truth,
                "fractional_error": frac,
                "replicate": replicate,
                **coverage,
            })

    if not parsed:
        raise ValueError(f"no {ARM_A} records found")
    groups = {
        "all_arm_a_records": parsed,
        "interior_truths_55_to_85": [r for r in parsed if 55 <= r["true_h0"] <= 85],
        "central_truth_h0_70": [r for r in parsed if r["true_h0"] == 70],
    }
    by_truth = {
        str(int(truth)): _summarize([r for r in parsed if r["true_h0"] == truth], 68)
        for truth in sorted({r["true_h0"] for r in parsed})
    }
    by_fraction_interior = {
        str(frac): _summarize(
            [r for r in parsed if r["fractional_error"] == frac and 55 <= r["true_h0"] <= 85], 68
        )
        for frac in sorted({r["fractional_error"] for r in parsed})
    }
    return {
        "status": "PASS",
        "scope": "independent fixed-truth interval-coverage comparison for synthetic Arm A",
        "not_real_gw_data": True,
        "source_file": path.name,
        "source_sha256": sha256,
        "arm": ARM_A,
        "arm_a_record_count": len(parsed),
        "reference_method": (
            "Invert central Gaussian distance interval using sigma recorded for each synthetic replicate; "
            "if d_obs - q*sigma <= 0, set the upper H0 endpoint to infinity."
        ),
        "reference_is_oracle_conditional": True,
        "nominal_levels": {"68": {"q": NormalDist().inv_cdf(0.84)}, "95": {"q": NormalDist().inv_cdf(0.975)}},
        "groups": {
            name: {str(level): _summarize(rows, level) for level in (68, 95)}
            for name, rows in groups.items()
        },
        "by_truth_68": by_truth,
        "by_fraction_interior_68": by_fraction_interior,
        "interpretation_warning": (
            "Bayesian credible coverage at fixed truths is not generally required to equal nominal frequentist "
            "coverage. Arm A's realized sigma depends on the simulation truth and is an oracle-style benchmark. "
            "This diagnostic does not establish astrophysical validity or universal posterior calibration."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.input_csv)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "arm_a_record_count": report["arm_a_record_count"],
        "source_sha256": report["source_sha256"],
        "interior_68": report["groups"]["interior_truths_55_to_85"]["68"],
        "interior_95": report["groups"]["interior_truths_55_to_85"]["95"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
