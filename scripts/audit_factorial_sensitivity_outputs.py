"""Independent structural audit for factorial sensitivity artifacts.

This validator checks experiment design, row counts, finite outputs, exact
catalogue masks, posterior interval ordering, and scenario/manifest agreement.
It does not independently validate the simulator's scientific assumptions.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

EXPECTED = {
    "complete_uniform_matched": (1.0, "complete", "uniform", "uniform"),
    "complete_luminosity_matched": (1.0, "complete", "luminosity", "luminosity"),
    "complete_luminosity_host_uniform_inference": (1.0, "complete", "luminosity", "uniform"),
    "random70_uniform_matched": (0.7, "random", "uniform", "uniform"),
    "random70_luminosity_matched": (0.7, "random", "luminosity", "luminosity"),
    "random70_luminosity_host_uniform_inference": (0.7, "random", "luminosity", "uniform"),
    "redshift70_uniform_matched": (0.7, "redshift", "uniform", "uniform"),
    "redshift70_luminosity_matched": (0.7, "redshift", "luminosity", "luminosity"),
    "redshift70_luminosity_host_uniform_inference": (0.7, "redshift", "luminosity", "uniform"),
}
POSTERIOR_FIELDS = ("median", "lower_68", "upper_68", "lower_95", "upper_95", "width_68", "width_95")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def finite(row: dict[str, str], key: str, label: str) -> float:
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"{label}: missing/non-numeric {key}") from exc
    if not math.isfinite(value):
        raise ValueError(f"{label}: non-finite {key}={value}")
    return value


def audit(root: Path) -> dict:
    summary_json = json.loads((root / "scenario_summary.json").read_text(encoding="utf-8"))
    summary = read_csv(root / "scenario_summary.csv")
    runs = read_csv(root / "per_run_results.csv")
    n = int(summary_json["requested_replicates_per_scenario"])
    if n < 2:
        raise ValueError("replicate count must be at least 2")
    try:
        true_h0 = float(summary_json["true_h0"])
        n_galaxies = int(summary_json["n_galaxies"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("JSON manifest must include true_h0 and n_galaxies") from exc
    if not math.isfinite(true_h0) or not math.isclose(true_h0, 70.0, abs_tol=1e-12):
        raise ValueError(f"unexpected true_h0; expected 70.0, got {true_h0}")
    if n_galaxies != 500:
        raise ValueError(f"unexpected n_galaxies; expected 500, got {n_galaxies}")
    expected_source_commit = "b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f"
    if summary_json.get("original_commit") != expected_source_commit:
        raise ValueError(
            "original simulator commit does not match the audited pinned source "
            f"({expected_source_commit})"
        )
    if int(summary_json["scenario_count"]) != len(EXPECTED) or len(summary) != len(EXPECTED):
        raise ValueError("expected exactly nine factorial scenarios")
    if len(runs) != len(EXPECTED) * n:
        raise ValueError(f"expected {len(EXPECTED)*n} per-run rows, found {len(runs)}")

    summary_by_name = {row["scenario"]: row for row in summary}
    if len(summary_by_name) != len(summary) or set(summary_by_name) != set(EXPECTED):
        raise ValueError("summary scenario names do not match prespecified design")

    # Check the JSON manifest independently against both the fixed design and
    # the CSV summary. A valid raw CSV must not conceal stale or edited metadata.
    manifest_scenarios = summary_json.get("scenarios")
    if not isinstance(manifest_scenarios, list) or len(manifest_scenarios) != len(EXPECTED):
        raise ValueError("JSON manifest must contain exactly nine scenario records")
    manifest_by_name = {item.get("scenario"): item for item in manifest_scenarios}
    if len(manifest_by_name) != len(manifest_scenarios) or set(manifest_by_name) != set(EXPECTED):
        raise ValueError("JSON manifest scenario names do not match prespecified design")
    try:
        base_seed = int(summary_json["base_seed"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("JSON manifest must contain an integer base_seed") from exc
    for scenario_index, (name, (target, missingness, host_generation, inference_weighting)) in enumerate(EXPECTED.items()):
        item = manifest_by_name[name]
        expected_seed = base_seed + scenario_index * 1_000_000
        try:
            manifest_seed = int(item["seed"])
            csv_seed = int(summary_by_name[name]["seed"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{name}: missing/non-integer scenario seed") from exc
        if manifest_seed != expected_seed or csv_seed != expected_seed:
            raise ValueError(
                f"{name}: seed mismatch; expected {expected_seed}, "
                f"JSON={manifest_seed}, CSV={csv_seed}"
            )
        for field, expected in (("target_completeness", target),
                                ("missingness", missingness),
                                ("host_generation", host_generation),
                                ("inference_weighting", inference_weighting)):
            if field not in item:
                raise ValueError(f"{name}: JSON manifest missing {field}")
            actual = item[field]
            if field == "target_completeness":
                if not math.isclose(float(actual), expected, abs_tol=1e-12):
                    raise ValueError(f"{name}: JSON manifest design mismatch for {field}")
            elif actual != expected:
                raise ValueError(f"{name}: JSON manifest design mismatch for {field}")
        if int(item.get("requested_replicates", -1)) != n:
            raise ValueError(f"{name}: JSON manifest replicate count mismatch")
        csv_item = summary_by_name[name]
        try:
            csv_requested_replicates = int(csv_item["requested_replicates"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{name}: CSV summary missing requested_replicates") from exc
        if csv_requested_replicates != n:
            raise ValueError(f"{name}: CSV summary replicate count mismatch")
        for field in ("target_completeness", "missingness", "host_generation", "inference_weighting"):
            if field in ("missingness", "host_generation", "inference_weighting"):
                if str(item[field]) != csv_item[field]:
                    raise ValueError(f"{name}: JSON/CSV mismatch for {field}")
            elif not math.isclose(float(item[field]), float(csv_item[field]), abs_tol=1e-12):
                raise ValueError(f"{name}: JSON/CSV mismatch for {field}")
        # The JSON embeds the same summary objects as the CSV export. Check
        # counts and all headline metrics so the two deliverables cannot drift.
        for field in ("n", "median_bias", "mean_bias", "rmse_median", "rmse_mean",
                      "mean_width_68", "mean_width_95", "coverage_68", "coverage_95",
                      "host_observed_fraction"):
            if field not in item or field not in csv_item:
                raise ValueError(f"{name}: JSON/CSV summary missing {field}")
            if not math.isclose(float(item[field]), float(csv_item[field]),
                                rel_tol=1e-9, abs_tol=1e-9):
                raise ValueError(f"{name}: JSON/CSV mismatch for summary metric {field}")

    grouped: dict[str, list[dict[str, str]]] = {name: [] for name in EXPECTED}
    for row in runs:
        name = row.get("scenario", "")
        if name not in EXPECTED:
            raise ValueError(f"unexpected scenario in raw output: {name!r}")
        grouped[name].append(row)

    checked = 0
    for name, design in EXPECTED.items():
        target, missingness, host_generation, inference_weighting = design
        srow = summary_by_name[name]
        for field, expected in (
            ("target_completeness", target),
            ("missingness", missingness),
            ("host_generation", host_generation),
            ("inference_weighting", inference_weighting),
        ):
            if field not in srow or (field in ("missingness", "host_generation", "inference_weighting")
                                     and srow[field] != expected):
                raise ValueError(f"{name}: summary design mismatch for {field}")
            if field == "target_completeness" and not math.isclose(float(srow[field]), expected, abs_tol=1e-12):
                raise ValueError(f"{name}: target completeness mismatch")
        if int(srow["n"]) != n or len(grouped[name]) != n:
            raise ValueError(f"{name}: replicate count mismatch")
        indices = set()
        for row in grouped[name]:
            label = f"{name} run={row.get('run', '?')}"
            indices.add(row.get("run"))
            for field, expected in (("missingness", missingness), ("host_generation", host_generation), ("inference_weighting", inference_weighting)):
                if row.get(field) != expected:
                    raise ValueError(f"{label}: raw scenario design mismatch for {field}")
            try:
                raw_scenario_seed = int(row["scenario_seed"])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(f"{label}: missing/non-integer scenario_seed") from exc
            expected_scenario_seed = int(srow["seed"])
            if raw_scenario_seed != expected_scenario_seed:
                raise ValueError(
                    f"{label}: raw scenario_seed mismatch; expected "
                    f"{expected_scenario_seed}, got {raw_scenario_seed}"
                )
            target_value = finite(row, "target_completeness", label)
            realized = finite(row, "realized_completeness", label)
            if not math.isclose(target_value, target, abs_tol=1e-12):
                raise ValueError(f"{label}: per-run target mismatch")
            if target == 0.7 and not math.isclose(realized, 0.7, abs_tol=1e-12):
                raise ValueError(f"{label}: expected exact 70% completeness, got {realized}")
            if target == 1.0 and not math.isclose(realized, 1.0, abs_tol=1e-12):
                raise ValueError(f"{label}: expected complete catalogue (100%), got {realized}")
            vals = {key: finite(row, key, label) for key in POSTERIOR_FIELDS}
            if not (vals["lower_95"] <= vals["lower_68"] <= vals["median"]
                    <= vals["upper_68"] <= vals["upper_95"]):
                raise ValueError(f"{label}: posterior interval ordering invalid")
            if not math.isclose(vals["width_68"], vals["upper_68"] - vals["lower_68"], rel_tol=1e-9, abs_tol=1e-9):
                raise ValueError(f"{label}: width_68 inconsistent with interval endpoints")
            if not math.isclose(vals["width_95"], vals["upper_95"] - vals["lower_95"], rel_tol=1e-9, abs_tol=1e-9):
                raise ValueError(f"{label}: width_95 inconsistent with interval endpoints")
            checked += 1
        expected_indices = {str(index) for index in range(n)}
        if len(indices) != n or indices != expected_indices:
            raise ValueError(f"{name}: run identifiers must be unique integers 0..{n-1}")
        # Recompute core scientific summaries from raw replicate rows rather
        # than merely checking their shape. median_bias is intentionally the
        # median of replicate-level median errors in this experiment runner.
        raw = grouped[name]
        medians = [finite(row, "median", name) for row in raw]
        means = [finite(row, "mean", name) for row in raw]
        lo68 = [finite(row, "lower_68", name) for row in raw]
        hi68 = [finite(row, "upper_68", name) for row in raw]
        lo95 = [finite(row, "lower_95", name) for row in raw]
        hi95 = [finite(row, "upper_95", name) for row in raw]
        host_observed = []
        for row in raw:
            value = row.get("host_observed", "").strip().lower()
            if value not in {"true", "false"}:
                raise ValueError(f"{name}: host_observed must be a boolean, got {value!r}")
            host_observed.append(value == "true")
        if target == 1.0 and not all(host_observed):
            raise ValueError(f"{name}: complete catalogue must have every simulated host observed")
        errors = [x - 70.0 for x in medians]
        mean_errors = [x - 70.0 for x in means]
        sorted_errors = sorted(errors)
        mid = len(sorted_errors) // 2
        median_bias = (sorted_errors[mid] if len(sorted_errors) % 2
                       else (sorted_errors[mid - 1] + sorted_errors[mid]) / 2.0)
        expected_metrics = {
            "median_bias": median_bias,
            "mean_bias": sum(mean_errors) / n,
            "rmse_median": math.sqrt(sum(x*x for x in errors) / n),
            "rmse_mean": math.sqrt(sum(x*x for x in mean_errors) / n),
            "mean_width_68": sum(b-a for a, b in zip(lo68, hi68)) / n,
            "mean_width_95": sum(b-a for a, b in zip(lo95, hi95)) / n,
            "coverage_68": sum(a <= 70.0 <= b for a, b in zip(lo68, hi68)) / n,
            "coverage_95": sum(a <= 70.0 <= b for a, b in zip(lo95, hi95)) / n,
            "host_observed_fraction": sum(host_observed) / n,
        }
        for metric, expected_value in expected_metrics.items():
            stored = finite(srow, metric, f"{name} summary")
            if not math.isclose(stored, expected_value, rel_tol=1e-9, abs_tol=1e-9):
                raise ValueError(
                    f"{name}: summary metric {metric} disagrees with raw rows "
                    f"(stored={stored}, recomputed={expected_value})"
                )
    return {"status": "PASS", "scenarios": len(EXPECTED), "replicates_per_scenario": n,
            "raw_rows_checked": checked, "summary_rows_checked": len(summary),
            "scientific_scope": "synthetic structural audit only"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir", nargs="?", default="results/factorial-sensitivity")
    args = parser.parse_args()
    print(json.dumps(audit(Path(args.results_dir)), indent=2))


if __name__ == "__main__":
    main()
