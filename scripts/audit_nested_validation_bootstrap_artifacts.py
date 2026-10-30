"""Independent structural and summary audit for nested validation-bootstrap artifacts.

This first-stage auditor checks artifact contracts and independently rebuilds
summary statistics from the saved calibration records. It does NOT regenerate
the parent population or replay posterior fits/RNG streams; those remain a
separate, explicitly required second-stage audit.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

METHODS = {"paired_bootstrap", "paired_bootstrap_fpc", "design_aware_srswor"}
KEY = ["scenario", "validation_size", "outer_draw", "method"]


def _require(condition: bool, message: str) -> None:
    if not bool(condition):
        raise ValueError(message)


def _same_number(a: float, b: float, *, atol: float = 1e-10) -> bool:
    if np.isnan(a) and np.isnan(b):
        return True
    return bool(np.isclose(a, b, rtol=1e-9, atol=atol, equal_nan=False))


def _independent_summaries(records: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Recompute calibration tables without importing the producer summarizer."""
    outer = records.drop_duplicates(["scenario", "validation_size", "outer_draw"])
    outer_rows = []
    for (scenario, size), group in outer.groupby(["scenario", "validation_size"], sort=True):
        values = group.loc[group.outer_status.eq("success"), "outer_posterior_mean"].to_numpy(float)
        outer_rows.append({
            "scenario": scenario, "validation_size": int(size),
            "n_outer_attempted": int(len(group)), "n_outer_successful": int(values.size),
            "outer_failure_fraction": float(1 - values.size / len(group)),
            "empirical_srswor_sd": float(np.std(values, ddof=1)) if values.size >= 2 else np.nan,
        })
    outer_summary = pd.DataFrame(outer_rows)
    scenario_rows = []
    for (scenario, size, method), group in records.groupby(
        ["scenario", "validation_size", "method"], sort=True
    ):
        boot = group.loc[group.bootstrap_status.eq("success")]
        oi = outer_summary.loc[
            outer_summary.scenario.eq(scenario) & outer_summary.validation_size.eq(size)
        ].iloc[0]
        sd = float(oi.empirical_srswor_sd)
        mean_se = float(boot.bootstrap_se.mean()) if len(boot) else np.nan
        ratio = mean_se / sd if np.isfinite(sd) and sd > 0 else np.nan
        covered = (
            (boot.interval_low.to_numpy(float) <= boot.finite_parent_target.to_numpy(float))
            & (boot.finite_parent_target.to_numpy(float) <= boot.interval_high.to_numpy(float))
        )
        n_outer = int(oi.n_outer_successful)
        scenario_rows.append({
            "scenario": scenario, "validation_size": int(size), "method": str(method),
            "n_outer_attempted": int(oi.n_outer_attempted),
            "n_outer_successful": n_outer,
            "outer_failure_fraction": float(oi.outer_failure_fraction),
            "n_bootstrap_attempted": n_outer,
            "n_bootstrap_successful": int(len(boot)),
            "bootstrap_failure_fraction": float(1 - len(boot) / n_outer) if n_outer else np.nan,
            "empirical_srswor_sd": sd, "mean_inner_bootstrap_se": mean_se,
            "mean_se_to_empirical_sd_ratio": ratio,
            "finite_parent_target_interval_coverage": float(np.mean(covered)) if len(boot) else np.nan,
            "mean_interval_width": float((boot.interval_high - boot.interval_low).mean()) if len(boot) else np.nan,
            "n_intervals_evaluated": int(len(boot)),
        })
    scenario = pd.DataFrame(scenario_rows)
    aggregate_rows = []
    for (size, method), group in scenario.groupby(["validation_size", "method"], sort=True):
        ratio = group.mean_se_to_empirical_sd_ratio.dropna()
        coverage = group.finite_parent_target_interval_coverage.dropna()
        widths = group.mean_interval_width.dropna()
        aggregate_rows.append({
            "validation_size": int(size), "method": str(method),
            "n_scenarios": int(group.scenario.nunique()),
            "n_scenarios_with_estimable_sd_ratio": int(len(ratio)),
            "scenario_balanced_mean_se_ratio": float(ratio.mean()) if len(ratio) else np.nan,
            "scenario_balanced_sd_se_ratio": float(ratio.std(ddof=1)) if len(ratio) > 1 else np.nan,
            "scenario_balanced_mean_coverage": float(coverage.mean()) if len(coverage) else np.nan,
            "scenario_balanced_sd_coverage": float(coverage.std(ddof=1)) if len(coverage) > 1 else np.nan,
            "scenario_balanced_mean_interval_width": float(widths.mean()) if len(widths) else np.nan,
            "scenario_balanced_outer_failure_fraction": float(group.outer_failure_fraction.mean()),
            "scenario_balanced_bootstrap_failure_fraction": float(group.bootstrap_failure_fraction.mean()),
        })
    return scenario, pd.DataFrame(aggregate_rows)


def audit(directory: Path) -> dict:
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    _require(manifest.get("schema") == "nested_validation_bootstrap_comparison_smoke_v1",
             "unsupported artifact schema")
    _require(manifest.get("status") == "SMOKE_ONLY_NOT_CONFIRMATORY",
             "artifact must remain explicitly smoke-only")
    _require(set(manifest.get("methods", [])) == METHODS, "declared method set mismatch")
    scenarios, draws = int(manifest["scenarios"]), int(manifest["outer_draws"])
    parent_n = int(manifest["parent_size"])
    sizes = [int(x) for x in manifest["validation_sizes"]]
    b = int(manifest["bootstrap_replicates_per_method_per_outer_draw"])
    _require(scenarios > 0 and draws >= 2 and parent_n >= 100 and b >= 2,
             "invalid manifest design counts")
    _require(len(sizes) > 0 and len(set(sizes)) == len(sizes)
             and all(10 <= n <= parent_n for n in sizes), "invalid validation-size grid")

    records = pd.read_csv(directory / "calibration_records.csv")
    attempts = pd.read_json(directory / "attempt_audit.jsonl", lines=True)
    expected = scenarios * draws * len(sizes) * len(METHODS)
    _require(len(records) == expected, f"expected {expected} calibration rows")
    _require(len(attempts) == expected, f"expected {expected} attempt rows")
    _require(not records.duplicated(KEY).any(), "duplicate calibration key")
    _require(not attempts.duplicated(KEY).any(), "duplicate attempt key")
    _require(set(map(tuple, records[KEY].itertuples(index=False, name=None))) ==
             set(map(tuple, attempts[KEY].itertuples(index=False, name=None))),
             "calibration and attempt keys differ")
    _require(set(records.method.unique()) == METHODS, "observed method set mismatch")
    _require(set(records.scenario.unique()) == set(range(scenarios)), "scenario IDs mismatch")
    _require(set(records.outer_draw.unique()) == set(range(draws)), "outer draw IDs mismatch")
    _require(set(records.validation_size.unique()) == set(sizes), "validation sizes mismatch")

    # Shared outer estimate and full-parent target must not vary by method/draw/size.
    for key, cell in records.groupby(["scenario", "validation_size", "outer_draw"]):
        _require(cell.outer_status.nunique(dropna=False) == 1, f"outer status differs across methods: {key}")
        _require(cell.outer_posterior_mean.nunique(dropna=False) <= 1, f"outer estimate differs across methods: {key}")
        _require(cell.finite_parent_target.nunique(dropna=False) == 1, f"target differs across methods: {key}")
    for scenario, cell in records.groupby("scenario"):
        _require(cell.finite_parent_target.nunique(dropna=False) == 1,
                 f"finite-parent target changes within scenario {scenario}")

    # Scenario-level hashes and deterministic seed declarations must be invariant.
    for scenario, cell in attempts.groupby("scenario"):
        for col in ["parent_frame_sha256", "catalogue_sha256", "gw_distance_sha256"]:
            _require(cell[col].nunique(dropna=False) == 1, f"{col} changes within scenario {scenario}")
        ss = np.random.SeedSequence([int(manifest["seed"]), int(scenario)])
        outer_seed = int(ss.spawn(3)[2].generate_state(1, dtype=np.uint32)[0])
        _require((cell.outer_seed.astype(int) == outer_seed).all(),
                 f"outer seed does not match deterministic scenario seed for {scenario}")
    for row in attempts.itertuples(index=False):
        if row.bootstrap_seed is not None and not pd.isna(row.bootstrap_seed):
            method_key = 1 if row.method == "paired_bootstrap" else 2 if row.method == "design_aware_srswor" else None
            if method_key is not None:
                expected_seed = int(np.random.SeedSequence([
                    int(manifest["seed"]), int(row.scenario), int(row.outer_draw),
                    int(row.validation_size), method_key
                ]).generate_state(1, dtype=np.uint32)[0])
                _require(int(row.bootstrap_seed) == expected_seed,
                         f"bootstrap seed mismatch for {row.method} / {row.scenario}")

    # Prefix construction must make outer samples nested across validation sizes.
    for (scenario, draw), cell in attempts.groupby(["scenario", "outer_draw"]):
        by_size = sorted(
            ((int(r.validation_size), set(r.validation_indices))
             for r in cell.drop_duplicates("validation_size").itertuples(index=False)),
            key=lambda pair: pair[0],
        )
        for (n0, a), (n1, z) in zip(by_size, by_size[1:]):
            _require(n0 < n1 and a.issubset(z),
                     f"outer SRSWOR samples are not nested for scenario={scenario}, draw={draw}")

    merged = attempts.merge(records[KEY + ["outer_status", "bootstrap_status",
        "outer_posterior_mean", "bootstrap_se", "interval_low", "interval_high",
        "finite_parent_target"]], on=KEY, suffixes=("_audit", "_record"), validate="one_to_one")
    for row in merged.itertuples(index=False):
        idx = row.validation_indices
        _require(isinstance(idx, list) and len(idx) == int(row.validation_size),
                 "outer validation index length mismatch")
        ix = np.asarray(idx)
        _require(np.issubdtype(ix.dtype, np.integer), "outer validation indices must be integers")
        _require(np.all((ix >= 0) & (ix < parent_n)) and len(np.unique(ix)) == len(ix),
                 "outer validation indices out of bounds or duplicated")
        _require(bool(row.parent_frame_sha256) and bool(row.catalogue_sha256)
                 and bool(row.gw_distance_sha256), "missing provenance fingerprints")
        if row.bootstrap_status == "success" and row.method != "paired_bootstrap_fpc":
            mask = np.asarray(row.bootstrap_success_mask, dtype=bool)
            means = np.asarray(row.bootstrap_attempt_posterior_mean, dtype=float)
            refs = np.asarray(row.bootstrap_reference_posterior_mean_by_attempt, dtype=float)
            _require(mask.shape == (b,) and means.shape == (b,) and refs.shape == (b,),
                     "attempt-aligned bootstrap arrays have wrong shape")
            _require(int(mask.sum()) == int(row.n_bootstrap_successful),
                     "success mask disagrees with successful count")
            _require(np.isfinite(means[mask]).all() and np.isfinite(refs[mask]).all(),
                     "successful attempts contain nonfinite posterior means")
            _require(np.isnan(means[~mask]).all() and np.isnan(refs[~mask]).all(),
                     "failed attempts must have NaN posterior means")
            _require(len(row.bootstrap_failure_records) == int((~mask).sum()),
                     "failure-record count disagrees with failure mask")
        elif row.bootstrap_status == "success":
            _require(pd.isna(row.bootstrap_success_mask),
                     "FPC diagnostic should not duplicate bootstrap attempt arrays")
        else:
            _require(row.bootstrap_status in {"failed", "not_run"},
                     "invalid bootstrap status")
        if row.method == "paired_bootstrap_fpc" and row.bootstrap_status == "success":
            paired = attempts.loc[
                attempts.scenario.eq(row.scenario)
                & attempts.validation_size.eq(row.validation_size)
                & attempts.outer_draw.eq(row.outer_draw)
                & attempts.method.eq("paired_bootstrap")
            ].iloc[0]
            expected_fpc = float(np.sqrt(max(0.0, 1.0-row.validation_size/parent_n)))
            _require(_same_number(float(row.fpc_factor), expected_fpc), "FPC factor mismatch")
            _require(row.bootstrap_seed == paired.bootstrap_seed, "FPC comparator seed mismatch")
            paired_record = records.loc[
                records.scenario.eq(row.scenario)
                & records.validation_size.eq(row.validation_size)
                & records.outer_draw.eq(row.outer_draw)
                & records.method.eq("paired_bootstrap")
            ].iloc[0]
            fpc_se = float(paired_record.bootstrap_se) * expected_fpc
            _require(_same_number(float(row.bootstrap_se), fpc_se),
                     "FPC standard error does not match paired SE times FPC")
            center = float(row.outer_posterior_mean)
            _require(_same_number(float(row.interval_low), center-1.96*fpc_se),
                     "FPC lower endpoint mismatch")
            _require(_same_number(float(row.interval_high), center+1.96*fpc_se),
                     "FPC upper endpoint mismatch")

    rebuilt_scenario, rebuilt_aggregate = _independent_summaries(records)
    for filename, rebuilt in [
        ("scenario_summary.csv", rebuilt_scenario),
        ("scenario_balanced_summary.csv", rebuilt_aggregate),
    ]:
        saved = pd.read_csv(directory / filename)
        _require(list(saved.columns) == list(rebuilt.columns), f"{filename}: column schema mismatch")
        _require(len(saved) == len(rebuilt), f"{filename}: row count mismatch")
        joined = saved.merge(rebuilt, on=[c for c in ["scenario", "validation_size", "method"] if c in saved.columns],
                             suffixes=("_saved", "_rebuilt"), validate="one_to_one")
        for col in rebuilt.columns:
            if col in {"scenario", "validation_size", "method"}:
                continue
            a = pd.to_numeric(joined[f"{col}_saved"], errors="coerce").to_numpy(float)
            z = pd.to_numeric(joined[f"{col}_rebuilt"], errors="coerce").to_numpy(float)
            _require(np.allclose(a, z, rtol=1e-9, atol=1e-10, equal_nan=True),
                     f"{filename}: independently reconstructed {col} differs")
    return {"status": "STRUCTURAL_AND_SUMMARY_AUDIT_PASS",
            "scope": "artifact structure, row grain, shared targets, attempt arrays, FPC arithmetic, summaries",
            "limitations": "does not regenerate parent frames or replay RNG/posterior calculations",
            "calibration_rows": len(records), "attempt_rows": len(attempts)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.artifact_directory), indent=2))


if __name__ == "__main__":
    main()
