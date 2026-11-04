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

from dark_siren_robustness.independent_calibration_audit import (
    reconstruct_calibration_summaries,
)

METHODS = {"paired_bootstrap", "paired_bootstrap_fpc", "design_aware_srswor"}
KEY = ["scenario", "validation_size", "outer_draw", "method"]


def _require(condition: bool, message: str) -> None:
    if not bool(condition):
        raise ValueError(message)


def _same_number(a: float, b: float, *, atol: float = 1e-10) -> bool:
    if np.isnan(a) and np.isnan(b):
        return True
    return bool(np.isclose(a, b, rtol=1e-9, atol=atol, equal_nan=False))


def audit(directory: Path) -> dict:
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    _require(
        manifest.get("schema") in {
            "nested_validation_bootstrap_comparison_smoke_v1",
            "nested_validation_bootstrap_comparison_smoke_v2",
        },
        "unsupported artifact schema",
    )
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
        if row.bootstrap_status in {"success", "non_evaluable"} and row.method != "paired_bootstrap_fpc":
            mask = np.asarray(row.bootstrap_success_mask, dtype=bool)
            means = np.asarray(row.bootstrap_attempt_posterior_mean, dtype=float)
            refs = np.asarray(row.bootstrap_reference_posterior_mean_by_attempt, dtype=float)
            _require(mask.shape == (b,) and means.shape == (b,) and refs.shape == (b,),
                     "attempt-aligned bootstrap arrays have wrong shape")
            _require(int(mask.sum()) == int(row.n_bootstrap_successful),
                     "success mask disagrees with successful count")
            _require(int(row.n_bootstrap_attempted) == b
                     and int(row.n_bootstrap_failed) == b - int(mask.sum()),
                     "inner attempt counts do not reconcile with success mask")
            _require((row.bootstrap_status == "success") == (int(mask.sum()) >= 2),
                     "success/non-evaluable status disagrees with successful count")
            _require(np.isfinite(means[mask]).all() and np.isfinite(refs[mask]).all(),
                     "successful attempts contain nonfinite posterior means")
            _require(np.isnan(means[~mask]).all() and np.isnan(refs[~mask]).all(),
                     "failed attempts must have NaN posterior means")
            _require(len(row.bootstrap_failure_records) == int((~mask).sum()),
                     "failure-record count disagrees with failure mask")
            _require(int(row.n_bootstrap_attempted) == b,
                     "bootstrap attempt count differs from manifest")
            if row.method == "paired_bootstrap":
                samples = row.bootstrap_sample_indices
                _require(isinstance(samples, list) and len(samples) == b,
                         "paired-bootstrap index trace has wrong attempt count")
                pseudo_trace = row.bootstrap_pseudo_population_indices
                _require(pseudo_trace is None
                         or (isinstance(pseudo_trace, float) and np.isnan(pseudo_trace))
                         or (isinstance(pseudo_trace, list) and len(pseudo_trace) == 0),
                         "paired bootstrap must not contain pseudo-population indices")
                for sample in samples:
                    ix = np.asarray(sample)
                    _require(ix.shape == (int(row.validation_size),)
                             and np.issubdtype(ix.dtype, np.integer)
                             and np.all((ix >= 0) & (ix < int(row.validation_size))),
                             "paired-bootstrap sample indices are malformed")
            elif row.method == "design_aware_srswor":
                pseudo = row.bootstrap_pseudo_population_indices
                samples = row.bootstrap_sample_indices
                positions = row.bootstrap_selected_pseudo_population_positions
                _require(all(isinstance(v, list) for v in (pseudo, samples, positions))
                         and len(pseudo) == len(samples) == len(positions) == b,
                         "pseudo-population index traces have wrong attempt count")
                for pidx, sidx, pos in zip(pseudo, samples, positions):
                    p = np.asarray(pidx)
                    s = np.asarray(sidx)
                    q = np.asarray(pos)
                    _require(p.shape == (parent_n,) and s.shape == (int(row.validation_size),)
                             and q.shape == (int(row.validation_size),),
                             "pseudo-population resampling index dimensions are malformed")
                    _require(np.issubdtype(p.dtype, np.integer)
                             and np.issubdtype(s.dtype, np.integer)
                             and np.issubdtype(q.dtype, np.integer),
                             "pseudo-population resampling indices must be integers")
                    _require(np.all((p >= 0) & (p < int(row.validation_size)))
                             and np.all((s >= 0) & (s < int(row.validation_size)))
                             and np.all((q >= 0) & (q < parent_n))
                             and len(np.unique(q)) == len(q),
                             "pseudo-population resampling indices violate their bounds/design")
                    _require(np.array_equal(s, p[q]),
                             "selected pseudo-population positions do not map to saved sample indices")
        elif row.bootstrap_status in {"success", "non_evaluable"}:
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

    rebuilt_scenario, rebuilt_aggregate = reconstruct_calibration_summaries(records)
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
