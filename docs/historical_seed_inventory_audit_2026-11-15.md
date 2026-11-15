# Historical seed inventory audit

Status: implementation added; **no historical inventory has yet been certified**. This is a fail-closed tool for auditing extracted artifacts, not an automatic claim that the supplied directories are complete.

## Purpose

Before the proposed density training and holdout experiments, compare every historical random seed actually recorded in retained development artifacts against both future seed blocks. Comparing only master-seed formulas is insufficient: old child configurations may have overridden seeds, and manifests may contain additional random streams.

## Input contract

The tool accepts one or more extracted artifact roots and a JSON schedule with non-empty `training` and `holdout` arrays. Each schedule row must explicitly provide integer `population_seed` and `validation_seed` values; any additional random stream must be represented by another integer field whose name contains `seed` (for example, `bootstrap_seed`). Seed values must be unique across all streams and replicates within a block, and the blocks must be mutually disjoint.

The scanner parses JSON, CSV, and NumPy NPZ files. It extracts scalar values from fields named exactly `seed`/`seeds` or ending in `_seed` (so `seed_stride` and aggregate fields such as `mean_pit_mean_across_seeds` are not treated as RNG seeds). NPZ loading uses `allow_pickle=False`; only arrays with seed-valued keys are inventoried. It hashes every supported source file with SHA-256 and records source-relative paths, byte counts, field names, and seed occurrences.

YAML files are deliberately unsupported and cause a non-evaluable result. Convert them to reviewed JSON/CSV without changing values, retain the originals in the audit package, and document the conversion. Malformed files and non-integer seed values fail closed rather than being skipped. NPZ object arrays are not loaded because pickle is disabled.

## Example schedule

```json
{
  "training": [
    {"population_seed": 2026111600, "validation_seed": 2026111601},
    {"population_seed": 2026111610, "validation_seed": 2026111611}
  ],
  "holdout": [
    {"population_seed": 2026112600, "validation_seed": 2026112601},
    {"population_seed": 2026112610, "validation_seed": 2026112611}
  ]
}
```

## Invocation

After extracting all relevant historical run artifacts into a reviewed directory:

```bash
python scripts/audit_historical_seed_inventory.py \
  --artifact-root artifacts/historical_density_runs \
  --artifact-root artifacts/historical_curvature_runs \
  --schedule configs/density_training_holdout_seed_schedule_candidate_v1.json \
  --output reports/historical_seed_inventory_audit.json
```

Exit code 0 means the supplied inventory passed the mechanical checks. Exit code 2 means it is not evaluable or has an overlap. The report contains the hashes and exact locations of every extracted seed occurrence.

## Critical completeness limitation

A clean scan proves only that no collision was found in the files supplied. It does **not** prove that every historical run or random stream was supplied. Before treating the report as a gate, maintain a separate source register covering all relevant workflows, run IDs, artifact IDs, extraction dates, expired/missing artifacts, and the reason each source is included or excluded. The audit is not complete while any relevant historical artifact is unavailable or unenumerated.

The repository currently does not contain the full set of extracted historical Actions artifacts needed to certify the proposed split. Therefore this implementation must not be run against an empty or partial directory and represented as a successful historical independence audit.

## What the audit does not establish

- Distinct integer seeds do not mathematically prove statistical independence.
- It does not validate the PRNG implementation or random-stream derivation.
- It does not infer seeds hidden in binary files or unrecognized field names.
- It does not download GitHub Actions artifacts or certify source-register completeness.
- It does not authorize the training run or open the holdout.

**Current decision: seed partition NOT YET VERIFIED; development pilot remains blocked.**


## Partial retained-artifact scan (2026-11-15; not a completeness certification)

A manual byte-level scan was performed on the following retained GitHub Actions artifacts. JSON/CSV fields whose names are exactly `seed`/`seeds` or end in `_seed`/`_seeds` were inventoried. ZIP SHA-256 hashes are included so the exact downloaded sources can be identified.

| Workflow run | Artifact ID | Artifact | ZIP SHA-256 |
|---|---:|---|---|
| 43617052462 | 14123634445 | expanded grid, floor 0.25 | `f7eb0576f0b9dc5badafb2a3b3b136bc125106856345df4749e173f94669c627` |
| 43617052462 | 14123867382 | expanded grid, floor 0.10 | `e83f944b22cc153c1b6354515f42ed25978db2cf7b94b2e66c13f16e04768d1a` |
| 43617052462 | 14123894369 | expanded grid, no floor | `26f9ccc6d711cff45ea7ac4a130eba48dcde7b243185d2f6c2fc2bb6c9a8af86` |
| 43671806461 | 14144815968 | seven-gamma engineering smoke | `290d77d9ff07bb949b10895d56fc510c84a8d6957d302fd3c18ab2f9f393a7e0` |
| 43674316806 | 14145680789 | compact paired curvature pilot | `c5972b12450fd70cc0e5041b934732d845b80eee97026c834a5035a66e09a77c` |
| 43674496952 | 14145857269 | numerical diagnostics, Python 3.10 | `b15717ef29cecf7e750181b6fde5ad1b5a196fae358c14cfa9a20bdff94a8e6b` |
| 43674496952 | 14146636130 | numerical diagnostics, Python 3.11 | `1f1ba6ba38c07eb67b2cea8e693e36c2da08080c1d3dccbcfa6e37cc9227375d` |
| 43674496952 | 14146830851 | numerical diagnostics, Python 3.12 | `a136dd553bee96fecd79de316e6386da9d3c9e4edeb9a4a2f81742e0b803b8ab` |
| 43674496952 | 14146747139 | five-arm KDE bandwidth pilot | `d1c2b65fb5f309a6d2fa43f3882315b2bf48f63a5fde22d74c9eb50b05c0f4d5` |

Across these nine retained ZIPs, the byte-level scan found **198 seed-bearing JSON/CSV files, 2,784,487 seed-field occurrences, and 260 distinct historical integer values**. No YAML files were present in these ZIPs. The three expanded-grid floor artifacts intentionally reuse the same replicate seed schedule; repeated occurrences across those three arms are not counted as additional distinct seeds.

For the proposed candidate schedule:
- training: 12 replicate pairs, 24 distinct stream seeds, range 2026111600–2026111711;
- holdout: 24 replicate pairs, 48 distinct stream seeds, range 2026112600–2026112831;
- training versus holdout overlap: **none**;
- training versus the historical values in these nine artifacts: **none**;
- holdout versus the historical values in these nine artifacts: **none**.

The seven-gamma development workflow currently specifies master seed 2026111500 with stride 10 and three replicates (population seeds 2026111500, 2026111510, 2026111520; validation seeds one greater). These are also disjoint from the candidate training and holdout sets. The full seven-gamma workflow has not yet been run, so these are planned rather than consumed seeds.

### Why this remains non-evaluable as a final gate

This was an independent manual scan of the listed retained artifacts, **not** a completed run of `audit_historical_seed_inventory.py` over a certified source register. It does not prove that all relevant historical Actions artifacts have been found; expired, missing, or older workflow artifacts may remain unenumerated. Binary files (including stored `.npz` reconstruction inputs) were not parsed for embedded seed metadata. The three numerical-diagnostics artifacts contain many bootstrap-seed records and were included in the overlap check, but the scan does not prove that every random stream in every historical experiment is represented in those artifacts.

Accordingly, the result is **PARTIAL_OVERLAP_SCAN_PASS / INVENTORY_COMPLETENESS_NOT_ESTABLISHED**. Keep the candidate schedule unapproved. Before training, finish the source register, identify and retrieve every relevant historical artifact (or document why it is excluded), inspect any additional seed-bearing binary/config formats, then run the repository inventory auditor over the extracted, hashed source roots. A no-overlap result from only the nine artifacts above is not permission to launch.


## Follow-up parser hardening and artifact-format review (2026-11-15)

The inventory implementation was hardened after a broader artifact-format review:

- Added NPZ support with `numpy.load(..., allow_pickle=False)` and scalar seed extraction from seed-named arrays.
- Added regression tests for integer seed arrays and rejection of object arrays without enabling pickle.
- Tightened the seed-field matcher to exclude plural suffixes such as `_across_seeds`. This prevents a numeric summary statistic such as `mean_pit_mean_across_seeds=0.5138` from being misread as a seed ID.
- Added a regression test for that exact class of summary-metric false positive.

The nine retained ZIPs listed above contain **81 NPZ files**. An independent inspection of their array keys found no seed-named arrays. This is now consistent with the scanner's NPZ support, but it does not cover other binary formats.

The strict follow-up scan also identified **four malformed JSON audit-report files** in the three numerical-diagnostics artifacts. Each contains a valid JSON object followed by a literal escaped newline token outside the object. The inventory auditor correctly fails closed on these files rather than silently skipping or repairing them. These are retained as source artifacts; any normalization must occur in a separate, hashed audit-staging copy with the original ZIPs and conversion provenance preserved.

A parseable-file-only pass over the nine archives found 201 seed-bearing JSON/CSV files, 2,784,484 seed-field occurrences and 258 distinct integer values; these are **provisional partial counts**, not a certified total, because the four malformed JSON reports prevent the strict inventory from returning PASS. The earlier 198 / 2,784,487 / 260 summary should likewise be treated as superseded and provisional until the corrected scanner is run against a complete, reviewed source register and all malformed reports are handled transparently.

Current gate remains **NOT EVALUABLE / INVENTORY COMPLETENESS NOT ESTABLISHED**. Do not authorize training or holdout from the no-overlap comparison of only the parseable subset.
