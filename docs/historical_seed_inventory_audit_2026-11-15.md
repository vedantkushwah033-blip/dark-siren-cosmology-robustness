# Historical seed inventory audit

Status: implementation added; **no historical inventory has yet been certified**. This is a fail-closed tool for auditing extracted artifacts, not an automatic claim that the supplied directories are complete.

## Purpose

Before the proposed density training and holdout experiments, compare every historical random seed actually recorded in retained development artifacts against both future seed blocks. Comparing only master-seed formulas is insufficient: old child configurations may have overridden seeds, and manifests may contain additional random streams.

## Input contract

The tool accepts one or more extracted artifact roots and a JSON schedule with non-empty `training` and `holdout` arrays. Each schedule row must explicitly provide integer `population_seed` and `validation_seed` values; any additional random stream must be represented by another integer field whose name contains `seed` (for example, `bootstrap_seed`). Seed values must be unique across all streams and replicates within a block, and the blocks must be mutually disjoint.

The scanner currently parses JSON and CSV files. It recursively extracts scalar values from field names containing `seed`, including nested JSON records and seed-named CSV columns. It hashes every supported source file with SHA-256 and records source-relative paths, byte counts, field names, and seed occurrences.

YAML files are deliberately unsupported and cause a non-evaluable result. Convert them to reviewed JSON/CSV without changing values, retain the originals in the audit package, and document the conversion. Malformed files and non-integer seed values fail closed rather than being skipped.

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
  --schedule configs/density_training_holdout_seed_schedule.json \
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
