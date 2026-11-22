# Global arm-cell coordinator

This runner schedules independently auditable arm × validation-size cells across
multiple pre-generated replicate/curvature bundles while persisting one global
four-state lifecycle checkpoint.

## 1. Prepare frozen input bundles

Use `scripts/generate_shared_inputs.py` once for each independent
replicate/curvature pair. Each bundle must remain inference-free and immutable.
Within a replicate, the parent, event, validation-frame, validation-sample
streams and H0 truth must match across curvature bundles to preserve paired
comparisons. Across independent replicates, the four seed-stream tuple must
differ. For prior-predictive studies, draw H0 truths independently by replicate;
repeated truths are appropriate only for an explicitly fixed-truth design.

Example bundle-plan JSON:

```json
{
  "bundles": [
    {"replicate": 1, "input_dir": "runs/shared/rep001_gamma_m1"},
    {"replicate": 1, "input_dir": "runs/shared/rep001_gamma_0"},
    {"replicate": 2, "input_dir": "runs/shared/rep002_gamma_m1"},
    {"replicate": 2, "input_dir": "runs/shared/rep002_gamma_0"}
  ]
}
```

Each bundle's `config.json` and frozen `shared_input_manifest.json` provenance
are checked against one another before the plan is accepted. The global
checkpoint also binds every cell to the exact config and shared-manifest hashes.

## 2. Execute all cells

```bash
python scripts/run_global_arm_cells.py \
  --bundle-plan runs/global_bundle_plan.json \
  --output-dir runs/global_arm_cells \
  --arms A_catalogue_only,B_parent_oracle,C_known_probability_ipw,D1_estimated_ipw,E1_integrated_ipw \
  --max-workers 4
```

One checkpoint contains every planned (replicate, gamma, validation size, arm)
cell. Each worker still writes a separate attempt directory and invokes the
independent per-cell auditor before a SUCCESS state can be committed.

Resume an interrupted run with the identical bundle plan and options:

```bash
python scripts/run_global_arm_cells.py \
  --bundle-plan runs/global_bundle_plan.json \
  --output-dir runs/global_arm_cells \
  --arms A_catalogue_only,B_parent_oracle,C_known_probability_ipw,D1_estimated_ipw,E1_integrated_ipw \
  --max-workers 4 --resume
```

The immutable checkpoint plan must match exactly. Terminal outcomes are not
silently rerun; interrupted NOT_RUN cells are recovered through the lifecycle
store. Preserve the original bundle directories for the entire run.

## 3. Independently audit

```bash
python scripts/audit_global_arm_cells.py --output-dir runs/global_arm_cells
```

The aggregate auditor reconstructs CSV exports from authoritative checkpoint
and cell-summary records, rechecks bundle/config/manifest provenance, validates
the run manifest, and repeats independent per-cell audits for every SUCCESS.

## Scope and limitations

This is engineering orchestration for a synthetic low-redshift Euclidean
methodology benchmark. It is not an observational H0 measurement, LVK/GWTC
analysis, real-survey validation, or evidence of estimator superiority.
The arm-specific NON_EVALUABLE policy remains unapproved; the global auditor
fails closed if such a row appears. A green smoke test is not confirmatory
calibration. Do not launch the 100/300-truth confirmatory design until
independent-truth seed planning, failure/evaluability policy, and all numerical
convergence/calibration gates are separately frozen and reviewed.
