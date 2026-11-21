# Confirmatory-run lifecycle ledger specification

**Status:** design gate; not yet implemented.  
**Applies to:** the synthetic five-arm prior-predictive H0 experiment only.  
**Purpose:** prevent a confirmatory run from silently converting execution, numerical, or inferential failures into successful posterior summaries.

## Why this gate is required

The current replicate orchestrator has one status for an entire child run and copies it to every arm/validation-size row. The child runner currently emits posterior summaries with `status=SUCCESS` for all five arms after the common likelihood calculation. The aggregate auditor requires a single status across every arm and validation size in a child. This is adequate for an all-or-nothing engineering smoke, but it cannot represent arm-specific inferential non-evaluability, a partially completed validation-size loop, or an interrupted run whose planned cells were never attempted.

Do not launch the 100-truth-per-cell calibration stage until the ledger contract below is implemented and independently tested.

## Status semantics (mutually exclusive)

| Status | Meaning | Required evidence | Included in success-only metrics? |
|---|---|---|---|
| `SUCCESS` | The planned arm × validation-size cell completed, passed its independent audit, and met every predeclared numerical/inferential evaluability rule. | Finite likelihood/posterior summaries; valid support and normalization; audit PASS; empty failure reason. | Yes |
| `NON_EVALUABLE` | Computation and audit completed, but the cell failed a predeclared inferential validity rule. This is a scientific outcome, not an execution error. | Raw outputs retained; audit PASS for the computation; one or more stable reason codes; diagnostic values where available. | No |
| `FAILED` | Execution, data-integrity, or independent-audit failure prevented a trustworthy result. | Exception/exit code or audit failure; error class and bounded message; audit status NOT_PASS where applicable. | No |
| `NOT_RUN` | The cell was planned but never started or was not completed because of interruption/skip. | Planned-cell record; explicit reason such as `RUN_INTERRUPTED`, `UPSTREAM_FAILURE`, or `NOT_SCHEDULED`; no fabricated output. | No |

A missing summary row must never be interpreted as `SUCCESS`. A cell can be `NON_EVALUABLE` only when its underlying computation is reproducible and the independent audit passes. An audit failure is `FAILED`, not `NON_EVALUABLE`.

## Unit of accounting

The canonical ledger key is:

`replicate × gamma × validation_size × arm`

Every planned key must appear exactly once in `planned_cells.csv` and exactly once in the final `cell_ledger.csv`. The final ledger must contain all four statuses, including `NOT_RUN` rows for planned cells not reached. Do not use a single child-level status as a substitute for cell-level status.

Required ledger fields:

- identity: `cell_id`, `replicate`, `gamma`, `validation_size`, `arm`;
- provenance: master seed and all five stream seeds (population, validation-frame, H0 truth, GW events, within-frame validation sample);
- truth: `h0_true`, `h0_truth_mode`;
- lifecycle: `status`, `reason_code`, `reason_detail`, `started_at`, `completed_at`;
- audit: `audit_status`, `audit_report_path`, `audit_error`;
- outputs: relative paths or nulls for posterior summary, likelihood curve, and diagnostics.

Reason codes must be enumerated and stable (for example `NONFINITE_LIKELIHOOD`, `INVALID_POSTERIOR_NORMALIZATION`, `SUPPORT_BOUNDARY_MODE`, `PREDECLARED_WEIGHT_INSTABILITY`, `CHILD_PROCESS_ERROR`, `CHILD_AUDIT_FAILED`, `RUN_INTERRUPTED`). Do not introduce a reason code after seeing confirmatory outcomes without versioning the protocol.

## Arm-specific evaluability

The preregistration must freeze exact thresholds before the confirmatory run. At minimum, check finite log-likelihood values, finite posterior quantiles, posterior normalization, monotone CDF/quantiles, and valid prior-support bounds. Weight-concentration diagnostics must be retained for every arm and validation size. If a weight-instability rule is used to label an IPW cell `NON_EVALUABLE`, its threshold and rationale must be predeclared; do not invent a universal ESS cutoff from the current pilot.

A stability warning alone is not automatically a failure. Keep continuous diagnostics (Kish ESS, entropy ESS, maximum normalized weight, top-1%/5% mass, support gaps) even when the cell remains evaluable. Separate numerical validity from scientific caution.

## Aggregation rules

1. Report planned, started, completed, audited, successful, non-evaluable, failed, and not-run counts for every design cell and arm.
2. Report success-only bias/RMSE/coverage/interval width with the denominator explicitly shown.
3. Report failure-inclusive bounds using the frozen support-derived error and interval-width bounds. Treat both `FAILED`, `NON_EVALUABLE`, and `NOT_RUN` as unresolved outcomes for conservative bounds; also show them separately rather than collapsing their counts.
4. Never label success-only coverage as unconditional calibration.
5. Paired contrasts use only replicate/design pairs where both compared cells are `SUCCESS`; report the paired denominator and the excluded-pair counts by status.
6. The independent auditor must reconstruct every aggregate from the cell ledger and retained cell outputs, and reject duplicate/missing planned keys, unknown statuses/reasons, invalid transitions, and summaries attached to non-success cells.

## Interruption-safe execution

The runner should write the immutable plan before starting any child. Each cell/child transition should be persisted atomically (write temporary file, flush, then replace) so an interrupted process leaves a recoverable ledger. On restart, verify existing artifacts against hashes and independent audits before reusing them; never silently overwrite or rerun a completed cell under the same run ID. At finalization, all planned keys not reached become `NOT_RUN` with `RUN_INTERRUPTED` or a more specific frozen reason. The manifest must include the plan, ledger, child artifacts, audit reports, and final aggregate tables.

The existing one-child-per-(replicate, gamma) architecture may remain for generation efficiency, but it must expose arm × validation-size completion records. A child-level process failure can mark all not-yet-trustworthy cells `FAILED`; a successfully audited child can still contain individual `NON_EVALUABLE` cells.

## Required tests before Stage 1

- all four statuses round-trip through plan, ledger, summaries, aggregates, and independent audit;
- each planned key occurs exactly once; duplicate, missing, and extra keys are rejected;
- `NON_EVALUABLE` requires audit PASS and a frozen reason; audit failure cannot be downgraded;
- `FAILED` and `NOT_RUN` never enter success-only metrics;
- failure-inclusive bounds use the correct unresolved count and are conservative;
- paired contrasts exclude non-success pairs and report their denominator;
- interruption after plan creation and after a subset of children leaves recoverable `NOT_RUN` rows;
- restart validates hashes/audits and does not overwrite verified outputs;
- mutation of status, reason, seed, truth, or output hash is detected by the independent auditor;
- full test matrix passes on supported Python versions.

## Gate decision

**Current status: BLOCKED.** The existing runner/auditor supports coarse `SUCCESS`/`FAILED` child-level accounting and success-only versus failure-inclusive summaries, but it does not yet implement this four-state, cell-level, interruption-safe contract. The ten-truth Stage 0 smoke remains engineering evidence only. No 100- or 300-truth confirmatory run should be launched until implementation, independent reconstruction, and preregistration are complete.
