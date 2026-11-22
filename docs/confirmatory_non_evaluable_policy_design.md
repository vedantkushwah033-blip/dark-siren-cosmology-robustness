# Confirmatory non-evaluable policy: design and freeze checklist

**Status: design note only — not yet a preregistered policy.**

This document defines the distinction between a scientifically non-evaluable cell and an engineering failure for the synthetic low-redshift Euclidean dark-siren inclusion benchmark. It deliberately does **not** invent numerical cutoffs for weight concentration, posterior boundary mass, or quadrature error. Those cutoffs must be justified from completed exploratory diagnostics, frozen before confirmatory outcomes are inspected, and recorded in a versioned configuration.

The benchmark asks when estimating galaxy-catalogue inclusion probabilities helps dark-siren inference, and when misspecification, finite validation samples, or unstable inverse-probability weights remove or reverse that benefit. It is not an observational H0 measurement, a real-survey validation, a public LVK/GWTC analysis, or a claim of universal estimator superiority.

## 1. Non-negotiable classification rule

A cell may be labelled `NON_EVALUABLE` only when all of the following hold:

1. Its frozen input bundle passes byte-level, configuration, seed, and provenance checks.
2. The declared model and estimator were executed as specified; there is no unresolved software exception, corrupted artifact, or audit disagreement.
3. The failure condition is a predeclared property of the scientific realization or estimand, not a post-hoc response to an inconvenient estimate or coverage result.
4. The exact trigger, diagnostic values, attempted numerical remedies, and affected arm are recorded.
5. An independent auditor reconstructs the trigger from raw artifacts and confirms the classification.

If any condition is not met, the cell is **not** eligible for `NON_EVALUABLE`.

A poor posterior median, large error, low coverage, an unfavorable arm comparison, or a result that contradicts expectations is never by itself a reason to exclude a cell.

## 2. Lifecycle states

| State | Meaning | Permitted use |
|---|---|---|
| `SUCCESS` | The planned computation completed and the independent cell audit passed. | Included in the relevant scientific summaries. |
| `NON_EVALUABLE` | The valid realization met a frozen scientific/numerical non-evaluability rule and the independent audit verified it. | Reported explicitly; never silently dropped. |
| `FAILED` | Execution, data integrity, or independent audit failed. | Engineering failure; not a scientific non-evaluable observation. |
| `NOT_RUN` | The cell was not completed, including not scheduled or interrupted states. | Incomplete run; never treated as a result. |

A run with any `FAILED`, `NON_EVALUABLE`, or `NOT_RUN` cells is not a fully successful confirmatory run. The planned denominator and state counts must always be reported.

## 3. Failure taxonomy and decision sequence

Apply the following sequence in order. Do not relabel a software or provenance failure as a scientific exclusion.

### 3.1 Input and provenance

- Missing or changed shared-input files, hash mismatch, config/manifest disagreement, invalid seed-stream pairing, or wrong cell-to-bundle binding: `FAILED / DATA_INTEGRITY_FAILURE`.
- A missing bundle, missing plan entry, or unstarted cell: `NOT_RUN` with the corresponding scheduling reason.
- These conditions are never `NON_EVALUABLE`.

### 3.2 Execution and audit

- Unexpected exception, child process failure, invalid output schema, or an unhandled numerical exception: `FAILED / CHILD_PROCESS_ERROR`.
- Independent reconstruction mismatch, absent audit report, or failed audit: `FAILED / CHILD_AUDIT_FAILED`.
- A cell cannot be `SUCCESS` or `NON_EVALUABLE` without an independent audit pass.
- An implementation defect discovered after a run requires fixing the implementation and rerunning under a new immutable run identity; it must not be converted into a scientific exclusion.

### 3.3 Posterior existence and numerical validity

The confirmatory implementation must define a deterministic, arm-independent sequence of numerical checks before a cell can be classified:

1. Verify finite inputs and finite log-likelihood values on the declared support.
2. Evaluate normalization using stable log-domain arithmetic.
3. Check that the normalization is finite and strictly positive.
4. Repeat with the predeclared finer H0 grid and higher quadrature order.
5. Compare posterior summaries and normalization diagnostics against frozen numerical tolerances.
6. Run the independent artifact auditor.

If the declared posterior cannot be normalized after the prescribed refinement sequence, and the auditor verifies that the issue is intrinsic to the valid realization rather than an implementation error, the cell may be classified as `NON_EVALUABLE / INVALID_POSTERIOR_NORMALIZATION`.

If non-finite values arise from a software exception, overflow that the specified stable implementation should have prevented, missing values, or corrupted inputs, classify as `FAILED`, not `NON_EVALUABLE`.

The current lifecycle schema also recognizes `NONFINITE_LIKELIHOOD`, but this policy does not yet define an active trigger for it. Keep that reason code reserved until a frozen rule distinguishes a mathematically undefined likelihood for an otherwise valid realization from floating-point instability, an implementation defect, or invalid input. Do not emit this classification merely because one grid point evaluates to a non-finite value.

The exact refinement sequence and numerical tolerances remain **unfrozen**.

### 3.4 H0 support-boundary behavior

A posterior mode near a grid endpoint is not sufficient evidence of non-evaluability. Before any boundary classification:

- Recompute on the predeclared expanded H0 support.
- Recompute at the predeclared finer grid resolution.
- Record posterior mass in explicitly defined boundary regions, the mode location, and changes in quantiles and normalization.
- Distinguish a genuine support-limited inference from an insufficiently wide or coarse numerical grid.

Only a frozen boundary diagnostic, applied after support expansion and independently reconstructed, may trigger `NON_EVALUABLE / SUPPORT_BOUNDARY_MODE`. A grid that is simply too narrow or too coarse is an engineering/design defect and must be repaired before confirmatory execution.

The expanded support, boundary-region width, posterior-mass cutoff, and allowed summary shifts remain **unfrozen**.

## 4. Weight-instability policy must be arm-specific

Do not impose one universal effective-sample-size threshold across all arms or claim that a pilot-derived threshold is a theorem.

For every weighted arm, preserve and report at least:

- Kish effective sample size and entropy effective sample size;
- maximum normalized weight;
- fraction of total normalized weight held by the largest 1% and 5% of weights;
- number of positive finite inclusion probabilities and minimum fitted inclusion probability;
- support gaps and host-density diagnostics;
- the realized KDE bandwidth and sensitivity to the frozen bandwidth/grid checks.

The policy must distinguish:

- **C — known-probability IPW:** instability of exact-probability weights under the generating design;
- **D — MAP plug-in estimated IPW:** instability introduced by finite validation and the fitted inclusion model;
- **E — uncertainty-integrated estimated IPW:** instability after integrating coefficient uncertainty, which must not be presumed to cure extreme weights.

A/B are unweighted reference arms and cannot be excluded for IPW weight instability.

A weight-instability exclusion may be enabled only after exploratory results support a fixed, predeclared rule for the relevant arm and design regime. The rule must specify the diagnostic, threshold, whether it is evaluated per event or over the full host sample, how zero/near-zero fitted probabilities are treated, and whether a single trigger or a conjunction is required. Freeze the rule before confirmatory H0 errors, coverage, or arm rankings are inspected.

Until then, `PREDECLARED_WEIGHT_INSTABILITY` is a reserved reason code, not an active exclusion rule. Record the diagnostics but retain the cell in the analysis.

## 5. Required arm-level accounting

The lifecycle unit is exactly one `(replicate, gamma, validation_size, arm)` cell. Do not mark all five arms non-evaluable because one arm is unstable.

For each planned cell, retain:

- immutable cell identity and input/config/manifest hashes;
- all independent seed-stream identifiers and H0 truth;
- worker attempt history and first-start timestamp;
- raw weight and posterior diagnostics;
- every numerical refinement attempt;
- lifecycle state, reason code, and reason detail;
- independent audit report and its hash.

For every design cell, publish the full planned count and counts of SUCCESS, NON_EVALUABLE, FAILED, and NOT_RUN. Report scientific metrics both with the complete lifecycle denominator visible and, where meaningful, conditional on evaluable cells. Never present a conditional-on-success metric without its exclusion/failure counts.

## 6. Sensitivity analyses required before freeze

Before activating any scientific exclusion rule:

1. Reconstruct diagnostics independently from raw artifacts.
2. Show how each candidate threshold changes the number and identity of excluded cells by arm, curvature, and validation size.
3. Check whether exclusions are associated with H0 truth, posterior error, or other outcomes. This is a bias-risk audit, not a reason to tune thresholds to improve results.
4. Compare conclusions with all valid cells retained and with the candidate predeclared exclusion rule.
5. Freeze one policy version and hash it into the confirmatory plan.
6. Do not change the policy after inspecting confirmatory outcomes. Any later policy change requires a new confirmatory run identity and a clearly labelled sensitivity analysis.

## 7. Freeze checklist — unresolved items

These items must be resolved and reviewed before confirmatory calibration:

- [ ] Exact H0 grid, expanded support, and refinement sequence.
- [ ] Numerical normalization and posterior-summary convergence tolerances.
- [ ] Operational trigger and independent audit rule for `NONFINITE_LIKELIHOOD`, or explicit decision to leave it inactive.
- [ ] Operational definition of boundary regions and boundary-mass trigger.
- [ ] Arm-specific weight-instability triggers, if any; no universal ESS cutoff assumed.
- [ ] Whether any probability floor is part of the estimator itself (not an after-the-fact diagnostic).
- [ ] Independent auditor coverage for every active `NON_EVALUABLE` reason.
- [ ] Immutable policy version/hash included in the global plan and every cell artifact.
- [ ] A dry-run proving that every state/reason combination is classified and audited correctly.
- [ ] Final preregistration reviewed before confirmatory outcomes are generated.

## 8. Current decision

**No confirmatory run is authorized by this document.** Existing smoke and pilot results can inform the freeze, but they do not establish universal diagnostic thresholds or calibration. Until the checklist is complete, keep all candidate weight-instability thresholds exploratory, do not activate `NON_EVALUABLE` for weight concentration, and do not launch the planned 100/300 independent-truth-per-cell confirmatory study.
