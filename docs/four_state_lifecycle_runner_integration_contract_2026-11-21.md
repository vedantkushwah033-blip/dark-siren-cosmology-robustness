# Four-state lifecycle integration design for the prior-predictive runner

**Status:** implementation contract; not evidence that the runner has been integrated.  
**Applies to:** synthetic quadratic-logit five-arm prior-predictive experiments only.

## Why the current runner cannot use the checkpoint safely yet

The current multi-replicate runner executes one child process for each replicate × curvature. That child computes all validation sizes and all five arms together. The parent runner then copies one child-level `SUCCESS` or `FAILED` status to every arm × validation-size row. Its output ledger is written only after the entire loop finishes.

The current aggregate auditor independently checks child artifacts and aggregate arithmetic, but requires all rows belonging to a child to share one status. This is incompatible with arm-level `NON_EVALUABLE`, per-cell restart, and the four-state checkpoint contract. Merely calling `mark_cell_started` around the existing child process would create false granularity: a single child failure would still affect many distinct planned cells.

## Frozen cell identity and plan

The immutable primary key is:

`(replicate, gamma, validation_size, arm)`

The plan must be fully materialized before the first simulation starts. Each row includes the cell key plus the exact replicate-level provenance needed to reproduce inputs: population/analysis seed, validation-frame seed, GW-event seed, validation-sample seed, H0-truth seed and truth, and the frozen design/configuration hash. The canonical plan hash covers the complete ordered plan. Duplicate keys, non-finite design values, invalid arms, and validation sizes outside the frozen frame must be rejected before checkpoint creation.

The checkpoint is initialized once, before any child work. The plan is immutable thereafter. A resumed process must verify the plan hash and all design/configuration fields against its requested invocation; it must not silently append, reorder, or reinterpret planned cells.

## Execution boundary required for genuine cell-level outcomes

Refactor the computational boundary so that each planned arm × validation-size cell can be evaluated and audited independently, while reusing common replicate-level generated inputs:

1. Generate and persist the replicate-level parent population, GW events, validation frame, and validation samples once per replicate × gamma input bundle (or once per replicate for quantities intentionally paired across gamma).
2. For each validation size, construct the five arm-specific analysis inputs from that immutable bundle.
3. Run one arm evaluator for one cell. It writes a cell-scoped result directory and a manifest binding the result to the plan hash, cell key, input hashes, code/config identifiers, and seed provenance.
4. Run an independent cell auditor. The auditor reconstructs the expected cell key and input provenance from the frozen plan and shared input bundle; it recomputes the posterior summary and diagnostics from the cell's saved inputs/outputs rather than trusting the runner's status field.
5. Only after the independent audit has completed may the parent transition that cell to a terminal state.

Shared input generation is a separate lifecycle stage. If it fails, all dependent planned cells receive `NOT_RUN/UPSTREAM_FAILURE` (or `FAILED/DATA_INTEGRITY_FAILURE` only where a cell was actually attempted and the failure is attributable to that cell). Do not label unattempted cells as attempted failures.

## Terminal-state mapping

| Status | Required evidence | Meaning |
|---|---|---|
| `SUCCESS` | finite valid posterior; independent audit PASS; result and audit paths; start/completion timestamps | Cell completed and passed all predeclared validity checks. |
| `NON_EVALUABLE` | independently audited reason; audit PASS; result and audit paths; timestamps | Cell ran and its result is scientifically non-evaluable under a rule frozen before outcomes are inspected. |
| `FAILED` | attempted cell; explicit execution/data-integrity error; completion timestamp; audit must not claim PASS | Cell execution or integrity failed. Preserve logs and partial artifacts. |
| `NOT_RUN` | no completion timestamp and no audit path; reason distinguishes not started, interrupted, upstream failure, or not scheduled | No valid terminal evaluation exists. |

A weight-instability rule may produce `NON_EVALUABLE` only if its diagnostic, threshold, scope (arm/cell), and handling were declared before the experiment. Do not invent a universal ESS threshold from the pilot results. In particular, do not retroactively exclude cells because their posterior is inconvenient.

The ledger's persisted `started_at` is the first attempt timestamp. A retry after interruption must preserve that timestamp and append attempt-level metadata to a separate immutable attempt log; it must not erase evidence of earlier attempts. Terminal rows are immutable. A retry policy must explicitly state whether an interrupted cell can be retried and how repeated attempts are represented.

## Aggregate and independent-audit contract

The aggregate auditor must stop enforcing one shared status per replicate × gamma child. Instead it must:

- reconcile exact planned keys against exactly one ledger row per key;
- validate lifecycle invariants and provenance for all four statuses;
- independently audit every SUCCESS and NON_EVALUABLE cell;
- require FAILED cells to have a reproducible failure record and forbid audit PASS;
- preserve NOT_RUN cells in denominators and report their reason codes;
- reconstruct success-only and failure-inclusive summaries from the cell ledger and cell artifacts;
- calculate paired arm contrasts only on explicitly matched evaluable pairs, while reporting the number and reasons for missing pairs;
- verify that pairing claims are supported by identical common-input hashes and the frozen seed schedule;
- hash every artifact in the run manifest and reject unplanned files that are presented as analysis inputs.

Do not collapse `NON_EVALUABLE`, `FAILED`, and `NOT_RUN` into one generic failure count. Report them separately by arm, curvature, and validation size. Success-only coverage/RMSE must always be accompanied by the planned denominator and failure-inclusive bounds.

## Restart and concurrency requirements

- On startup, validate the checkpoint and verify that no terminal row is overwritten.
- A cell is claimed atomically before execution; two workers must not both own the same cell.
- The attempt record must include worker/process identity, attempt number, start/end times, and output location.
- On restart, preserve terminal cells, identify active cells whose process died, and apply the frozen retry policy.
- A worker may not mark a cell SUCCESS or NON_EVALUABLE; only the coordinator, after independent audit, may commit the terminal transition.
- Test abrupt process termination between every transition and concurrent claims of the same key.
- Benchmark full-checkpoint JSON rewrites at the intended plan size before scaling to hundreds of truths.

## Required tests before any confirmatory run

1. Plan creation is complete before the first evaluator call; empty and duplicate plans are rejected.
2. A single cell can fail without changing sibling cells' states.
3. Each of the four states is independently audited and counted correctly.
4. NON_EVALUABLE reasons are rejected unless they belong to the frozen allowlist.
5. An audit failure cannot be represented as SUCCESS or NON_EVALUABLE.
6. A coordinator crash after evaluation but before audit leaves a recoverable nonterminal cell.
7. A coordinator crash after audit but before checkpoint commit does not create a false terminal state.
8. Restart preserves terminal cells and retries only eligible interrupted cells.
9. Two workers racing for one key cannot both claim it.
10. Independent aggregate reconstruction catches tampered plan hashes, input hashes, ledger states, cell outputs, and paired contrasts.
11. Full experiment output can be regenerated from the immutable plan and attempt records.
12. POSIX and Windows support policy is explicit. The current checkpoint module imports `fcntl` unconditionally and uses hard-link/directory-fsync semantics; this is not portable evidence for Windows support.

## Confirmatory gate

Do not launch the 100- or 300-truth calibration until the runner, cell auditor, aggregate auditor, restart logic, and tests above are integrated and pass exact-head CI. Freeze arm-specific evaluability rules, prior support, seeds, grids, quadrature, and all failure-inclusive metrics before the run. A successful infrastructure smoke is not calibration evidence.
