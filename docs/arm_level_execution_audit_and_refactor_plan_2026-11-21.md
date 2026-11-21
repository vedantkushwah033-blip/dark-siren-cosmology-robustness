# Arm-level execution audit and refactor plan (2026-11-21)

## Status and evidence boundary

This is a source-level review of PR #157 at inspected head
`c2ab97cc43ceaa7e1b4dee859941c11f4fdf913b`. It is a design/audit note, not
a claim that arm-level isolation has been implemented. At this head, the
repository-wide test workflow, four-state integration workflow, degree-two
pilot workflow, and H0 smoke workflow report success. The paired 20-replicate
degree pilot is still running; its independent artifact audit has not yet
started. These checks establish engineering behavior, not confirmatory
calibration.

## Confirmed architectural findings

### 1. The ledger is cell-granular; execution is not

The planned key is
`(replicate, gamma, validation_size, arm)`, but
`scripts/run_quadratic_five_arm_replicates.py` groups work by
`(replicate, gamma)`. It launches one
`run_quadratic_five_arm_smoke.py` child with every validation size and all
five arms. The parent marks all cells started together and records the whole
block in one `record_cell_outcomes` transaction.

Consequences:
- A process interruption replays the entire replicate × gamma block.
- A child-process or audit exception marks every validation-size × arm cell
  in that block FAILED.
- The current implementation does not provide arm-level retry or isolation,
  despite having one ledger row per arm and validation size.

### 2. NON_EVALUABLE is representable but not produced by the runner

The orchestrator explicitly requires every child summary row to have status
SUCCESS. Any other status raises an exception and the catch block records the
whole child block as FAILED. The child smoke currently writes SUCCESS for
every arm after a shared five-arm likelihood calculation.

The lifecycle schema permits NON_EVALUABLE, but there is no end-to-end
scientific classification path, evidence contract, or per-cell retry path.
This is an integration gap, not a failed CI test.

### 3. The five arms share a computational failure domain

In `scripts/run_quadratic_five_arm_smoke.py`, each validation size first
builds all five arm weights, then calls the matched five-arm likelihood
routine, then loops over the returned arm likelihoods to create summaries.
There is no per-arm exception boundary around weight construction, likelihood
evaluation, posterior normalization, or summary creation. An exception in
one shared stage aborts the child.

A safe refactor must preserve shared simulation inputs and paired design
while separating arm-specific calculations and outputs. Merely launching
one child per ledger row is not sufficient if each child still calculates
all five arms through the same shared failure path.

### 4. Aggregate lifecycle counts were improved, but execution is still block-level

Commit `c2ab97cc` adds separate
`non_evaluable_replicates`, `failed_replicates`, and
`not_run_replicates` fields to each gamma × validation-size × arm aggregate.
The independent aggregate auditor reconstructs these three counts from the
lifecycle ledger. A regression test covers SUCCESS, NON_EVALUABLE, FAILED,
and NOT_RUN.

The legacy `failed_or_non_evaluable` field remains for compatibility and
counts all non-success outcomes, including NOT_RUN. This accounting
improvement does not change the block-level execution boundary described
above.

Success-only metrics are conditional on SUCCESS; they must not be described
as unconditional calibration. Failure-inclusive bounds are useful only when
their finite error and interval-width bounds are justified by the frozen prior
support and estimand definition.

### 5. The aggregate auditor still assumes child-block outcome patterns

The aggregate auditor accepts a child block if all cells are SUCCESS or
NON_EVALUABLE with audit PASS, or if all cells are FAILED with no audit PASS.
It does not yet establish a fully general mixed-outcome model with
independently audited cell artifacts. This matches the current all-or-nothing
child layout.

## Recommended target architecture

Keep simulation and estimation separate:

1. **Replicate input stage** — generate and persist immutable parent
   population, GW events, validation frame, validation samples, H0 truth, and
   named RNG streams once per replicate × curvature. Record hashes and a
   versioned input manifest.
2. **Cell execution stage** — execute one
   `(replicate, gamma, validation_size, arm, degree)` estimation cell.
   Each cell reads immutable shared inputs and writes to a unique cell
   directory. Degree 1 and 2 must use identical shared inputs and truth.
3. **Cell audit stage** — independently reconstruct the cell's weights,
   likelihood, posterior normalization, quantiles, and diagnostics from the
   immutable inputs. A cell may become SUCCESS only after this audit passes.
4. **Aggregate stage** — consume only audited cell records; never infer a
   cell's state from a child-process return code or another arm's outcome.

The execution key should include degree when degrees are run in the same
checkpoint. If degree is instead frozen per run, record it in immutable plan
provenance and use a distinct run identity.

## Lifecycle classification contract

- **SUCCESS**: the predeclared estimand is defined; all required numerical
  outputs are finite and valid; the independent cell audit passes.
- **NON_EVALUABLE**: a predeclared, scientifically justified condition makes
  the target cell's estimand unavailable under the frozen protocol. Preserve
  raw inputs and diagnostic evidence, and independently audit the reason.
- **FAILED**: software exception, optimizer implementation failure, corrupted
  artifact, failed independent audit, or unexpected numerical defect. Do not
  relabel these as NON_EVALUABLE.
- **NOT_RUN**: no terminal attempt has been committed; retain interruption or
  not-scheduled reason.

Do not classify low ESS, high maximum weight, concentrated weights, or a
boundary mode as NON_EVALUABLE merely because the result looks undesirable.
Those are diagnostics unless an ex-ante scientific rule establishes that the
estimand is undefined. Thresholds must be frozen before confirmatory runs.

## Required invariants and tests

1. Inject a failure in exactly one arm and verify all other arms and sizes
   complete and remain auditable.
2. Inject a failure in one validation size and verify other sizes complete.
3. Kill a worker before and during atomic cell commit; resume without rerunning
   committed cells.
4. Replay the same cell concurrently; only one terminal transition may win.
5. Corrupt one cell output and refresh the ordinary manifest; the independent
   auditor must still detect disagreement with raw inputs/reconstruction.
6. Change an immutable input/configuration hash during resume; reject the run.
7. Verify every planned cell has exactly one terminal status at finalization.
8. Recompute status-specific counts, success-only metrics, and failure bounds
   independently from the ledger.
9. Report paired-comparison denominators and reasons for missing pairs; do not
   silently discard non-evaluable or failed pairs.
10. Test deterministic seed invariance across arm order, worker count, retries,
    and degree.

## Statistical reporting requirements

For every gamma × validation-size × arm × degree, report planned,
SUCCESS, NON_EVALUABLE, FAILED, and NOT_RUN counts separately. Report
success-only bias, RMSE, coverage, and interval width with explicit
denominators. Also report failure-inclusive bounds and the assumptions behind
them. Paired contrasts must include their eligible-pair denominator and
excluded-pair counts by lifecycle reason.

Keep numerical convergence, SBC/rank calibration, empirical prior-predictive
coverage, and estimator accuracy as distinct validation claims. A smoke test
or a small pilot cannot establish confirmatory calibration or estimator
superiority.

## Gate recommendation

- Checkpoint persistence, restart integration, and status-specific aggregate
  counts: useful engineering foundation; exact-head CI success does not remove
  the execution-granularity limitation.
- Arm-level execution, mixed lifecycle handling, and independently audited
  per-cell artifacts: required before confirmatory work.
- Confirmatory calibration: **NO-GO** until the architecture, frozen
  NON_EVALUABLE rules, independent per-cell auditor, and preregistered
  statistical protocol pass their tests.
