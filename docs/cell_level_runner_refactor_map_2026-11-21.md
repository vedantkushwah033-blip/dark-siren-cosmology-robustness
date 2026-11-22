# Cell-level runner refactor map

**Date:** 2026-11-21  
**Status:** source-inspection plan; no production runner refactor is claimed here.  
**Reviewed revision:** `146eb7354a02396ec73c3aec313a40e54743e8e6`

## Exact current execution boundary

The multi-replicate coordinator is `scripts/run_quadratic_five_arm_replicates.py`. It currently launches `scripts/run_quadratic_five_arm_smoke.py` once for each replicate × gamma. The child receives the full ordered list of validation sizes. Its `run()` function generates the parent population, GW events, catalogue, and validation frame, persists the common raw inputs, then loops over validation sizes. Within each size it calls `build_five_arm_weights(...)` and `matched_five_arm_kde_selected_log_likelihood(...)`, which jointly produce all five arms. It writes one `posterior_summaries.csv`, one `weight_diagnostics.csv`, and one likelihood archive for the entire child.

The current independent child auditor is `scripts/audit_quadratic_five_arm_smoke.py`. It reconstructs the complete child, including all validation sizes and all five arms. The replicate-level auditor then consumes the child report and copies one child outcome to every arm × validation-size ledger row. Therefore the current status is not a genuine cell-level status.

## Refactor target: three explicit layers

### Layer 1 — immutable shared-input bundle

Extract the pre-analysis portion of the smoke runner into a generator that accepts the frozen replicate/gamma seed provenance and writes a content-addressed input bundle:

- parent population and true inclusion probabilities;
- observed catalogue;
- GW events and simulator truth/trace;
- complete validation-frame truth and observed frame;
- configuration, seed schedule, and hashes.

The common-input bundle is generated once per replicate × gamma. It must preserve the existing cross-gamma pairing guarantees: population/GW streams are shared where specified, while the catalogue/validation streams follow their dedicated deterministic seeds. The plan must state which hashes are expected to match across gamma and which are expected to differ.

### Layer 2 — one arm × validation-size evaluator

Create a pure evaluator with inputs limited to:
- a verified shared-input bundle;
- one validation size;
- one arm identifier;
- frozen H0/redshift grids, quadrature, bandwidth, and likelihood settings.

The existing `build_five_arm_weights` and `matched_five_arm_kde_selected_log_likelihood` are joint five-arm APIs. The first safe implementation can reuse their joint calculation internally but must extract and persist exactly one arm's density, likelihood, posterior summary, and diagnostics under a cell-scoped directory. This is computationally redundant but gives honest independent lifecycle and audit boundaries. Optimize to arm-specific kernels only after the cell semantics and independent reconstruction are tested.

The evaluator must not write terminal lifecycle state. It writes immutable cell inputs/outputs plus a manifest containing the cell key, plan hash, common-input hashes, settings, and code identifier.

### Layer 3 — independent cell auditor and coordinator

Create an independent cell auditor that accepts only the frozen plan row, shared-input bundle, and cell artifacts. It must independently rebuild the selected arm's likelihood and posterior summary. It must not import the evaluator's likelihood or weight-building implementation. It emits a cell audit report.

Only the coordinator may transition the checkpoint after a cell audit:
- audited finite valid posterior → `SUCCESS`;
- audited, predeclared scientific non-evaluability → `NON_EVALUABLE`;
- execution/integrity failure → `FAILED`;
- not attempted/interrupted/upstream-blocked → `NOT_RUN`.

The coordinator must preserve completed cells on restart and record each attempt. A single cell failure must not alter sibling cell states.

## Migration order

1. Keep the existing smoke runner and auditor unchanged as regression references.
2. Add shared-input generation and a round-trip/hash test.
3. Add one-cell evaluation for each arm and at least two validation sizes; compare every output with the corresponding slice of the legacy five-arm child.
4. Add independent one-cell reconstruction. Require agreement with the legacy independent child auditor, while ensuring the new auditor has no production-likelihood imports.
5. Add coordinator/checkpoint orchestration with one-cell failure isolation, restart, and same-cell concurrency tests.
6. Rewrite aggregate auditing around exact `(replicate, gamma, validation_size, arm)` keys and all four lifecycle states. Reconstruct failure-inclusive bounds and paired contrasts from audited cell artifacts.
7. Benchmark checkpoint serialization, shared-input reuse, and redundant joint-kernel execution at the intended plan sizes before optimizing or launching calibration.

## Verification gate

The latest PR #162 exact-head CI is green on Python 3.10, 3.11, and 3.12, and its matched quadratic-logit H0 smoke plus independent child audit pass. Those checks validate the current repository and checkpoint primitives; they do **not** validate the proposed cell-level architecture. No confirmatory calibration should run until steps 1–7 are implemented, exact-head CI is green, and the frozen protocol and failure rules are reviewed.

The research remains a synthetic low-redshift methodology benchmark. It is not an observational H0 measurement, an LVK/GWTC analysis, or evidence that one estimator is universally superior.

## Validation update — 2026-11-22

**Reviewed revision:** `446249c45b91588dec8a7144ae8a30f73404c052`  
**Scope:** repository source inspection plus GitHub Actions evidence; no local test execution is claimed.

At this revision, the production boundary described above remains unchanged: the multi-replicate coordinator still launches one five-arm child per replicate × gamma, and one child outcome is copied to every arm × validation-size ledger row. The shared-input bundle and checkpoint primitives therefore remain preparatory infrastructure, not cell-level production execution.

Exact-head workflow evidence observed on this revision:

- The degree-matched quadrature sensitivity pilot completed successfully. Its workflow ran the targeted orchestration/quadrature test subset (**20 passed**), then evaluated 6 paired child runs (3 curvature values × 2 validation sizes) with 6/6 child audits passing. It recorded 35 input-pairing checks, 50 likelihood-invariance checks, 12 quadrature-convergence rows, and a maximum order-9 versus order-7 density L1 difference of approximately (1.25\times10^{-5}). This is a small engineering/numerical sensitivity pilot (400 parent objects, 4 detected events, 41 H0 grid points), not a calibration experiment and not evidence of estimator ranking.
- The matched quadratic-logit H0 smoke and its independent artifact audit completed successfully. The audit reports arithmetic/artifact integrity only; its own claim limit explicitly excludes calibration, SBC, coverage, and real-data validation.
- The full repository test-matrix workflow and paired-degree multi-replicate workflow were still pending when this note was written. The prior-predictive Stage0 workflow was skipped.

The numerical pilot is useful evidence that the paired-degree likelihood calculation and selected quadrature orders are behaving consistently at this small scale. It does not resolve the production lifecycle gap. The next implementation gate remains a true one-cell evaluator, independent one-cell auditor, checkpoint-owning coordinator, and aggregate auditor that reconciles all four states at exact cell-key granularity. Do not launch confirmatory calibration before that gate and exact-head CI pass.
