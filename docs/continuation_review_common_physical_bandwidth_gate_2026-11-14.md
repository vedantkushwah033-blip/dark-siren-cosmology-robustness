# Continuation review: common physical bandwidth gate — 2026-11-14

## Review scope and provenance

This review was prepared from the live GitHub state of PR #54, branch
`research/five-arm-kde-kernel-reuse-v2-20261114`, at head
`bb065b1a3574932b28e09e63a1ba7359b3a5da3f`. It is a source/design review,
not a claim that an unobserved workflow or artifact passed.

The separate five-arm implementation PR #44 remains open and draft at
`2abb4c5275235fd97f87c0818f59b4dcc1e3223c`. Its older branch contains the
five-arm likelihood and posterior-quadrature work; PR #54 is based on that
research branch and carries the later KDE/bandwidth development. Do not treat
the two heads as interchangeable.

## Findings confirmed in source

1. **Common physical bandwidth is implemented.** The five-arm likelihood
   accepts exactly one of a shared Scott multiplier or an absolute redshift
   bandwidth SD. In absolute-width mode, every A–D density and every E
   coefficient-node density receives the same requested physical kernel SD.
2. **The base-population measure is propagated once.** The five-arm likelihood
   passes the complete arm weights to the KDE and uses an identity post-KDE
   merger weight. The training/holdout preregistration explicitly forbids
   applying the base factor a second time.
3. **Arm E's H0 likelihood marginalization is joint.** Node-specific complete
   selected likelihoods are combined with the shared coefficient quadrature
   weights using log-sum-exp. This is not an average of separately normalized
   H0 posteriors.
4. **Density-stage Arm E estimands are distinguished.** The preregistration
   names the posterior-mean density IAE as primary and posterior expected
   node-level IAE as secondary. They must not be substituted for one another.
5. **The current multi-replicate pilot is explicitly exploratory.** Its runner
   writes a complete planned-replicate ledger and per-replicate artifacts. The
   outer auditor checks hashes, seeds, ledger reconciliation, child audits,
   and exact aggregate reconstruction. A structural outer-audit PASS does not
   imply every replicate succeeded; failures remain in the ledger.

## Blocking evidence and interpretation

The latest recorded review of the 20-replicate, floor-0.10 common-physical-width
pilot reports 19/20 successful replicates. Replicate 17 failed during the
validation logistic posterior mode fit. The workflow stopped before the full
child/outer mathematical audit. Manifest integrity and selected posterior/KDE
reconstructions were checked, but these do not convert the run into a complete
PASS.

The latest review also reports that 0.012, the smallest candidate width, had
the lowest descriptive mean IAE among the seven tested widths in the 19
successful replicates. Because the run is incomplete and the full audit did
not pass, this is only a grid-design signal. It is not a selected bandwidth.
The candidate grid may need to extend below 0.012, but any amendment must be
made before generating training or holdout outcomes.

The latest PR #54 head contains a narrow posterior-fit tolerance adjustment
and an independent auditor regression test. The available commit-status query
did not return status records for this head, so this review does not assert
that exact-head CI has passed. The corrected 20/20 rerun and complete
independent child plus outer audits have not been established by the evidence
reviewed here.

## Decision: NO-GO for training/holdout and H0

Do not launch the preregistered 600-replicate training/holdout experiment yet.
Do not freeze a physical bandwidth. Do not start the confirmatory five-arm H0
experiment.

The next permitted sequence is:

1. Obtain exact-head CI for the posterior fitter, independent posterior
   reconstruction, KDE auditor, and regression fixture.
2. Rerun the same exploratory 20-replicate seed schedule without changing the
   DGP, metrics, or candidate grid; require 20/20 successful attempts.
3. Complete every child raw-input audit and the outer aggregate audit.
4. If the smallest width remains best, amend and freeze a wider candidate grid
   before any training/holdout data are generated.
5. Only then run the preregistered disjoint training and holdout seed ranges.
6. Freeze a width only if all 45 primary cells meet the preregistered
   evaluability threshold, the training minimum is interior, and all audit
   gates pass. Holdout cannot reselect it.
7. Keep H0 evaluation on separate replicate seeds and preserve H0-grid,
   coefficient-quadrature, calibration, coverage, and failure-rate audits as
   separate gates.

## Scientific scope

All statements concern synthetic low-redshift Euclidean methodology. This
review provides no LVK/GWTC/public-data analysis, no observational H0
measurement, no real-survey validation, no H0 calibration or coverage result,
and no estimator-superiority claim.
