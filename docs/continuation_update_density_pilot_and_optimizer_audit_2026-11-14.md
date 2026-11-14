# Continuation update: density-pilot branch and optimizer audit — 2026-11-14

## Live state checked

This update is based on the repository's separate branch
`research/absolute-width-distribution-pilot-20261114` and workflow run
`43615384270`, code head
`3998718b12d33fa05debd13dec5b94ef459be054`. The branch is intentionally
separate from PR #54's KDE implementation branch; no merge or main-branch
change is made here.

At the time of this check, the three-floor workflow jobs were:
- floor 0.00: completed successfully;
- floor 0.10: completed successfully;
- floor 0.25: still in progress.

Artifacts were available for the completed floor-0.00 and floor-0.10 jobs.
The floor-0.25 artifact had not yet been produced. Therefore, the full
three-floor rerun and its complete artifact audit are still pending.

## Optimizer implementation reviewed

The production `fit_inclusion_posterior_quadrature` now has a deterministic
exact-Hessian Newton refinement after the `trust-exact` mode optimizer.
The final posterior score threshold remains strict at infinity-norm
`1e-6`; the refinement does not simply relabel an optimizer stop as success.
It uses a positive Newton decrement, a damped line search, and accepts a step
only when the objective decreases by the Armijo condition or the score norm
strictly improves. The posterior information matrix must remain positive
definite before the Laplace proposal and quadrature nodes are constructed.

The independent artifact auditor uses a different primary optimizer (BFGS),
then its own independently coded damped-Newton polish. It reconstructs the
Gaussian-prior logistic posterior, Laplace covariance, tensor Gauss-Hermite
nodes and normalized importance weights from observable validation records.
It does not import the production inclusion-posterior fitter. The independent
reference still enforces its own final score tolerance and checks posterior
nodes and weights against persisted outputs.

The regression suite includes a seeded reproduction of the previously
near-stagnating validation sample, plus checks for posterior-mean density
metrics, posterior expected node loss, missing ESS semantics for aggregate
Arm E densities, and rejection of metric/posterior corruption after manifest
refresh.

## Workflow evidence

Run `43615384270` had the focused density-pilot test step pass on the
three-floor workflow's code head. At the latest live poll:
- floor 0.00: job completed successfully; artifact
  `14121687969` available;
- floor 0.10: job completed successfully; artifact
  `14122536928` available;
- floor 0.25: job remained in progress; no artifact yet.

This is not a completed three-floor experiment. It is not yet evidence that
all 60 planned attempts succeeded or that all three floor-cell audits passed.
The two available artifacts still require direct download and inspection; a
workflow job's green status alone is not treated as the full independent
mathematical audit.

## Next actions and gates

1. Wait for the floor-0.25 job to finish and inspect its conclusion and logs.
2. Download all three artifacts and verify their outer manifests, every child
   manifest, raw-input reconstruction, seed schedule, and aggregate tables.
3. Report planned, successful, failed and non-evaluable replicate counts
   separately for each floor. An outer audit PASS must never hide a failed
   planned replicate.
4. Reproduce any failed seed exactly before modifying the optimizer again.
5. Only after the three-floor expanded-grid pilot is fully audited, decide
   whether the candidate physical-width range is adequate. This pilot has no
   training/holdout split and cannot select a final bandwidth.
6. Freeze a preregistered candidate grid and disjoint training/holdout seeds
   before running any selection study.
7. Keep the confirmatory five-arm H0 experiment blocked until density
   estimation, inclusion posterior quadrature, and all relevant numerical
   convergence gates are satisfied.

## Scientific scope

This is synthetic low-redshift Euclidean density-recovery methodology only.
It provides no H0 posterior, calibration, PIT, coverage, observational
measurement, public LVK/GWTC analysis, or estimator-superiority evidence.
