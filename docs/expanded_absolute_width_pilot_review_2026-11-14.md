# Expanded absolute-bandwidth pilot: audit review and required rerun

Date: 2026-11-14

## Existing expanded-grid run

Workflow run `43574613209` completed all three floor jobs and the independent child/outer auditors on the same paired-seed 20-replicate schedule used by the first three-floor pilot. It expanded the physical bandwidth grid from 0.012–0.040 to:

`0.006, 0.008, 0.010, 0.012, 0.016, 0.024, 0.032, 0.040, 0.050, 0.060, 0.080`.

The run used N=4,000 parent galaxies, validation sizes 100/300/1,000, quadrature order 5, and a 401-point redshift grid.

| Catalogue floor | Successful | Failed | Independent audit |
|---:|---:|---:|---|
| 0.00 | 20/20 | 0 | PASS |
| 0.10 | 20/20 | 0 | PASS |
| 0.25 | 19/20 | 1 | PASS, with failure preserved |

The floor-0.25 failure is replicate 5 (population seed 20261265, validation seed 20261266), a degree-one inclusion-posterior mode fit at validation size 1,000. The trust-exact optimizer stopped with score infinity norm (1.6759010191\times10^{-6}), above the then-required (10^{-6}). The outer audit PASS means the successful children, explicit failed-replicate ledger, manifests, and aggregate rows reconcile; it does not mean 20/20 completed.

Artifacts and independent audit reports:
- No floor: artifact 14106387924; 20 child audits; 19,800 metric rows and 1,740 Scott-width rows reconstructed.
- Floor 0.10: artifact 14107380789; 20 child audits; 19,800 metric rows and 1,740 Scott-width rows reconstructed.
- Floor 0.25: artifact 14105479927; 19 child audits; 18,810 metric rows and 1,653 Scott-width rows reconstructed.

## A–D-only descriptive pattern from the earlier expanded grid

The existing artifact can still describe the fixed arms A–D, subject to the incomplete floor-0.25 replicate count:

- For floors 0.00 and 0.10, the catalogue-only arm A had its smallest observed mean IAE at the upper candidate edge 0.080.
- Across all three floors, the parent-oracle B and known-probability IPW C arms had their smallest observed mean IAE around 0.016.
- The MAP plug-in D arm had its descriptive minimum around 0.016–0.024 for floors 0.00 and 0.10, and at 0.012 for floor 0.25.
- For floor 0.25, the single failed replicate was included in the planned denominator by assigning the normalized-density IAE upper bound 2. This changes the numerical means but not these descriptive minimizing grid points.

This is a warning against assuming a universal common bandwidth is optimal for every arm. It is not a final conclusion: the primary posterior-mean Arm E metric was absent from these artifacts, and floor 0.25 is 19/20. The predeclared equal-cell objective across arms and validation sizes must remain the basis for any later bandwidth selection; do not tune each arm separately after seeing the data.

## Important estimand defect discovered during review

The earlier runner persisted the posterior-mean Arm E density curve but did not write its own IAE/ISE/boundary-IAE metric row. It wrote:
1. node-specific Arm E losses; and
2. posterior-weighted expected node-specific losses.

These are not the primary posterior-mean density loss. By convexity, loss evaluated on the posterior-mean density is generally no greater than posterior expected node-specific loss. Treating the second quantity as a substitute would answer a different question.

Therefore, the existing expanded-grid artifacts are useful for:
- the A–D fixed-arm density summaries;
- the Arm E node-specific diagnostic summaries;
- the posterior expected node-loss diagnostic;
- realized Scott-width distributions.

They must **not** be used to report the primary posterior-mean Arm E density performance or to select a common bandwidth across all five primary arms.

The runner now writes a separate `E_posterior_mean_density` row at node = -1, computed from the quadrature-weighted mean density curve. The `E_posterior_weighted_node_metric` remains a separate secondary estimand. Kish ESS is missing for both aggregate E loss summaries because neither is a single weighted sample; node-specific E rows retain node-specific ESS.

The independent child auditor now independently reconstructs the posterior-mean density and its loss row. The outer auditor checks the exact missingness pattern for aggregate ESS values. A regression test verifies the convexity inequalities and corrupts the posterior-mean IAE while refreshing the manifest; the independent auditor must still reject it.

## Inclusion-posterior optimizer follow-up

The earlier floor-0.25 failure was reproduced and localized to a narrow trust-exact stopping condition. The production fitter now performs deterministic exact-Hessian Newton refinement if trust-exact stops above a tighter intermediate score target. A candidate is accepted only when the objective satisfies an Armijo decrease or the actual score norm improves. The final accepted score must still be at most (10^{-6}); the scientific acceptance threshold has not been relaxed.

A dedicated regression test regenerates the exact failed validation sample (seed 20261266, floor 0.25, n=1,000) and asserts the final posterior score and positive-definite covariance.

The independent raw-input auditor now separately reconstructs the degree-one Gaussian-prior logistic posterior, Laplace proposal, tensor Gauss-Hermite nodes, and normalized importance weights from validation records. The coefficient prior is explicitly persisted as independent Gaussian coefficients with mean (0,0) and standard deviations (5,2.5).

## Required next execution

The expanded-grid results above are superseded for primary Arm E evaluation. The updated standalone branch workflow should rerun all three floor cells on the same paired seeds and expanded 11-width grid. The rerun must pass:
- focused tests, including the exact floor-0.25 optimizer regression;
- independent validation-posterior reconstruction for every validation size;
- independent raw density reconstruction for every child;
- independent outer manifest, seed, lifecycle, aggregate metric, and Scott-width audits.

The floor-0.25 ledger must show 20/20 success for a complete pilot. An audit PASS with a preserved failure is still an incomplete pilot.

Even after a successful rerun, the result is an exploratory density-only sensitivity pilot, not bandwidth selection, H0 inference, calibration, PIT, coverage, or estimator superiority. Any bandwidth training/holdout protocol must be amended and frozen before generating its separate training/holdout data, with the posterior-mean Arm E loss as the primary Arm E estimand.

## Scope

Synthetic low-redshift Euclidean methodology only. No LVK/GWTC/public-data analysis, observational H0 measurement, real-survey validation, or claim of universal estimator superiority.
