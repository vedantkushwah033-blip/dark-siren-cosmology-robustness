# Exploratory common-physical-bandwidth pilot review — 2026-11-14

## Status

This note records a diagnostic review of the 20-replicate, floor-0.10, density-only pilot artifact produced by workflow run `43586031115` (commit `4fd2d85038d4cbac6410ee06216f9279616e6679`). It is **not** a confirmatory result, bandwidth selection, or H0 inference result.

The run attempted 20 planned replicates. Nineteen completed; replicate 17 failed during the degree-one Gaussian-prior logistic inclusion-posterior mode fit. The failure was a numerical stopping condition: the trust-exact optimizer reported “A bad approximation caused failure to predict improvement” with infinity-norm score (1.1534884\times10^{-6}), only slightly above the former (10^{-6}) acceptance threshold. The workflow correctly stopped before the independent aggregate audit. Therefore the pilot status remains incomplete/failed (19/20), not PASS.

## Integrity and focused independent checks

- All 238 files listed by the outer SHA-256 manifest were present and matched their recorded byte counts and hashes.
- All 19 successful child manifests were independently checked; all 12 files in each child matched its recorded size and SHA-256 digest.
- The degree-one Gaussian-prior logistic posterior was independently reconstructed from observable validation records for all 19 successful replicates at validation sizes 100, 300, and 1,000 (57 reconstructions). The independent reference used BFGS followed by a separately coded damped-Newton polish, then rebuilt the Laplace covariance, tensor Gauss-Hermite nodes, and normalized importance weights.
- Maximum absolute differences between independent and persisted posterior objects were (5.1\times10^{-8}) for the mode, (7.0\times10^{-8}) for quadrature nodes, and (2.7\times10^{-8}) for normalized quadrature weights.
- One representative full density cell (replicate 0, validation size 1,000, physical bandwidth 0.012) was independently reconstructed using reflected weighted Gaussian KDEs. A-D and all 25 Arm E node density curves agreed with persisted curves to better than (10^{-10}) maximum absolute density difference; the posterior-mean Arm E density agreed to (1.1\times10^{-13}).

These checks establish selected artifact-integrity and numerical-equivalence facts only. They do not replace the complete per-child density auditor and outer aggregate auditor, which did not run to completion for this artifact.

## Preliminary bandwidth-grid signal

A descriptive summary was computed from the 19 successful replicates. For every primary arm (A catalogue-only, B parent oracle, C known-probability IPW, D MAP plug-in IPW, and E posterior-mean density) and each nested validation size (100, 300, 1,000), the smallest tested physical bandwidth, 0.012, had the lowest mean IAE among the seven tested widths.

For this preliminary failure-inclusive upper-bound summary, the one failed replicate was assigned IAE = 2, the upper bound for the L1 distance between two normalized densities. The 19-success-only estimates and Monte Carlo standard errors were also retained in the working analysis. Because the run is incomplete and the full artifact auditor has not passed, these numbers are a grid-design diagnostic, not evidence that 0.012 is optimal.

The endpoint pattern indicates that the current exploratory grid may not extend far enough toward smaller physical bandwidths. It is consistent with the preregistered principle that a boundary minimum must not be frozen as a selected bandwidth. Before any training/holdout study begins, the candidate grid should be reconsidered through a documented preregistration amendment, with the amendment made before examining any training or holdout results. The pilot's floor-0.10 result must not be generalized to the floor-0.00 or floor-0.25 cells.

## Numerical follow-up

The following changes were made on the PR #54 research branch:

1. The proper Gaussian-prior logistic posterior fitter now accepts a score infinity norm up to (2\times10^{-6}), rather than (10^{-6}), while retaining finite-mode and positive-definite observed-information checks. The change addresses a narrow trust-exact termination case where the objective improvement is below floating-point resolution.
2. The independent BFGS/Newton reference uses the same explicit score bound and compares independently reconstructed modes, nodes, and weights at (10^{-7}) absolute tolerance. The tolerance is still materially tighter than the observed posterior object scales and was checked against all 57 successful persisted cases.
3. A regression test exercises the validation sample associated with the failed replicate and asserts the independently computed posterior score bound.

The latest branch commit when this note was added is recorded by GitHub on PR #54. A new exact-head CI result and a fresh 20/20 pilot with all child and outer audits passing have **not** yet been observed. Do not treat the numerical fix as CI-validated until that happens.

## Decision

- Do not freeze a physical bandwidth.
- Do not launch the 600-replicate training/holdout study yet.
- Do not begin H0 propagation based on this density-only pilot.
- First obtain exact-head CI evidence for the fitter, independent auditor, and regression test; then rerun the same 20-replicate seed schedule and require 20/20 success plus complete independent child and outer audits.
- If the corrected pilot confirms the lower-boundary signal, amend and freeze a wider candidate grid before the separate training/holdout data are generated.
- Keep the floor-based logistic misspecification experiment distinct from the separate quadratic-logit curvature family. This pilot tests only the floor-0.10 DGP.

## Scientific scope

Synthetic low-redshift Euclidean methodology only. No LVK/GWTC/public-data analysis, no observational H0 measurement, no real-survey validation, no H0 calibration or coverage result, and no estimator-superiority claim.
