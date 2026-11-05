# Posterior-level validation bootstrap: solver-integration audit

**Date:** 2026-11-05  
**Inspected branch:** `research/logistic-inclusion-v2-newton-2026-11-05`  
**Inspected HEAD:** `a60f73f115be3b98ac2d2250659f3e6ea2f0a0ab`  
**Status:** design gate; no solver-integration code changed in this commit.

## Finding

The posterior-level producer in
`src/dark_siren_robustness/posterior_validation_bootstrap.py` calls
`estimate_host_density_from_validation` from both its full-validation
reference path and each bootstrap attempt without specifying
`solver_version`. The estimator's default is `v1`. The producer therefore
uses V1 even though the lower-level density estimator and density-bootstrap
interfaces support V2.

The independent replay implementation in
`posterior_validation_bootstrap_audit.py` is also explicitly a two-parameter
linear logistic refit using bounded L-BFGS-B. It cannot currently verify V2
Newton fits, nor can it verify quadratic-model fits. Simply adding a
`solver_version` argument to the producer would create an unaudited path and
would be an unacceptable provenance claim.

## Required integration contract

A safe extension must add both `solver_version` and `model_type` as explicit,
validated arguments, retaining defaults `v1` and `linear` for backwards
compatibility. The selected pair must be passed consistently to:

1. the full-validation reference posterior (paired-bootstrap method);
2. every pseudo-population reference posterior (design-aware method);
3. every resampled-validation posterior.

The artifact must record both selections at top level and in its experiment
configuration/provenance. The auditor must reject missing, malformed, or
unsupported values rather than silently falling back to V1/linear.

## Independent replay requirements

The auditor must not call the producer's estimator or import its V2 fitter.
It must independently reconstruct:

- the exact linear or quadratic design matrix and standardized redshift basis;
- the selected solver's objective and optimization constraints;
- the fitted probabilities on the observed catalogue support;
- inverse weights, weighted moments, ESS and adaptive bandwidth;
- reflected KDE, merger weighting and normalized host density;
- selected likelihood, posterior grid and summary statistics.

For V2, replay must independently implement the damped-Newton update, projected
Armijo backtracking, fallback behavior, coefficient bounds, projected-gradient
termination criterion and boundary-limited classification. The replay should
compare fitted probabilities and downstream quantities, not require equality
of optimizer coefficient vectors where multiple parameter vectors can yield
equivalent fitted values (for example under rank deficiency).

For quadratic models, the auditor must use the same declared basis and bounds
as the selected fit. It must not silently refit a linear model or treat a
quadratic artifact as a linear one.

## Test plan before enabling V2

1. Preserve all existing V1/linear tests and exact artifact replay.
2. Add V1/linear default-equivalence regression tests.
3. Add producer routing tests proving that the requested solver/model pair is
   used for the full reference and all replicate branches.
4. Add V2/linear and V2/quadratic deterministic fixtures, including ordinary
   overlap, coefficient-boundary, rank-deficient, and separation-prone samples.
5. Add tamper tests for solver/model metadata and for replayed fitted
   probabilities/posteriors.
6. Verify that a deliberately mislabeled V1 artifact cannot pass as V2 and
   vice versa.
7. Run focused tests, then the full Python 3.10–3.12 CI matrix on the exact
   resulting HEAD.

## Scientific interpretation constraints

This bootstrap remains a conditional frequentist resampling sensitivity
analysis: the observed catalogue and GW-distance sample are fixed while the
validation records are resampled and the inclusion fit is refitted. Its
percentile summaries describe resampling variability of posterior means; they
are not H0 posterior credible intervals and do not establish Bayesian
marginalization over inclusion-model parameters. The design-aware
pseudo-population branch remains exploratory, and its scalar FPC adjustment
remains diagnostic only.

No confirmatory simulation campaign or manuscript claim should depend on
V2 posterior-bootstrap results until the producer and independent replay
auditor satisfy the contract above.
