# Independent posterior quadrature and invariance audit

**Date:** 28 November 2026  
**Branch:** `research/asa-model-audit-20261128`  
**Status:** independent reference calculations completed locally; CI execution of the newly added regression tests is pending.  
**Scope:** deterministic numerical validation of the positive-truncated discrete host-mixture posterior. This is not a new astrophysical model or a full grid-refinement audit of the archived 278,400 fits.

## 1. Why this audit was added

The 28 November source audit identified known-answer posterior checks and a full grid-refinement audit as the next numerical gates. The existing tests checked the positive-truncated likelihood against a direct formula, but did not yet cover the complete requested invariance set or compare posterior summaries against an independent continuous-integration implementation across several stress cases.

## 2. Independent method

A separate reference calculation was implemented using:
- direct component log densities `norm.logpdf(d; c z/H0, sigma) - norm.logcdf(c z/(H0 sigma))`;
- log-sum-exp mixture evaluation;
- adaptive Gauss–Kronrod integration through `scipy.integrate.quad` for posterior normalization and posterior mean;
- root finding through `scipy.optimize.brentq` for the 2.5%, 16%, 50%, 84%, and 97.5% quantiles.

The reference calculation does not call `posterior_h0` or `posterior_summary`. The production-style calculation uses the 4,001-point H0 grid and trapezoidal integration. The prior was fixed at [50, 90] km s^-1 Mpc^-1.

## 3. Numerical comparison

| Case | Data configuration | Maximum absolute discrepancy across mean and five quantiles (H0 units) |
|---|---|---:|
| Single host | z=0.030, d=130 Mpc, sigma=18 Mpc | 0.0000039 |
| Two listed hosts | z=(0.020, 0.040), weights=(0.7, 0.3), d=160 Mpc, sigma=20 Mpc | 0.0000039 |
| Hidden mixture | listed weights=(0.7, 0.3), hidden fraction=0.30 at z=0.030, d=180 Mpc, sigma=18 Mpc | 0.0000035 |
| Lower-prior-boundary stress | hidden fraction=0.45 at z=0.035, d=230 Mpc, sigma=9 Mpc | 0.0000181 |
| Upper-prior-boundary stress | hidden fraction=0.30 at z=0.025, d=75 Mpc, sigma=8 Mpc | 0.0000032 |

Across these five cases, the largest absolute discrepancy between the 4,001-point grid summaries and the independent adaptive-quadrature summaries was approximately **1.81 × 10^-5 H0 units**. The largest difference between 4,001- and 8,001-point grid summaries was approximately **1.57 × 10^-5 H0 units**. These are deterministic spot checks, not evidence that every archived replicate has this error bound.

## 4. Known-answer and invariance gates

The following regression tests have been added in `tests/test_independent_posterior_quadrature.py`:

1. One-host posterior compared with an independent continuous quadrature reference.
2. Zero hidden mass reduces exactly to the catalogue-only posterior.
3. Splitting a host into duplicate redshift entries with weights summing to the original weight leaves the posterior unchanged.
4. Permuting listed and hidden mixture components together with their weights leaves the posterior unchanged.
5. A one-host model represented as catalogue mass plus hidden mass at the same redshift is identical to the single-component model.
6. Four posterior-summary cases are compared with independent adaptive quadrature, including a posterior concentrated near the lower prior boundary.

The source-level mathematical identity behind the invariance checks is simply that a mixture density depends on the total weight assigned to each distinct component density, not on how that weight is partitioned across duplicate labels.

## 5. Interpretation

**Numerical disposition: encouraging, but not yet CI-confirmed.** The independent calculations agree closely with the 4,001-point grid on the tested configurations. The tests also make the principal mixture reductions explicit and guard against common bookkeeping errors.

The audit does not establish:
- grid stability for all 278,400 archived posterior fits;
- calibration or frequentist coverage;
- correctness of the scientific data-generating assumptions;
- realistic gravitational-wave detection selection;
- uncertainty propagation for estimated hidden-population parameters;
- external validity for LVK/GWTC observations.

The current hidden-host factorials still use fixed, known noise within each cell and fixed assumed hidden-population parameters in each inference arm. The oracle arm remains a simulation-only benchmark.

## 6. Next gate

Run the repository CI on this branch and inspect the exact-head workflow result. If the new tests pass, proceed to a prespecified full grid-refinement audit on a frozen confirmatory subset, recording changes in posterior means, all interval endpoints, and coverage classifications. Replication expansion should follow only after that numerical gate and after primary estimands and MCSE targets are frozen.
