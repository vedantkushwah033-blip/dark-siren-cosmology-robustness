# 100-replicate fixed-truth H0 calibration precision follow-up

**Status:** preregistered synthetic follow-up; the 600-cell GitHub Actions computation is running in workflow 44411667480. No calibration results are available yet.

## Motivation

The preceding 90-cell pilot used 10 replicates per curvature-by-validation-size cell. It exposed a substantial positive median bias even at the correctly specified reference curvature, but its coverage estimates have very large Monte Carlo uncertainty. The next experiment increases the number of paired replicates to 100 per cell while retaining the same low-redshift Euclidean toy model and fixed-truth target. It is a precision follow-up, not a change in the scientific scope.

## Fixed design

- 100 replicate seeds, crossed with catalogue-inclusion curvature \(\gamma\in\{-1,0,+1\}\) and validation sample size \(n_{\rm val}\in\{100,500\}\): 600 planned cells.
- Fixed generating \(H_0=70\;\mathrm{km\,s^{-1}\,Mpc^{-1}}\), 1,200 parent galaxies, and exactly four detected events per replicate.
- Validation frame size 500. The 100-unit sample is a prefix of the 500-unit simple random sample without replacement.
- The analysis seed initializes the parent-population, temporary catalogue-thinning, and GW-event generation sequence; this combined analysis sequence is reused across curvature conditions within replicate. The validation frame uses a separate `validation_seed`. This runner does not currently expose a third, independent event seed. Catalogue retention outcomes vary with curvature as intended.
- The DGP uses the quadratic-logit inclusion model; inference fits a degree-one logistic inclusion model. Thus \(\gamma=0\) is the correctly specified reference within the fitted family, while \(\gamma=\pm1\) are misspecified conditions.
- The fitted inclusion-coefficient posterior uses the declared independent Gaussian prior, adaptive Laplace-proposal Gauss-Hermite quadrature at orders 3, 5, 7, and 9, and a shared coefficient node across all four events.
- \(H_0\) prior support is uniform on [20, 300], represented on a 281-point grid (spacing 1 km s\(^{-1}\) Mpc\(^{-1}\)). A paired grid-refinement audit has now completed on 27 fixed-truth pilot cells: 281→561 passed both predeclared tolerances in 27/27 comparisons, while 141→281 passed only 1/27. Thus 281 points is supported relative to 561 on this subset; this remains a numerical sensitivity audit, not a calibration result.
- The runner writes a complete planned-cell ledger before executing children. Every cell remains visible as SUCCESS, NON_EVALUABLE, FAILED, or NOT_RUN; failed and non-evaluable cells are not silently dropped.

## Endpoints and analysis

### Primary descriptive endpoints

For each \((\gamma,n_{\rm val})\) cell, report:

1. Mean posterior-median bias relative to 70.
2. RMSE of posterior medians.
3. Empirical 68%, 90%, and 95% equal-tailed interval coverage.
4. Mean 68% interval width.
5. Median minimum catalogue Kish ESS and median maximum normalized catalogue weight.

Coverage is accompanied by a binomial Monte Carlo standard error and a 95% Wilson interval. These quantify simulation precision; they do not correct model misspecification or establish nominal frequentist coverage.

### Paired contrasts

- At each validation size, compare the posterior-median bias under \(\gamma=+1\) with \(\gamma=-1\), paired by replicate.
- Within each curvature condition, compare \(n_{\rm val}=500\) with \(n_{\rm val}=100\), paired by replicate.
- Report paired means, standard deviations, and Monte Carlo standard errors. Do not interpret the paired contrasts as universal monotonic effects.

### Integrity gates

Before interpreting numerical results:

- Run the focused auditor regression tests.
- Independently verify every child manifest and child audit report.
- Reconcile all 600 planned cell IDs and all status counts against the replicate ledger.
- Verify nested validation-ID prefixes, shared parent/event/frame identities, and cross-curvature validation ID/redshift pairing.
- Recompute all aggregate summary metrics from successful ledger rows.
- Preserve the complete artifact, including failed/non-evaluable attempts, precision summaries, and the independent orchestration audit.

## Interpretation limits

This is **fixed-truth repeated-sampling coverage**, not joint simulation-based calibration. The true inclusion coefficients are fixed by the DGP rather than jointly drawn from the fitted coefficient prior. It conditions on exactly four detected events and therefore does not validate event-rate or detection-count inference. The model remains a synthetic low-redshift Euclidean benchmark; it is not LVK/GWTC analysis, a real galaxy-survey validation, or an observational \(H_0\) measurement.

A 100-replicate cell still has non-negligible Monte Carlo uncertainty (for example, the standard error of a coverage proportion near 0.68 is about 0.047). Treat the results as a materially improved precision diagnostic, not definitive calibration evidence. Any apparent baseline bias at \(\gamma=0\) must be investigated separately rather than attributed automatically to curvature misspecification.

## Reproducibility

The branch-scoped workflow is `.github/workflows/h0-100rep-calibration-precision.yml`. It uses seed 2026111927, uploads the raw 600-cell artifact, and runs the independent orchestration auditor. No production likelihood code or `main` branch is changed by this follow-up.
