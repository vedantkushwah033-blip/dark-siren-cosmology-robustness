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


## Updated numerical-grid evidence

The second-stage paired grid experiment has now also completed on 27 fixed-truth pilot cells (run 44417364702). All 54 child runs at 561 and 1,121 points passed their independent raw-data audits; the seven observable/truth input tables were byte-identical within each pair. The 561→1,121 comparison passed both predeclared tolerances in all 27 cells (maximum posterior-density difference 5.93×10^-5; maximum reported-quantile change 0.01072 km s^-1 Mpc^-1). The 141-point baseline failed the quantile tolerance in all 27 cells. This supports the 281-point grid used in the planned 600-cell follow-up relative to finer grids on this 27-cell subset; it is not a calibration result and does not establish global convergence.


## Completed 100-replicate readout (workflow 44411667480)

**Execution and integrity:** the Python 3.12 workflow completed successfully. The planned 600 cells reconciled exactly: 599 SUCCESS, 1 NON_EVALUABLE, 0 FAILED, and 0 NOT_RUN. All 600 child audit reports were verified; the independent orchestration auditor checked 300 nested validation blocks, 200 cross-curvature paired-validation blocks, 14,106 aggregate-manifest entries, and recomputed the summary metrics. This is an integrity PASS, not a calibration PASS.

The 100-replicate results are:

| Curvature gamma | Validation n | Success / planned | Bias (km s^-1 Mpc^-1) | RMSE | 68% coverage | 90% coverage | 95% coverage | Mean 68% width |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -1 | 100 | 99 / 100 | +2.020 | 17.455 | 54.5% | 72.7% | 79.8% | 23.389 |
| -1 | 500 | 100 / 100 | +2.255 | 17.215 | 51.0% | 74.0% | 78.0% | 22.725 |
| 0 | 100 | 100 / 100 | +6.428 | 19.930 | 51.0% | 79.0% | 86.0% | 24.376 |
| 0 | 500 | 100 / 100 | +6.768 | 19.843 | 51.0% | 77.0% | 85.0% | 23.582 |
| +1 | 100 | 100 / 100 | +10.396 | 21.587 | 56.0% | 82.0% | 88.0% | 24.517 |
| +1 | 500 | 100 / 100 | +11.114 | 22.071 | 54.0% | 72.0% | 88.0% | 23.493 |

The 68% coverage Monte Carlo standard errors are approximately 0.050 in every cell. For the correctly specified reference gamma=0, 68% coverage is 51% at both validation sizes; the mean posterior-median bias is approximately +6.4 to +6.8. Increasing validation size from 100 to 500 did not resolve this baseline discrepancy in this experiment. The gamma=+1 cells show a larger positive bias (+10.4 to +11.1), while gamma=-1 cells show smaller positive bias (+2.0 to +2.3).

Paired contrasts reinforce the directional curvature effect: gamma=+1 minus gamma=-1 median bias was +8.273 (MCSE 0.512; 99 paired replicates) at n=100 and +8.859 (MCSE 0.466; 100 paired replicates) at n=500. The paired n=500 minus n=100 median-bias changes were +0.152 (MCSE 0.229) at gamma=-1, +0.340 (MCSE 0.266) at gamma=0, and +0.718 (MCSE 0.258) at gamma=+1.

**Scientific interpretation:** the reference gamma=0 bias and undercoverage are the main warning. They show that the current selected-H0 pipeline is not calibrated in this four-detected-event fixed-truth design even when the inclusion model is correctly specified. Therefore, the curvature contrasts cannot be interpreted as an isolated estimate of misspecification damage, and no arm should be declared superior. The effect of the number of detected events, finite-parent discreteness, selection normalization, and the fixed-truth versus prior-predictive design must be separated in a targeted diagnostic.

Weight concentration remains substantial: median minimum catalogue Kish ESS ranges from about 13.3 to 60.9, with median maximum normalized weights from 0.078 to 0.198 across these cells. This is consistent with an unstable-weight regime and should be reported alongside posterior summaries.

This experiment uses exactly four detected events per replicate, a 1,200-galaxy parent, fixed H0=70, and a uniform [20,300] prior. It is not joint SBC, not the 30-event five-arm smoke, and not observational inference. The single NON_EVALUABLE cell is retained in the denominator; the 68% coverage denominator is successful cells only and must always be read with the status count.
