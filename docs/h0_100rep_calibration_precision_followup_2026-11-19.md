# 100-replicate fixed-truth H0 calibration precision follow-up

**Status:** completed fixed-truth precision pilot; workflow 44411667480 succeeded at runner SHA `652021943b31490302d586b5b9c77c7ccdfcc8f7`. The 600-cell artifact passed the independent orchestration audit. Oracle-arm decomposition is being run separately in PR #97.

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


## Fixed-truth results (600-cell run)

The 600-cell run completed with 599 SUCCESS and one NON_EVALUABLE cell; there were no FAILED or NOT_RUN cells. The non-evaluable cell was retained in the ledger and excluded only from fitted-arm posterior summaries. The independent orchestration audit passed all 600 child audits, 300 within-curvature nested-validation blocks, 200 paired-validation blocks, and 14,106 aggregate manifest entries; all summary metrics were recomputed from the ledger.

| Curvature γ | Validation n | Successful / planned | 68% coverage | 90% coverage | 95% coverage | Mean median bias | RMSE | Mean 68% width | Median min Kish ESS | Median max normalized weight |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| -1 | 100 | 99 / 100 | 0.545 | 0.727 | 0.798 | +2.020 | 17.455 | 23.389 | 19.70 | 0.1552 |
| -1 | 500 | 100 / 100 | 0.510 | 0.740 | 0.780 | +2.255 | 17.215 | 22.725 | 60.89 | 0.0776 |
| 0 | 100 | 100 / 100 | 0.510 | 0.790 | 0.860 | +6.428 | 19.930 | 24.376 | 13.91 | 0.1976 |
| 0 | 500 | 100 / 100 | 0.510 | 0.770 | 0.850 | +6.768 | 19.843 | 23.582 | 42.41 | 0.0954 |
| +1 | 100 | 100 / 100 | 0.560 | 0.820 | 0.880 | +10.396 | 21.587 | 24.517 | 13.30 | 0.1623 |
| +1 | 500 | 100 / 100 | 0.540 | 0.720 | 0.880 | +11.114 | 22.071 | 23.493 | 34.53 | 0.0912 |

Bias and RMSE are in km s\(^{-1}\) Mpc\(^{-1}\). Coverage Monte Carlo standard errors are approximately 0.05 in every cell, so differences of a few percentage points are not precisely resolved. The 68% coverage values (0.51–0.56) are materially below nominal in this fixed-truth pilot; 90% coverage is 0.72–0.82 and 95% coverage is 0.78–0.88. These are descriptive repeated-sampling results, not joint SBC.

### Paired contrasts

- Mean posterior-median bias difference, γ=+1 minus γ=−1: +8.273 (MCSE 0.512; 99 paired replicates) at validation n=100, and +8.859 (MCSE 0.466; 100 paired replicates) at n=500.
- Mean bias difference, validation n=500 minus n=100: +0.152 (MCSE 0.229) for γ=−1, +0.340 (MCSE 0.266) for γ=0, and +0.718 (MCSE 0.258) for γ=+1.

### Interpretation

Larger validation samples clearly improve weight stability: median minimum catalogue Kish ESS rises in every curvature condition, while the median maximum normalized weight falls. However, this did not produce a clear reduction in posterior-median bias or an improvement in coverage. In particular, the γ=0 reference condition still has mean bias around +6.4 to +6.8 and 68% coverage of 0.51. Therefore, the undercoverage and baseline bias cannot be attributed solely to the omitted quadratic curvature in the fitted inclusion model.

The γ=+1 versus γ=−1 paired contrast is large relative to its Monte Carlo SE and supports a curvature-associated difference in this synthetic setting. It does not establish a universal causal or estimator-superiority claim. The full-parent and known-π oracle decomposition on these exact 100 event realizations is the next required analysis before assigning the baseline bias to catalogue estimation, finite-catalogue noise, or the selected-event likelihood.

The single NON_EVALUABLE cell is γ=−1, validation n=100; its fitted-arm result is excluded from fitted-arm means and paired fitted-arm contrasts. Oracle posteriors can still be reconstructed from its archived parent/catalogue/event inputs and should retain the full 100-replicate oracle denominator.

Artifact: workflow 44411667480, artifact 14485703115, ZIP SHA-256 `0d6141e8efda49aacba713f1be77548ae61002680461f8585c8d13cb96aad864`. The full raw artifact contains 14,708 files and is retained for 90 days.
