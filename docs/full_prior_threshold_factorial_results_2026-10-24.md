# Full prior-threshold factorial: first artifact audit (2026-10-24)

## Provenance and scope

- GitHub Actions run: [#4](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40098731627)
- Commit: `4839e3d6a4ec23b70f182bdb760a7c1dc0dfd405`
- Run conclusion: success; duration reported by GitHub: 44 seconds.
- Artifact: `full-prior-threshold-factorial-r300-g4001-seed20261028`
- Artifact SHA-256 reported by GitHub: `5f33a36b5dd38f918b6d9958697c818a83ddd685a4a505b3034a62f7e1a29ee0`
- Reproduced design manifest: 21,600 raw rows, 72 summary rows, 396 paired-contrast rows; 300 replicates per truth; 4,001 H0 grid points; seed 20261028.
- Scope: synthetic low-redshift Gaussian distance-threshold model only. This is not an LVK observation or a realistic detector-selection pipeline.

## Integrity checks performed

The downloaded artifact contained `per_run.csv`, `summary.csv`, `paired_contrasts.csv`, and `experiment_manifest.json`. The row counts and manifest matched the prespecified factorial design. Summary metrics were recomputed from the raw rows during the workflow's independent validator; the job and validator steps completed successfully. This establishes artifact consistency under the implemented checks, not scientific realism.

## Selected reference-cell results

For prior [50,90] and assumed upper threshold 120 Mpc:

| True H0 | Method | Bias | MAE | RMSE | 68% coverage | 95% coverage |
|---:|---|---:|---:|---:|---:|---:|
| 52 | Corrected | +10.617 | 10.617 | 11.493 | 6.7% | 94.0% |
| 52 | Uncorrected | +12.720 | 12.720 | 13.358 | 0.0% | 85.7% |
| 70 | Corrected | +0.466 | 5.656 | 6.625 | 93.3% | 100.0% |
| 70 | Uncorrected | +1.926 | 5.333 | 6.325 | 90.3% | 100.0% |
| 88 | Corrected | -14.433 | 14.433 | 16.211 | 3.3% | 84.0% |
| 88 | Uncorrected | -13.268 | 13.268 | 14.881 | 5.3% | 97.7% |

Coverage is the empirical fraction of intervals containing the known simulated truth. Nominal 68% and 95% intervals should achieve approximately those coverages under repeated sampling if calibrated for the stated experiment. The values above are based on 300 replicates per cell and have Monte Carlo uncertainty; do not interpret small differences as definitive.

## Matched corrected-minus-uncorrected comparison

For the same reference cell [50,90] / 120 Mpc, paired differences across the same 300 observations were computed. The intervals below are percentile-bootstrap 95% intervals, resampling matched replicate indices (3,000 bootstrap draws; seed 9032026). These are exploratory contrasts, not multiplicity-adjusted confirmatory tests.

| True H0 | Metric | Corrected − uncorrected | Bootstrap 95% interval |
|---:|---|---:|---:|
| 52 | Posterior median | -2.103 | [-2.139, -2.064] |
| 52 | Absolute error | -2.103 | [-2.138, -2.065] |
| 52 | 68% interval width | -1.149 | [-1.315, -0.993] |
| 52 | 95% interval width | -0.080 | [-0.130, -0.030] |
| 52 | 68% coverage difference | +6.7 percentage points | [+4.0, +9.7] pp |
| 52 | 95% coverage difference | +8.3 percentage points | [+5.3, +11.7] pp |
| 70 | Posterior median | -1.460 | [-1.526, -1.394] |
| 70 | Absolute error | +0.323 | [+0.156, +0.492] |
| 70 | 68% interval width | +0.294 | [+0.182, +0.403] |
| 70 | 95% interval width | +0.654 | [+0.586, +0.723] |
| 70 | 68% coverage difference | +3.0 percentage points | [+1.3, +5.0] pp |
| 70 | 95% coverage difference | 0.0 percentage points | [0.0, 0.0] pp |
| 88 | Posterior median | -1.165 | [-1.238, -1.093] |
| 88 | Absolute error | +1.165 | [+1.095, +1.241] |
| 88 | 68% interval width | +0.312 | [+0.188, +0.427] |
| 88 | 95% interval width | +0.828 | [+0.760, +0.895] |
| 88 | 68% coverage difference | -2.0 percentage points | [-3.7, -0.7] pp |
| 88 | 95% coverage difference | -13.7 percentage points | [-17.7, -10.0] pp |

The correction is not uniformly better: it improves absolute error and coverage at truth 52 in this reference cell, has a mixed result at truth 70, and worsens error and coverage at truth 88. The truths 52 and 88 are close to the [50,90] prior boundaries, so bounded-prior effects are a major interpretation concern.

## Across the 12 prior-threshold configurations

A descriptive average of the 12 cell-level corrected-minus-uncorrected contrasts per truth (not a pooled estimator and not adjusted for multiple comparisons) showed:

| True H0 | Mean change in absolute error | Mean change in 68% coverage | Mean change in 95% coverage |
|---:|---:|---:|---:|
| 52 | -2.636 | +11.0 pp | +2.9 pp |
| 70 | +0.121 | +3.0 pp | +0.1 pp |
| 88 | +1.332 | -3.8 pp | -3.7 pp |

This averaging is descriptive only: it gives each of the 12 design configurations equal weight and should not be interpreted as a single general performance guarantee.

## Interpretation and next scientific checks

1. The experiment demonstrates strong truth/prior-boundary dependence. At truth 52, the corrected reference-cell 68% coverage is only 6.7%; at truth 88 it is 3.3%. This is a serious calibration warning, not a result to conceal.
2. The selection correction is implemented as the mixture likelihood divided by β(H0) = P(detected | H0), then renormalized on prior support. The detection efficiency uses the same Gaussian observation model and finite interval as the generator. This is internally consistent for the stated toy event-selection model, subject to independent review of the full model assumptions.
3. The experiment reuses observations across methods/configurations, so method contrasts must remain paired. The supplied `paired_contrasts.csv` primarily compares each configuration against a reference configuration within each method; it does not itself provide every direct corrected-versus-uncorrected contrast. The reference-cell paired contrasts above were independently calculated from `per_run.csv`.
4. Before publication, perform a denser fixed-truth coverage map, separate prior-boundary effects from threshold misspecification, check grid convergence for posterior interval endpoints, and add simulation-based calibration/prior-predictive controls appropriate to the question.
5. Do not claim the corrected method is generally superior. These are synthetic stress tests, with only two host redshifts, fixed mixture weights, Gaussian distance noise, and a hard distance threshold; they are not calibrated to a real detector or galaxy catalogue.

## CI annotations

GitHub reported a Node.js 20 deprecation warning for third-party workflow actions and a notice about the future `ubuntu-latest` image migration. The run nevertheless completed successfully. These are maintenance warnings, not evidence of a failed experiment.
