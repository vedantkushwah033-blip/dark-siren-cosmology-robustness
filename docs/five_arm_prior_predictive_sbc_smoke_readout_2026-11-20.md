# Five-arm prior-predictive SBC smoke readout — 2026-11-20

## Execution and independent audit

Workflow 44674807151 completed successfully. The downloaded artifact was inspected independently.

- 225/225 planned lifecycle cells SUCCESS; 0 failures/non-successful cells.
- 15/15 replicate-by-curvature child runs passed independent audits.
- 293 outer manifest entries verified by the independent aggregate auditor.
- 45 aggregate metric rows and 360 paired contrast rows independently reconstructed.
- 15 cross-curvature pairing blocks verified.
- Five distinct H0 truths were drawn from the declared Uniform(30,110) prior and shared across curvature and validation-size conditions within each replicate.
- The prior-predictive summary contained 45 rows and 450 rank-histogram rows.

The audit establishes artifact integrity, lifecycle completeness, pairing and arithmetic. It does not establish calibrated inference from five draws.

## Descriptive fixed smoke results

Selected aggregate metrics (five replicates per cell):

| gamma | validation n | Arm | Median bias | RMSE | 68% coverage |
|---:|---:|---|---:|---:|---:|
| -1 | 100 | Catalogue-only A | -16.70 | 18.00 | 0/5 |
| -1 | 100 | Known-π IPW C | +0.07 | 5.95 | 3/5 |
| -1 | 100 | Integrated IPW E | -3.03 | 5.98 | 2/5 |
| 0 | 500 | Catalogue-only A | -15.97 | 17.43 | 0/5 |
| 0 | 500 | Known-π IPW C | -0.03 | 5.35 | 2/5 |
| 0 | 500 | Integrated IPW E | +0.22 | 4.92 | 2/5 |
| +1 | 500 | Catalogue-only A | -15.01 | 16.45 | 0/5 |
| +1 | 500 | Known-π IPW C | +0.37 | 5.18 | 3/5 |
| +1 | 500 | Integrated IPW E | +2.73 | 5.88 | 1/5 |

Across the remaining cells, E does not uniformly beat C: its RMSE advantage depends on curvature and validation size. These five-replicate frequencies are highly unstable. No rank-uniformity test or coverage comparison should be treated as confirmatory at this sample size.

## Interpretation

1. The prior-predictive generator, five-arm replicate pipeline, independent child reconstruction, lifecycle ledger and summary workflow all completed.
2. The catalogue-only arm is strongly displaced in this tiny draw set, while the corrected arms are closer to the prior-drawn truths. This is descriptive only.
3. Integrated estimated-probability IPW is not uniformly better than known-probability IPW, even in this small smoke.
4. The experiment does not isolate the coefficient-prior parameterization issue between degree-one and degree-two models.
5. The next calibration gate requires hundreds or preferably thousands of prior-predictive replicates, with explicit SBC rank diagnostics, prior sensitivity, weight-instability strata, and full independent reconstruction.

## Claim limit

Synthetic low-redshift Euclidean benchmark only. Not an observational H0 result, not LVK/GWTC data analysis, not a real galaxy-survey validation, and not a universal estimator ranking.
