# Full-parent prior-predictive simulation-based calibration

## Objective

Test whether the common full-parent selected-event likelihood is calibrated when the true H0 is generated from the same prior used for inference. This is a prerequisite baseline, not the five-arm catalogue-inclusion experiment.

## Preregistered design

- 200 independent replicates.
- For each replicate, draw H0 uniformly from [20, 300] km s^-1 Mpc^-1.
- Generate a fresh parent population of 1,200 galaxies and 30 detected events using the repository's synthetic low-redshift Euclidean generator.
- The generator applies the one-sided distance threshold when producing detected events.
- Inference uses the complete finite parent population, known host-rate weights, the Gaussian distance likelihood, and the same one-sided detection-efficiency normalization.
- The H0 prior in inference is uniform on [20, 300], matching the prior used to generate truth.
- Evaluate posterior ranks at 561 equally spaced H0 grid points.
- Report rank histogram, mean rank, KS and histogram diagnostics, posterior-median bias/RMSE, and 50/68/80/90/95% coverage.
- Reconstruct every truth draw, parent population, event sequence, posterior rank, and quantile independently from archived seeds.

## Interpretation

A calibrated continuous posterior under prior-predictive simulation should produce approximately uniform posterior ranks and nominal interval coverage, up to finite Monte Carlo and grid error. With 200 replicates, diagnostics are exploratory and have limited power; p-values are not binary proof of calibration. This is fixed-model full-parent SBC, not a test of estimated inclusion correction, catalogue misspecification, real galaxy surveys, LVK/GWTC data, or an observational H0 measurement.

## Reproducibility

The runner stores the complete design, replicate-level prior draws and posterior ranks, summary metrics, and SHA-256 manifest. The independent auditor regenerates the truth draws and simulation outputs from separate seed streams and reconstructs the selected-event posterior using a separately coded likelihood implementation.

## Completed 200-replicate result (2026-11-19)

Workflow `44440150116` completed successfully. Its independent raw-seed audit regenerated all 200 truth draws, 200 parent populations, 200 detected-event sequences, and all 200 posterior-rank/quantile rows. The audit reported `PASS`, verified all 3 source manifest entries, and independently recomputed the summary metrics. Artifact: `full-parent-h0-sbc`, ID `14490552961`; uploaded ZIP SHA-256: `a26fa4cc1af2c568ea0d309a606eced37fc4288bb90648f4132bf8780e1ad1fc`.

| Diagnostic | Result |
|---|---:|
| Mean posterior rank | 0.5063 |
| Rank SD | 0.2905 |
| KS statistic / p-value | 0.02756 / 0.9971 |
| 10-bin rank chi-square / p-value | 5.50 / 0.7887 |
| 50% interval coverage | 48.0% (MCSE 3.53 pp) |
| 68% interval coverage | 67.5% (MCSE 3.31 pp) |
| 80% interval coverage | 80.0% (MCSE 2.83 pp) |
| 90% interval coverage | 91.0% (MCSE 2.02 pp) |
| 95% interval coverage | 96.0% (MCSE 1.39 pp) |
| Mean posterior-median bias | -1.475 km s^-1 Mpc^-1 (MCSE 1.047) |
| Posterior-median RMSE | 14.848 km s^-1 Mpc^-1 |

The rank histogram counts across ten equal-width bins were `[18, 25, 14, 24, 16, 22, 20, 18, 21, 22]`. These results are consistent with approximately uniform ranks and nominal interval coverage at this sample size; the large KS p-value is not evidence that calibration has been proven. The median-bias MCSE is reported separately from RMSE: the latter includes replicate-to-replicate posterior uncertainty and is not an estimator-bias standard error.

**Grid qualification:** the source run used 561 H0 grid points. A separate independent comparison against 1,121 points is queued, so the numerical grid-sensitivity check is not yet complete. The 1,000-replicate extension is also pending and is not represented by these results.
