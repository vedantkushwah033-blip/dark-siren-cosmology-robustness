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

The latest exact-head rerun, workflow `44441591070` (master seed `2026111907`), completed successfully. Its independent raw-seed audit regenerated all 200 truth draws, 200 parent populations, 200 detected-event sequences, and all 200 posterior-rank/quantile rows. The audit reported `PASS`, verified all 3 source manifest entries, and independently recomputed the summary metrics. Artifact: `full-parent-h0-sbc`, ID `14492958677`; uploaded ZIP SHA-256: `4e6665738269e54e4a31b3e3cbf8e190ea99c59033e65cd9b100299fe93f8cd7`.

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

**Grid qualification:** the source run used 561 H0 grid points. The separate 561→1,121 comparison (PR #109, workflow `44444768331`) independently reconstructed all 200 posteriors and verified the four source manifest entries, but failed the preregistered worst-case numerical tolerances: max quantile shift 0.1633 km s^-1 Mpc^-1 (limit 0.05), max absolute rank shift 0.0696 (limit 0.002), and max posterior-density L1 difference 0.0812 (limit 0.001). Median differences were much smaller (quantile 0.00599, rank 0.0000273, L1 0.0000988), indicating a small tail of poorly resolved posteriors rather than a uniform shift. Replicate 100 was the worst case. Therefore the 200-replicate SBC summary is not yet demonstrated grid-stable for every replicate. Do not relax tolerances; inspect the worst posterior shapes and refine with denser/adaptive grids. The 1,000-replicate extension (PR #111) is pending and uses 1,121 points, but that grid should not be assumed converged without its own resolution check.
