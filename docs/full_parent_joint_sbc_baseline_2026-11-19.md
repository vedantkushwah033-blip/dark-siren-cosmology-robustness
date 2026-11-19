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
