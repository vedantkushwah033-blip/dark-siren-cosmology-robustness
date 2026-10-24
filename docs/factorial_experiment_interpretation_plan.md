# Interpretation plan for the prior-threshold factorial experiment

## Purpose

This document predefines the interpretation steps for the synthetic prior-support × assumed-detection-threshold experiment. It is intended to reduce the risk of selectively emphasizing favourable outcomes after inspecting the results.

## Design and estimands

The experiment simulates detected distance observations under a fixed toy model with true detection interval 0–120 Mpc, two candidate host redshifts (0.02, 0.04), host weights (0.7, 0.3), and Gaussian distance noise with sigma 20 Mpc. It compares three prior intervals, four assumed upper thresholds, and two inference approaches. Each simulated observation is reused across inference configurations, allowing paired comparisons.

This is a simplified toy model. The assumed threshold in the corrected likelihood is not the actual threshold used to generate the data for all configurations; that mismatch is deliberately part of the sensitivity experiment.

## Primary reporting quantities

For every true H0 × prior × assumed threshold × method cell, report:
- Bias of the posterior median.
- Mean absolute error and root-mean-square error.
- Empirical 68% and 95% credible-interval coverage.
- Wilson intervals for the coverage proportions.
- Mean 68% and 95% credible-interval widths.

For every paired contrast, report:
- Which configuration is the baseline and which is compared.
- The direction of subtraction (comparison minus baseline).
- Mean paired difference.
- Percentile-bootstrap 95% interval.
- Replicate count and bootstrap settings.

Do not treat a confidence interval containing zero as proof of equivalence. If practical equivalence matters, define a scientifically motivated margin before interpretation.

## Interpretation rules

1. Separate results by true H0 and prior support before discussing pooled patterns. A pattern that holds at one truth may not hold at another.
2. Distinguish posterior bias from interval calibration. A method can have small median bias but poor interval coverage, or vice versa.
3. Interpret observed coverage with its Monte Carlo uncertainty. With 300 replicates per truth/cell, the standard error of a binomial coverage estimate near 0.68 is approximately sqrt(0.68 × 0.32 / 300) ≈ 0.027; near 0.95 it is approximately 0.013. These are rough binomial standard errors, not a substitute for the Wilson intervals exported by the experiment.
4. Avoid treating the nominal 68% or 95% target as a hard pass/fail threshold based on a single cell.
5. Use paired differences for comparisons because configurations share the same simulated observation. Do not replace paired analysis with unpaired standard errors.
6. The paired bootstrap interval measures replicate-sampling uncertainty conditional on this simulation design. It does not account for uncertainty about whether the toy model is realistic.
7. Report numerical grid limitations, especially tail-quantile sensitivity identified in the separate grid-refinement audit.
8. Do not claim that the corrected method is universally better merely because it improves one metric or configuration. Examine bias, MAE/RMSE, coverage, interval width, and prior/threshold sensitivity together.
9. Do not describe the results as a real observational H0 measurement, a validation of LVK published analyses, or evidence that the model captures the full gravitational-wave selection function.
10. Treat all conclusions as conditional on the specified toy data-generating process.

## Reporting structure

1. Executive summary: main supported finding and one key limitation.
2. Simulation model: observation generator, host mixture, noise, detection interval, and inference methods.
3. Validation: software tests, artifact validation, grid-resolution checks.
4. Results: accuracy, coverage, interval width, and paired contrasts.
5. Sensitivity: differences across true H0, prior intervals, and assumed thresholds.
6. Limitations: simplified likelihood, selection rule, host model, noise, and catalogue assumptions.
7. Reproducibility: commit SHA, workflow run URL, input seed, replicate count, grid size, artifact names, and validator result.

## No-result rule

If the completed artifact does not support a clear advantage for one inference method, report that finding directly. A mixed or null result is scientifically useful when it is well calibrated, reproducible, and clearly scoped.
