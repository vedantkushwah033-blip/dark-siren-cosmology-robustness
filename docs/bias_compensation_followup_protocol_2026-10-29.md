# Prespecified follow-up: truth-dependent bias compensation

**Date:** 29 October 2026  
**Status:** follow-up protocol and interpretation guide; the four-cell results are already available and are reported separately. This is not a preregistration and the follow-up must not be described as prospectively registered before data inspection.

## Motivation

In the selection-function simulations, a misspecified shallow transition sometimes produced lower H0 RMSE than weighting with the true inclusion curve. That is not evidence that misspecification is preferable: RMSE combines squared bias and variance, so an incorrect model can accidentally compensate for bias under one truth setting. This follow-up checks whether the apparent advantage persists when the true H0 changes.

## Design

Use the existing four follow-up cells: high inclusion floor and broad transition, each at true H0 = 60 and 80 km s^-1 Mpc^-1. Each cell contains 250 paired synthetic replicates and compares the correct-weighting and shallow-transition arms. Preserve the paired replicate structure. Report per-arm bias, variance, squared bias, RMSE, and a paired bootstrap interval for

`Delta_RMSE = RMSE(correct weighting) - RMSE(shallow transition)`.

Positive values mean the shallow-transition arm has lower RMSE in that cell; negative values mean correct weighting has lower RMSE. The interval is a descriptive uncertainty summary, not a binary test of a universal mechanism.

## Interpretation rules

- An interval entirely above zero is evidence that the shallow arm has lower RMSE in that particular cell.
- An interval spanning zero is **inconclusive**; it does not affirm either model or support the bias-compensation hypothesis by itself.
- An interval entirely below zero is evidence that correct weighting has lower RMSE in that cell.
- Do not declare a hypothesis supported merely because an interval includes zero.
- Compare the bias-squared and variance components to describe the decomposition. Because RMSE is the square root of bias-squared plus variance (up to the chosen finite-sample variance convention), describe the components as explanatory diagnostics rather than causal mediation.
- A pattern that changes with true H0 is consistent with truth-dependent bias compensation, but four cells from a simplified simulation cannot establish that this is the only mechanism or generalize to real gravitational-wave data.
- The four cells are synthetic and conditional on the specified toy generator. Do not pool them as if they were independent real-world replications, and do not infer real-data H0 performance.

## Decision rule

Report the four cell estimates and intervals individually. The follow-up's conclusion is mixed/truth-dependent if the direction or magnitude of the RMSE contrast changes across H0 settings. Any interval that crosses zero is labelled inconclusive. No global “winner” is declared unless a separate estimand and aggregation rule are specified before a new analysis.

## Reproducibility checklist

1. Record repository commit, input artifact paths, file hashes, Python and dependency versions, and analysis seed.
2. Validate one row per replicate-arm, six arms per replicate where applicable, consistent seed and truth within replicate, finite metrics, and normalized posterior grids.
3. Recompute arm-level bias, RMSE, variance and squared bias from raw replicate records.
4. Bootstrap paired replicate contrasts, not independently resampled arms.
5. Preserve the original results; any corrections or alternate bootstrap seeds must be identified as reanalysis rather than silently replacing prior outputs.
