# Held-out validation of IPW weight diagnostics — 2026-11-20

## Design and provenance

This is a single fixed train/holdout split of the synthetic 30-event, five-arm quadratic-logit experiment, using the expanded H0 prior [30,110]. The training stream is master seed 2026111929 (workflow 44484212311); the disjoint holdout stream is master seed 2026111940 (workflow 44483923518). Both source artifacts were downloaded and verified by their embedded aggregate manifests in workflow 44672498194.

The analysis was run on PR #127 head `c2aa9f2ac77cc1cc68225b0b649df76c685b5b50`. Workflow 44672498194 completed SUCCESS and uploaded artifact 14599589819 (SHA-256 `ff4ec498ffe03724fe6e1ce0bad2e69d3761981b5b5fa2fa43409726b99f5f5a`).

Within each gamma × validation-size × arm stratum, training data alone selected the sign of the Spearman association between the pre-posterior diagnostic and absolute posterior-median error. It then froze both the oriented diagnostic 75th percentile alert threshold and the training absolute-error 75th percentile label threshold. These thresholds were applied unchanged to the disjoint holdout stream. D uses Kish ESS fraction; E uses minimum node Kish ESS.

## Results

There were 18 planned strata: 15 evaluable and 3 explicitly NON_EVALUABLE because the training association was too weak (absolute Spearman rho below 0.10).

| Curvature | Validation n | Arm | Training rho | Frozen risk direction | Holdout AUROC | Sensitivity | Specificity | Balanced accuracy |
|---:|---:|---|---:|---|---:|---:|---:|---:|
| -1 | 100 | D MAP IPW | +0.398 | Higher diagnostic | 0.732 | 0.357 | 0.826 | 0.591 |
| -1 | 100 | E integrated IPW | +0.328 | Higher diagnostic | 0.680 | 0.318 | 0.808 | 0.563 |
| -1 | 300 | D MAP IPW | +0.450 | Higher diagnostic | 0.742 | 0.474 | 0.839 | 0.657 |
| -1 | 300 | E integrated IPW | +0.381 | Higher diagnostic | 0.747 | 0.385 | 0.892 | 0.638 |
| -1 | 500 | D MAP IPW | +0.412 | Higher diagnostic | 0.811 | 0.481 | 0.863 | 0.672 |
| -1 | 500 | E integrated IPW | +0.364 | Higher diagnostic | 0.777 | 0.500 | 0.889 | 0.694 |
| 0 | 100 | D MAP IPW | +0.053 | NON_EVALUABLE | — | — | — | — |
| 0 | 100 | E integrated IPW | -0.099 | NON_EVALUABLE | — | — | — | — |
| 0 | 300 | D MAP IPW | -0.102 | Lower diagnostic | 0.559 | 0.417 | 0.671 | 0.544 |
| 0 | 300 | E integrated IPW | -0.202 | Lower diagnostic | 0.533 | 0.400 | 0.680 | 0.540 |
| 0 | 500 | D MAP IPW | -0.146 | Lower diagnostic | 0.498 | 0.320 | 0.707 | 0.513 |
| 0 | 500 | E integrated IPW | -0.156 | Lower diagnostic | 0.553 | 0.440 | 0.707 | 0.573 |
| +1 | 100 | D MAP IPW | -0.346 | Lower diagnostic | 0.657 | 0.406 | 0.868 | 0.637 |
| +1 | 100 | E integrated IPW | -0.261 | Lower diagnostic | 0.578 | 0.308 | 0.803 | 0.556 |
| +1 | 300 | D MAP IPW | -0.298 | Lower diagnostic | 0.689 | 0.391 | 0.805 | 0.598 |
| +1 | 300 | E integrated IPW | -0.192 | Lower diagnostic | 0.701 | 0.333 | 0.816 | 0.575 |
| +1 | 500 | D MAP IPW | -0.121 | Lower diagnostic | 0.641 | 0.296 | 0.808 | 0.552 |
| +1 | 500 | E integrated IPW | -0.097 | NON_EVALUABLE | — | — | — | — |

Across the 15 evaluable strata, the unweighted mean holdout AUROC is 0.660 (range 0.498–0.811). This pooled descriptive mean is not a meta-analytic estimate: strata are different arms/curvatures/validation sizes and should not be treated as exchangeable observations.

## Interpretation

1. The direction of the diagnostic-risk relationship is not stable across curvature. At gamma=-1, higher ESS-related diagnostic values were associated with larger absolute median errors; at gamma=0 and +1, the learned direction generally reversed.
2. Holdout discrimination is modest and inconsistent. AUROC ranges from essentially chance (0.498) to 0.811; several strata are only weakly discriminative.
3. Sensitivity is generally low (roughly 0.30–0.50), so the frozen top-quartile alert does not reliably capture the training-defined high-error cases. Specificity is higher, but does not rescue the low sensitivity.
4. Three strata could not learn a defensible direction. They remain NON_EVALUABLE; they are not silently dropped from the planned denominator.
5. These data do not support a universal ESS threshold or a pre-posterior diagnostic that reliably predicts posterior calibration failure.

## Limitations

- One train/holdout split and one fixed H0 truth; only 100 replicates per seed stream.
- The top-quartile absolute-error label is an exploratory operational definition, not a scientifically established failure threshold.
- Stratum-wise thresholding produces 18 separate small-sample tests; no multiplicity-adjusted confirmatory inference is claimed.
- AUROC and threshold metrics have substantial finite-sample uncertainty; Wilson intervals are reported for sensitivity and specificity in the machine-readable artifact.
- This is synthetic methodology validation only. It does not establish observational H0 performance, real galaxy-catalogue validity, or full inclusion-correction SBC calibration.

## Decision

**Exploratory negative result.** The tested ESS summaries are not sufficiently stable across curvature to serve as universal pre-posterior warning rules. Retain them as descriptive weight-stability diagnostics, not as calibrated posterior-failure predictors. A future confirmatory study would need a pre-registered error definition, larger independent seed streams, repeated outer splits, and uncertainty intervals for discrimination metrics.
