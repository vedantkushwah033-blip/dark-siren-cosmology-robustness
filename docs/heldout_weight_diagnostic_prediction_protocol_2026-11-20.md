# Prospective holdout protocol: do IPW weight-tail diagnostics predict H0 error?

**Status:** frozen protocol and implementation; no holdout result is available until the third seed stream completes. This is a synthetic fixed-truth predictive-validation experiment, not an observational analysis.

## Question

Can a pre-posterior weight-concentration diagnostic predict the absolute posterior-median H0 error on new replicate seeds, beyond a training-set mean-error baseline?

The earlier 200-replicate screen found that the same diagnostic/error association changed direction with curvature. That screen was exploratory and in-sample. This follow-up reserves a third, previously unused master seed for a prospective holdout.

## Frozen data split

- Training stream 1: master seed 2026111929, 100 replicates.
- Training stream 2: master seed 2026111940, 100 replicates.
- Holdout stream: master seed 2026111950, 100 replicates.
- All streams: N_parent=1,200; 30 detected events; fixed H0=70; gamma=-1,0,+1; validation sizes 100,300,500; five arms; common H0 prior [30,110].
- The two training streams are used only to fit the diagnostic-to-error mapping. The holdout stream is not used for model selection, threshold selection, or tuning.

## Primary outcome and predictors

Outcome: absolute posterior-median error, |median(H0 posterior)-70|, in km s^-1 Mpc^-1.

Only estimated-inclusion arms are tested:
- D (degree-one MAP IPW): `top_5_percent_mass`.
- E (degree-one uncertainty-integrated IPW): `maximum_node_top_5_percent_mass`.

Each predictor is computed before evaluating the H0 posterior. One predictor is used per arm/cell; no post-hoc selection among ESS, maximum weight, or alternative tail fractions is allowed in this holdout test.

## Frozen prediction rule

For each of 18 strata (3 curvatures × 3 validation sizes × 2 estimated-IPW arms):
1. Pool the two training streams (200 independent parent/event replicates).
2. Fit ordinary least squares with an intercept: absolute median error on the selected top-5% mass diagnostic. Standardization uses training data only.
3. Compare its holdout prediction with a constant baseline equal to the training mean absolute error.
4. Evaluate on the 100 untouched holdout replicates.

Primary metric: holdout MAE difference = diagnostic-model MAE minus training-mean-baseline MAE. Negative values favor the diagnostic model. Secondary metrics: holdout RMSE difference and holdout Spearman correlation.

## Uncertainty and multiplicity

For each stratum, use 3,000 deterministic nested bootstrap resamples. Each resample independently resamples the 200 training replicates and 100 holdout replicates, refits the OLS model, and calculates the paired MAE difference. Report percentile 95% intervals and a one-sided bootstrap tail probability for improvement. Apply Holm correction across all 18 prespecified tests.

A cell meets the predeclared predictive-value criterion only if:
- holdout MAE difference is negative;
- the upper endpoint of its unadjusted 95% bootstrap interval is below zero; and
- Holm-adjusted one-sided p < 0.05.

The 18-cell family is exploratory; the criteria are deliberately stricter than an unadjusted correlation screen. Report every cell, including cells that do not meet the criterion.

## Failure and audit handling

Every source artifact must contain the full 4,500-cell planned denominator, all cells SUCCESS, all child audits PASS, and an independent orchestration audit with 300 successful child runs. No failed or non-evaluable replicate may be silently dropped. If a source artifact or audit is incomplete, the holdout comparison is NON-EVALUABLE and no substitute seed is selected after seeing results.

## Interpretation limits

- One fixed truth and one synthetic DGP family only.
- This tests prediction of absolute posterior-median error, not coverage calibration, causal effects, or transfer to observational catalogues.
- A successful result would support predictive value only for these frozen strata and simulation settings; it would not establish a universal ESS threshold.
- A null result would mean these single tail-mass predictors did not improve on the training-mean baseline in this holdout design; it would not prove all diagnostics useless.
- The bootstrap intervals include training and holdout replicate resampling, but do not cover uncertainty across alternative DGP families.
- No LVK/GWTC/public event data or real galaxy survey is used.
