# Paired host-weight misspecification sweep: prespecified protocol

**Status:** revised experiment now includes an explicit prior-support factor. The initial 9-cell narrow-prior sweep is exploratory only: it showed that weight reversal can appear to have no effect when the alternative host's implied $H_0$ mode falls outside the prior. The revised run tests that interaction directly; results are pending its successful workflow and artifact audit.

## Question

How does reversing catalogue weights change posterior-median error and 68% interval coverage as host-redshift separation and distance uncertainty vary, and how much of that effect depends on prior support retaining alternative-host solutions?

This is a synthetic, low-redshift, known-support stress test. It is not a real LVK event analysis or a complete catalogue-incompleteness correction.

## Fixed design

- Injected value: $H_0=70$ km s$^{-1}$ Mpc$^{-1}$.
- Host 1: $z_1=0.02$; host 2: $z_2\in\{0.025,0.04,0.06\}$.
- Data-generation weights: $(0.7,0.3)$; reversed inference weights: $(0.3,0.7)$.
- Fractional distance-noise settings: $\{0.03,0.10,0.20\}$, with identical generated observation and noise standard deviation for both weighting models and both priors in each replicate.
- Two prior-support arms: narrow $[50,90]$ and broad $[20,220]$ km s$^{-1}$ Mpc$^{-1}$.
- The grid spans $[10,250]$ at step 0.025, retaining the previous grid spacing while allowing alternate-host solutions implied by the redshift pairs to remain in support.
- Default: 1,000 replicates per redshift/noise cell, base seed 20261028. There are 9 generated-data cells, 9,000 shared observations, 18 prior-support × design summaries, and four posterior fits per observation (36,000 total).
- Each observation is generated once and then evaluated under correct and reversed weights in both prior arms.

## Why prior support is a factor

For a distance generated at one host redshift, an alternative host redshift implies a different candidate $H_0$ mode proportional to $z_{\rm assumed}/z_{\rm true}$. Under $[50,90]$, some alternatives implied by $z=0.02$ versus $z=0.04$ or $0.06$ lie outside the prior. Once that component is excluded, changing its weight can have little or no effect on the normalized posterior. A near-zero weighting contrast in such a cell is therefore not evidence that host weights are intrinsically unimportant; it can be a support-truncation artifact.

The broad-prior arm is a diagnostic to retain these alternate modes. Its wider prior is deliberately not presented as a recommended astrophysical prior. The contrast between prior arms quantifies how much the observed weight effect depends on prior support.

## Primary estimands

Within each redshift/noise/prior cell, report reversed-minus-correct paired differences for:
1. absolute posterior-median error;
2. 68% truth-inclusion probability;
3. RMSE of posterior-median error.

For absolute error and coverage, bootstrap replicate-level paired differences. For RMSE, resample each replicate's correct/reversed errors jointly. Use deterministic percentile-bootstrap intervals with 4,000 resamples and recorded cell-specific seeds.

Positive absolute-error or RMSE differences indicate worse point-estimation performance under reversed weights. Negative coverage differences indicate lower empirical coverage. Coverage estimates are Monte Carlo proportions; interpret them with the paired interval and do not use a threshold-only significance rule.

## Validation and reporting gates

- The same observation, true host, and sigma must feed both weighting models and both prior arms within a replicate.
- Exactly 9,000 raw observation rows and 18 summaries must be present for the default run.
- Every summary estimate and interval endpoint must be finite; lower interval endpoints must not exceed upper endpoints.
- The workflow must confirm 4,000 bootstrap resamples per metric and cell before artifact upload.
- The report must show narrow and broad prior arms side by side and must not imply a universal ranking of weighting methods.
- This is a fixed-$H_0$ experiment. Generalization across $H_0$ requires a separate prespecified multi-truth design.
