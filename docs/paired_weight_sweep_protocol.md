# Paired host-weight misspecification sweep: prespecified protocol

**Status:** implementation updated to report paired bootstrap uncertainty; results from the revised sweep are pending a successful execution on this revision. The older single-run summaries in `host_misspecification_audit_results.md` are exploratory and are not a substitute for this paired design.

## Question

When two candidate host redshifts are plausible, how does reversing their catalogue weights affect posterior-median error and 68% interval coverage as host-redshift separation and distance uncertainty vary?

This is a synthetic, low-redshift, known-support stress test. It is not a model of a real LVK event or a complete catalogue-incompleteness correction.

## Fixed design

- Injected value: $H_0=70$ km s$^{-1}$ Mpc$^{-1}$.
- Prior: uniform on $[50,90]$ km s$^{-1}$ Mpc$^{-1}$.
- Candidate host 1: $z_1=0.02$.
- Candidate host 2: $z_2\in\{0.025,0.04,0.06\}$.
- True host weights: $(0.7,0.3)$; misspecified inference weights: $(0.3,0.7)$.
- Fractional distance-noise settings: $\{0.03,0.10,0.20\}$, with the same observation and noise standard deviation used for both inference weightings in each replicate.
- Default: 1,000 replicates in each of 9 design cells, with base seed 20261028.
- Each observation is evaluated twice, yielding 9,000 paired observations and 18,000 posterior fits.
- The posterior grid and likelihood implementation are held fixed to isolate weighting sensitivity; grid accuracy is covered by separate audits.

## Primary estimands

Within each design cell, report the reversed-minus-correct paired difference for:
1. absolute posterior-median error;
2. 68% truth-inclusion probability;
3. RMSE of posterior-median error.

For the first two, bootstrap the replicate-level paired differences. For RMSE, resample replicate rows jointly so the correct and reversed errors remain paired. Use deterministic percentile bootstrap intervals with 4,000 resamples and cell-specific recorded seeds.

The sign convention is explicit: positive absolute-error and RMSE differences indicate worse point-estimation performance under reversed weights; negative coverage differences indicate lower empirical coverage under reversed weights. Coverage is a Bernoulli outcome, so intervals can be noisy at 1,000 replicates and must be interpreted alongside the estimate and interval, not by a threshold-only significance rule.

## Validation and reporting gates

- The same generated observation, true host, and sigma must feed both weighting models in each replicate.
- All 9 design cells and exactly 1,000 replicate records per cell must be present in the default workflow.
- All summary estimates and interval endpoints must be finite; each lower endpoint must not exceed its upper endpoint.
- The workflow must confirm the expected 4,000 bootstrap resamples per metric and cell before uploading results.
- The report must preserve the synthetic-model limitation and must not generalize from one injected $H_0$ to a population of cosmological parameters.
- A follow-up should test multiple injected $H_0$ values and explicit omitted-host scenarios; do not combine those questions into this weighting-only experiment after inspecting results.

## Interpretation boundary

The weight reversal is a controlled misspecification intervention, not a universal ranking of host-weighting methods. Its effect can be weak when likelihood components overlap heavily, and stronger when host redshifts are separated or observations are more precise. The sweep is designed to quantify that interaction rather than assert in advance that reversed weights must always increase error.
