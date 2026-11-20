# Weight-concentration diagnostics versus H0 posterior error: 200-replicate screen

**Date:** 2026-11-20  
**Source:** two independent 100-replicate five-arm runs under H0 prior [30,110], master seeds 2026111929 and 2026111940.  
**Design:** N_parent=1,200; 30 detected events; fixed H0=70; gamma in {-1,0,+1}; validation sizes 100/300/500; five arms.  
**Outcome:** absolute posterior-median error.  
**Method:** within-cell Spearman rank association with percentile bootstrap intervals (3,000 resamples), pooling the two disjoint seed streams.

## Audit and denominator

The paired-prior workflow [44602409251](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44602409251) completed successfully. The analysis included:
- 200 replicate seeds per curvature condition across the two streams.
- 9,000 posterior-summary rows and 9,000 weight-diagnostic rows, all from child-audited five-arm artifacts.
- 117 within-cell diagnostic/outcome associations.
- 3,000 bootstrap resamples per reported association.
- A/B/C use only n_validation=500 to avoid treating repeated, validation-invariant rows as independent; D/E retain all three validation sizes.

The association analysis is exploratory and descriptive. It is not a trained prediction rule, causal estimate, calibration test, or proof that any diagnostic can forecast failure on new simulations.

## Main estimated-IPW results

The table reports selected correlations at n_validation=500. ESS fractions are expected to move opposite to concentration metrics when weights become more unstable. Confidence intervals are percentile bootstrap intervals.

| Curvature | Arm | Diagnostic | Spearman rho | 95% bootstrap interval |
|---:|---|---|---:|---:|
| -1 | D MAP | Kish ESS fraction | +0.438 | [+0.320, +0.546] |
| -1 | D MAP | Top-5% weight mass | -0.430 | [-0.534, -0.316] |
| -1 | E integrated | Minimum node Kish ESS | +0.418 | [+0.294, +0.533] |
| -1 | E integrated | Maximum node top-5% mass | -0.419 | [-0.523, -0.302] |
| 0 | D MAP | Kish ESS fraction | -0.094 | [-0.228, +0.053] |
| 0 | D MAP | Top-5% weight mass | +0.102 | [-0.042, +0.238] |
| 0 | E integrated | Minimum node Kish ESS | -0.096 | [-0.239, +0.052] |
| 0 | E integrated | Maximum node top-5% mass | +0.110 | [-0.026, +0.253] |
| +1 | D MAP | Kish ESS fraction | -0.215 | [-0.340, -0.079] |
| +1 | D MAP | Top-5% weight mass | +0.228 | [+0.097, +0.352] |
| +1 | E integrated | Minimum node Kish ESS | -0.183 | [-0.316, -0.036] |
| +1 | E integrated | Maximum node top-5% mass | +0.219 | [+0.091, +0.342] |

The positive-curvature associations are stronger at smaller validation sizes:
- gamma=+1, n=100: D Kish ESS rho -0.370 [-0.487,-0.239], D top-5% mass +0.378 [+0.246,+0.495]; E minimum-node ESS -0.265 [-0.395,-0.123], E maximum-node top-5% mass +0.271 [+0.138,+0.399].
- gamma=+1, n=300: D Kish ESS -0.337 [-0.454,-0.209], D top-5% mass +0.356 [+0.229,+0.467]; E minimum-node ESS -0.262 [-0.388,-0.126], E maximum-node top-5% mass +0.301 [+0.175,+0.418].

At gamma=-1, the direction is reversed: for D and E, larger ESS is associated with larger absolute median error, while higher top-5% concentration is associated with smaller absolute error. At gamma=0, the estimated-arm intervals generally include zero. This sign reversal is important: a single global ESS or weight-tail threshold cannot be assumed to predict H0 error across curvature regimes.

## Interpretation

1. **Positive curvature:** greater weight concentration (higher top-5% mass, lower ESS) is associated with larger absolute H0 median error for D/E in this synthetic sample. This is consistent with weight instability contributing to error under this misspecified degree-one inclusion fit.
2. **Negative curvature:** the association reverses. Weight concentration is not uniformly harmful in this fixed design; it may interact with the direction of model misspecification and the sign of the posterior bias.
3. **Zero curvature:** no clear monotonic association appears for the estimated arms in this screen.
4. **Integrated E is not immune:** its node-level ESS and concentration diagnostics show similar curvature dependence. Integrating coefficient uncertainty does not eliminate the omitted-quadratic-model problem.
5. **No causal or predictive conclusion:** correlations are computed on the same 200 replicates used to estimate them. The bootstrap intervals quantify sampling variability within each cell but do not correct for the 117 simultaneous exploratory comparisons.

## Limitations and next test

- The outcome is absolute posterior-median error at one fixed H0 truth, not coverage or a proper scoring rule.
- The two streams are independent seed schedules, but the first five replicates in seed-2026111940 overlap an earlier engineering smoke; they are included once in the 100-replicate stream and are not added again.
- A/B/C are evaluated at n=500 only because their diagnostic rows repeat across validation-size labels. D/E are evaluated at all n.
- The analysis is a within-sample association screen, not out-of-sample prediction.

The next defensible test is to freeze diagnostic thresholds on a training set of replicate seeds and evaluate error/coverage discrimination on a fresh holdout seed stream, stratified by curvature. Report sensitivity/specificity or calibration of predicted failure risk, with uncertainty and multiplicity handled prospectively. Until then, ESS and top-tail metrics are warning indicators, not validated predictors.
