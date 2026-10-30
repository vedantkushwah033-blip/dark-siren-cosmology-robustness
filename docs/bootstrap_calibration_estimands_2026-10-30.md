# Bootstrap calibration estimands for repeated SRSWOR validation

**Status:** analysis contract only. This document defines how to summarize a future persisted experiment; it reports no simulation results.

## Experimental hierarchy

The top-level independent unit is a synthetic finite-parent scenario. Within each scenario and validation size, repeated outer validation samples are drawn by SRSWOR from the same fixed parent frame. Each outer sample is then passed through one or more inner resampling procedures. Inner replicates and outer draws are not additional independent parent scenarios.

## Primary uncertainty-calibration comparison

For scenario (s), validation size (n), and method (m):

1. Let (hat\theta_{sni}) be the posterior mean from successful outer SRSWOR draw (i).
2. Estimate the empirical repeated-SRSWOR standard deviation as the sample SD of the successful (hat\theta_{sni}) values.
3. Let (\widehat{SE}_{snim}) be the inner-bootstrap standard error reported for outer draw (i). Average this over successful outer draws.
4. Form the scenario-level calibration ratio

   [
   R_{snm}=\frac{\operatorname{mean}_i(\widehat{SE}_{snim})}
   {\operatorname{SD}_i(\hat\theta_{sni})}.
   ]

A ratio near one indicates agreement of average reported standard error with empirical repeated-SRSWOR variability for that scenario and validation size. It does not prove nominal interval coverage.

## Interval calibration

For each successful outer draw, evaluate whether its method-specific inner percentile interval contains the fixed finite-parent target posterior mean (T_{s,N}). Compute coverage within each scenario first, then average scenario-level coverage with equal scenario weights. Do not pool all outer draws across scenarios as if they were independent scenarios.

The target is the posterior mean obtained by fitting the declared pipeline to all (N) actual paired validation-frame records. It is not the generating (H_0), a pseudo-population target, or a conditional H0 posterior credible interval.

## Failure and missingness accounting

All attempted outer draws remain in the denominator for outer failure fractions. Inner-bootstrap failures must remain in their own attempted-resample denominators when their failure rates are summarized. Calibration ratios and interval coverage use only successful estimates/intervals and must be reported alongside these failure fractions. A scenario with fewer than two successful outer estimates has no estimable empirical SD ratio and must remain visible as non-estimable rather than silently dropped.

## Implemented summary contract

`dark_siren_robustness.bootstrap_calibration.summarize_bootstrap_calibration` implements scenario-level and scenario-balanced summaries. It requires one row per scenario × validation size × outer draw × method, validates duplicate grain, fixed finite-parent targets, status values, successful metric finiteness, nonnegative standard errors, and ordered intervals.

This is a post-processing component only. The persisted nested outer-SRSWOR/inner-bootstrap runner, independent artifact replay for the combined experiment, pilot, and confirmatory study remain to be implemented. No calibration claim or empirical result is made by the presence of this summary function.
