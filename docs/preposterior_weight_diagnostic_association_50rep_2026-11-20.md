# Pre-posterior weight diagnostics versus posterior error: 50-truth pilot

**Date:** 2026-11-20  
**Scope:** Exploratory synthetic benchmark diagnostic; not a validated failure predictor.

## Question

Do pre-posterior weight-concentration diagnostics (Kish ESS, maximum normalized weight, top-5% weight mass, and Arm E quadrature-node ESS) reliably predict posterior calibration or H0 error?

This is relevant to the project's proposed diagnostic question: can warning signs visible before posterior construction identify when estimated inclusion correction will fail?

## Frozen data

- Source: 50-replicate prior-predictive experiment, seed `2026112022`, artifact ID `14613021797`.
- Exact refined-grid posterior ranks: 50-replicate exact-grid audit artifact ID `14629769931`, 2,250 curve rows.
- Selected sample size for this screen: validation size 500.
- Each curvature stratum contains 50 independent prior-predictive truths. Within a replicate, GW event truth and observed GW event files are identical across the three curvature settings; this was verified by file hashes.
- Posterior median errors below are taken from the archived replicate summary (561-point original H0 grid). The independent exact-grid audit found maximum 561→1,121 median change only 0.0031 km/s/Mpc, but the original 561-grid density-L1 criterion failed. Thus the median-error association screen is not a replacement for a full posterior re-computation on the refined grid.

## Main result: no universal monotonic predictor

Across the 150 arm/curvature observations per method, pooled Spearman associations between the weight diagnostics and absolute posterior-median error were weak (absolute rho ≤ 0.143). However, the 150 rows reuse the same 50 truths across curvature, so pooled p-values are not treated as inferential.

Within each curvature stratum (50 independent truths), associations were direction-dependent:

| Arm | Diagnostic | gamma −1 rho (p) | gamma 0 rho (p) | gamma +1 rho (p) |
|---|---|---:|---:|---:|
| C known-π IPW | Kish ESS vs absolute median error | +0.273 (0.055) | +0.120 (0.406) | −0.102 (0.482) |
| C known-π IPW | Top-5% mass vs absolute median error | −0.272 (0.056) | −0.135 (0.349) | +0.076 (0.598) |
| D MAP IPW | Kish ESS vs absolute median error | +0.452 (0.001) | +0.124 (0.392) | −0.203 (0.157) |
| D MAP IPW | Top-5% mass vs absolute median error | −0.444 (0.001) | −0.123 (0.396) | +0.302 (0.033) |
| E integrated IPW | Minimum node Kish ESS vs absolute median error | +0.471 (0.001) | +0.067 (0.646) | −0.154 (0.286) |
| E integrated IPW | Maximum node top-5% mass vs absolute median error | −0.486 (<0.001) | −0.112 (0.440) | +0.287 (0.043) |

The p-values are unadjusted exploratory screens across multiple diagnostics. They should not be read as validated significance claims.

## Interpretation

1. Weight concentration alone is not a stable monotonic predictor of posterior error across curvature settings.
2. For D and E, some diagnostics correlate with absolute error within gamma −1, become weak near gamma 0, and reverse sign at gamma +1.
3. The dominant risk is not merely a low ESS or large maximum weight. It is the interaction between weight instability and inclusion-model misspecification.
4. This suggests pre-posterior checks should be reported as a diagnostic panel (ESS, max weight, top-tail mass, support gaps, fitted inclusion curvature, and validation uncertainty), not compressed into a universal single threshold.
5. A stronger predictive-diagnostic study would freeze thresholds on a training set of independent truths and test them on held-out truths, stratifying by curvature and validation size.

## Limitations

- Only 50 independent truths per curvature.
- Multiple unadjusted correlation screens; no confirmatory claim.
- Validation sizes are nested and reuse each replicate's truth.
- The median-error response uses archived 561-grid posterior summaries; the exact-grid audit showed small median shifts but a density-L1 failure for the legacy grid.
- No observational H0 or real galaxy-catalogue inference is involved.
