# Monte Carlo precision plan for synthetic calibration and diagnostic-power studies

**Status:** planning document; no new simulations or runtime benchmarks were performed.  
**Date:** 2026-10-28  
**Scope:** fixed-truth coverage, prior-predictive coverage, diagnostic rejection/power, and paired multi-arm contrasts in the repository's synthetic models.

## 1. Separate the units of replication

Never report only the number of simulated events. Report at least:

- (N_{\rm event}): events within each synthetic dataset;
- (R_{\rm dataset}): independent simulated datasets/replicates;
- (K_{\rm fit}): posterior fits performed per dataset, including paired inference arms;
- number of independent seeds/configurations used for replication.

For coverage, the independent Monte Carlo unit is normally the simulated dataset/truth replicate, not each posterior-grid point and not each event within a dataset. For a diagnostic-power estimate, the unit is the independently generated dataset whose diagnostic outcome is recorded. If several arms analyse the same observations, preserve their pairing.

## 2. Single-method coverage or rejection-rate precision

For a Bernoulli outcome with probability (p) across (R) independent replicates, the planning approximation is

[
\operatorname{MCSE}(\hat p)\approx\sqrt{p(1-p)/R}.
]

Use this formula to plan, then report Wilson confidence intervals (or a prespecified exact/binomial interval) for the realized estimate. Do not substitute the planning approximation for the reported interval.

| Target quantity | Approximate MCSE at R=300 | At R=1,000 | At R=4,000 | R for approximate 95% half-width ≤ 0.02 | R for approximate 95% half-width ≤ 0.01 |
|---|---:|---:|---:|---:|---:|
| 68% coverage ((p=0.68)) | 2.69 pp | 1.47 pp | 0.74 pp | 2,092 | 8,366 |
| 95% coverage ((p=0.95)) | 1.26 pp | 0.69 pp | 0.34 pp | 457 | 1,825 |
| Power near 80% ((p=0.80)) | 2.31 pp | 1.26 pp | 0.63 pp | 1,537 | 6,147 |
| Type-I error near 1% ((p=0.01)) | 0.57 pp | 0.31 pp | 0.16 pp | 96 | 384 |

“pp” means percentage points. The final two columns use the approximate normal planning formula (R\approx p(1-p)(1.96/h)^2), where (h) is the desired absolute half-width; actual Wilson intervals are discrete and may differ, particularly for rare outcomes. At a 1% nominal test size, a small replicate count cannot estimate size precisely; do not infer reliable calibration from zero or one false rejection.

For a 68% coverage estimate, (R=300) has about 2.69 percentage points of MCSE, so it is a coarse diagnostic, not a high-precision estimate. If a scientific claim turns on a 1–2 percentage-point difference, hundreds of replicates are not enough.

## 3. Paired two-arm contrasts

For the same dataset (r) analysed under two prespecified arms, define

[
D_r=I_{A,r}-I_{B,r},\qquad \widehat\Delta=R^{-1}\sum_rD_r,
]

where (I_{A,r}) and (I_{B,r}) indicate whether the respective interval covers the same known truth. Estimate uncertainty from the replicate-level (D_r):

[
\operatorname{SE}(\widehat\Delta)=s_D/\sqrt R.
]

Do not calculate the paired contrast's uncertainty by treating the two coverage estimates as independent binomial proportions. Report the discordant-pair counts as a useful diagnostic and use a prespecified paired bootstrap or a suitable paired-proportion method for confidence intervals. Keep the bootstrap/resampling unit at the independent dataset level.

## 4. Four-arm interaction contrasts

For a (2\times2) host-error by distance-noise-error design, define

[
D_r=I_{11,r}-I_{10,r}-I_{01,r}+I_{00,r},
\qquad \widehat\Delta_{\rm int}=R^{-1}\sum_rD_r.
]

Compute the interaction uncertainty from the empirical distribution of this replicate-level contrast. Do not reuse the two-arm discordance formula, and do not combine four separate binomial standard errors as if the arms were independent. All four arms must analyse the same generated observations if the estimand is intended to isolate inference-model interaction. If a generator changes between arms, describe the contrast as a generator-plus-inference contrast instead.

For continuous outcomes such as posterior mean error, RMSE differences, or interval width, define the paired per-dataset contrast first and estimate uncertainty across independent datasets. If RMSE is the primary metric, resample whole datasets and recompute RMSE in each resample rather than treating event-level residuals as independent when they are not.

## 5. Diagnostic power and type-I error

The supplied nine-patch bundle reports a diagnostic-power study with 100,000 generated events divided into datasets of (N\in\{100,300,1000\}), corresponding to (R\in\{1000,333,100\}) independent datasets. The reported (N=1000) power estimates therefore use only 100 independent dataset outcomes; their uncertainty can be substantial (the bundle itself notes Wilson half-widths up to about 0.10).

For each test and condition:
- report the rejection count, denominator, estimated rejection rate, and 95% Wilson interval;
- report empirical type-I error under the matched control separately from power under each alternative;
- prespecify the test threshold and multiplicity treatment;
- do not infer high power from a small p-value in a single large pooled event sample;
- do not infer low power from a non-significant result without reporting its uncertainty interval;
- when comparing powers across conditions, preserve pairing if the same simulated datasets are analysed under multiple arms.

A power estimate near 0.5 from 100 independent datasets is very imprecise. To target a 95% interval half-width around 0.02 for power near 0.8 requires roughly 1,537 independent datasets under the normal approximation; for a 0.01 half-width, roughly 6,147. These are planning values, not a mandate to run a huge study before checking effect sizes and feasibility.

## 6. Prespecification and pilot separation

1. Use a pilot only to test code, estimate runtime, check numerical stability, and estimate plausible variance.
2. If pilot outcomes influence factor levels, diagnostics, thresholds, exclusions, or primary estimands, do not include those pilot replicates in confirmatory estimates.
3. Freeze the generator, inference arms, seed scheme, independent-replicate count, primary contrast, interval method, multiplicity handling, and stopping rule in version control before the confirmatory run.
4. Choose (R) from the desired precision for the primary estimand, not from a customary round number.
5. Report all failed fits and exclusions with prespecified reasons; do not silently drop numerical failures.
6. If multiple configurations or seeds are used, distinguish within-configuration Monte Carlo uncertainty from variation across configurations.

## 7. Current decision

The uploaded bundle already contains a diagnostic-power study and a coverage-heterogeneity decomposition, but the bundle is not the same thing as a verified production-repository run. Before commissioning another study, reconcile the nine-patch mbox against the claimed 13-commit rebased series and the live PR #3 head. Do not duplicate an existing experiment merely to increase project size.

No production runtime was measured in preparing this plan. No claim of statistical significance, calibrated power, or new scientific finding is made by this document.
