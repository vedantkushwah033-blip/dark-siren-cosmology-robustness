# Extreme-probability and IPW source audit

**Date:** 5 November 2026  
**Branch inspected:** `research/logistic-inclusion-v2-newton-2026-11-05`  
**Inspected branch head:** `6926c5fa981f1e32e1a0f4e4c76dde6242e7c124`  
**Scope:** source-level audit of V2 coefficient constraints, probability prediction, inverse-probability weights, weighted KDE, and the estimated-inclusion pilot likelihood path.  
**Status:** mathematical/source audit; no new simulation campaign or runtime test is claimed.

## Executive finding

The current coefficient boxes make logistic fitting numerically finite, but they do **not** impose a practically useful positivity floor. Across the full permitted predictor/coefficient region, the smallest possible fitted probability is approximately:

- Linear model: `expit(-75) = 2.6786e-33`, giving an inverse weight of about `3.733e32`.
- Quadratic model: `expit(-125) = 5.166e-55`, giving an inverse weight of about `1.936e54`.

These are mathematical worst cases over the declared domain and coefficient box, not observed pilot maxima. They demonstrate that finite coefficients and finite probabilities do not imply statistically usable IPW weights.

The production validation-to-density wrapper currently uses V1 only. V2 remains a standalone, opt-in fitter and is not yet integrated into host-density construction, pre-posterior diagnostics, the pilot runner, or the bootstrap path. Therefore no downstream V2 performance claim is currently supported.

## 1. Exact predictor range implied by the coefficient boxes

The declared standardization maps the fitted redshift interval to (x\in[-1,1]).

### Linear model

[
\eta(x)=a+bx,\qquad |a|\leq25,\quad |b|\leq50.
]

For every (x\in[-1,1]),

[
|\eta(x)|\leq |a|+|b||x|\leq75.
]

The endpoints attain both extremes, so the full box permits (\eta=-75) and (\eta=75). Since (\pi=\operatorname{expit}(\eta)),

[
\pi_{\min}=\operatorname{expit}(-75)\approx2.6786\times10^{-33},
\qquad
1/\pi_{\min}\approx3.733\times10^{32}.
]

### Quadratic model

[
\eta(x)=a+bx+qx^2,
\quad |a|\leq25,\quad |b|\leq50,\quad |q|\leq50.
]

The triangle inequality gives

[
|\eta(x)|\leq25+50|x|+50x^2\leq125.
]

The extremes are attainable at (x=1) or (x=-1) with corresponding coefficient signs. Thus

[
\pi_{\min}=\operatorname{expit}(-125)\approx5.166\times10^{-55},
\qquad
1/\pi_{\min}\approx1.936\times10^{54}.
]

This is not a claim that an ordinary fit will reach the corners of the coefficient box. It is a proof that the current box constraints alone cannot certify useful overlap or bound the inverse weights at a scientifically acceptable scale.


## 1A. Benchmark-truth probability range (important qualification)

The theoretical coefficient-box extrema above should not be confused with the true inclusion probabilities used by the planned curvature screen. The inspected pilot defaults are intercept = 0.25, slope = -1.2, and curvature q in {-2, 0, +2}; the quadratic term is centered under the parent intensity.

Recomputing that center from the repository's exact unnormalized intensity,

`lambda(z) = 1 + 2.2 exp(-0.5 ((z - 0.105)/0.025)^2) + 1.6 exp(-0.5 ((z - 0.275)/0.040)^2)`,

on the same 20,001-point grid gives `c = 0.28450056`. Over the declared domain z in [0.01, 0.40], the true probability ranges are:

| True curvature q | Minimum true inclusion probability | Maximum true inclusion probability | Maximum true inverse probability |
|---:|---:|---:|---:|
| -2 | 0.08463 | 0.73086 | 11.82 |
| 0 | 0.27888 | 0.81000 | 3.59 |
| +2 | 0.37778 | 0.94690 | 2.65 |

Thus the three planned true DGPs do **not** deliberately create near-zero inclusion probabilities: their continuous-domain minimum is approximately 0.0846. The coefficient-box extrema are much more extreme than the generating law. If an estimated fit produces probabilities far below this range, the issue would arise from finite-validation estimation, model misspecification, constrained fitting, or support effects—not from the intended true DGP having vanishing inclusion probabilities.

This is a design-specific calculation, not a universal positivity guarantee. It strengthens the case for recording true-versus-fitted probability ranges in synthetic experiments, while keeping truth-only quantities out of pre-posterior diagnostics.

## 2. Source trace: fitted probabilities to host density

The inspected path is:

1. `estimate_host_density_from_validation` in `estimated_inclusion.py` fits a V1 linear or quadratic model.
2. It predicts probabilities at observed catalogue redshifts.
3. It computes `inverse_weights = 1.0 / probability`.
4. It reports Kish ESS, maximum inverse weight, and upper weight quantiles.
5. It passes the probabilities to `catalogue_density.estimate_host_density`.
6. The density estimator computes inverse weights again, scales all weights by their maximum, computes weighted moments and ESS, constructs a weighted reflected Gaussian KDE, multiplies by the merger-rate weight once, and normalizes the final density to unit integral.
7. The pilot passes that density to `selection_aware_toy_likelihood`, which uses the same supplied density in each event integral and in (\beta(H_0)).

The wrapper docstring explicitly says no inverse-weight clipping is applied. This is appropriate for an unclipped baseline, but the project must not describe the coefficient box as an IPW stability safeguard.

## 3. What weight scaling does and does not do

The KDE implementation uses

[
\tilde w_i=\frac{w_i}{\max_j w_j}
]

before weighted moments, ESS and normalized KDE weights are computed. This common rescaling is numerically useful: it avoids overflow from a large common weight scale and leaves normalized weighted averages unchanged.

It does **not** reduce relative concentration. If one observation carries most of the original weight, it still carries most of the scaled and normalized weight. It does not repair positivity, reduce the variance induced by extreme relative weights, or correct inclusion-model misspecification.

The pre-posterior diagnostics also compute raw weight sums and squared sums before deriving normalized weights. Under the present theoretical maxima, individual squared weights are on the order of (10^{65}) (linear) or (10^{108}) (quadratic), which are representable in IEEE float64. This does not make them statistically benign. Larger future coefficient ranges or other transformations would require renewed overflow analysis.

## 4. Existing diagnostics: strengths and limits

The current diagnostics include:

- Fitted inclusion-probability minimum, selected quantiles and maximum.
- Kish ESS and ESS fraction.
- Entropy-based ESS fraction.
- Maximum normalized weight share.
- Top 1% and top 5% weight mass.
- Maximum inverse weight and upper-tail inverse-weight quantiles.
- Validation min/max support and catalogue fraction outside that range.
- Catalogue proximity to validation boundaries.
- Coefficient-boundary status.

These are useful complementary summaries. None is a universal pass/fail threshold without calibration under the declared simulation design.

The support-gap feature is a min/max-range warning proxy, not a local density or formal covariate-shift diagnostic. It can report no observations outside the range even when validation support is extremely sparse in a scientifically important region.

The diagnostics API currently accepts V1 model identifiers only. V2 output identifiers are not accepted. This is consistent with the fact that V2 is not integrated into the downstream estimation path.

## 5. Host-density and selection-normalization audit

In the inspected estimated-inclusion pilot:

- Catalogue thinning is performed using the synthetic inclusion law.
- The estimated-IPW arm reconstructs a weighted KDE from retained catalogue redshifts.
- The merger-rate weight is applied once after the KDE is constructed.
- The resulting normalized density is passed directly to the selection-aware toy likelihood.
- The same density is used in the event-distance integrals and in the detection efficiency (\beta(H_0)).
- No separate hidden-host mixture fraction appears in this inspected pilot path.

Therefore, **no IPW-plus-hidden-host-mixture double correction is present in this specific inspected path**. This is not a global statement about every likelihood or branch in the repository.

The likelihood conditions on a fixed number of detected events under the declared Gaussian-distance, hard-threshold toy selection. It is not an event-rate likelihood and is not a realistic detector selection model.

## 6. V2 integration status

V2 is currently a separately versioned fitter. The inspected `estimate_host_density_from_validation` wrapper selects only V1 fitters. The pilot imports V1 prediction and density functions. The pre-posterior diagnostic schema accepts V1 model identifiers. The bootstrap functions refit through the V1 wrapper.

Consequently:

- V2 unit tests do not establish downstream V2 host-density performance.
- V1 pilot results cannot be attributed to V2.
- A valid V2 comparison requires an explicit solver selection, preserved model-version metadata, matching diagnostic support, bootstrap handling, and independent artifact-audit updates.
- Such integration should be done only after fitter-level parity and adversarial validation pass.

## 7. Required next work

### Gate A — probability-domain policy

Before confirmatory experiments, specify the target probability domain and a scientifically justified positivity policy. Candidate approaches include:

- Keep the uncapped estimator as the primary reference and report all instability.
- Define a prespecified minimum fitted probability or maximum inverse-weight sensitivity arm.
- Compare prespecified truncation/stabilization alternatives as distinct estimators, not silent modifications.
- Report how often the rule activates and how it changes density error and H0 calibration.

A probability floor must not be chosen after inspecting H0 results. It changes the estimator and may introduce bias; it is not a harmless numerical patch.

### Gate B — V2 parity and adversarial testing

Extend identical-data comparisons to boundary-limited and ill-conditioned cases. Compare objective values, predictions, projected KKT residuals, Hessian rank/conditioning and termination reasons. Do not require identical coefficient vectors where rank deficiency makes coefficients non-identifiable.

### Gate C — downstream integration

Only after Gate B, integrate V2 behind an explicit solver option. Update schemas and auditors to preserve the solver version and failure state. Add end-to-end tests proving that the selected fitter is actually used.

### Gate D — confirmatory design

Freeze estimands, failure denominators, prior/domain, validation sampling design, estimator arms and Monte Carlo precision before a confirmatory campaign.

## 8. Current disposition

| Item | Finding |
|---|---|
| Finite logistic coefficients | Enforced by coefficient boxes |
| Practical positivity floor | **Not enforced** |
| IPW clipping | **Not applied** in the inspected estimated-inclusion wrapper |
| KDE common rescaling | Applied; numerical scaling only |
| Relative weight concentration | Diagnosed, not controlled |
| V2 downstream integration | **Not implemented** |
| IPW/mixture double correction | Not found in the inspected pilot path |
| Same density in event likelihood and beta | Yes, in the inspected toy likelihood |
| Confirmatory simulation authorization | **NO-GO pending frozen design and validated probability/weight policy** |

## Reproducibility and evidence boundary

This document records source-level inspection of the files named above at the stated branch head. It does not report a new test run, a new simulation, a historical D-01 replay, or a full-repository audit. The D-01 historical root cause remains unresolved unless its original records are recovered and replayed.


## 8. Addendum: subsequent V2 pipeline integration (5 November 2026)

This section supersedes the earlier implementation-status sentence in the executive finding that described V2 as standalone and not integrated. That sentence accurately described the earlier inspected snapshot, but is stale relative to the current research branch.

The current branch now has an opt-in `solver_version="v1"|"v2"` path through validation-to-host-density estimation, paired validation bootstrap, pre-posterior diagnostics, the estimated-inclusion pilot runner, and the independent pilot artifact auditor. V1 remains the default. The pilot CLI exposes `--solver-version v1|v2`. This is a software integration fact, not evidence that V2 improves H0 inference.

The workflow's paired smoke step runs the V1 and V2 pilots with the same seed and compares persisted inference payloads for exact equality of the redshift grid, catalogue redshifts, validation redshifts and outcomes, H0 grid, and observed GW distances; it also compares scoring payload parent redshifts and true inclusion probabilities. The smoke uses only two replicates, 1,200 parent objects, eight detected events, validation sizes 100 and 300, and coarse 41/201 H0/redshift grids. It is an integration/provenance check only, not a performance or calibration study.

The workflow subsequently runs the V2 pilot again into the same `artifacts/estimated_inclusion_pilot_v2_smoke` output directory used by the paired step. This is redundant generation and may overwrite the paired-step V2 files; the workflow should use one canonical generation or distinct output directories. Until an exact-head run and its artifacts are inspected, no successful paired artifact or current cross-version CI result is claimed here.

The original mathematical conclusions remain unchanged:
- coefficient bounds do not imply a useful positivity floor;
- uncapped IPW remains the baseline;
- empirical diagnostics do not prove population positivity;
- no universal ESS or weight-concentration threshold has been validated;
- no downstream V2 superiority claim is supported;
- the historical D-01 root cause remains OPEN / UNKNOWN.
