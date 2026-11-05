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
