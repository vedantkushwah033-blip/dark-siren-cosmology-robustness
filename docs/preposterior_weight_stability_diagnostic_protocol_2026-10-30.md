# Pre-posterior weight-stability diagnostic protocol

**Date:** 30 October 2026  
**Status:** prespecified follow-on design; not executed; no predictive performance is claimed.  
**Scope:** controlled synthetic estimated-inclusion dark-siren study.

## 1. Motivation and central distinction

The existing pilot records *oracle* population diagnostics (the true inclusion-probability distribution) and *fitted* diagnostics (estimated inverse weights). Oracle quantities are useful for explaining a simulation's data-generating mechanism, but they are unavailable in an actual survey analysis. A practical warning system must therefore be evaluated using only information available before the cosmological posterior is inspected.

This protocol asks whether pre-posterior diagnostics can flag conditions associated with unstable weighting and unreliable downstream inference. It does not assume that any single ESS or maximum-weight threshold is universally valid.

## 2. Prediction time and information firewall

Define the diagnostic time point as immediately after fitting the inclusion model and constructing host weights, but before evaluating the GW likelihood or producing an H0 posterior.

Candidate predictors may include:
- validation sample size and observed inclusion count;
- fitted-model convergence status and coefficient-boundary flags;
- fitted inclusion-probability quantiles among retained catalogue objects;
- inverse-weight entropy effective sample size, maximum normalized weight share, and the fraction of total weight carried by the largest 1% and 5% of catalogue objects;
- the fraction of catalogue redshifts near the observed validation-frame support boundaries;
- inverse-weight ESS and ESS fraction;
- maximum, 95th and 99th percentile estimated inverse weights;
- bootstrap variability of fitted coefficients and predicted inclusion probabilities, if computed using only the validation sample;
- validation-frame redshift support and the fraction of retained-catalogue redshifts outside or near its supported range.

Do **not** use true inclusion probabilities, true selection coefficients, posterior means, posterior intervals, coverage, absolute H0 error, or any other downstream outcome as predictors. True probabilities may be retained separately for simulation-only explanation and oracle scoring, never as inputs to the warning rule.

The pilot now records observable-only fitted-weight and validation-support diagnostics before posterior construction. It also records oracle population diagnostics separately for simulation explanation. Oracle quantities must remain excluded from any operational prediction model. The added concentration features are descriptive predictors, not universal warning thresholds.

## 3. Outcomes

Evaluate separate outcomes rather than collapsing all problems into one opaque label:

1. **Weight instability:** prespecified severe estimated-weight tail or very low estimated ESS fraction. This is a diagnostic target, not proof of posterior failure.
2. **Estimator instability:** fit failure, non-convergence, or coefficient-boundary limitation.
3. **Inference failure:** downstream posterior construction failure.
4. **Inferential degradation:** large absolute posterior-mean error and/or interval noncoverage relative to the injected truth. These are simulation outcomes and cannot be known for a real event.

Set thresholds for any binary outcome before opening held-out results. Also report continuous outcomes (ESS fraction, maximum weight, absolute error, interval width) to avoid threshold-only conclusions.

## 4. Evaluation design

- Generate independent simulation replicates and split by replicate seed into development and held-out test sets. Never split rows from the same replicate across sets.
- Keep all arms from a replicate in the same partition.
- Tune any diagnostic score and thresholds using development replicates only.
- Freeze the score, features, thresholds, and code before evaluating held-out replicates.
- Report held-out sensitivity, specificity, false-alert rate, precision, AUROC and area under the precision-recall curve where appropriate; include uncertainty intervals and prevalence.
- Report calibration (Brier score and reliability summaries) if a probabilistic risk score is used.
- Compare against simple baselines: validation size alone, ESS fraction alone, and a prespecified rule based on weight quantiles. A complex score is useful only if it improves held-out performance and calibration.
- Repeat the evaluation across curvature and validation-size conditions. Report every cell; do not pool away severe positivity conditions or differential failures.
- Use independent seeds for the held-out evaluation. A final untouched seed set is preferred if development choices are extensive.

## 5. Required ablations and stress tests

Evaluate the warning rule under:
- correctly specified linear selection (quadratic coefficient zero);
- positive and negative nonlinear curvature;
- small and large validation samples;
- low-retention and moderate-retention regimes;
- coefficient-boundary hits and fit failures;
- cases where weights are extreme but posterior error is not large, and cases where posterior error is large without obvious weight instability.

The last two cases are essential: they measure whether the diagnostic is informative rather than merely correlated with a convenient proxy.

## 6. Reporting requirements

For every condition report:
- attempted, successful, and paired-success replicate counts;
- warning prevalence and confusion matrix;
- sensitivity, specificity, false-positive and false-negative rates;
- precision-recall and ROC summaries with uncertainty;
- calibration and baseline comparisons;
- failure reasons and boundary-hit frequencies;
- downstream error, coverage, interval width and density error, stratified by warning status;
- a clear separation between oracle truth used for explanation and observable predictors used by the warning system.

Do not describe a warning rule as predicting H0 coverage or posterior failure unless held-out results support that exact claim. A warning about extreme weights is narrower and more defensible than a claim of reliable posterior-risk prediction.

## 7. Stop/go gates

**Gate 1 — information audit:** every predictor is available at the declared pre-posterior time point; no truth or downstream posterior field leaks into the feature matrix.

**Gate 2 — reproducibility:** deterministic replicate-level partitions, frozen seeds, and a machine-readable feature/target manifest.

**Gate 3 — held-out evaluation:** the untouched test set is evaluated once after the rule is frozen; all metrics and all prespecified conditions are reported.

**Gate 4 — external relevance:** before any real-survey recommendation, document which predictors can actually be computed from that survey's validation design and catalogue. Synthetic predictive performance alone does not establish transportability.

## 8. Interpretation limits

A model can only flag risk patterns represented in its simulated design. This protocol does not establish universal ESS thresholds, guarantee coverage, solve positivity violations, or validate the estimator on LVK or a real galaxy survey. The independent parent-frame validation assumption remains strong and must be relaxed in separate studies (for example, validation nonresponse, covariate shift, redshift error, and spatially structured missingness).

**Current status:** observable-only feature extraction is implemented in `src/dark_siren_robustness/preposterior_diagnostics.py` and integrated into the pilot. Independent artifact checks are being extended. The replicate-level predictive evaluation remains unimplemented; no predictive performance or held-out result is claimed.
