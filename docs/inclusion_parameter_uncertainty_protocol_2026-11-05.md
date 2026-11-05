# Uncertainty propagation for estimated catalogue-inclusion probabilities

**Date:** 5 November 2026  
**Status:** methodological protocol; implementation and empirical validation pending  
**Scope:** synthetic catalogue-inclusion correction in the dark-siren robustness study

## 1. Why this protocol is needed

The current estimated-inclusion pipeline fits logistic inclusion probabilities on an external validation sample, predicts probabilities for retained catalogue objects, forms inverse-probability weights, constructs a weighted KDE, and then evaluates an H0 likelihood. The current plug-in pipeline treats fitted inclusion-model coefficients as fixed after estimation.

That plug-in procedure does not propagate finite-validation-sample uncertainty in the inclusion coefficients into the host-density estimate or downstream H0 posterior. A narrow H0 posterior conditional on one fitted inclusion model can therefore understate uncertainty in the full two-stage procedure.

This is a statistical limitation, not a software defect that can be fixed by increasing the number of GW events or refining the H0 grid.

## 2. Data and notation

Let the complete parent frame be indexed by i=1,...,N, with redshift z_i and catalogue-retention indicator R_i. The validation design selects a subset V of size m by simple random sampling without replacement (SRSWOR) from the parent frame and observes (z_i,R_i) for i in V. The retained catalogue supplies redshifts z_i for objects with R_i=1, but does not supply the parent-frame redshifts for non-retained objects.

Let theta denote the inclusion-model coefficients. For the quadratic model,

$$
\pi_i(\theta)=P(R_i=1\mid z_i,\theta)
=\operatorname{logit}^{-1}(\alpha+\beta x_i+q x_i^2).
$$

The linear model is the nested restriction q=0. The centering convention for the true data-generating quadratic term must continue to use the parent-law moment documented in the statistical model; the fitted intercept absorbs a fixed shift in the quadratic basis.

Let C denote the observed retained catalogue, V the external validation data, and E the GW-event observations.

## 3. Target inferential construction

A principled Bayesian propagation target, conditional on the current validation sample and observed catalogue, is

$$
p(H_0\mid E,C,V)
\propto p(H_0)
\int p(E\mid H_0,C,\theta)\,
p(\theta\mid V)\,d\theta,
$$

provided the conditional likelihood p(E|H0,C,theta) is derived from a coherent generative model and does not count catalogue retention or GW detection selection twice.

This equation defines the intended integration order. The inclusion-parameter posterior is learned from validation data; uncertainty is then propagated through the nonlinear map theta → inverse weights → KDE → event likelihood. The event likelihood must be averaged over theta before normalizing over H0.

**Do not** average already-normalized H0 posteriors with equal weights and call the result Bayesian marginalization. Posterior draws may have different marginal likelihood/evidence contributions. Likewise, averaging log-likelihoods is not the required integral.

## 4. Candidate implementation strategy

### A. Bayesian inclusion-parameter integration (primary candidate)

1. Fit a proper Bayesian logistic inclusion model to the validation records (z_i,R_i), using a declared, weakly informative prior on coefficients.
2. Generate posterior draws theta^(s) from p(theta|V). Record convergence diagnostics, effective sample sizes, divergences (if using HMC), and prior sensitivity.
3. For every draw, predict pi_i(theta^(s)) for retained catalogue objects.
4. Apply the declared untruncated inverse-probability rule for the primary estimator. Any floor, cap, or stabilization rule must be a separately named sensitivity arm.
5. Recompute the weighted KDE for each theta draw, with bandwidth handling prespecified. If bandwidth is re-estimated per draw, that is part of the propagation map and must be recorded.
6. Compute log L_s(H0)=log p(E|H0,C,theta^(s)) on the same H0 grid and under the same event-selection convention.
7. Approximate the marginal likelihood by a stable log-mean-exp:
   
   $$
   \log \bar L(H_0)=\operatorname{logsumexp}_s\{\log L_s(H_0)\}-\log S.
   $$
   
8. Normalize p(H0) exp(log Lbar(H0)) over the declared prior domain using the repository's numerical quadrature convention.
9. Report Monte Carlo error in the integration over theta separately from posterior uncertainty in H0 and from simulation-replicate uncertainty.

This procedure is valid only if p(E|H0,C,theta) is a coherent conditional likelihood. It is not automatically valid merely because the numerical integral is computed accurately.

### B. Design-respecting repeated-sampling analysis (separate frequentist study)

For each fixed synthetic parent population, retention realization, catalogue and GW-event dataset, repeatedly draw validation samples under the actual SRSWOR design. Refit the inclusion model and rerun the complete downstream pipeline for every validation draw.

This estimates the repeated-validation-sample variability of the estimator/posterior summaries under that finite-frame design. It is a design-based uncertainty assessment, not Bayesian marginalization over theta. Report the finite-population frame, sampling fraction, replacement status, number of successful refits, and all failures.

Do not substitute an IID bootstrap for SRSWOR without a separate justification. Do not assume a scalar finite-population correction repairs a nonlinear logistic → inverse-weight → KDE → H0 functional.

### C. Comparison

The study should compare:
- fixed plug-in inclusion coefficients;
- Bayesian integration over inclusion coefficients;
- repeated SRSWOR validation-sample refitting.

The latter two answer different uncertainty questions and should not be collapsed into a single estimator ranking.

## 5. Required validation before implementation is considered scientifically usable

### Known-answer and simulation checks
- Under a correctly specified inclusion model, verify recovery and calibration of inclusion probabilities across validation sizes.
- Check prior sensitivity for intercept, slope and quadratic coefficient.
- Include low-overlap and near-separation cases; preserve optimizer or sampler failures.
- Verify deterministic seed lineage and exact reuse of parent/catalogue/GW data across uncertainty methods within replicate.
- Check the S=1 Bayesian integration limit against the corresponding plug-in likelihood calculation.
- Check increasing S convergence of log-mean-exp likelihood curves and posterior summaries.
- Compare numerical quadrature on the H0 grid and KDE grid independently.
- Verify the marginal likelihood is computed from likelihoods, not an unweighted average of normalized posteriors.
- Independently reconstruct a subset of posterior calculations from persisted draws and raw inputs.

### Calibration and coverage
- Distinguish posterior credible-interval calibration under the joint Bayesian model from fixed-truth repeated-sampling coverage under misspecification.
- Evaluate fixed-truth coverage across outer simulation replicates.
- Report interval width, posterior mean error, and inclusion-parameter uncertainty separately.
- Keep conditional-on-success results beside unconditional operational failure rates.
- Use a prespecified primary contrast and multiplicity strategy before any confirmatory campaign.

### Selection-model coherence
Before this method is connected to the dark-siren likelihood, derive the complete joint model and establish:
- whether the catalogue is conditioned on its observed size or includes a count likelihood;
- whether the host-population density is pre-detection or post-detection;
- whether the GW event numerator is unconditional or already conditioned on detection;
- whether beta(H0) uses exactly the same host mixture and selection mechanism as the numerator;
- whether the catalogue inclusion probability is applied exactly once;
- whether uncertainty in the inclusion model changes the population-level selection efficiency.

If these cannot be established from a single generative model, keep the uncertainty-propagation experiment at the host-density stage and do not describe it as a complete dark-siren posterior.

## 6. Release gates

**Gate U1 — model derivation:** write and independently review the joint probability factorization before coding a Bayesian sampler.

**Gate U2 — inclusion posterior:** validate coefficient posterior computation, prior sensitivity, and posterior predictive checks on synthetic validation data.

**Gate U3 — downstream integration:** establish the likelihood contract and implement log-mean-exp integration with an S=1 regression test.

**Gate U4 — repeated-sampling comparison:** implement SRSWOR outer validation resampling with complete failure accounting.

**Gate U5 — independent audit:** independently reconstruct selected posterior curves and calibration summaries from persisted raw draws and inputs.

No confirmatory simulation or manuscript claim about propagated inclusion uncertainty should be made before all gates pass.

## 7. Explicit non-claims

This protocol does not establish that Bayesian inclusion-parameter integration is already implemented, that the current plug-in posterior is calibrated, or that the toy likelihood is a realistic LVK analysis. The validation sample in the current pilot is drawn from the same synthetic parent frame; it is not a transportability study. All results remain synthetic.
