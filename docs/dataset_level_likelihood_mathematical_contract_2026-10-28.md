# Dataset-level likelihood contract for the diagnostic-identifiability gate

**Status:** mathematical contract / design review; not an implementation and not a claim of a validated hierarchical dark-siren likelihood.  
**Date:** 2026-10-28  
**Branch:** `audit/research-readiness-2026-10-28`

## Why this contract exists

The public API `posterior_h0` in `src/dark_siren_robustness/hardened_inference.py` accepts one observed distance and one noise scale per call, together with a set of candidate host redshifts and optional missing-host support. That is an event-level posterior calculation. A study about whether whole datasets can distinguish host-population misspecification from distance-noise misspecification needs an explicit dataset-level probability model before a simulator is appropriate.

This document fixes notation and identifies the extra assumptions required. It does not silently retrofit the event-level function into a population model.

## 1. Event-level mixture

Let event (i) have observed distance (d_i), reported uncertainty (s_i>0), catalogue host redshifts (z_{ij}), and normalized catalogue weights (w_{ij}). In the repository's low-redshift toy model,

[
mu_{ij}(H_0)=c z_{ij}/H_0.
]

With a Gaussian distance model and no event-selection conditioning, its catalogue-mixture likelihood is

[
L_i^{m cat}(H_0)
 = sum_j w_{ij},mathcal N(d_imidmu_{ij}(H_0),s_i^2).
]

If a discrete missing-host component is specified with mass (f_i), normalized missing-host weights (v_{ik}), and redshifts (z^{m miss}_{ik}), the finite-mixture likelihood is

[
L_i(H_0)=(1-f_i)L_i^{m cat}(H_0)
+f_isum_k v_{ik}mathcal N(d_imid c z^{m miss}_{ik}/H_0,s_i^2).
]

This expression assumes (f_iin[0,1]), weights sum to one within each component, and both components share the same event-level distance-noise family. It is a finite synthetic mixture, not a general model of the unobserved galaxy population. If the implementation uses a different convention for (f_i), the code and documentation must be reconciled before proceeding.

## 2. Independent-event composition without selection

Conditional on (H_0), the event observations are independent **only under an explicit assumption** that their event-level likelihoods share no additional latent variables or unmodeled dependencies. Under that assumption,

[
L(Dmid H_0)=prod_{i=1}^{N} L_i(H_0), qquad
log L(Dmid H_0)=sum_{i=1}^{N}log L_i(H_0).
]

For a common prior density (pi(H_0)), the dataset posterior is proportional to (pi(H_0)L(Dmid H_0)). The current scalar API returns a normalized posterior density, not a reusable log likelihood. Multiplying its normalized posterior outputs would multiply the prior once per event and therefore generally be wrong. A safe dataset-level implementation should expose event log likelihoods (or another explicitly prior-free quantity), apply the shared prior once, and normalize the final product once.

If host-population parameters, a hidden-host fraction, a noise scale, or other nuisance parameters are shared across events, they must be represented as shared parameters and integrated/marginalized once at the dataset level. Plugging event-specific posterior summaries into a product is not equivalent to a joint hierarchical model.

## 3. Selection conditioning is a separate probability model

For a declared detection event (S_i), the selected event likelihood is

[
L_i^{m sel}(H_0)=
rac{p(d_imid H_0,S_i=1)}{1}
=rac{p(d_imid H_0),mathbf 1[d_iinmathcal S]}{eta_i(H_0)},
quad
eta_i(H_0)=P(S_i=1mid H_0),
]

where (mathcal S) is the exact detection region and the indicator is one for observations admitted by that rule. This is the conditional-on-detection likelihood when the observed value is within the declared region. With independent selected events and a fixed observed event count, the product includes one (eta_i(H_0)^{-1}) factor per event, using the same event-specific population, weights, noise, and selection rule as the numerator.

This conditional likelihood is not interchangeable with a likelihood for the event count/rate. A rate-based population analysis may require a Poisson point-process likelihood and its own selection integral. The repository's finite distance-window selection function is a toy selection rule; it must not be presented as a real gravitational-wave detection model.

For the missing-host mixture, the numerator and (eta_i(H_0)) must integrate over the same complete mixture distribution. A selection denominator built only from catalogue hosts while the numerator includes hidden hosts is a model mismatch. Conversely, multiplying an extra catalogue-mass factor into an already complete mixture would double-count mass. Both are high-priority regression targets.

## 4. Positive-distance handling must match the generator

Two coherent toy models are possible, but they are different:

1. **Untruncated Gaussian observation:** (d_isimmathcal N(mu_i,s_i^2)), with the generator retaining non-positive draws. The likelihood must retain the full Gaussian density.
2. **Gaussian conditioned on positive observation:** the generator redraws until (d_i>0), so the density is
   [
   p(d_imid d_i>0,H_0)=
   rac{mathcal N(d_imidmu_i(H_0),s_i^2)}
   {Phi(mu_i(H_0)/s_i)},quad d_i>0.
   ]

The code's `positive_truncated` option corresponds to the second convention for an event-level Gaussian term. The generator, inference arm, and selection-efficiency calculation must all use mutually compatible conventions. Redrawing non-positive observations but fitting an untruncated Gaussian is a generator/likelihood mismatch, not a harmless implementation detail.

## 5. Required limiting-case and metamorphic tests

Before a dataset-level experiment, a separate module must pass these tests. They are requirements, not claims that the current code already passes them.

1. **Single event:** dataset log likelihood/posterior agrees with the existing event-level routine after matching prior support, grid, mixture weights, truncation convention, and selection conditioning.
2. **Two events:** the joint prior-free log likelihood equals the sum of the two event log likelihoods at every (H_0) grid point, to stated numerical tolerance.
3. **Permutation invariance:** reordering events leaves the dataset posterior unchanged.
4. **Duplication:** duplicating one event adds its log likelihood a second time; it must not duplicate the prior.
5. **Zero hidden mass:** (f=0) makes the missing-host redshifts/weights irrelevant and reproduces the catalogue-only likelihood.
6. **Zero catalogue mass limit:** if supported by the API, (f	o1) converges to the missing-host-only likelihood without a zero-mass logarithm.
7. **Mixture normalization:** each component's weights sum to one and the mixture masses sum to one; rescaling all weights within one component does not change its normalized shape.
8. **Selection consistency:** simulated detection frequencies agree with (eta_i(H_0)) within Monte Carlo error for fixed (H_0), and selected-likelihood tests use the identical detection interval.
9. **Truncation consistency:** for fixed parameters, generated and evaluated distance densities use the same positive-distance convention.
10. **Prior applied once:** a direct hand-computed two-event example catches multiplying normalized event posteriors rather than multiplying prior-free likelihoods.
11. **Grid convergence:** posterior summaries and primary diagnostics converge under a prespecified refinement rule; endpoint truths are analyzed separately from interior truths.
12. **Failure accounting:** zero/underflowing selection efficiency, invalid posterior support, and rejected simulated datasets are counted and reported, not silently discarded.

## 6. Diagnostic-identifiability study: minimum defensible estimand

Even a correct product likelihood does not make the error source identifiable. The next experiment should be limited to a prespecified synthetic model and a diagnostic vector whose statistics have mathematical definitions. It should compare a matched control, host-only mismatch, noise-only mismatch, and joint mismatch using paired datasets. Calibration and held-out evaluation seeds must be disjoint. The primary result should be a held-out estimate of the separation between host-only and noise-only diagnostic distributions, with uncertainty and a predeclared practical-equivalence criterion; coverage or RMSE curves alone are not source classifiers.

If the two misspecification mechanisms are not separable at the chosen sample size, that is a valid possible outcome. Do not tune thresholds on the held-out set, choose factor levels after looking at outcomes, or interpret a failed classifier as proof that the mechanisms are universally indistinguishable.

## 7. Go/no-go status

**Current status: NO-GO for a confirmatory dataset-level diagnostic-identifiability run.** The event-level implementation and existing synthetic experiments do not yet, by themselves, establish that a coherent shared-parameter dataset likelihood and prespecified diagnostic vector have been implemented and validated together. This document narrows the design requirements; it does not close that gate.

The appropriate next engineering step is a separate, minimal dataset-likelihood prototype plus limiting-case tests, after confirming the exact current branch head and test layout. Production inference code and existing headline artifacts should remain unchanged until this isolated prototype passes review.

## Scope limitations

This is a simplified low-redshift relation (d_L\approx cz/H_0) with finite candidate-host mixtures and toy Gaussian distance observations. It does not include a full gravitational-wave strain likelihood, real detector selection, a realistic galaxy-population model, a rate model, or external CHIMERA validation. No real-data (H_0) estimate or novelty claim follows from this design contract.
