# Dataset-level likelihood contract for the diagnostic-identifiability gate

**Status:** mathematical contract plus isolated prototype; not a validated hierarchical dark-siren likelihood.  
**Date:** 2026-10-28  
**Branch:** `audit/research-readiness-2026-10-28`

## Why this contract exists

The public API `posterior_h0` in `src/dark_siren_robustness/hardened_inference.py` accepts one observed distance and one noise scale per call, with candidate host redshifts and optional missing-host support. It is an event-level posterior calculation. A study asking whether whole datasets can distinguish host-population misspecification from distance-noise misspecification needs an explicit dataset-level probability model before simulation.

This document fixes notation and identifies assumptions. It does not silently retrofit the event-level function into a population model.

## 1. Event-level mixture

For event i, let observed distance be d_i, reported uncertainty s_i > 0, catalogue host redshifts z_ij, and normalized catalogue weights w_ij. In the repository's low-redshift toy model:

`mu_ij(H0) = c * z_ij / H0`

Without selection conditioning, the catalogue-mixture likelihood is:

`Lcat_i(H0) = sum_j w_ij * Normal(d_i | mu_ij(H0), s_i^2)`

If a discrete missing-host component has mass f_i, normalized weights v_ik, and redshifts zmiss_ik, the finite-mixture likelihood is:

`L_i(H0) = (1 - f_i) * Lcat_i(H0) + f_i * sum_k v_ik * Normal(d_i | c*zmiss_ik/H0, s_i^2)`

This assumes f_i in [0,1], weights sum to one within each component, and both components share the same event-level noise family. It is a finite synthetic mixture, not a general model of the unobserved galaxy population. If the implementation uses a different convention for f_i, code and documentation must be reconciled first.

## 2. Independent-event composition without selection

Conditional on H0, event observations are independent only if no additional shared latent variables or unmodeled dependencies remain. Under that explicit assumption:

`L(D | H0) = product_i L_i(H0)`

and:

`log L(D | H0) = sum_i log L_i(H0)`

For a common prior density pi(H0), the dataset posterior is proportional to pi(H0) * L(D | H0). The current scalar API returns a normalized posterior density, not a reusable prior-free likelihood. Multiplying normalized event posterior outputs would generally multiply the prior once per event and is therefore wrong. A safe dataset-level implementation should expose event log likelihoods (or another explicitly prior-free quantity), apply the shared prior once, and normalize the final product once.

If host-population parameters, hidden-host fraction, noise scale, or other nuisance parameters are shared across events, they must be represented as shared parameters and marginalized once at dataset level. Multiplying event-specific posterior summaries is not equivalent to a joint hierarchical model.

## 3. Selection conditioning is a separate probability model

For a declared detection region S, define:

`beta_i(H0) = P(d_i is detected | H0)`

For an observed value admitted by S, the conditional-on-detection likelihood is proportional to:

`Lsel_i(H0) = p(d_i | H0) / beta_i(H0)`

With independent selected events and a fixed observed event count, the dataset product includes one beta_i(H0)^(-1) factor per event, using the same population mixture, weights, noise convention, and detection rule as the numerator.

This conditional likelihood is not interchangeable with a likelihood for the event count or rate. A rate-based population analysis may require a Poisson point-process likelihood and its own selection integral. The repository's finite distance-window selection function is a toy rule, not a real gravitational-wave detection model.

For a missing-host mixture, the numerator and beta_i(H0) must integrate over the same complete mixture distribution. A selection denominator built only from catalogue hosts while the numerator includes hidden hosts is a model mismatch. Conversely, multiplying an extra catalogue-mass factor into an already complete mixture would double-count mass. Both deserve regression tests.

## 4. Positive-distance handling must match the generator

Two coherent toy models are possible:

1. **Untruncated Gaussian:** d_i ~ Normal(mu_i(H0), s_i^2), retaining non-positive draws. Fit the full Gaussian density.
2. **Gaussian conditioned on positive observation:** if the generator redraws until d_i > 0, then for d_i > 0 use:

   `p(d_i | d_i > 0, H0) = Normal(d_i | mu_i(H0), s_i^2) / Phi(mu_i(H0)/s_i)`

The code's `positive_truncated` option corresponds to the second convention for an event-level Gaussian term. Generator, inference arm, and selection-efficiency calculation must use compatible conventions. Redrawing non-positive observations but fitting an untruncated Gaussian is a generator/likelihood mismatch.

## 5. Required limiting-case and metamorphic tests

Before a dataset-level experiment, a separate module must pass these tests. These are requirements, not claims that the current code already passes them.

1. **Single event:** dataset posterior agrees with the event-level routine under matched prior, grid, mixture, truncation, and selection assumptions.
2. **Two events:** joint prior-free log likelihood equals the sum of event log likelihoods at every H0 grid point.
3. **Permutation invariance:** reordering events leaves the dataset posterior unchanged.
4. **Duplication:** duplicating one event adds its log likelihood a second time but does not duplicate the prior.
5. **Zero hidden mass:** f=0 makes hidden-host support irrelevant and reproduces the catalogue-only likelihood.
6. **Zero catalogue-mass limit:** if supported, f approaching 1 converges to the missing-host-only likelihood without zero-mass logarithm errors.
7. **Mixture normalization:** component weights and mixture masses sum to one; rescaling all weights within a component does not change its normalized shape.
8. **Selection consistency:** simulated detection frequencies agree with beta_i(H0) within Monte Carlo error, and inference uses the identical detection interval.
9. **Truncation consistency:** generated and evaluated distance densities use the same positive-distance convention.
10. **Prior applied once:** a hand-computed two-event example catches multiplying normalized event posteriors instead of prior-free likelihoods.
11. **Grid convergence:** summaries and diagnostics converge under a prespecified refinement rule; endpoint truths are analyzed separately from interior truths.
12. **Failure accounting:** zero/underflowing selection efficiency, invalid posterior support, and rejected datasets are counted and reported, not silently discarded.

## 6. Diagnostic-identifiability study: minimum defensible estimand

Even a correct product likelihood does not make the source of misspecification identifiable. Limit the next experiment to a prespecified synthetic model and a diagnostic vector whose statistics are mathematically defined. Compare a matched control, host-only mismatch, noise-only mismatch, and joint mismatch using paired datasets. Calibration and held-out evaluation seeds must be disjoint.

The primary result should be a held-out estimate of separation between host-only and noise-only diagnostic distributions, with uncertainty and a predeclared practical-equivalence criterion. Coverage or RMSE curves alone are not source classifiers.

If the mechanisms are not separable at the chosen sample size, that is a valid possible outcome. Do not tune thresholds on held-out data, choose factor levels after inspecting outcomes, or interpret a failed classifier as proof that the mechanisms are universally indistinguishable.

## 7. Go/no-go status

**Current status: NO-GO for a confirmatory dataset-level diagnostic-identifiability run.** Existing event-level inference and synthetic experiments do not by themselves establish that a coherent shared-parameter dataset likelihood and prespecified diagnostic vector have been implemented and validated together. This document narrows the design requirements; it does not close that gate.

The first isolated prototype now lives in `src/dark_siren_robustness/dataset_likelihood.py`, with 13 limiting-case and metamorphic tests in `tests/test_dataset_likelihood.py`. It composes prior-free event log likelihoods and applies the uniform H0 prior once. Tests cover one-event equivalence (including missing-host mixtures and positive-distance truncation), log-likelihood additivity, event-order invariance, event duplication, zero hidden-host mass, prior-support normalization, host-weight scale invariance, missing-weight boundary rejection, and direct three-event normalization. The workflow previously reported success on commit `97c7588e800659d7be9998d40647d0df008088ad`, but inspection showed its YAML did not invoke this new test file. Commit `8587d3d85ac97db9e17c8d73bbd1b526d8f18970` adds the dedicated pytest step; validation of that updated workflow is pending. This is only a prototype: selection conditioning, event-rate inference, and shared nuisance parameters remain out of scope. Production inference code and existing headline artifacts remain unchanged.

## Scope limitations

This is a simplified low-redshift relation, d_L approximately c*z/H0, with finite candidate-host mixtures and toy Gaussian distance observations. It does not include a full gravitational-wave strain likelihood, real detector selection, a realistic galaxy-population model, a rate model, or external CHIMERA validation. No real-data H0 estimate or novelty claim follows from this design contract.
