# Integration audit: inclusion-posterior quadrature and shared H0 likelihood

**Status:** source-level integration review; no multi-event calibration or coverage result is claimed.  
**Target branch:** `research/inclusion-posterior-quadrature-prototype`  
**Scope:** synthetic low-redshift Euclidean methodology only.

## 1. Review update

The first version of this note recorded the pre-integration base-branch API. The target branch has since added:

- `src/dark_siren_robustness/shared_inclusion_likelihood.py`
- `src/dark_siren_robustness/h0_posterior.py`
- corresponding tests.

The earlier statement that the target branch had no shared multi-event likelihood is therefore superseded. This revision reviews the new implementation and its remaining validation gap.

## 2. Implemented shared-node likelihood

`shared_inclusion_selected_log_likelihood` now:

1. Takes one validation-posterior quadrature object, a fixed discrete catalogue, and multiple event records.
2. Builds node-specific normalized host weights from the same inclusion-probability nodes:
   `w[j,k] ∝ base_weight[j] / pi(z[j]; theta[k])`.
3. Computes each event's raw Gaussian host-mixture log likelihood conditional on each shared coefficient node.
4. Sums event log likelihoods at fixed node k.
5. Computes a node-specific one-sided detection efficiency from those same host weights.
6. Subtracts `N_events * log(beta)` at each node.
7. Marginalizes the joint selected likelihood once over the validation-posterior node weights using log-sum-exp.

This is the correct shared-nuisance factorization for its declared common-noise, common-threshold, fixed-size toy sample. It does not independently marginalize each event. The optional `population_base_weights` are part of the population-mixture assumption and must be declared and matched to the generator.

`normalize_h0_posterior` then applies the uniform H0 prior and normalizes the already marginalized log likelihood. Its boundary interpolation is a useful improvement over silently truncating the prior interval to interior grid points.

## 3. Selection-rule compatibility: the critical remaining gate

The new shared-likelihood function explicitly implements a **one-sided** selection rule:

```
d_obs >= detection_threshold_mpc
beta(H0; theta_k) =
    sum_j w[j,k] * Phi((c*z[j]/H0 - threshold)/sigma).
```

This agrees algebraically with the one-sided `dataset_likelihood.event_detection_efficiency` helper, which also uses `P(d >= threshold)`.

It does **not** agree with `scripts/run_selected_likelihood_sbc.py`, whose generator and inference use the finite open interval

```
DETECTION_MIN < d_obs < DETECTION_MAX
```

and whose current experiment contains one event per replicate. The finite-interval `selection.detection_efficiency_h0` and `selection.posterior_h0_selected` are a distinct model. Their results cannot be used as calibration evidence for the new one-sided multi-event likelihood.

The current branch tests reconstruct the one-sided beta and event likelihood algebra independently at fixed nodes. In the source paths inspected for this review, no multi-event prior-predictive generator was found that generates several detected events under the exact one-sided rule and shared population model. Therefore the one-sided generator/inference match is **not yet established** by the existing one-event finite-interval SBC runner.

Before calibration, implement or identify a multi-event generator that:
- draws one H0 once per replicate;
- draws multiple hosts and noisy distances conditional on that same H0;
- applies exactly the one-sided `d_obs >= threshold` rule to each event;
- uses the same host population and declared base weights assumed by inference;
- records all event-level raw values, seed provenance and rejection counts;
- leaves validation data independent and uses the true inclusion law only for generation/scoring, never as an inference input.

Alternatively, change the likelihood to the finite interval and match a generator to that rule. Do not mix the two conventions.

## 4. Existing test coverage and missing tests

The target branch currently has tests for:
- node-specific host-weight normalization and diagnostics;
- independent reconstruction of per-node multi-event log likelihoods;
- independent reconstruction of one-sided beta from the same node-specific weights;
- shared-node log-sum-exp marginalization;
- normalized H0 posterior construction;
- rejection of heterogeneous event noise under the common-beta implementation;
- rejection of events below the one-sided detection threshold.

Still required before a scientific result:
1. **Full shared-node factorization adversary:** explicitly compare the implemented shared-node marginal likelihood with the incorrect product of event-wise marginalized likelihoods on a small two-node/two-event fixture.
2. **Event-order and node-order invariance:** permute events, and separately permute nodes together with their posterior weights.
3. **Single-event reduction:** compare the one-node shared likelihood and resulting posterior with the existing one-event one-sided selected-likelihood calculation.
4. **No-selection limit:** increase the threshold's acceptance probability toward one and verify beta tends to one and selected/unselected results coincide.
5. **Independent finite-host enumeration:** verify event numerator and beta on a tiny catalogue without importing the production helper.
6. **Matched multi-event SBC:** run a small smoke generator first; independently audit the raw data and PITs before scaling.
7. **H0 grid refinement after node marginalization:** compare final posterior densities and quantiles on nested grids.
8. **Quadrature refinement after likelihood propagation:** compare the final H0 log marginal likelihood across orders, not just coefficient moments/evidence/predictions.
9. **Prior-scale sensitivity:** halve/double slope prior scales and quantify effects on host weights and H0 posterior.
10. **Failure accounting:** preserve non-converged quadrature, invalid beta, underflow, and optimizer failures as explicit non-success attempts.

The existing unit tests establish algebraic consistency of the implemented primitive; they do not establish calibration, coverage, or robustness across the planned design.

## 5. Source-level cautions

- The new function rejects heterogeneous event noise because its selection denominator uses a single common sigma. This is explicit and preferable to silently applying one beta to heterogeneous events.
- The event distance is checked against the one-sided threshold, but the inference function alone cannot establish that the data were generated by that exact rule.
- The node-specific weights used for event likelihood and beta are derived from the same log-weight matrix, which avoids a separate normalization mismatch.
- The optional base weights must be justified as population mass/rate weights, not tuned after looking at H0 results.
- A high quadrature ESS is not proof of quadrature convergence; final H0 likelihood stability is a separate requirement.
- The modular/cut interpretation remains: validation data determine the inclusion-coefficient posterior; GW events do not update those coefficients.
- The construction is not an event-rate likelihood and does not establish a real detector selection function or real survey completeness model.

## 6. Decision

**Software direction: conditional GO.** The target branch now contains a mathematically coherent shared-node multi-event likelihood primitive for a one-sided, common-noise toy selection rule, with node-specific IPW weights reused in both the event mixture and beta.

**Calibration direction: NO-GO for confirmatory claims.** First establish a matched multi-event generator, complete the missing adversarial/limiting-case tests, add downstream H0 quadrature convergence, and independently audit a smoke artifact. Until then, describe this as a source-level likelihood prototype with unit-level algebraic tests—not as validated H0 posterior propagation.

No LVK/GWTC analysis, real observational H0 measurement, real detector validation, or universal estimator-superiority claim follows from this implementation.
