# Existing event-likelihood interface audit for inclusion-posterior propagation

**Status:** source-level interface audit; no end-to-end integration or calibration result is claimed.  
**Target branch:** `research/inclusion-posterior-quadrature-prototype`  
**Scope:** synthetic low-redshift Euclidean model only.

## 1. Files inspected

The audit inspected the current base-branch implementations of:

- `src/dark_siren_robustness/inference.py`
- `src/dark_siren_robustness/hardened_inference.py`
- `src/dark_siren_robustness/selection.py`
- `scripts/run_selected_likelihood_sbc.py`
- `docs/joint_selection_incompleteness_likelihood_derivation_2026-10-25.md`

The same interfaces must be rechecked on the target branch before implementation; this note does not assume that a later branch has identical source.

## 2. What the existing API actually provides

### Unselected single-event likelihood/posterior

`hardened_inference.posterior_h0` accepts one observed GW distance, a host-redshift vector, and optional host weights. Internally it forms the host-mixture log likelihood

```
log L_e(H0; theta) =
    logsumexp_j[
        log(w_j(theta)) +
        NormalLogPDF(d_e; c z_j / H0, sigma)
    ].
```

It then applies a uniform H0 prior on the requested support and normalizes the posterior density numerically. The function is a single-event posterior routine, not a multi-event joint log-likelihood API.

`inference.posterior_h0_low_z` is an older/simpler single-event routine. Its docstring explicitly says it omits GW detection-selection normalization. It is not the integration target for a selection-aware analysis.

### Event-selection normalization

`selection.detection_efficiency_h0` computes the finite-interval detection efficiency for the toy rule

```
d_min < d_obs < d_max,
d_obs | (H0, host j) ~ Normal(c z_j / H0, sigma^2).
```

For a normalized host mixture, it calculates

```
beta(H0; theta) =
    sum_j w_j(theta) * [
        Phi((d_max - c*z_j/H0)/sigma)
      - Phi((d_min - c*z_j/H0)/sigma)
    ].
```

This is the appropriate structural beta for that finite-interval generator, provided the exact same redshifts, normalized population weights, noise scale and open interval are used. It is not interchangeable with a one-sided threshold.

`selection.posterior_h0_selected` takes one observed event, calls the single-event posterior routine, divides the resulting density by beta, and renormalizes it. It therefore provides a single-event selected posterior, not the joint likelihood needed when many events share one uncertain inclusion-coefficient vector.

### Generator

`scripts/run_selected_likelihood_sbc.py` draws H0 once per replicate, then repeatedly draws a host and noisy distance conditional on that fixed H0 until the distance lies inside the finite detection interval. This is the matched one-event prior-predictive calibration design. A multi-event extension must keep H0 fixed across all events in a replicate and apply the same declared selection rule to each event.

## 3. Integration gap

The current API does **not** expose a raw event log likelihood for a vector of events and does not provide a joint selected log likelihood over a shared coefficient node. Calling `posterior_h0_selected` separately for each event and multiplying its outputs would be wrong for the intended shared-`theta` modular posterior: it would normalize each event posterior separately and, if each event were independently marginalized over quadrature nodes, would effectively give each event a different inclusion model.

The correct implementation boundary is a new raw-likelihood layer, not a wrapper around already normalized posterior densities.

## 4. Required shared-node calculation

For quadrature node k, first build the normalized host weights from that same node's fitted inclusion probabilities:

```
w_jk = (1 / pi(z_j; theta_k)) /
       sum_l (1 / pi(z_l; theta_k)).
```

For event e, H0 grid point h and node k, compute the **raw** host-mixture event log likelihood:

```
ell[e,k,h] = logsumexp_j(
    log(w[j,k]) +
    NormalLogPDF(d[e]; c*z[j]/H0[h], sigma[e])
).
```

For the current common-noise toy model, if all events share sigma and the same open detection interval, calculate beta[k,h] with `selection.detection_efficiency_h0` using the identical node-specific weights. The conditional joint selected log likelihood is then

```
joint[k,h] = sum_e ell[e,k,h] - N_events * log(beta[k,h]).
```

Finally marginalize **once** over the validation-posterior quadrature nodes:

```
log_marginal[h] = logsumexp_k(log(a[k]) + joint[k,h]).
```

If event-specific noise or event-specific detection regions are introduced, beta becomes event-specific and the denominator must be derived as the corresponding sum of log efficiencies. Do not retain the `N_events * log(beta)` shortcut without the common-selection assumptions.

The current `log_marginalize_shared_parameter` helper can perform the final node integration, but the raw event terms and beta arrays must be built first. The node-specific host-weight matrix must be the common source for both event likelihood and beta; independently normalizing or reconstructing weights in the two paths risks a silent mismatch.

## 5. Required numerical and compatibility tests before a simulation

1. **Single-event reduction:** one event and one quadrature node reproduce the existing `posterior_h0_selected` density after applying the same prior and grid normalization.
2. **Fixed-node multi-event product:** at a fixed theta node, the joint numerator equals the product of the event-level host-mixture likelihoods, evaluated in log space.
3. **Shared-node distinction:** a two-node/two-event fixture must show that one marginalization of the joint product differs from multiplying event-wise marginalized likelihoods.
4. **Selection match:** beta computed from the same node-specific weights agrees with direct finite-host summation of Gaussian interval probabilities.
5. **Complete-catalogue / no-IPW limit:** constant inclusion probabilities produce uniform normalized host weights and recover the catalogue-only selected likelihood.
6. **No-selection limit:** as the detection interval captures essentially all probability, beta tends to one and the selected result approaches the unselected joint result.
7. **Event-order invariance:** permuting events leaves the joint likelihood unchanged.
8. **Node-order invariance:** permuting quadrature nodes and their weights leaves the marginalized likelihood unchanged.
9. **Numerical stability:** compare direct log-sum-exp calculations with independent small-host enumeration; reject non-positive/non-finite beta rather than replacing it by machine tiny.
10. **Prior and grid handling:** apply the H0 prior once, after joint likelihood marginalization, and normalize once on the actual supported grid.

The one-event SBC runner is not a substitute for a multi-event shared-nuisance calibration study. Once the raw joint likelihood is implemented and independently checked, the next calibration design should draw one H0 per replicate, one shared validation dataset and inclusion-posterior node distribution per replicate, then multiple conditionally independent detected events under that fixed H0. Report PIT/SBC separately from fixed-truth repeated-sampling coverage.

## 6. Stop conditions and claim boundary

- Do not connect the current quadrature prototype to `posterior_h0_selected` by multiplying normalized per-event posterior densities.
- Do not claim shared-coefficient uncertainty propagation until the raw joint likelihood and beta use the same node-specific population weights.
- Do not run a large calibration matrix before the single-event reduction, multi-event factorization, selection-match and independent-enumeration tests pass.
- A successful implementation would remain a modular/cut posterior conditional on the validation data and the declared synthetic finite-interval selection model. It would not be a full joint posterior, real detector model, LVK/GWTC analysis, or observational H0 measurement.
