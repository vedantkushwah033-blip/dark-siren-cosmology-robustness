# Canonical estimands and uncertainty semantics for the five-arm study

**Status:** source-review methods note; not an experimental result or preregistration amendment.  
**Code reference:** PR #55 head 95e8164bd1c17756e5416566ec3a28475b15f1b9.  
**Scope:** synthetic low-redshift benchmark only.

## 1. Probability model and target

Let a parent galaxy redshift be Z distributed according to g(z) on Z-domain [0.01, 0.40]. Let R in {0,1} indicate inclusion in the observed catalogue, with P(R=1|Z=z)=pi_0(z), where 0 < pi_0(z) <= 1. Let b(z)>0 be the declared host/merger-rate weight. The intended intrinsic host-redshift density is

q_0(z) = g(z)b(z) / integral[g(u)b(u) du].

In the current synthetic generator, b(z)=1+1.5(z-z_min)/(z_max-z_min). This is an intrinsic host prior before GW detection conditioning, not the redshift distribution of detected GW events.

The catalogue redshift density is g_C(z)=g(z)pi_0(z)/bar_pi, where bar_pi=integral[g(u)pi_0(u)du]. The population identity behind known-probability IPW is

[g_C(z)b(z)/pi_0(z)] / integral[g_C(u)b(u)/pi_0(u)du] = q_0(z).

Equivalently, E_C[f(Z)b(Z)/pi_0(Z)] / E_C[b(Z)/pi_0(Z)] = E_g[f(Z)b(Z)] / E_g[b(Z)].

This requires positivity on target support and correct inclusion probabilities. It does not establish finite-sample stability, KDE accuracy, or cosmological posterior calibration.

## 2. Five arms and their intended roles

For catalogue observations z_i^C and parent observations z_i^P:

- A: catalogue-only weights b(z_i^C).
- B: full-parent oracle weights b(z_i^P), unavailable in a real incomplete-catalogue analysis.
- C: known-probability IPW weights b(z_i^C)/pi_0(z_i^C).
- D: plug-in IPW weights b(z_i^C)/pi(z_i^C; theta_MAP).
- E: node-specific weights b(z_i^C)/pi(z_i^C; theta^(k)), with validation-posterior node masses w_k.

For the target q_0(z) proportional to g(z)b(z), the weighted KDE must apply b exactly once, either through the sample weights or as a post-KDE multiplier—not both.

## 3. Validation-sample likelihood

The quadratic-logit generator draws an independent validation frame from the same parent redshift law, then selects a simple random sample without replacement from that frame. Conditional on sampled redshifts, the retention labels have Bernoulli likelihood

L_V(theta)=product_j pi(z_j^V;theta)^{r_j^V}[1-pi(z_j^V;theta)]^{1-r_j^V}.

With a proper Gaussian coefficient prior p(theta), p(theta|V) is proportional to L_V(theta)p(theta).

The validation sample is independent of the analysis parent/catalogue sample by construction. The sizes 100, 300 and 1000 are nested prefixes of one random permutation, so comparisons across validation sizes are paired and nested, not independent replications.

The fitted model in D/E is degree-one logistic. At nonzero gamma, the true DGP includes a quadratic term. Therefore D/E are structurally misspecified in those cells; their posterior is conditional on the degree-one model, not a posterior over the true quadratic model.

## 4. Arm E: distinct uncertainty summaries

For posterior nodes theta^(k) and masses w_k, let qhat_k(z) be the node-specific density estimate. Two density-stage summaries answer different questions:

Posterior-mean density loss: L(sum_k w_k qhat_k, q_0).

Posterior expected node loss: sum_k w_k L(qhat_k, q_0).

They are not generally equal. Report both with exact labels and retain node-level losses. Neither is cosmological likelihood marginalization.

For a future analysis with one shared inclusion-model coefficient vector across all GW events, the validation-posterior marginal likelihood must be

p(d_vector|H0,V) = integral p(d_vector|H0,theta)p(theta|V)dtheta
                 ~= sum_k w_k p(d_vector|H0,theta^(k)).

If event likelihoods factor conditionally on the shared coefficient vector, this is sum_k w_k product_e p(d_e|H0,theta^(k)), not product_e sum_k w_k p(d_e|H0,theta^(k)). The latter incorrectly redraws the common catalogue-inclusion model independently for each event.

The existing log_marginalize_shared_parameter helper encodes the shared-node summation pattern, but this review has not verified an end-to-end H0 implementation. Averaging weights, averaging normalized densities, and marginalizing the final H0 likelihood are different operations.

## 5. Keep selection mechanisms separate

Catalogue retention R and GW detection S_e are distinct stochastic mechanisms. Catalogue IPW reconstructs the intrinsic galaxy/host distribution. GW selection normalization belongs to the event-population likelihood and depends on the GW detection model and cosmological/population parameters.

A future H0 likelihood must state whether it conditions on the observed number of detections or models the count, and include the corresponding GW selection normalization. It must not apply an additional missing-host mixture correction on top of an already IPW-corrected catalogue term unless a derivation shows the two factors address non-overlapping components.

## 6. Positivity and estimability

The population identity requires pi_0(z)>0 wherever q_0(z)>0. In finite samples, very small pi_0(z) can make IPW practically unusable. Report parent-measure and catalogue-measure summaries of pi_0, target-mass fractions below prespecified probability thresholds, theoretical and realized maximum inverse weights, Kish and entropy ESS, top 1% and 5% normalized weight mass, validation support in low-inclusion regions, and all failure/non-evaluable rates.

Do not silently clip weights. Any stabilization or truncation analysis changes the weighting rule and potentially the target; define and interpret it separately.

## 7. Immediate source-level blocker: possible double application of host weight

Source trace at PR #55 head:

1. run_absolute_bandwidth_density_smoke.py defines base_p, base_c and merger using the same linear function b(z).
2. It passes base_p and base_c into build_five_arm_weights.
3. build_five_arm_weights places those base weights into arms A-E.
4. The smoke runner passes the resulting arm weights as sample_weights to estimate_host_density and also passes merger_weight=merger.
5. estimate_host_density multiplies the weighted KDE by merger_weight before normalization.

For arm B, for example, the KDE is already based on b(z_i^P), and the final density is multiplied by b(z) again. In the large-sample limit this targets g(z)b(z)^2, not g(z)b(z). For arms C/D/E the same duplication occurs around the inclusion correction.

**This is a potentially consequential estimand bug. It blocks interpretation of current five-arm density metrics as recovery of the stated target until resolved.** The repair depends on the intended estimand:

- If the target is g(z)b(z), use base weights b(z_i) with no second post-KDE multiplier; or use unweighted/selection-corrected KDE followed by one multiplier b(z).
- If another target is intended, define it explicitly and change the truth density, oracle and every arm consistently.

No production code was changed in this review. Implement the repair on a dedicated branch after tests demonstrate the population target for all arms.

## 8. Verification status

- Source-level tracing of the weight construction and KDE multiplication: reviewed.
- Double-application issue: identified, not fixed.
- Unit tests proving the target density for each arm: not verified.
- Seven-gamma density experiment: not run.
- H0 shared-node marginal likelihood: helper exists; end-to-end implementation not verified.
- No experimental result or publication claim is made here.
