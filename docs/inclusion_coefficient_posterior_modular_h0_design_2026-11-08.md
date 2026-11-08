# Inclusion-coefficient posterior and modular H0 propagation: design contract

**Status:** design specification; implementation and numerical validation not yet complete.  
**Branch:** `research/inclusion-posterior-modular-design`  
**Scope:** synthetic low-redshift Euclidean dark-siren methodology only.

## 1. Why this stage is necessary

The fixed-analysis-catalogue experiment established a finite-validation effect on synthetic host-redshift density recovery. It used plug-in logistic inclusion probabilities. The current `InclusionFit.covariance` is the inverse observed information at an unpenalized MLE: it is an asymptotic covariance diagnostic, not a Bayesian posterior distribution.

The next stage must propagate uncertainty in the inclusion coefficients through the *same* host-mixture likelihood and the GW detection-selection normalization. It must not draw coefficients from the inverse-information Gaussian and call those draws a posterior.

## 2. Validation-data model and proper prior

For validation observations \((z_i,R_i)\), with \(R_i\in\{0,1\}\), define standardization using the validation covariates only:

\[
x_i=(z_i-\bar z_V)/s_V,\qquad
X_i=(1,x_i,\ldots,x_i^d),\qquad d\in\{1,2\}.
\]

The inclusion model is

\[
R_i\mid z_i,\theta\sim\operatorname{Bernoulli}
\left[\operatorname{logit}^{-1}(X_i\theta)\right].
\]

Use a proper, declared weakly informative Gaussian prior for the first benchmark:

\[
\theta_0\sim N(0,5^2),\qquad
\theta_j\sim N(0,2.5^2),\quad j\geq1.
\]

The prior is on coefficients of the standardized-redshift model. It is not informed by the synthetic truth, target retention, analysis catalogue, GW events, or H0. Prior-scale sensitivity must be reported (at minimum, halve and double the slope scale). The fixed prior is a benchmark convention, not a claim that these scales are universally optimal.

The validation posterior is

\[
p(\theta\mid V)=\frac{1}{Z_V}
\exp\left[
\sum_i R_i X_i\theta-
\sum_i\log(1+e^{X_i\theta})
-\frac12(\theta-\mu_0)^T\Sigma_0^{-1}(\theta-\mu_0)
\right].
\]

The posterior is proper because the Gaussian prior is proper. Separation therefore does not make the posterior undefined, although it can make the data weakly informative and the posterior prior-dominated. Report that condition rather than silently regularizing or discarding it.

## 3. Preferred deterministic integration for this low dimension

The logistic model has only two or three coefficients. Before adding an MCMC dependency, use **adaptive Gauss-Hermite importance quadrature** around the posterior mode as the first numerical implementation:

1. Find the posterior mode \(\hat\theta\) by minimizing the negative log posterior.
2. Compute the positive-definite Hessian \(H\) of the negative log posterior at the mode, including prior precision; set \(C=H^{-1}\).
3. Construct a multivariate Gaussian proposal \(q(\theta)=N(\hat\theta,C)\).
4. For tensor Gauss-Hermite nodes \(u_k\) and normalized standard-normal quadrature weights \(a_k\), transform \(\theta_k=\hat\theta+L u_k\), where \(LL^T=C\).
5. Compute log importance weights
   \[
   \ell_k=\log a_k+\log p(V,\theta_k)-\log q(\theta_k),
   \]
   normalize with log-sum-exp, and retain the normalized quadrature weights \(\tilde a_k\).

The Gaussian proposal is a numerical integration device only; the target remains the exact logistic likelihood times the declared proper prior. Do not report the Laplace proposal itself as the posterior.

Run an order ladder (initially 7, 11, 15, 19, 23 in each dimension). The largest tensor has 12,167 nodes for degree 2. Compare posterior means, standard deviations, selected tail probabilities, posterior-predictive inclusion probabilities, and the final H0 log marginal likelihood across consecutive orders. Stop only when predeclared tolerances are met; otherwise return a non-converged status and expand the ladder. Monitor normalized quadrature-weight ESS, maximum normalized weight, and tail-node influence. A high ESS alone is not a convergence proof.

The implementation must explicitly fail or return NON_EVALUABLE for rank-deficient validation covariates, invalid binary labels, non-finite prior settings, non-positive Hessian eigenvalues, failed optimization, non-finite quadrature weights, or lack of quadrature convergence. Do not clip fitted inclusion probabilities. If numerical probabilities underflow or produce unusable inverse weights, record the failure and its cause.

This deterministic quadrature approach is attractive because the dimension is at most three and it avoids MCMC warm-up, chain mixing, and convergence-diagnostic ambiguity. If quadrature fails in stress cases, add a separately validated MCMC fallback rather than silently switching methods.

## 4. Host-mixture construction for each coefficient node

Let \(C=\{z_j^{\rm cat}\}\) be the observed analysis catalogue, which is fixed across all coefficient nodes and is disjoint from the independent validation frame in the FAC design.

For each quadrature node \(\theta_k\), compute

\[
\pi_{jk}=\operatorname{logit}^{-1}
\left(X(z_j^{\rm cat})\theta_k\right),\qquad
w_{jk}=1/\pi_{jk}.
\]

Normalize weights stably (scale by their maximum before summing). Record the minimum inclusion probability, maximum raw inverse weight, Kish ESS, entropy ESS, maximum normalized weight, top-1%/5% weight mass, and any numerical failure. Do not truncate or winsorize weights in the primary analysis; any stabilized/truncated alternative must be a separately declared sensitivity arm.

For a discrete galaxy catalogue, the event host-mixture likelihood can use the normalized discrete weights directly. The existing KDE density estimator remains useful for density-recovery diagnostics, but the event likelihood should not introduce an arbitrary second smoothing layer unless that is explicitly the model being studied.

## 5. Correct shared-nuisance marginalization

For detected event \(e\), let \(\mathcal L_e(H_0;\theta_k)\) be the **unconditional** event likelihood formed from the observed catalogue redshifts and node-specific normalized IPW weights. Let \(\beta(H_0;\theta_k)\) be the detection efficiency computed from the same declared population mixture, the same node-specific weights, and the same toy detector/noise model.

For a fixed-size sample of \(N_{\rm det}\) conditionally independent detected events, conditional on one shared catalogue-inclusion coefficient vector, the joint selected log likelihood is

\[
\log \mathcal L_{\rm sel}(H_0;\theta_k)
=\sum_{e=1}^{N_{\rm det}}\log\mathcal L_e(H_0;\theta_k)
-N_{\rm det}\log\beta(H_0;\theta_k).
\]

The coefficient vector is shared across events. Therefore the modular/cut marginal likelihood is

\[
\log \mathcal L_{\rm marg}(H_0)
=\operatorname{logsumexp}_k
\left[
\log\tilde a_k+
\log\mathcal L_{\rm sel}(H_0;\theta_k)
\right].
\]

**Do not** integrate each event separately and multiply the separately marginalized event likelihoods. That would incorrectly give each event an independent inclusion model. **Do not** average log likelihoods. Use log-sum-exp of the joint conditional likelihood over the validation posterior quadrature weights.

The fixed-size conditional-on-detection likelihood is not an event-rate likelihood. Rate information and Poisson counts remain out of scope unless a separate generative rate model is derived.

## 6. Catalogue inclusion and GW detection are separate mechanisms

- \(\pi(z;\theta)\) is the probability that a galaxy enters the observed catalogue.
- \(\beta(H_0;\theta)\) is the probability that a GW event is detected under the declared distance-noise/threshold model and the population mixture.
- The same node-specific IPW population weights may enter the population mixture used in \(\beta\), but catalogue inclusion is not itself the GW detection probability.
- If the catalogue is already IPW reconstructed, do not add a second missing-host mixture or completeness correction unless a separate joint likelihood derives it.
- The current toy \(\beta\) implementation must be matched to the exact data-generating detection rule; a one-sided threshold and a finite interval are not interchangeable.

The resulting construction is explicitly a **modular/cut posterior**: validation data determine \(p(\theta\mid V)\), and GW data do not feed back into \(\theta\). It is not a fully joint posterior. That choice must appear in all output metadata and scientific descriptions.

## 7. Required numerical and statistical validation gates

### Unit and property tests
- Recover the analytic logistic posterior in a near-Gaussian large-validation regime against an independent high-accuracy reference.
- Check the log posterior against direct Bernoulli probabilities on hand-computable cases.
- Verify quadrature weights sum to one, are finite/non-negative, and log-sum-exp is invariant to a common log-likelihood offset.
- Verify a one-node/constant-likelihood integration returns the expected result.
- Verify shared-theta marginalization differs from product-of-eventwise marginalizations in a constructed example.
- Verify the same node-specific weights enter event likelihood and \(\beta\).
- Reject invalid priors, rank deficiency, non-finite values, non-positive inclusion probabilities and failed convergence.

### End-to-end synthetic checks
- Simulate validation and analysis frames independently, with validation labels generated under the declared inclusion law.
- Generate GW events under the declared toy population and selection mechanism.
- Keep the true \(\theta\) for scoring only; never pass it to the inference routine.
- Compare catalogue-only, plug-in estimated-IPW, quadrature-marginalized estimated-IPW, and oracle-inclusion diagnostics.
- Check the numerical integration over \(\theta\) at multiple quadrature orders and the H0 grid at multiple resolutions.
- Store per-node weights/diagnostics or a reproducible compressed representation; store raw validation records and event records.
- Independently reconstruct the posterior target, quadrature weights, host weights, event likelihood, selection efficiency, joint log likelihood and final H0 posterior from raw artifacts.

### Calibration and coverage
Density recovery, numerical quadrature convergence, H0-grid convergence, posterior calibration, and frequentist coverage are distinct claims. Run enough independent simulation replicates to estimate 68% and 95% coverage with binomial uncertainty. Report PIT/SBC diagnostics separately from frequentist fixed-truth coverage. Stratify by curvature, validation size, retention, and prior scale. Report failure rates without imputing posterior outcomes.

## 8. Pre-implementation review questions

1. Confirm whether the event likelihood uses a discrete host mixture or a KDE-integrated redshift density. The existing dataset-level likelihood naturally supports the discrete mixture; the FAC KDE result is a separate density target.
2. Match the selection normalization to the exact finite or one-sided detection rule in the generator.
3. Confirm that the finite validation frame is representative of the parent redshift distribution and has complete retention labels. If not, model validation sampling design explicitly.
4. Treat prior-scale sensitivity as a first-class result, especially in weakly identified validation cells.
5. Do not use the old unpenalized MLE covariance as a posterior approximation in the primary analysis.

## 9. Explicit claim boundary

A successful implementation would establish only that a declared proper validation posterior can be numerically propagated through this synthetic, low-redshift, Euclidean toy likelihood under stated assumptions. It would not establish performance on LVK/GWTC data, a real galaxy survey, a real detector-selection function, or a measured H0 value. No novelty or estimator-superiority claim follows automatically.


## 10. Prototype implementation snapshot (2026-11-08)

The current prototype branch now contains the following implementation pieces:

- `inclusion_posterior.py`: a proper Gaussian-prior logistic validation posterior, adaptive tensor Gauss-Hermite importance quadrature, node-weight diagnostics, and a single canonical order ladder. The ladder checks posterior means, covariance, posterior-predictive inclusion probabilities, and the log evidence normalizer. It requires two consecutive passing order transitions by default; quadrature ESS and maximum node mass are diagnostics, not convergence criteria.
- `shared_inclusion_likelihood.py`: a discrete catalogue host-mixture likelihood with one shared inclusion-coefficient node across all events. Event likelihoods and the one-sided detection efficiency use the same node-specific log-normalized IPW host weights. The current selection implementation deliberately requires a common Gaussian distance-noise scale across events and uses the one-sided rule (d_{\rm obs}\geq d_{\rm threshold}).
- `h0_posterior.py`: normalization of the marginalized log likelihood under a declared uniform H0 prior, explicit inclusion of prior boundaries in the integration grid, posterior quantiles, and a coarse/fine density-and-quantile comparison utility.
- Tests include finite-difference checks for the logistic target derivatives, shared-parameter marginalization checks, independent reconstruction of the nodewise selected likelihood and its posterior-node marginalization, and constant-likelihood recovery of a uniform H0 prior.

The duplicate quadrature-ladder module was removed in favor of the canonical implementation inside `inclusion_posterior.py). The canonical ladder now includes log-evidence stability in its convergence decision.

**Validation status:** these are implemented prototype components and authored tests, not yet a validated end-to-end inference pipeline. The latest branch head has not reported a CI run or combined status check at the time of this update. The previously successful workflow runs correspond to earlier commits and do not validate the current branch head. No full synthetic H0 calibration, SBC, frequentist coverage, prior-scale sensitivity, or independent end-to-end raw-artifact audit has been completed. The code remains limited to the declared synthetic toy model.
