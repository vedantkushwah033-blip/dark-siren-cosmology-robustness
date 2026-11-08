# Design contract: propagate finite-validation inclusion uncertainty into toy dark-siren inference

**Status:** specification only. No H0 implementation or result is claimed here.

**Scope:** synthetic low-redshift Euclidean toy model only. No LVK/GWTC analysis, real detector model, real-survey completeness validation, or observational H0 measurement.

## 1. Purpose and separation of stages

The finite-validation density-recovery benchmark evaluates reconstruction of a declared host-redshift density. It does not evaluate downstream H0 posterior calibration. Density IAE alone cannot guarantee posterior calibration because the event likelihood is nonlinear in the density and the GW detection normalization also depends on that density.

The next stage must propagate inclusion-model uncertainty through both (a) the host-redshift mixture in each event likelihood and (b) the GW detection-efficiency normalization conditional on the fixed number of detected events. Galaxy-catalogue inclusion and GW-event detection remain distinct processes.

## 2. Validation model and modular posterior

For validation observations (z_j,R_j), with validation units sampled independently of retention outcomes and both labels observed, use a declared degree-d standardized-redshift logistic model:

$$
R_j|z_j,\theta \sim Bernoulli(\pi(z_j;\theta)),\quad
logit\,\pi(z;\theta)=\sum_{k=0}^{d}\theta_k x(z)^k.
$$

The validation likelihood is

$$
L_V(\theta)=\prod_{j\in V}\pi(z_j;\theta)^{R_j}
[1-\pi(z_j;\theta)]^{1-R_j}.
$$

A proper prior p(theta) must be declared before simulation. The modular validation posterior is p(theta|V) proportional to L_V(theta)p(theta). This is a cut construction: GW data do not update inclusion coefficients. It is not automatically the full joint posterior of a generative model.

**Fixed-basis prior requirement:** the current MLE fitter standardizes redshift using the realized validation sample's mean and standard deviation. That is a numerical MLE parameterization, but a fixed prior placed naively on those sample-standardized coefficients would induce a validation-sample-dependent prior on the inclusion curve. Define the Bayesian prior in a fixed domain coordinate such as u=(z-z_min)/(z_max-z_min) in [0,1], using basis [1,u] or [1,u,u^2], or transform the MLE coefficients into that fixed basis before applying the prior. Declare prior location and scale before confirmatory simulation.

The inverse observed-information covariance from the current MLE fitter is an asymptotic covariance estimate, not a Bayesian posterior. It must not be silently relabelled as posterior uncertainty or treated as exact posterior draws. A downstream implementation must use a declared posterior sampler/approximation, assess convergence and tails, and retain failed/non-evaluable draws. A Laplace approximation may be a separately labelled sensitivity arm only after comparison against a more reliable reference on selected cases.

## 3. Density for each coefficient draw

For posterior draw theta_b, compute probabilities only for analysis-catalogue objects and construct

$$
\widehat g_b(z)=\sum_{i\in C}K_h(z-z_i)/\pi(z_i;\theta_b),
\qquad
\widehat f_b(z)=\frac{\widehat g_b(z)r(z)}
{\int \widehat g_b(u)r(u)du}.
$$

Freeze the bandwidth rule, boundary treatment, redshift grid and merger-rate weight r(z) across draws and comparison arms. Use the same validation posterior draw across the entire H0 grid for a given Monte Carlo index. Do not independently resample coefficients at each H0 value.

The analysis catalogue is fixed conditional on the simulated dataset. Oracle inclusion probabilities remain simulation-only diagnostics and are never supplied to the estimated-inclusion fitter.

## 4. Event likelihood and GW selection

For the existing toy model, d_k|z_k,H0 is Normal(c z_k/H0,sigma^2), and an event is detected iff d_k >= d_th.

For each coefficient draw b, calculate the unconditional event-mixture log likelihood

$$
\ell_{uncond}(H_0;\theta_b)=
\sum_{k=1}^{N_{det}}\log\left[
\int \widehat f_b(z)p(d_k|z,H_0)dz\right].
$$

For the same draw, calculate

$$
\beta_{GW}(H_0;\theta_b)=
\int \widehat f_b(z)P(d\ge d_{th}|z,H_0)dz.
$$

Conditional on the observed number of detected events, the per-draw selected log likelihood is

$$
\ell_b(H_0)=\ell_{uncond}(H_0;\theta_b)
-N_{det}\log\beta_{GW}(H_0;\theta_b).
$$

**Critical ordering:** compute beta_GW separately for every coefficient draw and every H0 value, subtract the selection term inside each draw's log likelihood, and only then marginalize over theta.

## 5. Marginalization

Given B approximately independent draws theta_b from p(theta|V), estimate

$$
\log L_{marg}(H_0)=logsumexp_b\{\ell_b(H_0)\}-\log B.
$$

Then p(H0|d,V,C) is proportional to p(H0) exp(log L_marg(H0)). Declare the H0 prior explicitly; a flat prior in H0 is not flat in a transformed parameter. Check that posterior mass is not truncated at either grid boundary.

Do not replace this with a plug-in likelihood at posterior-mean coefficients, an MLE plus plug-in covariance, an average beta followed by one selection correction, an average unconditional likelihood followed by a selection correction, or an average of normalized per-draw posteriors when the target is the marginal likelihood above.

Assess Monte Carlo error by increasing B, using independent posterior-draw batches, and checking log-likelihood differences and posterior summaries. H0-grid convergence is a separate check.

## 6. Required comparison arms

1. Parent/oracle density diagnostic (not an operational estimator).
2. Catalogue-only unweighted density.
3. Known-inclusion IPW (oracle simulation diagnostic).
4. Estimated-inclusion plug-in point estimate (comparator only).
5. Estimated-inclusion modular marginalization over p(theta|V).

The fixed-analysis-catalogue and disjoint finite-parent designs answer different questions and must not be pooled. In the former, validation information changes while the analysis catalogue is fixed. In the latter, increasing validation size also removes units from the available analysis frame.

## 7. Deterministic and statistical gates

Deterministic tests:
- A degenerate coefficient distribution agrees with the existing plug-in selected likelihood.
- B=1 equals the sole draw's selected log likelihood.
- A hand-computable two-draw case verifies log-mean-exp.
- An adversarial test rejects average-beta selection correction.
- Every density is finite, nonnegative and integrates to one.
- Every beta_GW lies in (0,1]; invalid draws are not silently dropped.
- Reordering draws does not change results beyond floating-point tolerance.
- Repeated identical draws reproduce the plug-in likelihood.
- An independent auditor reconstructs posterior normalization, quantiles and boundary mass from raw H0 grids and density vectors.

Statistical validation:
- Establish redshift-grid and H0-grid convergence first.
- Separately establish Monte Carlo convergence in coefficient draws.
- Use simulation-based calibration for parameters generated from the declared prior/model, with rank diagnostics and sampler convergence checks.
- Use repeated fixed-truth datasets for frequentist interval coverage and report Monte Carlo uncertainty.
- Distinguish SBC, fixed-truth coverage, PIT diagnostics, posterior predictive checks and numerical grid stability.
- Report all failed, non-evaluable and not-run replicates/draws.

## 8. Implementation order

1. Finish and inspect both 100-replicate density-recovery confirmatory artifacts and independent audits.
2. Reconcile disjoint-frame and fixed-analysis-catalogue estimands without pooling.
3. Specify and test the coefficient posterior sampler independently of GW likelihood code.
4. Add a pure function mapping one coefficient draw to one normalized host density and selected H0 log-likelihood vector.
5. Add stable log-mean-exp marginalization and adversarial selection-order tests.
6. Add raw draw-level artifacts, immutable config, manifest hashes and an independent auditor.
7. Run and audit a smoke benchmark, then freeze the confirmatory design.
8. Run SBC and fixed-truth coverage only after sampler, density, likelihood, selection and numerical-convergence gates pass.

**State at creation:** finite-validation probability and density-recovery stages exist. The finite-validation H0 marginalization described here is not yet implemented or validated.
