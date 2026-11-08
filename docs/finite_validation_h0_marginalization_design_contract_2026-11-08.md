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

Given proposal draws theta_b from a declared proposal q(theta), with unnormalized importance weights w_b proportional to L_V(theta_b)p(theta_b)/q(theta_b), estimate the validation-posterior expectation by self-normalized importance sampling:

$$
\log L_{marg}(H_0)\approx
logsumexp_b\{\ell_b(H_0)+\log w_b\}
-logsumexp_b\{\log w_b\}.
$$

The current implementation candidate uses a multivariate Student-t proposal centered at the fixed-basis posterior mode, with shape based on the inverse posterior Hessian. Its coefficient proposals are **not equally weighted posterior draws**; normalized importance weights must be retained for every downstream expectation. Report importance ESS, maximum normalized weight, and Monte Carlo error, and mark low-ESS cases NON_EVALUABLE rather than silently proceeding. Exact/effectively exact posterior draws are a separate possible implementation and recover the unweighted log-mean-exp formula as the special case of equal weights.
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


## Candidate fixed-basis prior and prior-predictive check (exploratory)

The DGP's true inclusion model, expressed in the fixed coordinate u=(z-0.01)/0.39, has coefficients

| Curvature gamma | Intercept in fixed basis | Linear coefficient | Quadratic coefficient |
|---:|---:|---:|---:|
| -1.2 | -0.3201 | -1.0 | -1.2 |
| 0.0 | -0.6259 | -1.0 | 0.0 |
| +1.2 | -0.9841 | -1.0 | +1.2 |

A candidate weakly regularizing prior for the Bayesian stage is independent Normal(0, s^2) coefficients in the fixed [1,u] or [1,u,u^2] basis, with s=1.0 as the candidate reference and s in {0.75,1.5} as prior-sensitivity arms. This is **not yet a frozen confirmatory prior**.

A reproducible prior-predictive simulation used 500,000 coefficient draws per degree/scale and evaluated probabilities on 101 equally spaced u values. For the degree-2 model:
- s=0.75: fraction of prior-predictive grid probabilities outside [0.01,0.99] = 0.0000273.
- s=1.00: fraction = 0.0009129; central 95% probability range at u=0 is [0.1242,0.8766] and at u=1 is [0.0325,0.9676].
- s=1.50: fraction = 0.0176122; central 95% ranges are [0.0505,0.9501] at u=0 and [0.0062,0.9941] at u=1.

For degree 1, the corresponding extreme-grid fractions are 0.0000010, 0.0001792, and 0.0095059 for s=0.75, 1.00, and 1.50. These are prior-predictive diagnostics, not posterior calibration results. The candidate s=1.0 prior covers the DGP coefficient magnitudes in this benchmark without placing much prior mass on near-deterministic inclusion across the domain. Prior sensitivity remains mandatory, and the candidate must be evaluated for separation regimes and coefficient-posterior calibration before confirmatory use.

The sampler implementation should use the fixed basis directly. It should not assign an identical Normal prior to the current fitter's validation-sample-standardized coefficients.


## Implementation update (2026-11-08; focused CI verified)

The initial implementation-order list above records the plan at document creation. The following software components have since been added on the research branch:

1. `inclusion_posterior.py`: four-chain proper-prior logistic posterior sampler in the fixed domain basis, with split-Rhat and effective-sample-size diagnostics.
2. `bayesian_inclusion.py`: independent multivariate-Student-t importance sampler around the fixed-basis posterior mode. It retains normalized importance weights and reports importance ESS, maximum weight and relative evidence Monte Carlo error. It is not treated as an equal-weight posterior sample.
3. `finite_validation_h0_marginalization.py`: high-level propagation from coefficient proposals through IPW host-density reconstruction, per-draw unconditional event likelihood, per-draw GW detection efficiency, draw-specific selection correction, and weighted marginalization.
4. `finite_validation_h0.py`: batched/vectorized likelihood evaluation, explicit H0-prior posterior normalization, 68% and 95% intervals, and a boundary-mass diagnostic. The vectorized likelihood is checked against the existing single-density selection-likelihood implementation.

The fixed-basis importance sampler was tested against two-dimensional numerical quadrature and the existing multi-chain MCMC posterior. The H0 pipeline tests include an actual weighted inclusion-posterior output passed through host-density reconstruction and selected-likelihood marginalization, and an explicit test showing that the average-beta shortcut differs from the correct per-draw selection correction.

**Focused CI:** inclusion-posterior workflow [#42479890034](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/42479890034) completed successfully on Python 3.10, 3.11 and 3.12 at PR head `c47ab95e26adaede691f0d71306eaf3b517fa892`. The focused sampler, vectorized H0 likelihood, posterior normalization, and finite-validation propagation tests passed in all three jobs.

The general numerical-diagnostics workflow for that exact head (#42479890015) was cancelled, so this is **not** a claim that the entire repository test matrix passed at that head. Earlier, the fixed-catalogue independent-analysis-auditor tests passed on Python 3.10, 3.11 and 3.12 in workflow #42478908404. The full confirmatory artifact's saved contrast JSON was independently recomputed locally from the raw attempt CSV and matched to floating-point precision.

### Still not established

- No confirmatory H0 simulation has been run through the new pipeline.
- No SBC, fixed-truth H0 coverage, PIT assessment, or posterior predictive validation has been performed for finite-validation inclusion uncertainty.
- The candidate Gaussian prior scale s=1.0 and sensitivity values 0.75/1.5 have not yet been compared in an H0 experiment; the prior remains provisional.
- The importance-sampling approximation has not yet been stress-tested across the full planned curvature/validation-size/separation design.
- No claim of improved H0 calibration or coverage is warranted.

### Next gates

1. Create a small reproducible end-to-end H0 pilot runner with a raw artifact, immutable config, manifest and independent auditor.
2. Compare the MCMC and Student-t importance posterior summaries and the marginalized H0 likelihood on the same validation datasets; report Monte Carlo uncertainty and ESS.
3. Freeze the prior only after the prior-scale sensitivity and prior-predictive checks are reviewed.
4. Establish redshift-grid, H0-grid, coefficient-draw and KDE-bandwidth convergence separately.
5. Only then predeclare and run SBC and repeated fixed-truth coverage experiments with complete SUCCESS/NON_EVALUABLE/FAILED/NOT_RUN accounting.
