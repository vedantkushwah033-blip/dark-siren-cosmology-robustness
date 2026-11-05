# First-order perturbation of normalized IPW host reconstruction

**Date:** 2026-11-05  
**Status:** Mathematical clarification for the synthetic benchmark; not an empirical result.

## Executive correction

The familiar expression
```
(f_hat(z) - f(z)) / f(z) = -[delta(z) - E_f delta]
```
is the first-order perturbation for a *pointwise, unsmoothed population density ratio* under a small relative inclusion-probability error. It is not, without qualification, the exact first-order expression for the finite-sample, reflected, weighted KDE implemented in this repository.

For a fixed-bandwidth population KDE, the local error is a kernel average of `delta`, not generally `delta(z)`. The current implementation also chooses bandwidth from weighted moments and Kish effective sample size, reflects observations at both boundaries, multiplies by a merger-rate weight, and renormalizes on a finite grid. Those operations add terms or discretization effects not represented by the elementary pointwise formula.

A second correction is essential downstream: the perturbation of the dark-siren log likelihood includes both the event-integral terms and the selection-normalization term `-n log beta(H0)`. A formula containing only an event-kernel average is incomplete for the likelihood used by this project.

## 1. Population ratio derivation (unsmoothed target)

Let `f(x)` be the normalized parent host density and let the true inclusion probability be `pi(x)>0`. The observed-catalogue density is proportional to `f(x) pi(x)`. Suppose the fitted inclusion probability is
```
pi_hat(x) = pi(x) [1 + delta(x)]
```
with `|delta(x)|` uniformly small on the target support and `1+delta(x)>0`.

The population Hájek/IPW density is
```
f_hat(x) = { f(x) pi(x)/pi_hat(x) } /
           { integral f(u) pi(u)/pi_hat(u) du }.
```
Since `pi/pi_hat = (1+delta)^(-1) = 1-delta+O(delta^2)`, and writing
`mu_delta = integral f(u) delta(u) du`,
```
f_hat(x) = f(x) [1-delta(x)] /
           [1-mu_delta] + O(delta^2)
         = f(x) {1-delta(x)+mu_delta} + O(delta^2).
```
Therefore
```
[f_hat(x)-f(x)]/f(x)
    = -{delta(x)-E_f[delta(X)]} + O(delta^2).
```

This is a population-level first-order result. The centering term removes a common multiplicative error in the fitted inclusion probabilities; only relative, covariate-dependent error changes the normalized density shape at this order. Positivity and uniformly small *relative* error are needed. Merely having small absolute error in `pi_hat-pi` is not sufficient when `pi` is close to zero.

## 2. Fixed-bandwidth KDE: the perturbation is kernel-averaged

Let `K_h(z-x)` be a fixed-bandwidth kernel and define the parent-smoothed density
```
g_h(z) = integral f(x) K_h(z-x) dx.
```
The corresponding normalized population IPW KDE is
```
g_hat_h(z) =
  { integral f(x) [1-delta(x)] K_h(z-x) dx } /
  { 1-E_f[delta(X)] } + O(delta^2).
```
Define the kernel-local average
```
delta_h(z) =
  integral f(x) delta(x) K_h(z-x) dx / g_h(z).
```
Then
```
[g_hat_h(z)-g_h(z)]/g_h(z)
  = -[delta_h(z)-E_f delta] + O(delta^2).
```

Thus the pointwise expression in Section 1 is recovered only in a suitable small-bandwidth/localization limit, with smooth `f` and `delta`, away from boundary complications. At finite bandwidth, the KDE responds to a local average of the inclusion-model error. Reflection changes the effective boundary kernel; the same ratio derivation applies with the actual reflected kernel operator, but `delta_h(z)` must be defined using that operator rather than an unreflected Gaussian kernel.

## 3. What the repository estimator adds

The current `estimate_host_density` implementation in `src/dark_siren_robustness/catalogue_density.py`:

1. forms inverse-probability weights `1/pi_hat_i`;
2. rescales them by their maximum for arithmetic stability;
3. computes Kish ESS and weighted mean/variance;
4. sets the default bandwidth to weighted SD times `ESS^(-1/5)` (or uses an explicit bandwidth factor);
5. reflects the observed redshifts at both grid boundaries;
6. evaluates a weighted Gaussian KDE;
7. multiplies by the merger-rate weight and normalizes by trapezoidal integration on the finite grid.

The max-rescaling does not change the normalized weighted KDE, but it does not control relative weight concentration or statistical instability.

The fixed-bandwidth formula above isolates the direct effect of `delta`. In the actual default estimator, changing `pi_hat` also changes the weighted moments and ESS, and therefore can change the bandwidth. A full first-order derivative of the implemented estimator must include the derivative of this data-adaptive bandwidth, the reflected-kernel boundary operator, merger weighting, and finite-grid normalization. It is not correct to label the simple pointwise formula as the complete first-order expansion of the code path.

The usual KDE bias order `O(h^2)` and variance order `O(1/(N h))` are generic one-dimensional, smooth-density, fixed-bandwidth heuristics under standard regularity conditions. They are not automatically the complete error rates here. In particular:
- `N` in the KDE sampling term refers to the catalogue sample size/effective weighted sample under the relevant assumptions, not automatically the validation size `n_val`;
- finite validation creates an additional random perturbation through `pi_hat`, whose covariance depends on the validation design and fitted model;
- estimated bandwidth and weight concentration alter the simple textbook variance expression;
- reflection reduces some boundary bias but does not eliminate all boundary or finite-grid error.

Consequently, do not write the KDE variance as `O(1/(n_val h))` unless a separate derivation establishes that rate for the chosen estimator and sampling design.

## 4. Merger-rate weighting and normalization

Let `m(z)>0` be the merger-rate weight. If the pre-merger host density is `f(z)`, the normalized target density entering the event likelihood is
```
p(z) = f(z)m(z) / integral f(u)m(u)du.
```
Under the unsmoothed first-order perturbation from Section 1, the perturbation after merger weighting and its own normalization is
```
[p_hat(z)-p(z)]/p(z)
  = -[delta(z)-E_p delta] + O(delta^2),
```
where `E_p` is expectation under the merger-weighted density, not under the unweighted parent density `f`.

For a smoothed density, replace `delta(z)` with the appropriate local kernel average before merger weighting, and then recenter under the normalized merger-weighted target. This distinction matters because the repository applies the merger-rate weight after KDE construction.

## 5. First-order propagation into the dark-siren likelihood

For event `j`, define
```
q_j(z;H) = Normal(d_j | c z/H, sigma_j)
I_j(H) = integral p(z) q_j(z;H) dz.
```
Let the GW detection probability at redshift `z` be `s(z;H)`, and define
```
beta(H) = integral p(z) s(z;H) dz.
```
The fixed-detected-count log likelihood has the form
```
ell(H) = sum_j log I_j(H) - n log beta(H).
```
Write the first-order relative perturbation of the normalized host density as
`p_hat(z)=p(z)[1-epsilon(z)]+O(epsilon^2)`, with `E_p epsilon=0`. Then
```
Delta log I_j(H) = - E_{p_j,H}[epsilon(Z)] + O(epsilon^2),
```
where `p_j,H(z) = p(z)q_j(z;H)/I_j(H)`, and
```
Delta log beta(H) = - E_{p_sel,H}[epsilon(Z)] + O(epsilon^2),
```
where `p_sel,H(z) = p(z)s(z;H)/beta(H)`.

Therefore
```
Delta ell(H) =
  - sum_j E_{p_j,H}[epsilon(Z)]
  + n E_{p_sel,H}[epsilon(Z)]
  + O(epsilon^2).
```

This expression shows that the H0 likelihood responds to the covariance/overlap of the host-density error with both the event kernels and the GW-selection kernel. Omitting the `+ n E_sel epsilon` term would ignore the change in the selection normalizer and would not describe the project's fixed-count likelihood.

This is a perturbation of the log likelihood at fixed `H`; translating it into a posterior-mean or credible-interval shift requires another approximation involving posterior curvature and prior/boundary effects. The full simulation remains necessary.

## 6. Finite-validation randomness is not removed by the expansion

Conditional on a realized validation sample, `delta(x)` is a realized function and the expressions above describe its first-order effect. Across repeated validation samples, `delta` is random. Its covariance depends on:
- the validation sampling design (e.g. SRSWOR versus independent Bernoulli sampling);
- the finite parent frame and any finite-population correction;
- the inclusion-model estimator and its constraints;
- separation, near-separation, rank deficiency, and boundary solutions;
- model misspecification and support gaps.

The conditional expansion does not itself provide a variance estimator for `delta`, nor does it propagate that variance through KDE and H0 inference. Use design-respecting repeated-sampling experiments or a fully specified joint model for that purpose. Do not insert an FPC or a Wald covariance formula without deriving it for the actual validation design and estimator.

## 7. Recommended manuscript wording

> For intuition, consider a population Hájek reconstruction with fitted inclusion probability `pi_hat(z)=pi(z)[1+delta(z)]`. Under positivity and uniformly small relative error, the normalized unsmoothed reconstruction has first-order relative distortion `-[delta(z)-E_f delta]`. With a fixed-bandwidth KDE, the pointwise error is replaced by a kernel-local average of `delta`; the actual implementation further uses reflection, data-adaptive bandwidth selection, merger-rate weighting, and finite-grid normalization. The expansion is conditional on the fitted inclusion function and does not quantify finite-validation sampling variability or the nonlinear effect on posterior summaries. In the fixed-detected-count dark-siren likelihood, first-order perturbations enter both the event integrals and the GW-selection normalization. We therefore use this calculation only as an interpretive guide and evaluate operating characteristics with end-to-end simulation.

## 8. Scope and claim guardrails

- This derivation does not establish that estimated IPW improves H0 inference.
- It does not validate a bootstrap, FPC, Bayesian nuisance-parameter integration, or a universal ESS threshold.
- It does not prove the missing-host mixture is redundant; that depends on the precise generative target and likelihood.
- It does not establish real-survey or LVK external validity.
- The numerical optimizer's finite constrained solution remains distinct from existence of an unconstrained logistic MLE under separation.
- The contribution remains a synthetic application/benchmark connecting established two-phase sampling, IPW, KDE, and dark-siren inference components.

## Source-code basis

The implementation details above were checked against:
- `src/dark_siren_robustness/catalogue_density.py` (`estimate_host_density`);
- the documented fixed-count likelihood specification and prior project audits.

A future code-level derivation should additionally inspect the exact production likelihood function and test the finite-grid perturbation numerically against small controlled perturbations of fitted inclusion probabilities.
