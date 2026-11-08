# Finite-validation inclusion posterior sampler: implementation note

**Status:** experimental sampler implemented; downstream host-density/H0 propagation is not implemented by this change.

## Fixed-coordinate target distribution

For validation records \((z_j,R_j)\), use the fixed domain coordinate
\[
u_j=(z_j-z_{\min})/(z_{\max}-z_{\min}),
\qquad
\operatorname{logit}\pi(z_j;\theta)=\sum_{k=0}^{d}\theta_k u_j^k.
\]
For this synthetic project, the fixed domain is \([z_{\min},z_{\max}]=[0.01,0.40]\). The sampler uses the proper independent Gaussian prior
\[
\theta_k\sim\mathcal N(0,1.0^2).
\]
The prior is defined on the fixed-domain polynomial coefficients and therefore does not change when the realized validation sample's mean or standard deviation changes. This is important: applying the same numerical prior to sample-standardized coefficients would induce a validation-sample-dependent prior on the inclusion curve.

The prior scale 1.0 is a candidate reference informed by the separate prior-predictive audit, not a frozen confirmatory choice. Sensitivity analysis at 0.75 and 1.5 is required before confirmatory posterior claims. The target is the validation-only modular posterior; GW data do not update \(\theta\).

## Algorithm and diagnostics

- Find the posterior mode by BFGS using the analytic gradient.
- Compute the posterior negative Hessian at the mode. Its inverse is used only to precondition proposals; it is **not** treated as the posterior distribution.
- Run four independently seeded random-walk Metropolis chains with multivariate Gaussian proposals based on that local curvature.
- Adapt the proposal scale during warmup toward acceptance 0.234; freeze it during retained sampling.
- Return all chain draws, per-chain post-warmup acceptance rates, classical split-Rhat, and an initial-positive-sequence approximate ESS.
- Mark rank-deficient/ill-conditioned validation designs and one-class validation samples NON_EVALUABLE. Mark failed MAP/Hessian checks FAILED. Retain draws but mark MCMC diagnostic failures NON_EVALUABLE.

Default run: 2,000 warmup and 2,000 retained draws per chain, with gates max split-Rhat 1.05 and min approximate ESS 100 per coefficient. These are minimum software gates, not a guarantee of tail accuracy or independent draws. Draw-level Monte Carlo error must be checked by longer runs and independent batches before downstream marginalization.

## Implementation boundary

New module: `src/dark_siren_robustness/inclusion_posterior.py`.

This implementation does not yet:
- build an IPW density for each posterior draw;
- calculate event-mixture likelihoods or GW detection efficiency per draw;
- marginalize the selected likelihood over coefficient draws;
- produce (H_0) posteriors, SBC, or fixed-truth coverage.

The existing MLE fitter and its inverse observed-information covariance remain separate. The sampler does not use that covariance as posterior draws. The provisional prior must be varied and justified before a confirmatory H0 study.


## Relationship to the weighted importance sampler

The production candidate for downstream marginalization is `bayesian_inclusion.py`, which uses a fixed-basis multivariate Student-t importance proposal and retains normalized/log importance weights. This file's `inclusion_posterior.py` random-walk Metropolis implementation is an **independent equal-weight reference sampler**, not the proposal distribution used by the H0 propagation pipeline. It exists to cross-check posterior coefficient means and selected-likelihood integrals on selected synthetic cases. Its chain draws must not be passed to the weighted importance-sampling API as if they carried importance weights.

The dedicated test suite compares both samplers under the same fixed-domain Gaussian prior and checks a small toy likelihood integral. Agreement on those tests is a numerical cross-check, not proof of global sampler accuracy, tail adequacy, or posterior calibration.


## Single-replicate prior-sensitivity check (exploratory, not confirmatory)

I ran the independent MCMC reference sampler on one actual validation subset from the audited disjoint density-recovery artifact:

- Source: density-recovery workflow run #42476772981, artifact #13616821671.
- Selected cell: replicate 0, curvature \(\gamma=+1.2\), target retention 0.25, validation size 300.
- Validation labels: 70 retained and 230 not retained.
- Held-out parent units: 3,700 (analysis-catalogue and non-retained units combined).
- MCMC: four chains, 2,000 warmup and 2,000 retained draws per chain, fixed-basis Gaussian prior.

The DGP's exact fixed-basis coefficients, reconstructed from the stored true inclusion probabilities, were approximately \((-0.9841,-1.0,+1.2)\).

| Degree | Prior SD | Posterior coefficient means | Minimum ESS | Maximum split-Rhat | Held-out Brier | Held-out log loss |
|---:|---:|---|---:|---:|---:|---:|
| 1 | 0.75 | (-1.445, 0.553) | 955 | 1.006 | 0.18948 | 0.56697 |
| 1 | 1.00 | (-1.568, 0.740) | 918 | 1.004 | 0.18990 | 0.56835 |
| 1 | 1.50 | (-1.696, 0.952) | 861 | 1.006 | 0.19049 | 0.57033 |
| 2 | 0.75 | (-1.412, 0.095, 0.604) | 845 | 1.009 | 0.18950 | 0.56709 |
| 2 | 1.00 | (-1.504, 0.108, 0.760) | 714 | 1.005 | 0.19001 | 0.56869 |
| 2 | 1.50 | (-1.546, 0.026, 0.973) | 608 | 1.006 | 0.19042 | 0.56995 |

All six MCMC runs passed the declared Rhat/ESS software gates. For the correctly specified degree-2 model, the true coefficient vector lies inside the componentwise 95% intervals for all three prior scales. For degree 1, the model is misspecified and there is no true degree-1 coefficient vector to cover.

**Interpretation restriction:** this is one validation sample from one simulation replicate. It shows modest prior sensitivity in posterior coefficient means and small changes in held-out probability scores; it does not establish prior robustness, posterior calibration, or a general predictive advantage. The prior-scale comparison must be repeated across independent validation replicates and checked against the weighted importance sampler before any confirmatory H0 use.
