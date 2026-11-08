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
\theta_k\sim\mathcal N(0,2.5^2).
\]
The prior is defined on the fixed-domain polynomial coefficients and therefore does not change when the realized validation sample's mean or standard deviation changes. This is important: applying the same numerical prior to sample-standardized coefficients would induce a validation-sample-dependent prior on the inclusion curve.

The prior scale is provisional. Sensitivity analysis at 1.5, 2.5, and 5.0 is required before confirmatory posterior claims. The target is the validation-only modular posterior; GW data do not update \(\theta\).

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
