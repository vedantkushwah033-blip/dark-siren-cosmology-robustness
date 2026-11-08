# Finite-validation inclusion posterior sampler: implementation note

**Status:** experimental sampler implemented; downstream host-density/H0 propagation is not implemented by this change.

## Target distribution

For validation records \((z_j,R_j)\), use the standardized polynomial logistic model
\[
\operatorname{logit}\pi(z_j;\theta)=\sum_{k=0}^{d}\theta_k x_j^k,
\qquad x_j=(z_j-\bar z_V)/s_V.
\]
The sampler uses the proper independent Gaussian prior
\[
\theta_k\sim\mathcal N(0,2.5^2)
\]
for each coefficient. The prior scale is a provisional modelling choice and must receive sensitivity analysis (at minimum 1.5, 2.5, and 5.0) before confirmatory posterior claims. Because the basis is standardized separately for each validation sample, this prior is defined on standardized-basis coefficients.

The target is the validation-only modular posterior
\[
p(\theta\mid V)\propto
\left[\prod_j \pi(z_j;\theta)^{R_j}
(1-\pi(z_j;\theta))^{1-R_j}\right]
\exp\{-\tfrac12\theta^T\theta/2.5^2\}.
\]
GW data do not update \(\theta\). This is a cut/modular posterior, not automatically a full joint posterior.

## Algorithm

- Find the posterior mode by BFGS using the analytic gradient.
- Compute the posterior negative-Hessian at the mode. Its inverse is used only to precondition proposals; it is **not** treated as the posterior distribution.
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
