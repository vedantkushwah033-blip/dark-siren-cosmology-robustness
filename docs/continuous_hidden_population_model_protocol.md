# Continuous hidden-population model: implementation contract

**Status:** first production implementation on research branch; not yet independently run in CI.  
**Date:** 28 November 2026.

## Scientific model

For each independent event, define (u_i=\log d_i). Conditional on host redshift (z_i) and Hubble parameter (H_0),

\[
u_i\mid z_i,H_0\sim\mathcal N\left(\log\frac{cz_i}{H_0},\sigma_u^2\right).
\]

The listed component has redshifts (z_j) and normalized weights (q_j). The hidden component has fraction (f\in[0,1]) and redshift distribution

\[
z\mid\mu_h,s_h\sim\mathcal N(\mu_h,s_h^2)\text{ truncated to }[z_{\min},z_{\max}].
\]

The event density is

\[
p(u_i\mid\theta)=(1-f)\sum_jq_j\,\varphi\left(u_i;\log(cz_j/H_0),\sigma_u\right)
+f\int_{z_{\min}}^{z_{\max}}
\varphi\left(u_i;\log(cz/H_0),\sigma_u\right)
\,p_{TN}(z\mid\mu_h,s_h)\,dz.
\]

The event log likelihood is the sum of log event densities. This is a density in (u=\log d), not in (d). If the data are represented as distances rather than log-distances, the Jacobian (1/d) must be included when converting the density; for fixed observed distances it is parameter-independent and does not change the posterior over these parameters.

## Numerical method

The hidden-redshift integral is evaluated by Gauss-Legendre quadrature over the fixed truncation interval. Mixture calculations are performed in log space. The four-dimensional posterior is evaluated on a Cartesian grid with a uniform prior over the supplied grid bounds and normalized using tensor-product trapezoidal weights.

The implementation lives in `src/dark_siren_robustness/hidden_population.py`. It is intentionally separate from the earlier Gaussian-distance positive-truncation likelihood.

## Prior and grid contract for the initial reproduction

The initial reproduction should use the analogue's stated bounds, where applicable:
- (H_0\in[50,90]);
- (f\in[0,1]);
- (mu_h\in[0.03,0.095]);
- (s_h) fixed at 0.012 for the first parameter-recovery reproduction, then varied or estimated only in a separately declared experiment;
- hidden redshift support ([0.005,0.10]);
- listed redshifts ((0.02,0.04,0.06)) with equal weights if reproducing the supplied analogue, not the existing discrete-host paper experiment;
- (sigma_u=0.1) for the supplied analogue baseline.

Do not use a coarse four-dimensional grid as final evidence. First verify quadrature and grid convergence, then increase grid density or use a validated sampler. Report posterior boundary mass and sensitivity to grid bounds.

## Critical identifiability caution

The hidden-redshift truncation bounds are fixed. Therefore the transformation that scales (H_0,\mu_h,s_h) together does not generally preserve the same truncated-normal family. Exact scale invariance at (f=1) must not be claimed without a proof under a transformed-support model. Assess local weak directions numerically using standardized parameters and do not equate a Fisher matrix with global identifiability.

## Validation gates

1. Truncated-normal density integrates to one.
2. Gauss-Legendre quadrature converges as order increases.
3. (f=0) exactly reduces to the listed-host likelihood.
4. (f=1) exactly reduces to the hidden-component likelihood.
5. Identical-component limits are fraction invariant.
6. Posterior tensor-grid normalization is one.
7. Marginal credible intervals are ordered and lie on their grids.
8. Parameter recovery uses data generated from the same model and fixed seeds.
9. Oracle/fixed-wrong/plug-in/integrated comparisons use paired simulated datasets.
10. Finite-sample coverage is reported with Monte Carlo uncertainty and is kept separate from prior-predictive SBC.

Passing unit tests establishes implementation identities, not astrophysical validity or inferential calibration.
