# Mathematical framework for finite-validation catalogue inclusion

**Status:** derivation and design contract; numerical results remain pending.  
**Scope:** synthetic parent-frame experiments with low-redshift Euclidean GW distances. No real-survey or detector claim is made.

## 1. Parent population and catalogue selection

Let a parent-frame galaxy have redshift (Z\sim f(z)) on ([z_{\min},z_{\max}]), with predictors (X=x(Z)). Conditional on its predictors, catalogue retention is

\[
R\mid Z=z\sim\operatorname{Bernoulli}(\pi_0(z)),\qquad
\pi_0(z)=\operatorname{expit}\{a_0+b_0x+\gamma_0(x^2-c)\}.
\]

The parent-frame mean inclusion probability is
\[
\bar\pi_0=\int f(z)\pi_0(z)\,dz.
\]
The redshift density among retained galaxies is therefore
\[
f_{\rm cat}(z)=\frac{f(z)\pi_0(z)}{\bar\pi_0}.
\]
This is catalogue selection. It is not gravitational-wave detection selection.

For the benchmark, (x=(z-z_{\min})/(z_{\max}-z_{\min})), and
\[
c=E_f[x(Z)^2]
=\frac{\int x(z)^2\lambda(z)\,dz}{\int\lambda(z)\,dz},
\]
where the declared parent intensity is (\lambda(z)). The implementation evaluates these integrals on the same fine grid used to sample the parent redshift distribution. Thus the curvature regressor is centered under the declared parent law, not assigned an arbitrary constant.

## 2. What inverse-inclusion weighting targets

For a fixed positive candidate inclusion function (\tilde\pi(z)), the large-sample weighted catalogue measure is proportional to
\[
f_{\rm cat}(z)\frac{1}{\tilde\pi(z)}
\propto f(z)\frac{\pi_0(z)}{\tilde\pi(z)}.
\]
After normalization, the reconstructed density is
\[
f_{\rm IPW}(z;\tilde\pi)=
\frac{f(z)\pi_0(z)/\tilde\pi(z)}
{\int f(u)\pi_0(u)/\tilde\pi(u)\,du}.
\]
It equals (f(z)) if (\tilde\pi(z)=\pi_0(z)) almost everywhere on the parent support and positivity holds. Under misspecification, the multiplicative distortion is (\pi_0(z)/\tilde\pi(z)), followed by global renormalization. In particular, a model that underpredicts inclusion in a redshift region overweights that region; a model that overpredicts inclusion underweights it.

For any integrable test function (g), the limiting weighted expectation is
\[
E_{\rm IPW}[g]=
\frac{E_f[g(Z)\pi_0(Z)/\tilde\pi(Z)]}
{E_f[\pi_0(Z)/\tilde\pi(Z)]}.
\]
This identity is the primary analytic benchmark for density-reconstruction bias. It also shows why low fitted probabilities can create large variance even when the normalized weighted expectation is asymptotically correct.

## 3. Linear logistic pseudo-true parameter under quadratic truth

Let (h(z)=(1,x(z))^T) be the linear fitted design and
\[
\pi_\beta(z)=\operatorname{expit}(h(z)^T\beta).
\]
The validation sample is a simple random sample without replacement from the full finite parent frame, before conditioning on (R). Under the superpopulation formulation, the expected Bernoulli negative log likelihood is minimized at the pseudo-true coefficient (\beta^\dagger) satisfying
\[
E_f\left[h(Z)\{\pi_0(Z)-\pi_{\beta^\dagger}(Z)\}\right]=0.
\]
This is not generally the coefficient obtained by dropping the quadratic term from the true coefficient vector: the fitted intercept and slope adjust jointly to the nonlinear probability curve. The pseudo-true model matches two weighted score moments, not the entire function (\pi_0(z)). The resulting IPW density bias follows by substituting (\tilde\pi=\pi_{\beta^\dagger}\) into the preceding ratio.

For a finite validation sample, the MLE solves the empirical score equation
\[
\sum_{i\in V}h(z_i)\{R_i-\pi_\beta(z_i)\}=0.
\]
The current implementation estimates this MLE only. It does not yet implement Firth, Bayesian, profile-likelihood, sandwich, or design-based finite-population uncertainty methods.

## 4. Validation sampling and disjoint analysis catalogue

The primary design draws (V\subset\{1,\ldots,N\}) uniformly without replacement, independently of retention outcomes. The validation records contain ((z_i,R_i)), and validation units are excluded from the inference catalogue:
\[
C=\{i\notin V:R_i=1\}.
\]
The validation and analysis IDs are disjoint by construction. The fitted inclusion model receives only validation redshifts and retention labels. True inclusion probabilities are retained in scoring-only artifacts and are never supplied to the fitter.

The complement (V^c) is held out for probability scoring. Scores on (V^c) are distinct from the IPW diagnostics calculated on (C). The benchmark records the full parent partition so an independent auditor can reconstruct both.

## 5. Finite-population correction

A simple random validation sample without replacement has a finite-population correction for design-based sample means: conditional on a fixed finite frame, the variance of a sample mean contains the factor (1-n_V/N). It is not valid to multiply an ordinary model-based logistic inverse-Hessian covariance by this factor automatically. For logistic regression, a design-based covariance requires the finite-frame score contributions and a design-based covariance for their sample total/mean, followed by a linearization using the information matrix. That calculation is not implemented in the current fitter. The reported inverse observed-information covariance is therefore a model-based local approximation, not an FPC-adjusted design-based covariance.

## 6. Catalogue inclusion versus GW detection selection

Let (f_{\rm host}(z)) be the parent host-redshift density after applying any declared merger-rate weight and normalizing. In the toy distance model,
\[
d\mid z,H_0\sim N(c_{\rm light}z/H_0,\sigma^2),\qquad
D=\{d\ge d_{\rm th}\}.
\]
The GW detection efficiency is
\[
\beta(H_0)=\int f_{\rm host}(z)
P(D\mid z,H_0)\,dz
=\int f_{\rm host}(z)
\Phi\!\left(\frac{c_{\rm light}z/H_0-d_{\rm th}}{\sigma}\right)dz.
\]
For a fixed number (N_D) of independent detected events, the selected event-data likelihood is
\[
L_{\rm sel}(H_0)=
\prod_{j=1}^{N_D}
\frac{\int f_{\rm host}(z)p(d_j\mid z,H_0)\,dz}
{\beta(H_0)}.
\]
Equivalently, its log likelihood is the sum of unconditional event log likelihoods minus (N_D\log\beta(H_0)). The host density in both numerator and denominator must be the same parent host density. Catalogue IPW is one possible upstream estimator of the galaxy component of that density; it does not replace or duplicate (\beta(H_0)).

This is a fixed-detected-count conditional likelihood. It is not an event-rate likelihood and does not include realistic detector response, measurement pipelines, or correlated population nuisance parameters.

## 7. Statistical interpretation and limits

- Correct inclusion probabilities plus positivity are required for standard IPW identification.
- Estimated probabilities introduce additional uncertainty; replacing them by a plug-in estimate generally understates uncertainty unless the estimation step is propagated.
- Stabilization or truncation changes the estimating procedure and generally introduces bias; it must be reported as a separate sensitivity arm.
- A validation sample selected conditional on retention does not identify (P(R=1\mid Z)) without additional design information. The primary design avoids that selection.
- Probability calibration, parent-density reconstruction, host-KDE error, posterior calibration, and credible-interval coverage are separate outcomes.
- Numerical grid convergence is not evidence of statistical calibration.
- The current finite-validation stage does not yet connect estimated probabilities to the host KDE or H0 posterior. No downstream claim is made.
