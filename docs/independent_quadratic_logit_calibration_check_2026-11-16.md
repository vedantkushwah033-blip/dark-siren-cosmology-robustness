# Independent quadratic-logit calibration check — 2026-11-16

## Purpose

This note records a separate numerical reconstruction of the deterministic parent-measure centering and parent-average inclusion-rate calibration used by the seven-cell quadratic-logit density pilot.

This is a **calibration arithmetic check only**. It is not a run of the seven-gamma simulation, not a validation of the production generator/auditor, and not evidence about density-recovery performance.

## Independently reconstructed specification

The parent redshift density is proportional to

\[
g(z) \propto 1 + 2.2\exp\{-\tfrac12[(z-0.105)/0.025]^2\}
       + 1.6\exp\{-\tfrac12[(z-0.275)/0.040]^2\},
\]

on \([0.01,0.40]\). The standardized coordinate is

\[
x=(z-0.205)/0.195,
\]

and the inclusion model is

\[
\pi_\gamma(z)=\operatorname{logit}^{-1}
\{\alpha_\gamma - 3.9x + \gamma(x^2-c_g)\}.
\]

I independently evaluated the parent integrals on a 20,001-point equally spaced redshift grid using trapezoidal integration, with a scalar root solve for each nonzero-gamma intercept. The centering constant is

\[
c_g=E_g[x^2]=0.2845005579191241.
\]

The reference parent-average inclusion rate at \(\gamma=0\), \(\alpha_0=0\), and \(\beta=-3.9\) is

\[
r_0=E_g[\operatorname{logit}^{-1}(-3.9x)]
   =0.5101919962291073.
\]

For each gamma, the intercept was solved so that \(E_g[\pi_\gamma(z)]=r_0\).

## Reconstructed calibration values

| gamma | alpha_gamma | achieved parent-average inclusion |
|---:|---:|---:|
| -2.0 | -0.24047847028208044 | 0.5101919962291073 |
| -1.0 | -0.12472169347234588 | 0.5101919962291073 |
| -0.5 | -0.06297071445790560 | 0.5101919962291072 |
| 0.0 | 0.00000000000000000 | 0.5101919962291073 |
| 0.5 | 0.06312269530994806 | 0.5101919962291072 |
| 1.0 | 0.12530550566919094 | 0.5101919962291070 |
| 2.0 | 0.24244138342969593 | 0.5101919962291060 |

The reconstructed values agree with the constants documented in the cross-gamma auditor to the stated rounding precision.


## Analytic inclusion and weight-stress profile

Using the same parent-measure grid and calibrated intercepts, I also evaluated the range of the true inclusion probability and the second-moment implication for ideal inverse-probability weights. Let r0 be the common parent-average retention rate. Under the retained-catalogue redshift law proportional to g(z) pi_gamma(z), the population coefficient of variation of the ideal weight 1/pi_gamma is

\[
\operatorname{CV}_{\mathrm{catalogue}}(1/\pi_\gamma)
=\sqrt{r_0 E_g[1/\pi_\gamma]-1}.
\]

| gamma | min pi_gamma | max pi_gamma | E_g[1/pi_gamma] | catalogue-law CV of 1/pi_gamma |
|---:|---:|---:|---:|---:|
| -2.0 | 0.003790 | 0.902891 | 13.1461 | 2.3889 |
| -1.0 | 0.008661 | 0.955203 | 8.4510 | 1.8198 |
| -0.5 | 0.013116 | 0.970092 | 6.9865 | 1.6014 |
| 0.0 | 0.019840 | 0.980160 | 5.9077 | 1.4192 |
| 0.5 | 0.029912 | 0.986885 | 5.1124 | 1.2682 |
| 1.0 | 0.044822 | 0.991344 | 4.5265 | 1.1443 |
| 2.0 | 0.097281 | 0.996217 | 3.7850 | 0.9649 |

These are deterministic population-level quantities, not empirical finite-sample ESS values or observed pilot results. They show that the symmetric gamma grid is **not symmetric in practical positivity/weight stress** under the chosen negative linear coefficient: negative curvature creates substantially smaller minimum inclusion probabilities and more variable ideal catalogue weights than positive curvature. This is scientifically relevant when interpreting gamma contrasts. It also means the seven gamma cells vary both curvature direction and the severity of inverse-weight instability; analyses should report the actual inclusion/weight profile rather than describe the stress grid as symmetric in difficulty.

## Interpretation and limitations

1. Parent-measure centering is important: using the uniform-coordinate value or an arbitrary constant would change the meaning of gamma and could shift the linear component.
2. Intercept recalibration holds the *marginal parent inclusion rate* fixed across gamma. It does not make the inclusion functions equal, nor does it hold their variance, local slopes, minimum probabilities, or inverse-weight distributions fixed.
3. Consequently, comparisons across gamma are controlled for the overall retention rate but deliberately vary the redshift pattern of retention. That is the intended curvature stress.
4. The degree-one fitted model in Arms D/E remains misspecified when gamma is nonzero. Its posterior integrates coefficient uncertainty conditional on the degree-one model; it does not automatically integrate uncertainty about the omitted quadratic term.
5. The known-probability Arm C is a simulation oracle benchmark. It is not available to a real analyst without knowing the true inclusion mechanism.

## Execution status at note creation

- Independent deterministic centering/intercept arithmetic: **SUCCESS**.
- Exact-head CI for PR #55 head `6a7b79bd0ffc9edbdd50f8de169bda9029fdf87b`: **NOT VERIFIED**; no workflow run was returned for that commit.
- Seven-gamma workflow dispatch: **NOT RUN / NOT VERIFIED**.
- Seven gamma cells with 20 replicates each: **NOT RUN / NOT VERIFIED**.
- Cross-gamma raw-artifact audit: **NOT RUN / NOT VERIFIED**.
- Density-performance conclusions: **NOT EVALUABLE until the complete run and independent audit exist**.

The manual workflow is `.github/workflows/quadratic-logit-curvature-density-pilot.yml`. The workflow is explicitly density-only and must not be treated as authorising H0 propagation, bandwidth selection, or an estimator-superiority claim.

## Reproduction outline

The independent check used SciPy's `trapezoid`, `expit`, and `brentq` with the parent intensity and coordinate definition stated above. For each nonzero gamma, it solved

\[
\frac{\int_{0.01}^{0.40}g(z)
\operatorname{logit}^{-1}\{\alpha-3.9x(z)+\gamma[x(z)^2-c_g]\}\,dz}
{\int_{0.01}^{0.40}g(z)\,dz}
-r_0=0.
\]

This note is intentionally limited to the deterministic calibration layer. It does not replace the workflow's independent raw-input audits or its paired cross-gamma artifact audit.
