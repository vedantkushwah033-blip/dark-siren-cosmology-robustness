# Adaptive-quadrature validation of selection-posterior normalization

**Status:** the first three-version run completed successfully on 24 October 2026. This report describes a synthetic toy-model numerical audit, not a real gravitational-wave analysis.

## Method

For three fixed observed distances, the script compares the normalizing integral of the unnormalized selected likelihood
\[
\mathcal L_{\rm sel}(d\mid H_0)=
\frac{\sum_j w_j\,\mathcal N(d\mid cz_j/H_0,\sigma^2)}{\beta(H_0)}
\]
over the prior support \(H_0\in[50,90]\). It compares composite trapezoidal integration on a 4,001-point grid with SciPy adaptive quadrature. The test tolerance is a relative absolute difference below \(10^{-4}\).

The model uses redshifts 0.02 and 0.04, weights 0.7 and 0.3, \(\sigma=20\) Mpc, and the toy detection interval \(0<d<180\) Mpc.

## Results

| Observed distance (Mpc) | Grid integral | Adaptive quadrature | Estimated quadrature absolute error | Relative absolute difference |
|---:|---:|---:|---:|---:|
| 50 | 0.15279860323260594 | 0.15279858857095670 | 5.17e-14 | 9.60e-8 |
| 90 | 0.50719054090332390 | 0.50719061281032940 | 1.63e-12 | 1.42e-7 |
| 140 | 0.18155628366986507 | 0.18155621509578350 | 1.49e-14 | 3.78e-7 |

All three relative differences are far below the predeclared \(10^{-4}\) tolerance. The artifact CSV schema and row count were also checked by the script.

## CI evidence and version caveat

The matrix workflow passed on Python 3.10, 3.11, and 3.12 in run [40073404985](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40073404985), at commit `50bdcdaf24b5d49729a2eae87016f708559b6c01`. The subsequent grid-validation-only adjustment is included in commit `7720a35ce2a41d6d2e08753a03764496c8051e4e); its corresponding run was still in progress at the time this report was prepared. Thus the numerical result is validated by the completed earlier workflow, while the newest revision's CI status must be checked separately.

## Interpretation limits

This checks numerical integration of the stated toy likelihood; it does not establish that the likelihood or detection function is physically realistic. It is also not a direct comparison of the normalized posterior density at every grid point, nor a calibration test. Those are separate validation targets. Small integration differences do not resolve the fixed-truth coverage issues near prior boundaries.
