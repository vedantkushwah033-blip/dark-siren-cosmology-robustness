# Randomized posterior numerical stress audit

**Status:** Passed with the prior-boundary correction.  
**Workflow:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40067659869  
**Design:** 120 seeded synthetic cases, two grid sizes per case, 240 posterior summaries total.  
**Artifact:** `randomized-numerical-stress-results` (CSV).

## Why this test matters

The audit deliberately varies injected H0, three candidate host redshifts, the data-generating host weights, the inference weights, and the distance-noise scale. It includes injected values near the prior boundaries. For each case, grid summaries are compared with the adaptive-quadrature reference calculation.

This audit also served as a regression check after fixing a numerical defect: trapezoidal integration had previously counted artificial half-bin mass where the posterior density jumped from a nonzero value at the prior boundary to zero outside the prior. The corrected code normalizes and summarizes on the prior-supported grid.

## Absolute quantile error

Errors are in km s^-1 Mpc^-1. P95 is the empirical 95th percentile across 120 cases; maximum is the largest error in this fixed seeded sample.

| Grid points | Quantile | Median absolute error | P95 absolute error | Maximum absolute error |
|---:|---|---:|---:|---:|
| 1,001 | Posterior median | 0.00009780 | 0.00101453 | 0.00342452 |
| 1,001 | Lower 68% endpoint | 0.00025711 | 0.00161912 | 0.00487366 |
| 1,001 | Upper 68% endpoint | 0.00019247 | 0.00169365 | 0.00304955 |
| 1,001 | Lower 95% endpoint | 0.00045373 | 0.00309862 | 0.02469987 |
| 1,001 | Upper 95% endpoint | 0.00020786 | 0.00295202 | 0.00670155 |
| 4,001 | Posterior median | 0.00000657 | 0.00005859 | 0.00023066 |
| 4,001 | Lower 68% endpoint | 0.00001505 | 0.00009360 | 0.00034513 |
| 4,001 | Upper 68% endpoint | 0.00001328 | 0.00010089 | 0.00019588 |
| 4,001 | Lower 95% endpoint | 0.00002524 | 0.00020608 | 0.00145402 |
| 4,001 | Upper 95% endpoint | 0.00001247 | 0.00020073 | 0.00039009 |

## Interpretation and limitations

1. **Grid refinement substantially improves agreement in this seeded test.** The median absolute error for the posterior median drops from approximately 9.78e-5 at 1,001 points to 6.57e-6 at 4,001 points.
2. **Tail quantiles remain the most sensitive.** The largest error is in the lower 95% endpoint, especially at 1,001 points. At 4,001 points, its maximum absolute error is approximately 0.00145 km s^-1 Mpc^-1 in this sample.
3. **This is not a universal error guarantee.** It is one seeded, 120-case stress test. Additional adversarial cases, narrower likelihood components, more complicated mixtures, and boundary configurations should be tested.
4. **Numerical agreement is not astrophysical validation.** The test uses a simplified low-redshift distance relation and a Gaussian-mixture likelihood. It does not establish realism of selection effects or catalogue incompleteness and is not a real LVK measurement.

## Reproduce

```bash
python scripts/run_randomized_numerical_stress.py --cases 120 --seed 20261024 --output-dir results/randomized-numerical-stress
pytest -q tests/test_randomized_numerical_stress.py
```
