# Independent posterior reference audit

**Status:** Passed on 2026-10-24  
**Workflow:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40067150477  
**Artifact:** `independent-reference-audit-results` (42 row CSV)

## Purpose and scope

This audit checks numerical integration and quantile extraction in the simplified synthetic two-host dark-siren model. It compares the repository's grid-based posterior summaries with an independently calculated reference: adaptive quadrature of the weighted Gaussian-mixture likelihood over the uniform prior interval [50, 90] km s^-1 Mpc^-1, followed by root-finding for posterior quantiles.

This is a numerical implementation check only. It is not observational validation, does not establish realistic selection modelling, and does not imply a real-data H0 measurement.

## Design

- Seven representative configurations spanning injected H0 values 52, 70, and 88; two-host redshift choices; two noise levels; and correct or reversed host weights.
- Two deterministic host realizations per configuration, giving 14 observations.
- Each observation summarized on 1,001-, 4,001-, and 16,001-point grids: 42 grid summaries.
- Compared posterior median and central 68% and 95% interval endpoints.
- Adaptive quadrature's reported absolute error estimate was at most 6.77e-12 for the normalization integral in these cases.

## Absolute difference from adaptive-quadrature reference

All errors below are in km s^-1 Mpc^-1. “Median” is the median absolute error across the 14 observations; “P95” is the empirical 95th percentile; “Maximum” is the largest observed absolute error.

| Grid points | Summary | Median absolute error | P95 absolute error | Maximum absolute error |
|---:|---|---:|---:|---:|
| 1,001 | Posterior median | 0.012031 | 0.044006 | 0.053647 |
| 1,001 | Lower 68% endpoint | 0.014129 | 0.036180 | 0.043934 |
| 1,001 | Upper 68% endpoint | 0.015383 | 0.040617 | 0.041980 |
| 1,001 | Lower 95% endpoint | 0.009998 | 0.047010 | 0.048893 |
| 1,001 | Upper 95% endpoint | 0.025375 | 0.048226 | 0.048648 |
| 4,001 | Posterior median | 0.003033 | 0.011201 | 0.013814 |
| 4,001 | Lower 68% endpoint | 0.003546 | 0.009071 | 0.011006 |
| 4,001 | Upper 68% endpoint | 0.003860 | 0.010178 | 0.010508 |
| 4,001 | Lower 95% endpoint | 0.002580 | 0.011660 | 0.012230 |
| 4,001 | Upper 95% endpoint | 0.006338 | 0.012027 | 0.012159 |
| 16,001 | Posterior median | 0.000760 | 0.002813 | 0.003475 |
| 16,001 | Lower 68% endpoint | 0.000887 | 0.002270 | 0.002753 |
| 16,001 | Upper 68% endpoint | 0.000966 | 0.002546 | 0.002628 |
| 16,001 | Lower 95% endpoint | 0.000640 | 0.002910 | 0.003058 |
| 16,001 | Upper 95% endpoint | 0.001586 | 0.003005 | 0.003040 |

## Interpretation

1. **Grid refinement reduces numerical error in these test cases.** The median absolute error in the posterior median fell from about 0.0120 at 1,001 points to 0.00303 at 4,001 points and 0.000760 at 16,001 points.
2. **Tail endpoints deserve attention.** The 95% upper endpoint has larger typical error than the median, especially on the coarsest grid. This is expected to be an important diagnostic for mixture-shaped or boundary-adjacent posteriors.
3. **The 4,001-point grid is numerically close in this tested set, not universally guaranteed.** Its maximum observed error among these 14 cases was about 0.0139 for the median and 0.0123 for the interval endpoints. This small deterministic set cannot rule out larger errors in other configurations.
4. **The reference is independent in integration and quantile extraction, not in the scientific likelihood.** Both methods intentionally use the same weighted Gaussian-mixture likelihood and prior interval, so this audit does not test whether that model is physically adequate.

## Reproducibility

Run locally after installing the package:

```bash
python scripts/run_independent_reference_audit.py --output-dir results/independent-reference-audit
pytest -q tests/test_independent_reference_audit.py
```

The workflow validates that all 42 comparisons are written and uploads the CSV artifact. The main test suite for the commit containing this audit also passed on Python 3.10, 3.11, and 3.12.

## Next recommended validation

Add analytic limiting-case tests and a wider randomized numerical stress test, including posteriors near prior boundaries and strongly overlapping or separated host likelihood components. Only after numerical behavior is characterized should the project use these routines inside a more realistic distance–redshift and selection model.
