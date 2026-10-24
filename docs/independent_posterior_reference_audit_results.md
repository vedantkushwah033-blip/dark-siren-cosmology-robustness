# Independent posterior reference audit

**Status:** Passed after fixing a prior-boundary integration defect.  
**Workflow:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40067663955  
**Artifact:** `independent-reference-audit-results` (42-row CSV)

## Why the audit changed

The initial comparison exposed a numerical bookkeeping problem in the posterior implementation. The prior density was set to zero outside the prior interval, but trapezoidal integration across the discontinuous transition from a nonzero boundary value to an exterior zero introduced artificial half-bin area. The posterior summary then integrated that representation across the entire computational grid. In some multimodal cases, this materially distorted tail quantiles.

The implementation now normalizes the posterior only over the prior-supported grid and summarizes only that support. Regression tests cover this boundary behavior. The earlier numerical table in this report was superseded; the values below are from the corrected implementation.

## Method and scope

The reference calculation integrates the weighted Gaussian-mixture likelihood adaptively over the uniform prior interval [50, 90] km s^-1 Mpc^-1 and uses root-finding to obtain posterior quantiles. This is independent of the grid-based CDF calculation, but intentionally uses the same likelihood and prior.

Seven representative configurations × two deterministic host realizations = 14 observations. Each observation was summarized on 1,001-, 4,001-, and 16,001-point computational grids, yielding 42 comparisons. All errors below are absolute differences from adaptive-quadrature quantiles, in km s^-1 Mpc^-1.

| Grid points | Summary | Median absolute error | P95 absolute error | Maximum absolute error |
|---:|---|---:|---:|---:|
| 1,001 | Posterior median | 0.00009186 | 0.00087759 | 0.00108800 |
| 1,001 | Lower 68% endpoint | 0.00016820 | 0.00085554 | 0.00120046 |
| 1,001 | Upper 68% endpoint | 0.00019448 | 0.00095555 | 0.00124826 |
| 1,001 | Lower 95% endpoint | 0.00037893 | 0.00157012 | 0.00198855 |
| 1,001 | Upper 95% endpoint | 0.00020666 | 0.00191249 | 0.00196482 |
| 4,001 | Posterior median | 0.00000737 | 0.00005510 | 0.00007536 |
| 4,001 | Lower 68% endpoint | 0.00001288 | 0.00004843 | 0.00007053 |
| 4,001 | Upper 68% endpoint | 0.00000834 | 0.00007588 | 0.00007968 |
| 4,001 | Lower 95% endpoint | 0.00003692 | 0.00010070 | 0.00014941 |
| 4,001 | Upper 95% endpoint | 0.00001619 | 0.00009775 | 0.00012269 |
| 16,001 | Posterior median | 0.00000038 | 0.00000243 | 0.00000358 |
| 16,001 | Lower 68% endpoint | 0.00000087 | 0.00000269 | 0.00000285 |
| 16,001 | Upper 68% endpoint | 0.00000056 | 0.00000363 | 0.00000458 |
| 16,001 | Lower 95% endpoint | 0.00000268 | 0.00000533 | 0.00000633 |
| 16,001 | Upper 95% endpoint | 0.00000068 | 0.00000581 | 0.00000710 |

## Interpretation

- The corrected grid summaries agree closely with adaptive quadrature in these selected cases, and the discrepancy decreases with grid refinement.
- The 95% endpoints show larger errors than the median, as expected for tail quantiles, but the largest error in this 14-observation set at 4,001 points is about 0.00015 km s^-1 Mpc^-1.
- These results establish numerical agreement for the tested likelihood and cases only. They do not validate the physical model, realistic survey selection, catalogue completeness assumptions, or any real-data H0 inference.

## Reproducibility

```bash
python scripts/run_independent_reference_audit.py --output-dir results/independent-reference-audit
pytest -q tests/test_independent_reference_audit.py tests/test_hardened_inference.py
```

The GitHub Actions workflow validates all 42 comparisons and uploads the CSV artifact. The workflow now reruns when the hardened inference implementation changes.
