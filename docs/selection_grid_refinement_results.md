# Posterior grid-refinement audit: results

**Workflow:** [Toy selection posterior grid refinement](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40072015994)  
**Commit tested:** `d1da2e382f15c5c1857815f4dc4ad5328d55deee)  
**Status:** workflow succeeded; focused tests passed (2 passed). This is a deterministic numerical audit of a synthetic model.

## Setup

The workflow compared posterior summaries at 1,001, 4,001, and 8,001 grid points for observed distances of 70, 100, and 140 Mpc. The prior support was \(H_0\in[50,90]\), with host redshifts 0.02 and 0.04, weights 0.7 and 0.3, Gaussian distance uncertainty 20 Mpc, and upper detection threshold 180 Mpc. The grid points spanned 25 to 125, with posterior support restricted to the prior.

The reference for differences was the 8,001-point result. This audit compares posterior summaries across grids; unlike the separate adaptive-quadrature audit, it does not directly compare the normalization integral to quadrature.

## Results

| Observed distance (Mpc) | Grid points | Median | Absolute median difference vs 8,001 | Maximum absolute difference across 68%/95% interval endpoints |
|---:|---:|---:|---:|---:|
| 70 | 1,001 | 75.9423621 | 2.7233e-5 | 2.9763e-4 |
| 70 | 4,001 | 75.9423882 | 1.1407e-6 | 1.7523e-5 |
| 70 | 8,001 | 75.9423893 | 0 | 0 |
| 100 | 1,001 | 63.7521021 | 1.2151e-4 | 1.1501e-4 |
| 100 | 4,001 | 63.7519844 | 3.8248e-6 | 4.2098e-6 |
| 100 | 8,001 | 63.7519806 | 0 | 0 |
| 140 | 1,001 | 72.8150900 | 3.5928e-4 | 7.6775e-5 |
| 140 | 4,001 | 72.8154290 | 2.0299e-5 | 9.9332e-6 |
| 140 | 8,001 | 72.8154493 | 0 | 0 |

## Interpretation

- The 4,001-point medians differ from the 8,001-point medians by at most approximately \(2.03\times10^{-5}\) across these three test observations.
- The maximum 4,001-vs-8,001 difference among the reported 68% and 95% interval endpoints is approximately \(1.75\times10^{-5}\).
- The 1,001-point grid shows larger discrepancies, especially for the 140 Mpc median.
- These are deterministic comparisons for three selected observations, not a universal error bound over all possible distances or model parameters.

Together with the independent adaptive-quadrature audit, this supports the use of 4,001 points for these particular toy-model checks. It does not establish physical realism, calibration, or acceptable frequentist coverage.

## Reproduction

```bash
python scripts/run_selection_grid_refinement.py --output-dir results/selection-grid-refinement
python -m pytest tests/test_selection_grid_refinement.py -q
```
