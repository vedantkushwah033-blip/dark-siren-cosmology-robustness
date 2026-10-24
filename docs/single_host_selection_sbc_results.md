# Single-host selection-conditioned SBC: results

**Status:** completed synthetic control; workflow passed.  
**Workflow:** [Toy single-host selection SBC](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40072622191)  
**Commit tested:** `18436c037761209ddc4225abef98554c249b04e7`

## Design

This control uses one host at redshift 0.03, Gaussian distance uncertainty of 20 Mpc, the toy detection interval (0<d_{obs}<180) Mpc, a uniform (H_0) prior on [50, 90], and a 4,001-point posterior grid. It runs 2,000 prior-predictive replicates with seed 20261031 and compares the selection-normalized likelihood against a deliberate negative control that omits the selection normalization.

For prior-predictive SBC, (H_0) is sampled from the prior and the observation is generated conditional on detection for that (H_0). The parameter is not rejected or resampled based on whether a detection occurred; doing so would change the prior distribution.

## Results

| Diagnostic | Corrected | Uncorrected |
|---|---:|---:|
| Replicates | 2,000 | 2,000 |
| Mean posterior-CDF rank | 0.49834 | 0.47818 |
| Rank variance | 0.08586 | 0.08851 |
| KS statistic | 0.01746 | 0.05062 |
| KS p-value | 0.56946 | 0.0000680 |
| Decile chi-square | 9.79 | 34.68 |
| Decile chi-square p-value | 0.36775 | 0.0000678 |
| Minimum (eta(H_0)) over prior grid | 0.50248 | 0.50248 |
| Maximum (eta(H_0)) over prior grid | 0.99997 | 0.99997 |

The corrected decile counts were [217, 209, 181, 210, 180, 191, 187, 214, 216, 195]. The uncorrected counts were [256, 235, 170, 213, 173, 182, 193, 200, 198, 180].

## Interpretation

In this particular one-host synthetic model, the corrected ranks are compatible with uniformity under the reported KS and decile chi-square diagnostics, while the uncorrected negative control shows clear rank non-uniformity. This supports the claim that omitting the (H_0)-dependent selection normalization can impair prior-predictive calibration in this toy setup.

It does **not** prove that the likelihood captures realistic gravitational-wave detection selection, that every aspect of the code is correct, or that the same quantitative effect holds in a real analysis. A non-significant p-value is not proof of calibration; the rank summaries, design assumptions, and independent numerical checks should also be reviewed.

## Relation to fixed-truth coverage

These are complementary diagnostics, not contradictory results. SBC averages over parameters drawn from the prior and tests the joint prior-predictive inference procedure. The fixed-truth coverage map conditions on particular values of (H_0) and assesses repeated-sampling coverage at those points. A method can have satisfactory prior-predictive ranks while showing poor pointwise coverage near prior boundaries.
