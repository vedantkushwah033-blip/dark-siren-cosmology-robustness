# Paired host-weight sensitivity sweep

**Status:** completed; workflow passed 9,000 synthetic observations evaluated under two weightings each (18,000 posterior fits).  
**Workflow:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40066226236  
**Scope:** synthetic low-redshift model only.

## Design

This experiment addresses a weakness of the prior weight stress test: it evaluates correct and reversed weights on the **same observed distance for each replicate**. The true $H_0$ is fixed at 70 km s$^{-1}$ Mpc$^{-1}$, with host redshifts $z_1=0.02$ and $z_2 \in \{0.025,0.04,0.06\}$. Data are generated with host weights (70%, 30%), then inferred using either the correct (70%, 30%) or reversed (30%, 70%) weights.

The distance-noise scale is varied across 3%, 10%, and 20% of $cz_{ref}/H_0$, using $z_{ref}=0.03$. Each of the 9 design cells has 1,000 paired replicates. Metrics include posterior-median RMSE, 68% interval coverage and the paired change in absolute median error when weights are reversed.

## Results

| $z_2$ | Noise fraction | Correct-weight RMSE | Reversed-weight RMSE | Correct 68% coverage | Reversed 68% coverage | Mean paired change in absolute error (reversed − correct) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.025 | 3% | 7.261 | 10.677 | 72.4% | 56.2% | +3.899 |
| 0.025 | 10% | 6.921 | 7.899 | 81.6% | 73.9% | +1.105 |
| 0.025 | 20% | 6.444 | 6.850 | 93.2% | 90.6% | +0.514 |
| 0.040 | 3% | 2.758 | 2.758 | 70.2% | 70.2% | approximately 0 |
| 0.040 | 10% | 6.448 | 6.369 | 76.4% | 78.7% | −0.048 |
| 0.040 | 20% | 5.229 | 5.896 | 94.6% | 94.6% | +0.740 |
| 0.060 | 3% | 2.890 | 2.890 | 66.7% | 66.7% | approximately 0 |
| 0.060 | 10% | 6.294 | 6.294 | 72.4% | 72.4% | approximately 0 |
| 0.060 | 20% | 6.526 | 6.468 | 85.6% | 86.8% | −0.039 |

RMSE and error metrics are in km s$^{-1}$ Mpc$^{-1}$.

## Interpretation

- **Weight misspecification matters most in some close-host configurations.** When $z_2=0.025$ and noise is 3%, reversing weights increased median RMSE from 7.26 to 10.68 and reduced empirical 68% coverage from 72.4% to 56.2%.
- **There is no universal monotonic effect.** For larger redshift separation, the two mixture components can lead to distinct distance-supported regions, and in many observations the likelihood itself dominates the prior weights. In several cells, the posterior median was effectively unchanged by reversing the weights.
- **Higher noise does not automatically make weighting irrelevant.** The effect depends jointly on host separation, noise scale and the observed distance. Some cells show only small differences, while others retain a measurable contrast.
- These are cell-wise point estimates from 1,000 replicates each. Small differences should not be over-interpreted without uncertainty intervals and replication across seeds.

## Limitations

The model uses only two fixed host redshifts, a low-redshift distance approximation, a fixed injected $H_0$, and Gaussian distance noise. It does not include realistic GW detection selection, galaxy-catalogue construction, cosmological distance integration or real event data. These results should not be presented as a measured effect in LVK data.

## Reproducibility

- Script: `scripts/run_paired_weight_sweep.py`
- Tests: `tests/test_paired_weight_sweep.py`
- Workflow: `.github/workflows/paired-weight-sweep.yml`
- Artifact: `paired-weight-sweep-results`
