# Full-parent 1,000-replicate prior-predictive H0 SBC readout

**Workflow:** [44446994455](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44446994455)  
**Status:** COMPLETE; independent raw-seed audit PASS.  
**Scope:** synthetic full-parent baseline only; no catalogue-inclusion correction or observational inference.

## Frozen experiment

- 1,000 prior-predictive replicates.
- H0 drawn independently from Uniform(20,300) km s^-1 Mpc^-1.
- 1,200 parent galaxies and 30 detected events per replicate.
- 1,121-point H0 grid.
- Full finite-parent selected-event likelihood with matching one-sided detection-efficiency normalization.
- Independent seed streams and raw-seed reconstruction.

## Results

| Diagnostic | Result |
|---|---:|
| Mean posterior rank | 0.497308 |
| SD posterior rank | 0.284898 |
| 10-bin rank counts | 88, 109, 105, 101, 108, 103, 93, 98, 94, 101 |
| Rank histogram chi-square | 4.14 |
| Rank histogram p-value | 0.90196 |
| Rank KS statistic | 0.02081 |
| Rank KS p-value | 0.77127 |
| Mean posterior-median bias | +0.1965 km s^-1 Mpc^-1 |
| Bias Monte Carlo SE | 0.4688 km s^-1 Mpc^-1 |
| RMSE | 14.8189 km s^-1 Mpc^-1 |
| 50% interval coverage | 51.0% |
| 68% interval coverage | 69.9% |
| 80% interval coverage | 81.1% |
| 90% interval coverage | 90.0% |
| 95% interval coverage | 94.4% |

## Independent audit

The independent auditor completed successfully:

- 1,000 replicates regenerated from raw seed streams.
- 1,000 parent populations regenerated.
- 1,000 detected-event sequences regenerated.
- 1,000 posterior-rank and quantile rows reconstructed.
- Summary metrics independently recomputed.
- Three aggregate manifest entries verified.

The workflow's final status is successful and the independent audit reports PASS.

## Interpretation

The rank distribution and interval coverage are consistent with a well-calibrated **synthetic full-parent baseline** under this prior-predictive design. This is an important reference point for debugging the separate catalogue-inclusion pipeline.

It does not validate estimated-inclusion MAP/IPW or coefficient-integrated IPW. It also does not resolve the fixed-truth four-event pilot's positive median bias: that experiment has a different event count, H0 prior, and inference design. Those experiments must not be compared as if they were a single controlled contrast.

The next useful diagnostic is a controlled event-count bridge under the same generator, H0 prior, and likelihood family, followed by the paired five-arm 30-event repeated-seed experiment. No universal estimator-superiority claim follows from this baseline.
