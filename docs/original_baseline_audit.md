# Original baseline calibration: audit notes

## Source report

The original repository's `docs/baseline_calibration.md` reports the following for 1,000 complete-catalogue, unknown-host trials:

| Quantity | Reported value |
|---|---:|
| 68% coverage | 94.4% |
| 95% coverage | 99.9% |
| Mean 68% interval width | 22.372 km s^-1 Mpc^-1 |
| Mean 95% interval width | 35.820 km s^-1 Mpc^-1 |
| Median bias | +1.986 km s^-1 Mpc^-1 |
| Median RMSE | 6.190 km s^-1 Mpc^-1 |
| Median-pull standard deviation | 0.670 |

The same report gives a known-host control with 68% coverage of 74.9%, 95% coverage of 98.6%, and a mean 68% interval width of 17.045 km s^-1 Mpc^-1.

The report also records a separate 2,000-run result of 94.35% coverage at nominal 68%.

## Interpretation supported by the report

These values indicate that the current synthetic baseline has conservative credible intervals relative to the nominal levels, especially for the unknown-host mixture. The unknown-host intervals are wider on average than the known-host control. The reported median-pull standard deviation below one is also consistent with standardized residuals being less dispersed than a unit-scale reference under the report's definition.

This pattern is consistent with host-mixture uncertainty broadening the posterior, but it does not by itself prove that this is the only cause. Positive bias, distance likelihood assumptions, prior/grid effects, and the data-generating process should also be checked.

## Important reproducibility distinction

The 94.4% value comes from a documented 1,000-run baseline report; 94.35% comes from a separately reported 2,000-run calibration. Neither number should be described as independently reproduced by the new robustness repository until the manual workflow finishes, its artifact is inspected, and its source commit and dependencies are recorded.

The new single-host diagnostic is a separate simplified model and must not be presented as a reproduction of either original result.

## Next audit questions

1. Does the pinned source commit match the implementation used to generate the baseline report?
2. Does the workflow reproduce the 1,000-run baseline within Monte Carlo variation?
3. Is the nominal 68% interval exactly the equal-tail interval documented in the original inference code?
4. How much of the conservative coverage remains in a known-host control using the same data generator and likelihood?
5. Are the reported pull diagnostics consistent with the saved per-run values?
