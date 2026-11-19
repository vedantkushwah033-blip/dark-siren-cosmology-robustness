# Full-parent prior-predictive SBC: 1,000-replicate extension

**Status:** preregistered computational validation extension; no results until the complete run and independent audit pass.

## Design
- 1,000 independent prior-predictive replicates.
- Draw each true H0 independently from Uniform(20, 300) km s^-1 Mpc^-1.
- Generate a fresh parent population of 1,200 galaxies and 30 detected events per replicate.
- Use the complete finite parent catalogue and the matching selected-event likelihood with its detection-efficiency normalization.
- Evaluate on 1,121 equally spaced H0 points.
- Use independent deterministic seed streams for truth draws and population/event generation.
- Independently regenerate all 1,000 prior draws, parent populations, event sequences, posterior ranks, quantiles, and summary statistics.
- Retain the complete rank histogram, KS and histogram diagnostics, posterior-median bias/RMSE, and 50/68/80/90/95% coverage.

## Why extend the 200-replicate pilot?
The 200-replicate experiment is exploratory. A 1,000-replicate run provides a more informative empirical rank histogram and materially tighter Monte Carlo uncertainty, although it still does not prove exact calibration.

## Interpretation guardrails
- This tests only the full-parent baseline likelihood.
- It does not test catalogue-inclusion correction or finite-validation estimation.
- P-values are descriptive diagnostics, not binary proof of calibration.
- It remains a synthetic low-redshift Euclidean model, not LVK/GWTC data or an observational H0 measurement.
