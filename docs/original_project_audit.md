# Audit of the Version-1 dark-siren simulation

**Audit status:** source-level review of selected files and documentation on 2026-10-22. This is not an independent rerun of the workflows or a complete review of every file.

Reviewed files include `README.md`, `src/bayesian_inference.py`, `src/simulations_v2.py`, `src/experiments.py`, `src/completeness.py`, `src/metrics_v2.py`, selected tests, and `docs/baseline_calibration.md` / `docs/selection_conditioning.md`.

## Findings

### 1. The reported 94.35% is empirical coverage, not evidence of nominal 68% calibration

The baseline report describes 94.4% empirical coverage for central 68% credible intervals across 1,000 trials, with a later reported 94.35% result across 2,000 trials. The report also records 99.9% coverage for nominal 95% intervals in the 1,000-run control.

This means the intervals are substantially conservative under this simulated setup. It does **not** mean the model is calibrated at the nominal 68% level. The result is important and should be retained, not tuned away. The numerical report is documented in the source repository; I have not independently retrieved the workflow artifacts or rerun the complete experiments in this audit.

### 2. The distance model is intentionally low-redshift and simplified

The code uses the approximation (d_L \approx cz/H_0). This is suitable as a transparent first-stage toy model, but it is not a full cosmological distance-redshift relation. Results should remain explicitly conditional on this approximation.

### 3. The likelihood is a weighted mixture over candidate redshifts

In `src/bayesian_inference.py`, each candidate redshift contributes a Gaussian distance likelihood centred at (cz/H_0), and the contributions are combined using normalized host weights. The prior is uniform over the supplied grid range, by default 50–90 km/s/Mpc.

This is a useful controlled baseline, but it is not a full LVK population analysis or a complete event-level hierarchical likelihood.

### 4. Catalogue incompleteness and GW detection selection are not jointly modelled

The original project explicitly conditions on the existence of a simulated GW event and defers an H0-dependent GW-detection normalization. Catalogue masks are separate toy mechanisms. This is an acknowledged limitation, not a hidden feature; any later selection-aware model must derive the corresponding normalization consistently rather than insert an arbitrary correction.

### 5. The missing-host component is an oracle validation control

When enabled, the missing-host term uses the simulated underlying population, which is known to the generator. This can test whether a correction works under controlled conditions, but it is not directly available in a real survey. Results from this oracle condition must be reported separately from catalogue-only inference.

### 6. The simulation design has useful controls, but further validation is needed

The source includes repeatable seeds, multiple catalogue masks, multiple host-generation and inference-weighting choices, coverage/bias/RMSE/width/pull summaries, and tests for repeatability and mask sizes. However, the selected tests do not by themselves establish:
- simulation-based calibration under a clearly defined prior-generative model;
- correct nominal coverage;
- numerical accuracy against an independent cosmology implementation;
- full validation of all likelihood edge cases;
- realism of any survey selection function.

The current `posterior_h0` function also has limited input validation for finite/strictly increasing H0 grids, finite positive distance uncertainty, and finite non-negative redshifts/weights. These should be added and tested before extending the likelihood.

## Recommended next steps

1. Reproduce the small known-host control with a fixed seed and inspect the posterior numerically.
2. Add tests for input validation, posterior normalization, and comparison with an independent reference for the low-redshift model.
3. Investigate why nominal 68% intervals have around 94% empirical coverage. Check interval construction, posterior grid, priors/boundaries, and the repeated-sampling definition before changing any method.
4. Add simulation-based calibration under a precisely specified generative model, and report its limitations separately from misspecification experiments.
5. Only after these checks, design a more realistic selection-aware extension.

## Scope of this audit

This is a source-level audit of selected files, not a claim that the whole project has been independently reproduced. No reported scientific result was changed or re-labelled as independently verified.
