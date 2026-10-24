# Cosmological distance-model robustness protocol

## Status

Implementation and CI validation are in progress as of 2026-10-24. The initial paired synthetic experiment and focused tests have been added; report no result as final until the workflow succeeds and its artifact is inspected. See `scripts/run_cosmological_model_robustness.py`, `src/dark_siren_robustness/cosmological_inference.py`, and `tests/test_cosmological_inference.py`.

## Research question

Do the estimated effects of synthetic galaxy-catalogue incompleteness and host-population weighting persist when the luminosity-distance relation is changed from the low-redshift approximation to a flat Lambda-CDM relation?

This is a model-robustness experiment for a synthetic simulator, not an analysis of real LVK observations.

## Models

### A. Low-redshift baseline

Use the existing approximation

\[
d_L(z;H_0) \approx \frac{cz}{H_0}.
\]

Reproduce the existing baseline before comparing models. Preserve the original implementation and its pinned commit.

### B. Flat Lambda-CDM

For fixed matter-density parameter \(\Omega_m\), use

\[
d_L(z;H_0,\Omega_m) =
\frac{(1+z)c}{H_0}
\int_0^z \frac{dz'}{
\sqrt{\Omega_m(1+z')^3 + 1 - \Omega_m}
}.
\]

Initially fix \(\Omega_m=0.3\) and state this assumption explicitly. The initial extension must not imply that uncertainty in \(\Omega_m\) has been marginalized over.

### C. Deliberate distance-model mismatch

Generate synthetic distances under flat Lambda-CDM, then infer with the low-redshift approximation. Compare with the matched Lambda-CDM generation/inference case. This isolates the impact of this specified form of model misspecification; it does not establish that all other modelling assumptions are correct.

## Design

1. **Baseline reproduction:** rerun the original low-redshift setup and verify it against previously validated artifacts before interpreting any change.
2. **Matched Lambda-CDM:** generate and infer with the same flat Lambda-CDM relation.
3. **Mismatch control:** generate with flat Lambda-CDM and infer with the low-redshift relation.
4. **Catalogue mechanism comparison:** for each distance-model setting, compare the existing complete-catalogue and 70%-completeness scenarios, including random and inverse-redshift-weighted synthetic missingness where supported by the current experiment design.
5. **Host-weighting comparison:** retain the existing uniform and luminosity-weighted host-generation/inference conditions where those scenarios are defined and validated.
6. Use matched seeds and identical simulated populations for paired contrasts when scientifically valid. Record precisely which random draws are shared and which are independent.
7. Start with 500 replicates per scenario as a diagnostic run. Increase to 1,000 or more if Monte Carlo uncertainty obscures the contrasts. Do not treat a small p-value or a bootstrap interval excluding zero as a substitute for effect size, model checks, or multiplicity-aware interpretation.

## Required implementation tests

- At \(z=0\), the computed luminosity distance approaches zero.
- For small positive redshift, the flat Lambda-CDM distance approaches \(cz/H_0\) to the expected low-redshift order.
- For fixed \(z\) and \(\Omega_m\), luminosity distance scales inversely with \(H_0\).
- Distance is finite, positive for positive redshift, and monotonic over the experiment's supported redshift range.
- Numerical integration agrees with an independent high-accuracy reference at selected test points.
- The matched-generation and inference paths call the intended distance-model implementation; no accidental mixture of models is permitted.
- Posterior normalization and credible-interval bounds remain valid under both models.
- Output validation checks replicate counts, finite posterior summaries, interval ordering, model labels, and the complete configuration needed to reproduce each run.

## Outcomes to report

For each scenario, report at minimum:

- empirical coverage of the nominal 68% and 95% credible intervals, with Monte Carlo uncertainty;
- bias and RMSE of the posterior median relative to the generating \(H_0\);
- mean or median credible-interval width;
- number of valid replicates and any excluded or failed runs;
- paired contrasts between pre-specified scenarios, with uncertainty intervals;
- exact random seeds, parameter values, redshift support, noise model, completeness mechanism, host-generation weights, and inference weights.

Pre-specify the primary contrasts before inspecting the outputs. Treat additional contrasts as exploratory and state whether multiplicity is addressed.

## Interpretation limits

- This remains a synthetic simulation, not a real-data measurement or validation against an observational catalogue.
- The inverse-redshift selection mechanism is a designed stress test, not a calibrated survey selection function.
- Fixed-size completeness masks condition on the retained count and do not model uncertainty in the number of missing galaxies.
- Fixing \(\Omega_m\) does not capture its uncertainty or its degeneracy with \(H_0\).
- A successful test suite establishes software properties covered by those tests; it does not establish that the scientific model is realistic.
- Do not claim improved accuracy, calibrated coverage, or robust cosmological conclusions until validated run artifacts are inspected.

## Reproducibility checklist

- [ ] Record the exact source commit and dependency versions.
- [ ] Run existing tests and the baseline reproduction first.
- [ ] Save machine-readable configuration and per-replicate outputs.
- [ ] Validate all outputs before analysis.
- [ ] Preserve raw results and analysis scripts as workflow artifacts.
- [ ] Report Monte Carlo uncertainty and failed-run handling.
- [ ] Separate pre-specified confirmatory contrasts from exploratory analyses.
