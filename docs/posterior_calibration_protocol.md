# Posterior-level H₀ calibration: prespecified protocol

**Scope:** synthetic statistical validation only. This experiment is not a real gravitational-wave analysis and does not model the full LVK event likelihood, detector selection, catalogue incompleteness, redshift uncertainty, or hierarchical population inference.

## Question

In a deliberately small two-host low-redshift model, how do host-weight misspecification, prior support, and prior measure affect posterior median error, interval coverage, and interval width across repeated synthetic observations?

## Fixed design

- Injected H₀ values: 55, 65, 75, and 85 km s⁻¹ Mpc⁻¹.
- Fractional distance-noise settings: 3%, 10%, and 20%.
- Two possible host redshifts: z₁ = 0.02 and z₂ = 0.04.
- Data-generating host probabilities: (0.7, 0.3).
- Inference weight arms: correct (0.7, 0.3) and reversed (0.3, 0.7).
- Prior supports: narrow [50, 90] and broad [20, 220] km s⁻¹ Mpc⁻¹.
- Prior measures: uniform in H₀ and log-uniform in H₀ (density proportional to 1/H₀).
- Default 200 replicates per truth × noise cell: 12 cells, 2,400 raw observations, and 96 cell-by-inference-arm summaries.
- Fixed seed: 20261030 for generation; 20261031 is the base seed for deterministic bootstrap analyses.
- Posterior grid: 2,401 equally spaced points from 10 to 250 km s⁻¹ Mpc⁻¹, restricted to the selected prior support.

For each truth/noise cell, distance uncertainty is set to noise_fraction × c × z₂ / H₀_truth. A host is sampled from the fixed data-generating probabilities, and a Gaussian distance observation is drawn around c × z_host / H₀_truth; nonpositive simulated distances are redrawn. The exact same raw observation is evaluated under all eight inference arms. The inference likelihood is a two-component Gaussian mixture in observed distance.

## Primary outputs

For every truth/noise/support/measure/weighting cell, report:
- mean posterior-median error (signed bias);
- median absolute posterior-median error and RMSE;
- empirical coverage of central 68% and 95% posterior intervals;
- mean 68% and 95% interval widths.

Pooled summaries must retain equal representation of the 12 prespecified truth/noise cells. They are descriptive averages across this design, not estimates for a population of real events. Report truth/noise-specific uncertainty as well as pooled uncertainty; never let a pooled metric replace the cell-level view.

## Validation gates

1. **Raw provenance audit:** independently replay the declared seeded random-number generator and data-generating process, checking row order, truth/noise/replicate keys, host draws, noise scales, distances, and the raw-file hash.
2. **Posterior reconstruction audit:** independently recompute all posterior medians and 68%/95% interval endpoints from each raw observation and each inference arm; verify every raw record, summary row, and input hash.
3. **Uncertainty analysis:** run deterministic stratified bootstrap resampling within each fixed truth/noise cell. Do not treat observations as exchangeable across different truth/noise cells.
4. **Paired contrasts:** compare correct versus reversed host weights using the same observation in both arms; resample replicate IDs within each cell and report differences with uncertainty intervals.
5. **Reproducibility:** retain raw CSVs, summaries, metadata, audit outputs, uncertainty reports, paired contrasts, seeds, and workflow provenance as a downloadable artifact.

A passing unit test is not a substitute for the independent raw-data and posterior audits. If any gate fails, the experiment is not considered validated until the failure is investigated and corrected.

## Interpretation guardrails

- The narrow prior contains all four injected truths and is strongly informative relative to the broad support. Better errors under the narrow support do not imply it is generally preferable.
- Uniform and log-uniform priors encode different prior measures; neither should be called universally superior based on this experiment.
- Incorrect weights may change posterior modes and tails. The direction and magnitude of effects are conditional on this specific two-host model.
- Empirical coverage is a repeated-sampling property of this simulator, not the posterior probability that a particular real-event interval contains H₀.
- Bootstrap intervals quantify Monte Carlo uncertainty conditional on the simulation design. They do not quantify model discrepancy or astrophysical uncertainty.
- Do not describe synthetic values as observed or measured H₀, and do not claim astrophysical realism or a full dark-siren analysis from these controls.

## Reproduction

From the repository root in an environment with the project development dependencies installed:

```bash
python scripts/run_posterior_calibration.py --replicates 200 --seed 20261030 --output-dir results/posterior-calibration
python scripts/audit_posterior_calibration_generation.py results/posterior-calibration
python scripts/audit_posterior_calibration.py results/posterior-calibration
python scripts/analyze_posterior_calibration_uncertainty.py results/posterior-calibration --seed 20261031 --bootstrap-resamples 4000
python scripts/analyze_paired_host_weight_contrasts.py results/posterior-calibration --seed 20261031 --bootstrap-resamples 4000
```

The independent posterior audit changes metadata from PENDING to PASS; uncertainty and paired-contrast scripts refuse to run unless that audit and its hashes pass. A successful workflow establishes internal reproducibility under this declared toy model, not observational validity.
