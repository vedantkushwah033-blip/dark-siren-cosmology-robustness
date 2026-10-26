# Selected-likelihood calibration grid-sensitivity protocol

## Purpose

Quantify whether numerical posterior-grid resolution materially changes prior-predictive PIT diagnostics in the repository's finite synthetic selected-likelihood model. This is a numerical sensitivity study, not a test of astrophysical realism.

## Prespecified design

The CI screening experiment uses:

- Grid resolutions: 1,201, 2,401, and 4,801 equally spaced points on the same H0 support [40, 100].
- Replicates per resolution: 300.
- Random seed: 20261025.
- Identical prior, two-host mixture, host weights, Gaussian distance noise, and open detection interval across resolutions.
- Independent raw-record replay and PIT/summary reconstruction for every resolution.

The generator uses a single seeded RNG stream, and grid evaluation occurs after each replicate's truth, host, and observed distance are drawn. With the same seed and generation settings, each grid run should therefore reuse the same simulated truth/host/distance records. The posterior grid only changes the numerical evaluation of the posterior CDF and resulting PIT.

## Outputs

For each resolution, the runner writes:

- `pit_replicates.csv`: replicate-level truth, host, observation, rejection count, and PIT.
- `pit_summary.json`: model/protocol metadata, PIT mean and Monte Carlo standard error, variance, descriptive KS statistic and p-value, decile counts, and Wilson intervals.
- `audit_report.json`: independent replay/reconstruction outcome.

The top-level `grid_sensitivity.csv` contains one row per resolution with grid size, replicate count, seed, PIT summary diagnostics, and audit status. CI archives these files under `artifacts/calibration_grid_sensitivity/`.

## Interpretation rules

1. Compare PIT means relative to their Monte Carlo uncertainty, but do not treat overlapping intervals as a formal equivalence test.
2. Compare the PIT distribution and decile frequencies across resolutions; a mean near 0.5 alone is insufficient.
3. Because the generated records should be paired across resolutions, differences in PIT summaries are driven primarily by grid discretization rather than independent data realizations.
4. Treat 300 replicates as a CI screening run. If the screening indicates meaningful grid sensitivity, or for any publication-facing claim, repeat with at least 1,000 replicates and additional seeds.
5. A non-significant uniformity test is not proof of calibration. The synthetic model's internal consistency does not establish correctness for real gravitational-wave selection, event likelihoods, or galaxy catalogues.

## Follow-up

A stronger study should quantify per-replicate PIT changes between grids, predefine acceptable tolerances, repeat over multiple seeds and selection/noise regimes, and distinguish numerical discretization error from finite Monte Carlo variability. The current experiment deliberately isolates grid resolution while keeping the data-generation model fixed.
