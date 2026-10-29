# Cross-condition independent audit of the confirmatory selection experiment

**Date:** 2026-10-29  
**Scope:** independent reconstruction of exported posterior summaries and arm-level aggregates from all three condition-specific workflow artifacts. This checks saved computational outputs; it does not rerun the simulation generator from raw random seeds or validate a physical gravitational-wave model.

## Provenance and integrity

All three artifacts came from workflow run [#40819870231](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40819870231) on source commit `9aadda0ffde113ebb2584aa5425424dbe154dd2b`. GitHub artifact ZIP SHA-256 digests:

- Low inclusion floor (artifact ID 12873663060): `1919b6ffbf3aae489c8b08225a129d30e244d4f748126b77d937e4cd92bc0712`
- High inclusion floor (artifact ID 12873737617): `f88015f19f11940f7e5228e6b13d2b1c4ad657cbc3140df1991873dcff54e03a`
- Broad transition (artifact ID 12873982241): `3ce2d085cc6669e694a953190fe995b2fd75363efaad3164a55907ef5cf2a79b`

For each ZIP, the CSV and JSON byte lengths and SHA-256 digests match the included manifest. Each artifact has 1,500 rows, 250 replicates per six arms, no duplicate `(arm, replicate)` keys, one seed and one true H0 per replicate across arms, and exactly 30 detected events per row.

## Independent posterior and summary reconstruction

For all 4,500 rows across the three conditions, I parsed the saved 81-point posterior-density grid, integrated it on the declared H0 grid, and reconstructed posterior mean, median, and the 16th/84th percentiles. The largest absolute discrepancy in any condition was:

- Density integral versus one: **4.45e-16**
- Posterior mean: **5.69e-14 km s^-1 Mpc^-1**
- Posterior median and either 68% interval endpoint: **1.43e-14 km s^-1 Mpc^-1**

I also recomputed bias, RMSE, empirical 68% coverage, density integrated squared error, and mean effective-sample-size fraction from the raw table. The largest difference from the stored aggregate summaries across these fields and conditions was **3.56e-15**. This supports internal consistency of the exported posterior grids, row summaries and aggregate summaries to floating-point precision.

## Results across the three declared synthetic conditions

True H0 is 70 km s^-1 Mpc^-1. Values below are recomputed from raw replicate records; coverage is empirical 68% interval coverage.

| Catalogue selection setting (floor, midpoint, width) | Arm | Mean bias | RMSE | Coverage | Mean density ISE | Mean ESS fraction |
|---|---|---:|---:|---:|---:|---:|
| 0.05, 0.22, 0.05 | Catalogue only | -10.617 | 11.040 | 4.8% | 1.0742 | 1.000 |
| 0.05, 0.22, 0.05 | Correct weighting | +1.185 | 4.095 | 65.6% | 0.0428 | 0.560 |
| 0.05, 0.22, 0.05 | Shifted transition | -0.488 | 3.846 | 68.4% | 0.1774 | 0.745 |
| 0.05, 0.22, 0.05 | Shallow slope | -1.990 | 4.201 | 60.0% | 0.1423 | 0.772 |
| 0.25, 0.22, 0.05 | Catalogue only | -4.717 | 5.898 | 38.0% | 0.4533 | 1.000 |
| 0.25, 0.22, 0.05 | Correct weighting | +0.938 | 4.056 | 69.6% | 0.0290 | 0.828 |
| 0.25, 0.22, 0.05 | Shifted transition | +0.950 | 4.052 | 70.8% | 0.0666 | 0.882 |
| 0.25, 0.22, 0.05 | Shallow slope | -0.094 | 3.873 | 70.0% | 0.0522 | 0.896 |
| 0.10, 0.22, 0.08 | Catalogue only | -4.994 | 6.153 | 33.6% | 0.4724 | 1.000 |
| 0.10, 0.22, 0.08 | Correct weighting | +1.304 | 4.154 | 70.8% | 0.0329 | 0.808 |
| 0.10, 0.22, 0.08 | Shifted transition | +0.549 | 3.931 | 71.6% | 0.0550 | 0.878 |
| 0.10, 0.22, 0.08 | Shallow slope | -0.545 | 3.864 | 66.0% | 0.0692 | 0.907 |

### Paired absolute-error contrast

For each replicate, `D = |error_catalogue_only| - |error_correct_weighting|`; positive values favor correct weighting on absolute posterior-mean error. Using a fresh percentile bootstrap (20,000 resamples; seed 20261029) gives:

| Condition | Mean paired D | Fresh bootstrap 95% interval | Fraction of pairs favoring correct weighting |
|---|---:|---:|---:|
| Low floor | 7.4358 | [6.8050, 8.0406] | 90.4% |
| High floor | 2.2462 | [1.7217, 2.7476] | 74.4% |
| Broad transition | 2.3537 | [1.7767, 2.9286] | 72.4% |

These intervals are independently recomputed and therefore need not exactly equal the workflow summary's percentile intervals, which may use a different bootstrap seed or implementation convention. They are consistent in scale and direction.

## Interpretation and limitations

- In these configurations, catalogue-only inference has materially larger bias and RMSE than correct inclusion weighting.
- Correct weighting is not the best arm on every metric in every condition. Shifted-transition or shallow-slope alternatives sometimes have lower RMSE, which is a reminder that bias, variance, density error, and interval coverage need not rank methods identically.
- The low-floor condition has a mean ESS fraction of 0.560 under correct weighting, signalling substantial weight concentration and a practical variance/positivity concern. This should remain visible in the interpretation, not hidden behind the bias reduction.
- The three condition-specific jobs have overlapping seed ranges (150 seeds shared by low/high floor, 150 by high-floor/broad-transition, and 50 by low-floor/broad-transition, as recorded in the release-gate notes). Treat each within-condition contrast as paired, but do **not** describe the three conditions as independent replications.
- The audit reconstructs the exported posterior grid and summary table. It does not regenerate the original latent populations, rerun likelihood calculations from source, establish robustness outside these toy settings, or validate CHIMERA/LVK/GWTC inference.

**Outcome:** PASS for artifact file integrity, posterior-grid summary reconstruction, row pairing structure, and aggregate-summary consistency. Scientific conclusion remains limited to the declared synthetic experiment.
