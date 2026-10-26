# Calibration grid-sensitivity screening results — 2026-10-26

## Provenance and status

The values in this report come from the archived Python 3.11 artifact for [workflow run 40247354524](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40247354524), commit `f4046f892376cd4090c2d1c19a416eedf7a0687b`. The Python 3.11 job completed successfully, including the 300-replicate, three-resolution calibration sensitivity step, its independent audits, and artifact upload. The Python 3.10 and 3.12 jobs for that commit were still running at the time this report was prepared. Later commits extend the output with explicit paired-difference metrics and are being validated separately; the numerical results below remain tied to the cited archived artifact.

## Design

- Posterior-grid resolutions: 1,201, 2,401, and 4,801 equally spaced points on [40, 100].
- Replicates per resolution: 300.
- Seed: `20261025`.
- Identical prior, host redshifts/weights, noise model, and detection interval across resolutions.
- Each resolution's raw records and summary were independently audited.
- The same seeded truth, host, distance, and rejection-attempt records were confirmed across the three resolutions.

This paired design isolates posterior-grid discretization: differences in PIT values arise from the posterior CDF's numerical grid rather than different simulated observations.

## Results

| Grid points | Mean PIT | Monte Carlo SE | PIT variance | KS statistic | Descriptive KS p-value | Maximum decile-frequency deviation | Independent audit |
|---:|---:|---:|---:|---:|---:|---:|---|
| 1,201 | 0.48845956 | 0.01612666 | 0.07802077 | 0.04037250 | 0.69683395 | 0.02667 | PASS |
| 2,401 | 0.48845974 | 0.01612668 | 0.07802098 | 0.04037183 | 0.69685337 | 0.02667 | PASS |
| 4,801 | 0.48845978 | 0.01612669 | 0.07802103 | 0.04037169 | 0.69685753 | 0.02667 | PASS |

I independently compared replicate-level PITs from the downloaded CSVs against the finest (4,801-point) grid:

| Grid points compared with 4,801 | Mean paired PIT difference | RMS paired PIT difference | Maximum absolute paired difference |
|---:|---:|---:|---:|
| 1,201 | -2.21e-7 | 2.20e-6 | 1.48e-5 |
| 2,401 | -3.52e-8 | 4.41e-7 | 3.05e-6 |
| 4,801 | 0 | 0 | 0 |

All three grid resolutions passed the independent raw-record/PIT/summary audit. The main 1,000-replicate selected-likelihood calibration artifact in the same job also passed its audit.

## Interpretation

For this particular matched synthetic model and seed, the PIT is numerically very stable across the tested grid resolutions. The root-mean-square per-replicate PIT difference is about (2.2\times10^{-6}) for 1,201 points and (4.4\times10^{-7}) for 2,401 points, relative to 4,801 points. Those discrepancies are much smaller than the roughly 0.0161 Monte Carlo standard error of the mean PIT in this 300-replicate screening run.

The mean PIT is approximately 0.4885, about 0.71 Monte Carlo standard errors below 0.5. The descriptive KS results do not indicate an obvious departure from uniformity in this run, but neither that result nor numerical stability proves exact calibration.

## Limitations and follow-up

- This is one seed and one finite synthetic model configuration; the findings should not be generalized to all noise levels, selection thresholds, host mixtures, or real events.
- Three hundred replicates are suitable for CI screening, not a high-powered calibration claim.
- A stronger study should repeat the paired grid comparison with multiple independent seeds and additional noise/selection scenarios, define acceptable per-replicate PIT tolerances in advance, and preserve raw artifacts.
- Internal audit success establishes that this synthetic experiment is reproducible and internally consistent. It does not establish astrophysical validity or validate real detector-selection and galaxy-catalogue models.
