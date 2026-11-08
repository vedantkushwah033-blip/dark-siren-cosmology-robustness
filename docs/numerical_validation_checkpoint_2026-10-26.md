# Numerical validation checkpoint — 2026-10-26

## Provenance

This report summarizes the completed GitHub Actions numerical-diagnostics workflow:

- Workflow: [run 40271224634](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40271224634)
- Tested commit: `5d502507d744506087a94c5422f8e2d4e594cb19`
- Python versions: 3.10, 3.11, and 3.12
- Workflow result: **success on all three versions**.
- Each version completed the selected-likelihood PIT run, independent raw-data and summary audit, grid-sensitivity experiment and audit, multi-seed aggregation and audit, focused convergence tests, posterior-grid stability diagnostic, and artifact upload.

The workflow run is the source of truth for the artifacts and exact environment. The archived CI experiment uses small multi-seed settings to keep continuous integration tractable; those settings are not a publication-grade Monte Carlo design.

## 1. Selected-likelihood PIT calibration

The Python 3.10 log reports the 1,000-replicate run with seed `20261025`:

| Diagnostic | Result |
|---|---:|
| Replicates | 1,000 |
| Mean PIT | 0.4966 |
| Monte Carlo SE of mean PIT | 0.0089 |
| Uniform KS statistic | 0.0249 |
| Descriptive KS p-value | 0.5588 |
| Maximum absolute decile-frequency deviation from 0.1 | 0.0230 |

The independent auditor returned `PASS` and reconstructed both replicate-level PIT values and the summary. These values are consistent with uniform PIT behavior for this particular synthetic experiment; they do not prove exact calibration. The KS p-value is descriptive, not a certificate of model validity.

## 2. Paired grid-resolution sensitivity

The workflow evaluated 300 paired replicates on 1,201, 2,401, and 4,801 point grids, using the same simulated records across resolutions. The independent comparison audit returned `PASS`.

The archived Python 3.11 summary reported:

| Grid points | Mean PIT | Monte Carlo SE | KS statistic | Descriptive KS p-value |
|---:|---:|---:|---:|---:|
| 1,201 | 0.48845956 | 0.01612666 | 0.04037250 | 0.69683395 |
| 2,401 | 0.48845974 | 0.01612668 | 0.04037183 | 0.69685337 |
| 4,801 | 0.48845978 | 0.01612669 | 0.04037169 | 0.69685753 |

For the 1,201-point grid relative to 4,801 points, the mean paired PIT difference across replicates was approximately −2.21e-7, with RMS paired difference approximately 2.20e-6 and maximum absolute paired difference approximately 1.48e-5. For 2,401 points, the corresponding RMS difference was approximately 4.40e-7 and maximum absolute difference approximately 3.05e-6. These are small relative to the Monte Carlo SE of the mean PIT in this 300-replicate screening run.

This supports numerical stability for this finite model configuration and these grids only. The finest grid is a reference, not an exact solution.

## 3. Multi-seed grid sensitivity: CI smoke test

The workflow also ran two seeds (101 and 202), 100 replicates per seed, and 1,201- and 2,401-point grids. All per-seed audits and the aggregate audit passed.

Across the two seeds, the reported mean PIT mean was about 0.51384, with sample SD about 0.04446 and SE about 0.03144. The mean paired PIT difference for the 1,201-point grid relative to the 2,401-point reference was approximately −2.25e-7; its across-seed SD was approximately 3.84e-8. This indicates very small discretization effects for this smoke configuration, while the seed-level estimate is highly imprecise because it uses only two seeds.

**Do not interpret two seeds as a robust estimate of between-seed variability.** The protocol recommends a materially larger set of independent seeds for substantive numerical claims.

## 4. Targeted posterior-grid convergence

A separate deterministic synthetic case was evaluated at successively finer grids:

| Grid points | Posterior median | Upper 95% endpoint | Change in upper endpoint from preceding grid |
|---:|---:|---:|---:|
| 1,001 | 51.71391501 | 60.84837304 | — |
| 2,001 | 51.71359073 | 60.74224929 | 0.10612375 |
| 4,001 | 51.71348568 | 60.71791973 | 0.02432955 |
| 8,001 | 51.71345074 | 60.71190518 | 0.00601455 |
| 16,001 | 51.71344642 | 60.71042810 | 0.00147708 |
| 32,001 | 51.71344513 | 60.71005463 | 0.00037347 |
| 64,001 | 51.71344470 | 60.70996170 | 0.00009293 |

The endpoint changes decrease consistently with refinement, which is evidence of convergence in this one synthetic case. This is not a global error bound and should not be generalized to every posterior shape, boundary condition, or selection regime.

## 5. What this checkpoint establishes—and what it does not

### Established for the tested commit

- CI passed on Python 3.10, 3.11, and 3.12.
- The selected-likelihood PIT raw records and summary passed independent reconstruction audits.
- Grid-sensitivity raw records, paired comparisons, and multi-seed aggregates passed their audits.
- The targeted synthetic posterior summaries stabilized as the grid was refined.

### Not established

- Calibration under every possible prior, noise level, host mixture, or detection threshold.
- Correctness of a realistic gravitational-wave detector selection function.
- Completeness or representativeness of a real galaxy catalogue.
- Validity of a real-event cosmological inference or agreement with observational data.
- A publication-ready uncertainty estimate from the two-seed CI smoke test.

## 6. Next research step

The next meaningful step is a pre-specified multi-seed stress study, not simply increasing grid density. Use at least 20 independent seeds as an initial sensitivity design, and vary selection-window boundaries and distance-noise levels while retaining paired simulations across grid resolutions. Twenty seeds is a practical starting point, not a guarantee of adequate precision. Before interpreting the study, define acceptable numerical tolerances and report seed-level results, uncertainty, all failed audits, and the limits of the toy model.

The scientific question is whether conclusions persist when the stochastic realization and key assumptions change—not merely whether the code reproduces its own outputs.
