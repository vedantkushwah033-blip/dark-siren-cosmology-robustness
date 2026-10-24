# Matched host-structure control: results and interpretation (2026-10-24)

## Scope

This is a controlled synthetic simulation, not an analysis of real LVK event data and not a measurement of the observed Universe's Hubble constant. It tests whether replacing a two-host redshift mixture with a single effective redshift changes inference under a deliberately simplified distance/noise/selection model.

Workflow: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40100518091  
Code commit: `c0b076633136019c30222e4ecafac9d94b880b93`

## Design

- 300 detected observations at each of nine fixed truths: H0 = 52, 55, 60, 65, 70, 75, 80, 85, 88 km s^-1 Mpc^-1 (2,700 observations total).
- Two-host mixture: redshifts [0.02, 0.04], intrinsic host weights [0.7, 0.3].
- Collapsed approximation: one effective redshift z = 0.026.
- Distance noise sigma = 20 Mpc; detection interval (0, 120) Mpc.
- Inference prior [50, 90] km s^-1 Mpc^-1; 4,001-point grid spanning [25, 125].
- Each model was analysed with and without the toy selection-normalization term. Paired bootstrap intervals used 5,000 draws.
- Contrast is absolute error of collapsed model minus absolute error of two-host model. Positive values favor the two-host model for that truth; negative values favor the collapsed approximation on absolute error.

## Aggregate metrics

These are equal-weight averages across the nine selected truths, not estimates of performance under a population distribution of true H0 values.

| Analysis | Mean bias | Mean absolute error | Mean RMSE | Mean 68% coverage | Mean 95% coverage |
|---|---:|---:|---:|---:|---:|
| Two-host, uncorrected | +0.8238 | 8.7436 | 10.0977 | 60.37% | 96.59% |
| Two-host, corrected | -0.7248 | 8.6098 | 10.1659 | 61.96% | 96.81% |
| Collapsed, uncorrected | +8.7878 | 11.3474 | 12.1233 | 42.41% | 74.19% |
| Collapsed, corrected | +5.2927 | 9.7323 | 11.2748 | 55.22% | 90.70% |

The average metrics conceal strong truth-dependent behavior. In particular, mean coverage should not be interpreted as uniform coverage across the parameter range.

## Paired absolute-error contrasts

The table gives mean( |error_collapsed| - |error_two-host| ) and its paired-bootstrap 95% interval, in H0 units.

| True H0 | Uncorrected contrast (95% interval) | Corrected contrast (95% interval) |
|---:|---:|---:|
| 52 | +9.50 [9.35, 9.64] | +5.71 [5.45, 5.96] |
| 55 | +9.33 [9.17, 9.50] | +5.90 [5.65, 6.13] |
| 60 | +9.06 [8.88, 9.24] | +6.04 [5.72, 6.33] |
| 65 | +6.84 [6.48, 7.18] | +4.34 [3.87, 4.82] |
| 70 | +4.06 [3.48, 4.62] | +2.41 [1.83, 2.96] |
| 75 | +0.55 [-0.12, 1.22] | ~0.00 [-0.63, 0.60] |
| 80 | -3.11 [-3.74, -2.47] | -3.01 [-3.51, -2.50] |
| 85 | -6.42 [-6.73, -6.11] | -5.66 [-5.88, -5.42] |
| 88 | -6.39 [-6.68, -6.11] | -5.63 [-5.82, -5.42] |

## Detection changes the host mixture

The intrinsic weights 70/30 are not the mixture weights among detected events. Detection probability depends on each host's distance distribution and on H0. In this simulation, the observed fraction assigned to the z=0.04 host was 0%, 0%, 0%, 0%, 0%, 0.67%, 3.00%, 6.67%, and 5.67% at H0 = 52, 55, 60, 65, 70, 75, 80, 85, and 88 respectively. The expected detected fractions under the simulation's Gaussian-noise and truncation assumptions are approximately 0%, 0%, 0.002%, 0.029%, 0.230%, 1.00%, 2.84%, 5.93%, and 8.22%, respectively. The differences at 300 observations per truth are compatible with finite-sample variation; this is a simulation-specific check, not an empirical astrophysical estimate.

## Interpretation and limitations

1. **Host-structure compression matters in this toy setup.** The single-redshift approximation has larger absolute error at lower tested truths, while it has smaller absolute error at higher tested truths. There is no universal winner over this grid.
2. **The result is truth-dependent.** The sign reversal around H0 between 75 and 80 means a single pooled metric can hide meaningful variation.
3. **Selection correction is not uniformly beneficial.** The correction changes average bias and error differently for the two model structures; it should not be described as universally improving inference.
4. **Coverage remains imperfect.** The 68% interval coverage averages are below nominal for all four analyses. More repetitions alone would narrow Monte Carlo uncertainty but would not fix model misspecification, prior-boundary effects, or an inadequate likelihood.
5. **This is a deliberately simplified control.** It uses a two-host synthetic mixture, Gaussian distance noise, a hard detection interval, a fixed prior and grid, and a toy selection normalization. It does not model the full gravitational-wave likelihood, real galaxy-catalogue completeness, peculiar velocities, cosmological parameter covariance, or the real LVK selection function.

## Reproducibility

The experiment code, protocol, tests, and GitHub Actions workflow are tracked in the repository at the code commit above. The workflow completed successfully, and the repository's test-suite workflow also completed successfully:
- Experiment: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40100518091
- Test suite: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40100518110
- Additional tests: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40100518250

For publication, preserve the raw artifact and manifest with the final analysis and state the exact artifact digest. Do not present these results as observational constraints on H0.
