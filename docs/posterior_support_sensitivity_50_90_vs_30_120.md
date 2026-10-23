# Posterior-support sensitivity audit: 50–90 versus 30–120 km/s/Mpc

**Status:** both 1,400-fit audits completed successfully; this is a sensitivity analysis, not a calibration certificate.  
**Narrow-support workflow:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39962269634  
**Wide-support workflow:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39972965523  
**Wide-support test suite:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39972965547

## Design

The two audits use the same fixed truth grid (55, 60, 65, 70, 75, 80, 85 km/s/Mpc), 100 replicates per truth, two scenarios, and the same seed scheme. The posterior grid support changes from 50–90 to 30–120 km/s/Mpc. Each audit contains 1,400 posterior fits.

The scenarios are deliberately simple:
- `known_host`: inference uses only the simulated host's redshift.
- `unknown_host_complete_catalogue`: inference sums over the full simulated catalogue with uniform weights.

These are not realistic end-to-end LVK analyses and do not use real GW events or a real galaxy survey.

## Results

Coverage is the empirical fraction of central posterior intervals containing the simulated truth. Widths and errors are in km/s/Mpc. Each cell has n=100; the Monte Carlo standard error for coverage around 0.68 is about 4.7 percentage points.

| H0 truth | Scenario | Support | Mean signed median error | RMSE of median | 68% coverage | 95% coverage | Mean 68% width | Mean 95% width |
|---:|---|---|---:|---:|---:|---:|---:|---:|
| 55 | Known host | 50–90 | +5.97 | 8.29 | 68% | 98% | 15.99 | 29.37 |
| 55 | Known host | 30–120 | +4.02 | 10.04 | 64% | 97% | 20.15 | 41.65 |
| 55 | Unknown host | 50–90 | +9.36 | 11.79 | 56% | 100% | 19.05 | 32.56 |
| 55 | Unknown host | 30–120 | +7.63 | 16.66 | 69% | 99% | 28.95 | 56.90 |
| 60 | Known host | 50–90 | +4.24 | 7.98 | 70% | 97% | 17.02 | 30.68 |
| 60 | Known host | 30–120 | +4.48 | 11.25 | 64% | 92% | 21.51 | 43.51 |
| 60 | Unknown host | 50–90 | +8.01 | 11.09 | 60% | 100% | 20.56 | 33.88 |
| 60 | Unknown host | 30–120 | +10.73 | 19.90 | 54% | 98% | 32.26 | 61.87 |
| 65 | Known host | 50–90 | +2.10 | 7.19 | 78% | 100% | 17.38 | 31.06 |
| 65 | Known host | 30–120 | +3.63 | 10.76 | 71% | 96% | 22.20 | 44.26 |
| 65 | Unknown host | 50–90 | +4.61 | 8.36 | 94% | 100% | 21.44 | 34.97 |
| 65 | Unknown host | 30–120 | +9.60 | 19.39 | 62% | 100% | 33.34 | 63.62 |
| 70 | Known host | 50–90 | +1.26 | 7.24 | 75% | 99% | 16.99 | 29.97 |
| 70 | Known host | 30–120 | +4.74 | 12.34 | 67% | 92% | 23.44 | 45.33 |
| 70 | Unknown host | 50–90 | +1.87 | 6.75 | 92% | 98% | 22.20 | 35.51 |
| 70 | Unknown host | 30–120 | +10.55 | 19.71 | 56% | 98% | 35.56 | 66.47 |
| 75 | Known host | 50–90 | −1.39 | 7.15 | 74% | 98% | 16.15 | 28.39 |
| 75 | Known host | 30–120 | +3.08 | 11.90 | 59% | 97% | 23.46 | 44.75 |
| 75 | Unknown host | 50–90 | −1.41 | 4.79 | 94% | 100% | 23.05 | 36.38 |
| 75 | Unknown host | 30–120 | +7.71 | 15.30 | 81% | 100% | 36.62 | 68.91 |
| 80 | Known host | 50–90 | −3.27 | 6.89 | 80% | 99% | 15.07 | 26.42 |
| 80 | Known host | 30–120 | +2.99 | 11.77 | 63% | 98% | 23.80 | 44.43 |
| 80 | Unknown host | 50–90 | −6.06 | 7.64 | 85% | 100% | 22.83 | 36.16 |
| 80 | Unknown host | 30–120 | +4.18 | 14.46 | 84% | 100% | 36.37 | 67.97 |
| 85 | Known host | 50–90 | −4.80 | 6.76 | 78% | 97% | 13.31 | 23.19 |
| 85 | Known host | 30–120 | +4.26 | 11.55 | 70% | 96% | 24.10 | 43.81 |
| 85 | Unknown host | 50–90 | −9.25 | 9.88 | 78% | 98% | 23.20 | 36.47 |
| 85 | Unknown host | 30–120 | +5.05 | 12.67 | 92% | 98% | 37.71 | 69.66 |

## Interpretation

1. **Widening the support did not simply repair calibration.** It substantially increased posterior interval widths, but coverage changed in both directions across truth values and scenarios.
2. **The posterior median errors changed systematically.** In the wide-support audit, mean signed median errors are positive in all 14 cells, including at truths 75–85 where the narrow-support audit produced negative mean errors. This points to more than a simple edge-truncation explanation and requires source-level inspection of the likelihood, prior handling, and H0–distance/redshift relationship.
3. **RMSE increased in every cell** under the wide support in these matched-seed runs. For example, at truth 70 in the unknown-host scenario, RMSE increased from 6.75 to 19.71 km/s/Mpc. Wider credible intervals alone do not make the estimator better calibrated or more accurate.
4. The support sensitivity is large enough that the current posterior implementation should **not** be treated as validated. Do not select the range that yields the most attractive coverage. First inspect the posterior implementation and test it on a simple known analytic likelihood.
5. Coverage at n=100 is noisy. The Monte Carlo standard error depends on observed coverage; values near 0.68 have SE around 4.7 percentage points. These differences are diagnostic leads, not definitive hypothesis tests.
6. The two scenarios differ in what information is given to inference. The unknown-host scenario marginalizes across 500 catalogue redshifts, while the known-host scenario conditions on one selected host. The results do not by themselves isolate catalogue incompleteness or realistic host weighting.

## Next required steps

- Inspect `posterior_h0` implementation and how `prior_min`/`prior_max` interact with the supplied grid.
- Verify the distance-redshift relation and the generation of the simulated event distance/noise.
- Build a minimal unit test with a known Gaussian likelihood and flat prior; compare numerical posterior normalization, median, and credible intervals against analytic values.
- Add simulation-based calibration/rank diagnostics under a fully matched generative/inference model only after the posterior implementation passes those unit tests.
- Do not make claims about real LVK/GWTC events from these toy simulations.

## Reproducibility

Both workflows completed successfully and validated 14 summary rows and 1,400 per-run rows. The wide-support workflow also passed the repository test suite. The CSV artifacts are available from the linked GitHub Actions runs.


## Source-code inspection update (23 October 2026)

The pinned source at revision `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f` was inspected directly:

- `src/cosmology.py` generates distances with the low-redshift relation (d_L = cz/H_0).
- `src/simulations_v2.py` generates observed distance as a Gaussian draw around that same distance, with `sigma = 0.15 * true_distance`.
- `src/bayesian_inference.py` evaluates a Gaussian distance likelihood with mean (cz/H_0), normalizes candidate-galaxy weights, applies a flat bounded prior through the grid mask, and normalizes the posterior using trapezoidal integration.

The distance-redshift relation is matched between simulation and inference in this pinned revision; it is not the first mismatch to blame. The substantial support sensitivity remains unexplained and warrants a minimal single-host likelihood benchmark plus checks of posterior-grid and summary calculations. A targeted test has been added to check recovery near a known injected (H_0) for a narrow-noise single-host likelihood. This is a unit test, not a full simulation-based calibration test.
