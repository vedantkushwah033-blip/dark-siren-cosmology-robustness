> **Supersession notice — 2026-10-27:** These results were generated before correction of the direct missing-host mixture likelihood, which previously applied the catalogue population-mass factor twice. Treat all numerical results below as historical/provisional, not as the current validated result. Regenerate the experiment from the corrected implementation and pass the independent raw-artifact audit before using these estimates in a paper, abstract, or application. See [`missing_host_mixture_likelihood_impact_audit_2026-10-27.md`](missing_host_mixture_likelihood_impact_audit_2026-10-27.md).

# Paired host-misspecification sensitivity surface: audited results (2026-10-26)

## Executive summary

The full 108-cell sensitivity experiment completed successfully, and its independent auditor passed. It generated 21,600 unique synthetic observations, applied four inference arms to each observation (86,400 posterior fits), independently reconstructed all posterior outputs, and reconstructed all 1,296 cell-level paired-bootstrap contrast rows.

Across equal-weight averages of the 27 truth/noise/redshift cells at each hidden-host fraction, omitting a hidden host became progressively more damaging as its population fraction increased. The omitted-host arm's posterior-median RMSE rose from 6.278 at fraction 0 to 9.929 at fraction 0.45, while its 95% interval coverage fell from 98.3% to 86.6%. The modelled-hidden-host arm had RMSE 7.499 and 95% coverage 98.6% at fraction 0.45. The paired, equal-cell bootstrap estimates support a conditional conclusion that representing the hidden component improves performance in this specific simulator; the modelled arm's 68% intervals overcover, so this is not proof of general calibration.

## Reproducibility record

- Successful workflow and independent audit: [GitHub Actions run](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40280537559)
- Raw artifact: [21,600 observations, per-run posterior summaries, cell summaries, and paired contrasts](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40280537559/artifacts/12632197386)
- Artifact SHA-256: `452cc37bab97fb95605d38a9e1eff1643e8cb6d66a417ddb107af7beff56bd67`
- Runner/auditor/workflow revision: `cbcca36982d3e2b4aee96e7b30067361797eedb5`
- Data seed: `20261029`; cell bootstrap seed: `20261030`
- Full design: 108 cells, 200 observations per cell, four inference arms; 86,400 posterior fits
- Outputs: 432 cell-arm summaries and 1,296 paired contrast rows
- Independent audit: PASS; all 86,400 posterior outputs, summaries, and paired intervals reconstructed
- Protocol: [prespecified sensitivity design](host_misspecification_sensitivity_surface_protocol.md)
- Reproducible equal-cell summary/aggregate-bootstrap code: [summary script](../scripts/summarize_host_misspecification_sensitivity_surface.py); it checks the raw-to-summary reconstruction before calculating marginal effects.

The factorial crosses injected H0 = 60, 70, 80; fractional distance-noise levels 0.05, 0.10, 0.20; hidden-host fractions 0, 0.15, 0.30, 0.45; and hidden-host redshifts 0.025, 0.030, 0.035. Each component likelihood is conditioned on positive distance, matching the simulator's rejection of non-positive observations. The prior is uniform on [50, 90].

## Equal-cell mean performance by hidden-host fraction

Each entry is the equal-cell average over 27 cells (three injected truths × three noise levels × three hidden-host redshifts), with 200 observations per cell. Errors and interval widths are in km s^-1 Mpc^-1; coverage is a proportion.

| Hidden fraction | Inference arm | Mean median error | RMSE | 68% coverage | 95% coverage | Mean 68% width | Mean 95% width |
|---:|---|---:|---:|---:|---:|---:|---:|
| 0.00 | Base weights 0.7/0.3; hidden omitted | +0.589 | 6.278 | 0.790 | 0.983 | 15.790 | 26.709 |
| 0.00 | Reversed weights 0.3/0.7; hidden omitted | +1.858 | 6.579 | 0.775 | 0.983 | 15.799 | 26.595 |
| 0.00 | Equal weights; hidden omitted | +1.195 | 6.326 | 0.789 | 0.983 | 15.871 | 26.702 |
| 0.00 | Hidden host modelled | +0.589 | 6.278 | 0.790 | 0.983 | 15.790 | 26.709 |
| 0.15 | Base weights 0.7/0.3; hidden omitted | +0.258 | 7.916 | 0.748 | 0.944 | 16.099 | 27.032 |
| 0.15 | Reversed weights 0.3/0.7; hidden omitted | +1.988 | 8.144 | 0.729 | 0.943 | 15.988 | 26.879 |
| 0.15 | Equal weights; hidden omitted | +1.083 | 7.901 | 0.745 | 0.944 | 16.145 | 27.006 |
| 0.15 | Hidden host modelled | +0.417 | 7.145 | 0.828 | 0.987 | 18.479 | 29.802 |
| 0.30 | Base weights 0.7/0.3; hidden omitted | -0.014 | 9.011 | 0.710 | 0.904 | 16.367 | 27.490 |
| 0.30 | Reversed weights 0.3/0.7; hidden omitted | +1.981 | 9.209 | 0.692 | 0.899 | 16.170 | 27.193 |
| 0.30 | Equal weights; hidden omitted | +0.977 | 9.023 | 0.707 | 0.902 | 16.331 | 27.424 |
| 0.30 | Hidden host modelled | +0.096 | 7.163 | 0.844 | 0.987 | 19.414 | 30.484 |
| 0.45 | Base weights 0.7/0.3; hidden omitted | -0.411 | 9.929 | 0.661 | 0.866 | 16.814 | 27.914 |
| 0.45 | Reversed weights 0.3/0.7; hidden omitted | +2.080 | 10.192 | 0.629 | 0.865 | 16.326 | 27.544 |
| 0.45 | Equal weights; hidden omitted | +0.900 | 10.020 | 0.651 | 0.866 | 16.694 | 27.816 |
| 0.45 | Hidden host modelled | +0.272 | 7.499 | 0.821 | 0.986 | 19.456 | 30.834 |

At hidden fraction zero, the base and hidden-modelled arms are exactly equivalent in all six summary metrics across all 27 cells. This is the prespecified negative control.

## Paired, equal-cell contrasts across the hidden-host fraction

The following marginal contrasts are **equal-cell averages over 27 cells per fraction**. Their intervals were calculated with 1,000 stratified paired bootstrap resamples, using the same resampled replicate indices across arms within each cell (seed 20261031). These are separate marginal intervals, not a pooled all-design interval. Coverage differences are absolute proportions: -0.080 is -8.0 percentage points.

### Hidden host omitted minus hidden host modelled

| Hidden fraction | Metric | Difference | Paired bootstrap 95% interval |
|---:|---|---:|---:|
| 0.00 | RMSE | 0.000 | [0.000, 0.000] |
| 0.15 | RMSE | +0.771 | [+0.680, +0.867] |
| 0.30 | RMSE | +1.848 | [+1.680, +1.993] |
| 0.45 | RMSE | +2.430 | [+2.215, +2.635] |
| 0.00 | 68% coverage | 0.000 | [0.000, 0.000] |
| 0.15 | 68% coverage | -0.080 | [-0.088, -0.072] |
| 0.30 | 68% coverage | -0.134 | [-0.144, -0.123] |
| 0.45 | 68% coverage | -0.160 | [-0.173, -0.147] |
| 0.00 | 95% coverage | 0.000 | [0.000, 0.000] |
| 0.15 | 95% coverage | -0.042 | [-0.048, -0.037] |
| 0.30 | 95% coverage | -0.083 | [-0.089, -0.076] |
| 0.45 | 95% coverage | -0.120 | [-0.128, -0.112] |
| 0.15 | Mean 68% width | -2.380 | [-2.506, -2.256] |
| 0.30 | Mean 68% width | -3.047 | [-3.219, -2.866] |
| 0.45 | Mean 68% width | -2.642 | [-2.848, -2.446] |
| 0.15 | Mean 95% width | -2.770 | [-2.897, -2.643] |
| 0.30 | Mean 95% width | -2.994 | [-3.136, -2.852] |
| 0.45 | Mean 95% width | -2.920 | [-3.069, -2.774] |

Negative coverage contrasts mean the omitted-host arm covers less often than the modelled-host arm. Negative width contrasts mean the omitted-host arm reports narrower intervals. The omitted-host model has both higher RMSE and lower coverage as the hidden-host fraction grows. The modelled arm improves coverage but has 68% coverage of 82.1–84.4% for nonzero fractions, above the nominal 68% level. It therefore remains overconservative in these cells.

### Reversed two-host weights minus base weights

| Hidden fraction | Mean median-error difference | Paired 95% interval | RMSE difference | Paired 95% interval | 68% coverage difference | Paired 95% interval |
|---:|---:|---:|---:|---:|---:|---:|
| 0.00 | +1.269 | [+1.192, +1.348] | +0.302 | [+0.228, +0.384] | -0.015 | [-0.020, -0.011] |
| 0.15 | +1.730 | [+1.616, +1.845] | +0.228 | [+0.110, +0.336] | -0.019 | [-0.024, -0.013] |
| 0.30 | +1.995 | [+1.865, +2.123] | +0.198 | [+0.091, +0.308] | -0.019 | [-0.025, -0.011] |
| 0.45 | +2.491 | [+2.344, +2.644] | +0.263 | [+0.146, +0.375] | -0.033 | [-0.040, -0.025] |

Reversing the weights shifts the posterior-median error upward in every hidden-fraction stratum and modestly increases RMSE; its effect grows for the mean error at higher hidden-host fractions. These are equal-cell averages over this design, not universal effects of host weighting.

## Cell-level heterogeneity check

The uploaded paired-contrast table also permits a descriptive count of the 27 individual cells in each hidden-fraction stratum where the cell-level 95% bootstrap interval lies entirely on the expected side of zero. These are exploratory counts across many cellwise comparisons, with no multiple-comparison adjustment; they should not be treated as confirmatory hypothesis tests.

| Hidden fraction | RMSE: estimate positive / CI entirely positive | 68% coverage: estimate negative / CI entirely negative | 95% coverage: estimate negative / CI entirely negative |
|---:|---:|---:|---:|
| 0.00 | 0 / 0 | 0 / 0 | 0 / 0 |
| 0.15 | 25 / 21 | 19 / 18 | 13 / 12 |
| 0.30 | 22 / 19 | 20 / 18 | 15 / 12 |
| 0.45 | 22 / 19 | 20 / 16 | 15 / 13 |

At nonzero fractions, the direction is not identical in every cell. The macro-average pattern is supported by many cells, but variation with injected H0, noise level, and hidden-host redshift remains material and is a reason to avoid a universal claim.

## Interpretation and caveats

1. **The zero-fraction negative control passed.** With no hidden population, the modelled-hidden-host and base-weight arms give identical posteriors and summaries. This helps rule out an implementation artifact in the added mixture component.
2. **Omitted-host degradation is systematic across the tested grid.** Increasing hidden fraction worsens average RMSE and reduces 68% and 95% coverage for the base-weight arm. Paired contrasts against the correctly specified mixture grow in magnitude.
3. **Representing the hidden component broadens intervals.** That improves coverage here, but does not mean the resulting credible intervals are calibrated: 68% coverage is still too high for the modelled arm and 95% coverage averages about 98.6–98.7% at nonzero fractions.
4. **Weight reversal is a distinct, smaller effect in this design.** It shifts mean median error and slightly increases RMSE, but the magnitude depends on the hidden-host fraction and the other cell factors.
5. **Uncertainty is conditional on the design.** Bootstrap intervals quantify finite-simulation uncertainty within cells and equal-cell marginal summaries. They do not include model uncertainty, uncertainty over the chosen factor ranges, or realism of the data-generating model.

## Limits and next work

This remains a discrete-host, Gaussian-distance toy model with a fixed prior and a single seeded run. It omits realistic gravitational-wave likelihoods, detector selection, galaxy-catalogue completeness, peculiar velocities, and full hierarchical population inference. No real observational H0 estimate is made. These findings do not establish that the project is scientifically original or publication-ready.

The next useful research phase is to test prior-support sensitivity jointly with the hidden-host factorial, investigate the coverage behavior near prior boundaries, and validate the likelihood/selection assumptions against a realistic simulated or public gravitational-wave/catalogue benchmark before making astrophysical claims.
