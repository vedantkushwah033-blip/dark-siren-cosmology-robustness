> **Supersession notice — 2026-10-27:** These results were generated before correction of the direct missing-host mixture likelihood, which previously applied the catalogue population-mass factor twice. Treat all numerical results below as historical/provisional, not as the current validated result. Regenerate the experiment from the corrected implementation and pass the independent raw-artifact audit before using these estimates in a paper, abstract, or application. See [`missing_host_mixture_likelihood_impact_audit_2026-10-27.md`](missing_host_mixture_likelihood_impact_audit_2026-10-27.md).

# Paired host-misspecification experiment: audited results (2026-10-26)

## Executive summary

The paired experiment completed on GitHub Actions and passed an independent reconstruction audit. It generated 1,000 synthetic observations from a three-host population and analyzed every observation under all four inference assumptions (4,000 posterior fits). The auditor independently regenerated the seeded observations, reconstructed all 4,000 posteriors and summary statistics, and recomputed all 12 paired bootstrap contrast rows.

The most important result in this fixed toy design is that explicitly modelling the 30% hidden-host component lowered posterior-median RMSE from 8.216 to 4.902 km s^-1 Mpc^-1 and raised empirical 68% interval coverage from 68.7% to 87.4%, while producing wider intervals. This is evidence about this simulation setup only, not a general performance guarantee.

## Reproducibility record

- Workflow: [Paired host-misspecification experiment, successful run](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40280238225)
- Artifact: [raw records, summaries, and paired contrasts](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40280238225/artifacts/12631703296)
- Artifact SHA-256: `abc6f451311b921a0370e94356848a65734ecfb554f8443b43b740bf2b9f0d47`
- Runner/auditor/workflow revision: `49eb51371465bb0e819e5a8a7b5d0704b32dea95`
- Data seed: `20261027`
- Bootstrap seed: `20261028`
- Replicates: 1,000 shared observations; 4 inference arms per observation
- Paired bootstrap: 4,000 percentile-bootstrap resamples
- Audit: PASS; 4,000 posterior reconstructions and 12 paired contrast rows reconstructed
- Protocol: [paired host-misspecification protocol](paired_host_misspecification_protocol.md)

## Inference-arm results

All rows use the same 1,000 simulated observations, with injected H0 = 70 km s^-1 Mpc^-1 and a uniform prior on [50, 90].

| Inference arm | Mean median error | RMSE of posterior median | 68% coverage | 95% coverage | Mean 68% width | Mean 95% width |
|---|---:|---:|---:|---:|---:|---:|
| Two-host weights 0.7/0.3; hidden host omitted | +0.387 | 8.333 | 68.7% | 99.3% | 16.706 | 28.888 |
| Two-host weights reversed to 0.3/0.7; hidden host omitted | +3.065 | 8.346 | 68.3% | 99.3% | 16.255 | 28.391 |
| Equal two-host weights; hidden host omitted | +1.694 | 8.216 | 68.7% | 99.3% | 16.694 | 28.679 |
| Equal two-host weights plus 30% hidden-host component | +1.586 | 4.902 | 87.4% | 99.3% | 21.880 | 32.750 |

Widths and errors are in km s^-1 Mpc^-1. The 0.7/0.3 arm is only nominally “correct” for the relative weights of the two listed hosts; the generator also includes z = 0.03 with probability 0.30. It therefore omits a true population component and is not a fully specified oracle model.

## Paired contrasts

The contrasts are **left arm minus right arm**. Coverage contrasts are absolute proportion differences, so 0.187 means 18.7 percentage points. Intervals are percentile-bootstrap 95% intervals using the same replicate indices for both arms.

| Contrast | Metric | Difference | Paired bootstrap 95% interval |
|---|---|---:|---:|
| Reversed minus nominal two-host weights | Mean median error | +2.678 | [+2.239, +3.135] |
| Reversed minus nominal two-host weights | RMSE | +0.014 | [-0.142, +0.168] |
| Reversed minus nominal two-host weights | 68% coverage | -0.004 | [-0.026, +0.017] |
| Reversed minus nominal two-host weights | 95% coverage | 0.000 | [0.000, 0.000] |
| Reversed minus nominal two-host weights | Mean 68% width | -0.450 | [-0.932, +0.002] |
| Reversed minus nominal two-host weights | Mean 95% width | -0.497 | [-0.741, -0.264] |
| Hidden host omitted minus modelled | Mean median error | +0.108 | [-0.286, +0.514] |
| Hidden host omitted minus modelled | RMSE | +3.314 | [+3.016, +3.611] |
| Hidden host omitted minus modelled | 68% coverage | -0.187 | [-0.213, -0.160] |
| Hidden host omitted minus modelled | 95% coverage | 0.000 | [0.000, 0.000] |
| Hidden host omitted minus modelled | Mean 68% width | -5.186 | [-5.590, -4.776] |
| Hidden host omitted minus modelled | Mean 95% width | -4.071 | [-4.379, -3.757] |

The degenerate bootstrap intervals for 95% coverage occur because, in this particular finite sample, the 95% coverage indicators agree across arms replicate by replicate. They should not be interpreted as evidence that the arms have identical 95% coverage in the underlying data-generating population.

## Interpretation

1. **Pairing makes the comparison sharper.** All inference arms are applied to identical observations, so the paired contrasts remove between-arm differences caused merely by drawing different host identities and noise realizations.
2. **The omitted-host penalty is material in this design.** Modelling the hidden component improves RMSE and 68% coverage, but it also broadens posterior intervals. Better coverage is not automatically better calibration: 87.4% is above the nominal 68% level.
3. **Reversing the two listed host weights changes average error but barely changes RMSE here.** The mean-error contrast is clearly positive, while the RMSE and 68% coverage intervals include zero. The 95% widths are modestly narrower under reversed weights. These effects are not universal statements about catalogue weights.
4. **The 95% intervals overcover.** Every arm has 99.3% empirical coverage in this run. That persistent overcoverage remains a calibration warning and motivates varying the truth, prior support, noise level, and host-mixture configuration.
5. **The paired bootstrap measures finite-simulation uncertainty only.** It does not include uncertainty in the model specification or establish generality across other seeds and designs.

## Limitations and next step

This is a simplified synthetic model with one injected H0, one noise scale, one prior, one host-redshift configuration, and one main data seed. It does not include realistic gravitational-wave likelihoods, detector selection, galaxy-catalogue completeness, peculiar velocities, or full hierarchical population inference. It is not a real-data validation and does not establish scientific originality or publication readiness.

The next research step is a prespecified sensitivity surface that varies injected H0, hidden-host fraction/redshift, noise, and prior support, preserving paired observations within every cell and independently auditing the complete factorial output.
