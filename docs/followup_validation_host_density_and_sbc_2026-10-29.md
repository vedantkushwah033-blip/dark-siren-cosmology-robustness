# Follow-up validation note: host-density sensitivity and H0-dependent-noise calibration

**Date:** 2026-10-29  
**Branch:** `research/audit-posterior-grid-integrity`  
**Scope:** Synthetic methodology only. None of these experiments is a real LVK analysis or a new measurement of the Hubble constant.

## Executive result

The newest diagnostic runs completed successfully and the Python 3.10, 3.11, and 3.12 numerical-diagnostics jobs passed. Two findings deserve explicit treatment in the manuscript:

1. **Host-density bandwidth and selection-law misspecification matter in the toy model.** A 13-configuration paired screening campaign used 40 replicates per configuration, 6,000 parent galaxies, and 40 detected events per replicate. It is useful for identifying sensitivities, not for making precise coverage claims.
2. **Prior-averaged calibration does not establish conditional calibration at fixed H0.** In a separate 1,200-replicate selected-likelihood experiment, the prior-drawn mode had mean rank 0.514 and a descriptive rank-decile chi-square p-value of 0.450. Fixed-truth modes showed substantial nonuniform ranks and conditional bias in this discrete-host, three-event model. The full raw artifact was independently reconstructed to floating-point precision.

## Host-density campaign

Artifact: GitHub Actions run [40797592794](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40797592794), artifact `host-density-validation-campaign-5bbe5a434878af396e8c58e3ddfe7af0356e29ed`.

Design: paired seeds 20261029–20261068; 40 replicates per configuration; 6,000 parent galaxies; 40 detected events; fixed true H0=70 km/s/Mpc; Gaussian distance noise; fixed-threshold detection. Configurations covered seven KDE bandwidth choices (weighted Scott plus factors 0.05, 0.1, 0.2, 0.5, 0.8, 1.2), coarse/fine integration grids, and supplied selection-law midpoint/width misspecification. Results are marked exploratory in the campaign artifact.

| Configuration | Correct-completeness H0 RMSE | Correct-completeness mean bias | Catalogue-only H0 RMSE | Notes |
|---|---:|---:|---:|---|
| Weighted Scott | 3.896 | +1.180 | 9.165 | Baseline screening setting |
| Bandwidth factor 0.05 | 3.787 | +0.951 | 9.297 | Low bandwidth; coverage estimate is very noisy |
| Bandwidth factor 0.1 | 3.807 | +0.990 | 9.271 | Similar to 0.05 |
| Bandwidth factor 0.2 | 3.883 | +1.153 | 9.157 | Similar to Scott in H0 error |
| Bandwidth factor 0.5 | 4.440 | +2.137 | 8.122 | Smoothing worsens oracle/weighted H0 error |
| Bandwidth factor 0.8 | 4.990 | +2.954 | 6.202 | Strong oversmoothing in this setup |
| Bandwidth factor 1.2 | 5.330 | +3.446 | 4.238 | Catalogue-only appears less poor, but both density and oracle results are distorted |
| Coarse grid (81 H0 × 401 z) | 3.896 | +1.180 | 9.165 | Baseline grid |
| Fine grid (201 H0 × 1001 z) | 3.896 | +1.180 | 9.165 | Nearly identical posterior means/RMSE at displayed precision |
| Supplied midpoint 0.18 | 3.801 | +0.750 | 9.165 | Misspecified law |
| Supplied midpoint 0.26 | 3.719 | +0.709 | 9.165 | Lower H0 RMSE here does not imply the wrong selection law is better generally |
| Supplied width 0.025 | 4.669 | +2.578 | 9.165 | Stronger error under narrow transition misspecification |
| Supplied width 0.10 | 4.245 | −2.536 | 9.165 | Bias changes sign; RMSE and interval behavior also change |

**Interpretation and limits.** In this toy design, ignoring catalogue incompleteness yields severe H0 bias under the Scott baseline, while inverse-inclusion weighting brings the error closer to the oracle reference. This is not a universal guarantee for weighting: misspecified inclusion probabilities can materially change the result, and the deliberately misspecified midpoint 0.26 even has lower H0 RMSE than the correctly supplied midpoint in this finite sample. That is a reminder that RMSE is a realized, design-specific score—not a measure of whether the model is correctly specified.

The 40-replicate campaign is a screening exercise. Its empirical coverage values have wide uncertainty and must not be used as evidence of calibrated 68% intervals. The fine/coarse comparison suggests integration-grid effects are small for this particular summary, but does not replace a general convergence analysis.

**Audit distinction.** The archive manifest and byte counts were checked, and aggregate metrics were independently rederived from each raw CSV. The new auditor then reconstructed all 1,560 fit rows in separate per-configuration local runs, and all passed. The clean all-configurations CI reproduction/audit passed at [#40802415000](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40802415000). Artifact `12866753211` contains the full audit JSON; it records PASS, all 1,560 fits reconstructed, 403 aggregate fields checked with zero difference, and 14 manifest files verified. One boundary-sensitive coverage-flag mismatch for oracle/bandwidth-factor-0.5/replicate-25 is disclosed: the independent endpoint lies 9.37×10⁻⁵ from the truth, below the prespecified 0.10 endpoint tolerance. The older compact two-replicate smoke audit no longer defines the full scope of the local audit.

## H0-dependent-noise selected-likelihood calibration

Artifact: GitHub Actions run [40797963979](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40797963979), artifact `h0-dependent-noise-sbc-e9a512d46aeeb7070649ac8f851c1d86103f4379`.

Design: 1,200 rows total, 300 replicates in each mode (prior-drawn H0, fixed H0=60, fixed H0=70, fixed H0=80); three detected events per replicate; 12 discrete host redshifts; H0 grid with 161 points. This is a discrete-host distance-only model, not realistic gravitational-wave inference.

| Truth mode | Mean posterior-mean bias (km/s/Mpc) | 68% interval coverage | 90% interval coverage | Mean SBC rank | Rank-decile chi-square p (descriptive) |
|---|---:|---:|---:|---:|---:|
| H0 drawn from prior | −0.664 | 72.7% (218/300) | 93.3% (280/300) | 0.514 | 0.450 |
| Fixed H0=60 | +4.901 | 67.3% (202/300) | 92.3% (277/300) | 0.391 | 7.14e−7 |
| Fixed H0=70 | +0.949 | 74.3% (223/300) | 96.3% (289/300) | 0.470 | 2.34e−12 |
| Fixed H0=80 | −4.368 | 88.3% (265/300) | 97.0% (291/300) | 0.604 | 2.23e−42 |

The independent auditor reconstructed all 1,200 records, posterior summaries, and ranks; maximum absolute posterior-summary reconstruction error was about 4.3e−14, and the audit status was PASS. The rank-decile p-values are descriptive diagnostics, not a standalone proof or a multiple-testing-adjusted decision rule.

**Interpretation.** The prior-drawn mode looks broadly compatible with prior-averaged rank calibration, but the fixed-truth modes show clear conditional rank distortions and H0-dependent bias. Therefore, prior-predictive calibration must not be presented as proof of conditional calibration across the H0 range. The result is a useful warning about this specific likelihood/selection/data-generating combination, not evidence that all dark-siren likelihoods are miscalibrated.

This 1,200-row artifact is a distinct model and experiment from the separate exact-density oracle diagnostic with 300 replicates per mode. Keep their numbers and claims separate.

## Validation status

- Host-density sensitivity workflow: both campaign and exact-density oracle jobs succeeded.
- H0-dependent-noise calibration: simulation, independent reconstruction, and artifact upload succeeded; auditor status PASS.
- Numerical diagnostics matrix: Python 3.10, 3.11, and 3.12 jobs completed successfully on the current branch. Python-3.12-only full numerical experiments were run in the 3.12 job; corresponding steps are skipped on 3.10/3.11 by design.
- This establishes reproducibility and internal consistency for the declared synthetic experiments. It does **not** validate a real detector likelihood, real LVK event catalogues, real host-galaxy catalogues, or external scientific performance.

## Recommended manuscript language

Use the Stage B locked paired contrast as the primary confirmatory result, and present these newest results as supporting diagnostics/sensitivity analyses. Report the fixed-truth calibration failures rather than hiding them. Do not select a bandwidth or misspecified selection law merely because it gives a lower RMSE in these 40 replicates. For a stronger sensitivity claim, increase replicates under a preregistered plan and independently reconstruct the full campaign—not only the compact benchmark. Real-data validation remains a separate, currently unmet gate.

## Cross-check against the canonical campaign report

The canonical screening report is `docs/host_density_sensitivity_screening_results_2026-10-29.md`. It should be treated as the authoritative full table; the compact table above is a summary. To make the interpretation explicit, the default-bandwidth oracle reference has bias +1.123, RMSE 3.936, 68% coverage 52.5%, and density L1 error 0.058. The known-inclusion arm has bias +1.180, RMSE 3.896, coverage 55.0%, and density L1 error 0.082. The catalogue-only arm has bias −8.634, RMSE 9.165, coverage 10.0%, and density L1 error 0.537. The paired absolute-error contrast for known-inclusion weighting minus oracle is −0.057, with bootstrap 95% interval [−0.199, +0.083], which includes zero.

For full sensitivity details, including catalogue-only and weighted metrics at every bandwidth, grid differences, and all selection-curve variants, use the canonical report rather than extrapolating from the shortened table here. The sensitivity campaign's 40-replicate-per-setting design remains an evidence limit despite the full independent posterior reconstruction; nominal coverage is still imprecise at this replicate count.


## Exact-density oracle audit (latest revision)

The exact-density oracle job in [run #40802415000](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40802415000) completed successfully. Its raw CSV has 1,200 rows (300 per truth regime), 1,129,316 bytes, SHA-256 `a38c64a2b9be9afce46c8f7b88659fa9f34f08a146d68376142b2184ba14cbee`; summary JSON 3,881 bytes, SHA-256 `a26c3b62d2612f04f34d9c2fb416ca1ae6024a69845ef7622666f1fe2efe3d72`. The independent oracle auditor passed, replayed all seeds and detected distances exactly, reconstructed every posterior and summary aggregate, and verified the manifest. Maximum posterior summary error was below 1.9×10⁻¹³. Artifact ID `12867070639`; archive digest `sha256:4853794c2840980a8efc8bd6d82e1b82fc492aea3af4b291818db24464bd4134`. This continuous-density fixed-noise oracle is distinct from the H₀-dependent-noise discrete-host study.
