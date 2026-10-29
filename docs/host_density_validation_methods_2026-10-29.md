# Host-density validation protocol and reproducibility status

Date: 2026-10-29. Scope: synthetic low-redshift (H_0) toy only.

## What changed

The benchmark now accepts a physical KDE bandwidth factor, defaulting to weighted Scott's rule, and deliberately misspecified catalogue-inclusion midpoint/width values for the inference correction. The generator still uses its fixed true selection law (floor 0.1, midpoint 0.22, width 0.05). The scoring-only target density is computed from the known toy galaxy intensity multiplied by the declared merger weight.

Each replicate records L1 host-density error, maximum-to-mean inverse-weight ratio, and effective-sample-size fraction alongside posterior mean/median, equal-tail 68% interval, and coverage. The command-line runner writes a manifest with SHA-256 hashes for its CSV and summary. The independent audit re-simulates selected replicate seeds, uses a separately implemented reflected-kernel KDE and Simpson quadrature on finer grids, and checks posterior mean/median/interval endpoints against fixed tolerances.

## Paired sensitivity design

The expedited screening campaign uses 40 replicates per setting, initial seed 20261029, 6000 parent galaxies, 40 detected events, Gaussian distance noise (sigma=140) Mpc, threshold 250 Mpc, and fixed truth (H_0=70). Settings are Scott plus physical bandwidth factors 0.05, 0.1, 0.2, 0.5, 0.8, and 1.2; coarse and fine numerical grids; and selection-law assumptions that shift the midpoint or widen/narrow the transition. The script computes arm bias/SE, RMSE, Wilson intervals for empirical 68% coverage, mean interval width, host-density L1 error, weight diagnostics, and paired absolute-error differences versus the oracle arm.

The screening campaign completed in GitHub Actions run [#40797592794](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40797592794) at source commit `5bbe5a434878af396e8c58e3ddfe7af0356e29ed`; artifact ID `12865051001`. Its 13-configuration archive contains each configuration's per-replicate CSV, aggregate JSON, and a SHA-256 manifest. I independently rechecked SHA-256 hashes and byte counts for all 14 manifest-listed files, then independently re-aggregated the per-replicate CSV values to verify the headline screening metrics. The separate posterior reconstruction audit covered a compact two-replicate benchmark only; it did **not** independently recompute every posterior in all 40-replicate campaign configurations. Thus the full screen is checksum-verified and its aggregates re-derived, while record-by-record posterior reconstruction remains limited to the smoke audit.

## Interpretation guardrails

- A correctly specified inclusion law is a strong toy-model control; real surveys rarely provide the exact true inclusion law.
- A finite weight ratio and reasonable effective sample size in this toy do not guarantee stable inverse-probability weighting under other selection laws.
- Oversmoothing can flatten the inferred host density and change (H_0) error without correcting catalogue incompleteness.
- Prior-averaged Bayesian calibration and fixed-truth estimator bias answer different questions. A posterior mean can be biased at fixed truth while remaining calibrated when truth is drawn from the prior.
- The study does not model a realistic GW detector, event-rate uncertainty, sky/inclination marginalization, peculiar velocities, redshift uncertainty, or full galaxy-survey hierarchical selection.

## Scientific gate

GO for exploratory synthetic-methods reporting with explicit limitations. NO-GO for claims that this result validates real LVK inference, universally corrects incompleteness, or measures real (H_0). No human expert review or CHIMERA v2.2 benchmark-data recheck is claimed by this protocol.


## Expedited execution and decision rule

The first 100-per-setting CI run remained in progress longer than useful for rapid iteration. The current workflow therefore uses 40 paired replicates per configuration and 300 replicates per oracle truth regime on coarser 401-redshift/81-H0 grids. This is explicitly a screening campaign. Do not use its empirical coverage as a high-precision estimate or make a confirmatory claim from the exploratory sensitivity surface. For any configuration whose conclusion appears decision-relevant, extend it using a new, disjoint, preregistered seed range and a precision target based on paired variability. The independent fine-grid smoke audit is separate.


## Completed screening result and provenance (2026-10-29)

The screening used 40 paired replicates per configuration, 13 configurations, seeds 20261029–20261068, 6,000 parent galaxies, 40 detected events, fixed truth H0=70 km s⁻¹ Mpc⁻¹, Gaussian distance noise 140 Mpc, and detection threshold 250 Mpc. Campaign archive digest: `sha256:bfd9132783c354a929f0e633e0b0536361e8427cbb06ae1f94e451aae48470e8`. All 14 manifest-listed campaign files passed SHA-256 and byte-count checks after download. The separate two-replicate default-setting audit also passed its independent resimulation/posterior reconstruction and manifest checks.

### Baseline and bandwidth screen

| Configuration | Catalogue-only bias | RMSE | 68% coverage | Density L1 | Known-inclusion bias | RMSE | 68% coverage | Density L1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Weighted Scott default | −8.634 | 9.165 | 10.0% | 0.537 | +1.180 | 3.896 | 55.0% | 0.082 |
| Physical factor 0.05 | −8.786 | 9.297 | 7.5% | 0.537 | +0.951 | 3.787 | 57.5% | 0.075 |
| Physical factor 0.10 | −8.757 | 9.271 | 7.5% | 0.536 | +0.990 | 3.807 | 55.0% | 0.055 |
| Physical factor 0.20 | −8.625 | 9.157 | 10.0% | 0.537 | +1.153 | 3.883 | 55.0% | 0.076 |
| Physical factor 0.50 | −7.400 | 8.122 | 20.0% | 0.531 | +2.137 | 4.440 | 55.0% | 0.231 |
| Physical factor 0.80 | −4.983 | 6.202 | 40.0% | 0.488 | +2.954 | 4.990 | 55.0% | 0.297 |
| Physical factor 1.20 | −1.736 | 4.238 | 62.5% | 0.413 | +3.446 | 5.330 | 50.0% | 0.313 |

Bias and RMSE are in km s⁻¹ Mpc⁻¹. The “known-inclusion” arm uses the exact declared logistic inclusion curve for these bandwidth configurations. With 40 replicates, empirical coverage changes in 2.5 percentage-point increments and has broad uncertainty; these are descriptive screening numbers, not precise coverage estimates.

Heavy smoothing can make the catalogue-only posterior mean appear less biased at this one fixed truth, while host-density L1 error remains large. At factor 1.20, catalogue-only H0 bias is closer to zero but density L1 error is still 0.413; the known-inclusion arm's RMSE worsens to 5.330 compared with 3.896 at Scott default. This is bias compensation, not evidence that oversmoothing solves catalogue selection. It motivates reporting density diagnostics alongside H0 metrics and preregistering bandwidth rather than selecting it on a single H0 outcome.

### Grid and inclusion-law screen

Relative to the 161-point H0 / 801-point redshift default grid, posterior means changed by at most 0.000851 km s⁻¹ Mpc⁻¹ on the coarse 81 × 401 grid and 0.000103 on the fine 201 × 1001 grid, among paired saved rows. This is practical stability for this toy and tested range, not proof of asymptotic convergence.

| Assumed inclusion curve | Bias | RMSE | 68% coverage | Density L1 | Paired absolute-error difference vs oracle (95% interval) |
|---|---:|---:|---:|---:|---:|
| Midpoint 0.18 (truth 0.22) | +0.750 | 3.801 | 50.0% | 0.118 | −0.118 [−0.297, +0.064] |
| Midpoint 0.26 (truth 0.22) | +0.709 | 3.719 | 55.0% | 0.155 | −0.252 [−0.475, −0.043] |
| Width 0.10 (truth 0.05) | −2.536 | 4.245 | 60.0% | 0.259 | +0.293 [−0.714, +1.232] |
| Width 0.025 (truth 0.05) | +2.578 | 4.669 | 42.5% | 0.227 | +0.613 [+0.201, +1.021] |

A negative contrast means the assumed-selection arm had lower absolute H0 error than the oracle arm. The midpoint 0.26 condition looked better for that single error contrast, but that can arise from bias compensation; it is not evidence that the misspecified curve is preferable. At 40 replicates, all contrasts are screening estimates. The width 0.025 result is a candidate for a higher-replicate prospective follow-up, not a general conclusion.

### Decision

The screen supports three limited conclusions: (1) ignoring the imposed inclusion curve is severely harmful at the default bandwidth; (2) the true-inclusion control tracks the complete-parent oracle reasonably closely at the default and narrow bandwidths; and (3) bandwidth and selection-law assumptions change host-density and H0 summaries. It does **not** establish universal benefit from inverse weighting, realistic detector-selection validity, or real-H0 inference. A confirmatory follow-up must predeclare multiple H0 values and catalogue/noise regimes and set replicate count from paired variability.
