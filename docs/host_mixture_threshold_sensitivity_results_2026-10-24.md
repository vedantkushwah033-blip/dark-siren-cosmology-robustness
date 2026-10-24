# Host-mixture and detection-threshold sensitivity results

**Status: completed and independently audited (2026-10-24).** This report describes a controlled synthetic diagnostic, not an analysis of real gravitational-wave observations or a real galaxy catalogue.

## Reproducibility and audit

- Prespecified design: [protocol](host_mixture_threshold_sensitivity_protocol.md).
- Workflow: [successful GitHub Actions run #40101614384](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40101614384).
- Code revision: `f699f631d5a4184b227e26a8c7ca812a33feebd5`.
- Artifact: `host-mixture-threshold-sensitivity-results`; SHA-256 digest: `75f747b109ab6f321d58e1dc94efb2f155accffc79a30641f849480bc6c3f9a8`.
- Independent audit status: **PASS**. It rechecked 32,400 raw rows across 108 cells, 756 summary rows, 3,024 method/cell metric groups, 216 paired contrasts and their 216 bootstrap intervals. Detection bounds, finite values, error arithmetic, summary metrics, paired means and bootstrap endpoints were recomputed from the archived raw rows.
- Full test suite: successful on Python 3.10, 3.11 and 3.12 in [run #40101614170](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40101614170).

The manifest records seed 20261031, 300 replicates per cell, 5,000 paired bootstrap draws, Gaussian distance noise with standard deviation 20 Mpc, inference prior H0 in [50, 90] km s^-1 Mpc^-1 and a 4,001-point grid on [25, 125]. The crossed design varies four fixed truths (60, 70, 80, 88), three host-redshift pairs, three intrinsic-weight pairs and three upper detection thresholds (90, 120, 150 Mpc).

## Equal-weight design averages

The table averages cell-level metrics equally across the 108 prespecified design cells. These are descriptive averages over this deliberately balanced design, **not** population-weighted performance estimates or an estimate of real-world H0.

| Inference method | Mean bias | Mean absolute error | RMSE | 68% coverage | 95% coverage |
|---|---:|---:|---:|---:|---:|
| Two-host, uncorrected | -0.849 | 8.353 | 9.557 | 61.33% | 98.27% |
| Two-host, selection-corrected | -3.738 | 8.669 | 10.158 | 65.69% | 97.41% |
| Collapsed, uncorrected | +4.588 | 9.686 | 10.614 | 47.54% | 82.31% |
| Collapsed, selection-corrected | +0.871 | 8.994 | 10.607 | 58.56% | 92.09% |

The nominal 68% and 95% credible-interval coverage rates are empirical repeated-sampling fractions under this simulator. In particular, the 68% intervals under-cover on average; a selection-normalization term does not by itself guarantee calibrated intervals. Mean bias can also conceal strong dependence on the fixed truth and design cell.

## Primary paired contrast

For each generated observation, define
**Delta = |error(collapsed)| - |error(two-host)|**.
Positive values favor the two-host likelihood on absolute error; negative values favor the collapsed approximation. Identical observations are used for both models. The table averages cell-level paired contrasts across host redshifts, intrinsic weights and thresholds; interval endpoints are the averages of the corresponding cell-level percentile-bootstrap endpoints, so they summarize heterogeneity rather than form one pooled confidence interval.

| Fixed true H0 | Uncorrected mean Delta | Average cell 95% bootstrap limits | Corrected mean Delta | Average cell 95% bootstrap limits |
|---:|---:|---:|---:|---:|
| 60 | +6.364 | [+6.104, +6.617] | +4.919 | [+4.654, +5.174] |
| 70 | +4.488 | [+4.186, +4.781] | +2.518 | [+2.110, +2.914] |
| 80 | -0.919 | [-1.369, -0.473] | -1.984 | [-2.406, -1.561] |
| 88 | -4.603 | [-4.901, -4.304] | -4.150 | [-4.421, -3.874] |

Across all 108 cells, the equal-cell average Delta is +1.333 without selection correction and +0.326 with correction. These grand averages are not a sufficient model-selection rule: the contrast changes sign with truth. For H0 = 60 and 70, the two-host model has lower absolute error on average; for H0 = 80 and 88, the collapsed model has lower absolute error on average. Selection correction changes the magnitude and does not eliminate this truth dependence.

## Threshold dependence

Equal-cell averages within each threshold reinforce that the results depend on the detection setup. For example, two-host uncorrected MAE is 8.579, 8.217 and 8.264 for thresholds 90, 120 and 150 Mpc, respectively; corresponding collapsed-uncorrected MAE is 10.261, 9.255 and 9.542. The corrected methods also vary across thresholds. These summaries are averaged over different truths, redshift separations and intrinsic weights, so they should not be interpreted as an isolated causal effect of changing the threshold.

## Host fractions and numerical grid check

The simulation records both intrinsic host weights and host fractions among detected observations, and compares the latter with analytic expectations under its Gaussian-noise and hard-truncation assumptions. Across 108 design cells, the mean absolute difference between observed and expected host fractions was 0.00763 (about 0.76 percentage points); the largest cell-level difference was 0.04310 (about 4.31 percentage points), consistent with finite-sample variation at 300 replicates per cell. Mean absolute differences were 0.00369, 0.00789 and 0.01131 for upper thresholds 90, 120 and 150 Mpc, respectively. Intrinsic weights are not the same as detected-sample weights: detection changes the composition, sometimes substantially. The archived `summary.csv` contains the per-cell observed and expected fractions for inspection.

A prespecified subset (H0 = 70, redshifts 0.02/0.04, intrinsic weights 0.7/0.3, thresholds 90/120/150 Mpc) was recomputed on 2,001- and 4,001-point grids using the same observations. Across its 12 method/threshold rows, the largest absolute median difference was about 0.000092 in H0 units; the largest absolute difference among the reported 68%/95% interval endpoints was about 0.000107. This supports numerical stability for this subset only, not universal grid convergence across all possible data or posterior shapes.

## Interpretation and limitations

1. This is a synthetic stress test using a Gaussian distance proxy, low-redshift distance relation, two discrete hosts, hard detection truncation and fixed parameter values. It is not a detector-level likelihood, real survey catalogue or LVK measurement.
2. The collapsed effective redshift is a weighted mean and is not generally likelihood-equivalent to a host mixture.
3. The selection-normalization correction and catalogue incompleteness are distinct mechanisms. This experiment varies the threshold but does not simulate missing galaxies or a realistic detector selection function.
4. The bounded prior [50, 90] can affect fixed-truth coverage, especially near the boundaries. Average coverage across these selected truths is not interchangeable with prior-predictive simulation-based calibration.
5. Each cell has 300 replicates; bootstrap intervals quantify Monte Carlo uncertainty for paired contrasts under this simulator and do not establish astrophysical significance. Many cells are examined, so individual interval exclusion should not be read as a multiplicity-adjusted discovery.
6. The results support a limited methodological conclusion: **the relative performance of the two-host and collapsed models is truth- and design-dependent; neither is a universal winner in this setup.** Further work would need realistic selection and catalogue models, external validation and a full observational likelihood before making claims about real dark-siren inference.

## Artifact files

The workflow artifact contains `per_run.csv`, `summary.csv`, `experiment_manifest.json`, `independent_audit_report.json` and `grid_sensitivity.csv`. The artifact and run link above are the source of the numerical tables in this report.
