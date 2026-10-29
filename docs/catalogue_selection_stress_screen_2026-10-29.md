# Catalogue-selection stress screen — 29 October 2026

## Status and scope

**Status: all three stress-screen jobs and the exact-head Python 3.10/3.11/3.12 test matrix passed.** This report summarizes the three new 50-replicate stress configurations generated from code commit `4a95ff40b68849e0ff33c9cae8916a3e96eeaa7b`. These are small, exploratory stress screens; they do not replace the locked 500-replicate Stage B experiment and do not validate a real gravitational-wave analysis.

- [Exact-head test matrix #40811163394](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40811163394): Python 3.10, 3.11, and 3.12 all passed.
- [Selection-misspecification pilot #40811157211](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40811157211): passed.
- [Catalogue-selection stress screen #40811159599](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40811159599): all three jobs passed.
- Every downloaded artifact's raw CSV and summary JSON matched the manifest byte counts and SHA-256 digests. All three stress configurations used source commit `4a95ff40...`.

## Design

Each configuration contains 50 paired replicates and six inference arms (300 fit rows). Each replicate uses 4,000 parent galaxies, 30 detected events, fixed true H0 = 70 km s⁻¹ Mpc⁻¹, an 81-point H0 grid over [50, 90], a 401-point redshift grid over [0.01, 0.40], and the same distance/noise setup as the main selection-misspecification pilot. Seeds are disjoint across the three screens:

| Screen | True inclusion floor | Midpoint | Width | Seed range |
|---|---:|---:|---:|---:|
| Low inclusion floor | 0.05 | 0.22 | 0.05 | 20270000–20270049 |
| High inclusion floor | 0.25 | 0.22 | 0.05 | 20270100–20270149 |
| Broad transition | 0.10 | 0.22 | 0.08 | 20270200–20270249 |

The paired primary contrast was recomputed from raw replicate rows as
`D = |posterior_mean_catalogue_only - H0_true| - |posterior_mean_correct_weighting - H0_true|`.
Positive values favor correct weighting on absolute posterior-mean error. Percentile-bootstrap intervals use 10,000 paired-resample means and RNG seed 20271000. The analysis was independently recomputed from downloaded CSVs; this is a separate analysis script, not a claim that the existing per-row posterior auditor was rerun on these artifacts.

## Main results

| True inclusion setting | Catalogue-only bias | Catalogue-only RMSE | Catalogue-only 68% coverage | Correct-weighting bias | Correct-weighting RMSE | Correct-weighting 68% coverage | Mean paired D (km s⁻¹ Mpc⁻¹) | Bootstrap 95% interval | Weighting lower absolute error |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Floor 0.05, width 0.05 | −11.096 | 11.470 | 2% (1/50) | +0.830 | 4.301 | 66% (33/50) | +7.944 | [6.462, 9.240] | 44/50 (88%) |
| Floor 0.25, width 0.05 | −4.913 | 5.901 | 38% (19/50) | +0.761 | 3.677 | 72% (36/50) | +2.280 | [1.089, 3.363] | 37/50 (74%) |
| Floor 0.10, width 0.08 | −5.928 | 6.740 | 24% (12/50) | +0.155 | 3.641 | 74% (37/50) | +3.486 | [2.355, 4.491] | 42/50 (84%) |

The corresponding analytic-truth and oracle-KDE controls had RMSEs of 4.113 and 4.155 at floor 0.05; 3.581 and 3.666 at floor 0.25; and 3.726 and 3.710 at width 0.08. The shallow-slope arm's RMSEs were 4.655, 3.544, and 3.841, respectively. The shifted-transition arm's RMSEs were 4.199, 3.697, and 3.612. These small screens show that selection mismatch changes the size of catalogue-only error, while the correct-weighting arm substantially reduces it in these three chosen synthetic settings. They do **not** establish universal superiority: each cell has only 50 replicates, interval-coverage estimates are noisy, and a misspecified arm can sometimes perform similarly or better on a particular finite-sample metric due to bias compensation.

## Interpretation and limitations

1. The main pattern is directionally consistent with the locked Stage B experiment: omitting the known catalogue-inclusion correction can produce substantial downward H0 bias in this synthetic setup.
2. The size of the error depends on the inclusion curve. It is largest for the floor-0.05 screen and smaller for the floor-0.25 screen; the broad-transition screen remains materially biased when correction is omitted.
3. Correct weighting is not identical to the oracle/analytic truth controls. It estimates a population density from a finite catalogue and can have higher RMSE than the analytic truth control.
4. A 50-replicate screen is not a precision calibration study. Do not interpret a displayed 68% coverage rate as established calibration.
5. The fixed H0, simplified host population, toy selection/noise model, and lack of real detector response, event-rate marginalization, uncertain redshifts, and external likelihood remain major scope limits.

**Decision:** these stress screens strengthen the internal synthetic robustness evidence and are suitable as clearly labeled supplementary results. They do not change the project release gate: no real H0 measurement, no LVK validation, no claim of universal inverse-weighting performance, and no external CHIMERA benchmark until exact compatible v2.2 inputs are obtained and run with provenance.

## Artifact provenance

- Low-floor artifact: [12870297736](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40811159599)
- High-floor artifact: [12870541126](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40811159599)
- Broad-transition artifact: [12870855101](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40811159599)

The artifacts are available for the repository's configured retention period. Their manifests report 300 rows each, and each raw CSV and summary JSON passed byte-count and SHA-256 verification after download.
