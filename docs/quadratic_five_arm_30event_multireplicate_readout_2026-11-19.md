# Quadratic five-arm 30-event repeated-seed pilot: complete readout

**Workflow:** [44458318640](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44458318640)  
**Executed commit:** `86e4b0688983c06b603146fc7eac2711d26d4cb9`  
**Status:** workflow successful; independent orchestration audit PASS.  
**Scope:** synthetic low-z Euclidean benchmark only; exploratory smoke, not confirmatory calibration.

## Frozen design

- 10 independent parent/event replicate seeds.
- Curvature gamma = -1, 0, +1; beta = -3.9.
- 1,200 parent galaxies; exactly 30 detected events per replicate.
- Fixed generating H0 = 70 km s^-1 Mpc^-1; posterior prior [50,90].
- Validation frame 500; nested SRS sizes 100, 300, 500.
- Five arms: A catalogue-only, B full-parent oracle, C known-probability IPW, D degree-one MAP plug-in IPW, E coefficient-posterior-integrated IPW.
- 281 H0 grid points; 501 redshift grid points; five-node-per-dimension coefficient quadrature.
- 3 × 3 × 5 × 10 = 450 planned arm-cells.

## Execution and independent orchestration audit

| Check | Result |
|---|---:|
| Planned cells | 450 |
| SUCCESS | 450 |
| FAILED / NON_EVALUABLE / NOT_RUN | 0 |
| Successful child runs | 30 |
| Aggregate metric rows reconstructed | 45 |
| Paired contrast rows reconstructed | 720 |
| Cross-curvature paired validation blocks | 30 |
| Aggregate manifest entries verified | 578 |
| Independent orchestration audit | PASS |

All workflow steps passed, including compilation, focused regression tests, the 30 child experiments, independent orchestration audit, denominator assertions, and artifact upload. The audit explicitly limits its claim to lifecycle, pairing, child-audit attestations, and aggregate arithmetic; it is not itself an independent calibration reconstruction. The artifact is available from the workflow run above.

## Results

Bias is posterior median minus true H0 (km s^-1 Mpc^-1); RMSE is across the 10 replicates; coverage is the empirical fraction of nominal 68% equal-tail intervals containing H0=70; width is mean 68% interval width. Each row has 10/10 successful replicates.

| gamma | n | Arm | Bias | RMSE | 68% coverage | Width |
|---:|---:|---|---:|---:|---:|---:|
| -1 | 100 | A catalogue-only | -15.220 | 15.326 | 0.0 | 6.548 |
| -1 | 100 | B parent oracle | +2.097 | 4.390 | 0.7 | 7.592 |
| -1 | 100 | C known-probability IPW | +0.607 | 6.174 | 0.4 | 7.150 |
| -1 | 100 | D MAP plug-in IPW | -1.975 | 6.978 | 0.2 | 7.387 |
| -1 | 100 | E uncertainty-integrated IPW | -1.804 | 6.507 | 0.3 | 8.417 |
| -1 | 300 | A catalogue-only | -15.220 | 15.326 | 0.0 | 6.548 |
| -1 | 300 | B parent oracle | +2.097 | 4.390 | 0.7 | 7.592 |
| -1 | 300 | C known-probability IPW | +0.607 | 6.174 | 0.4 | 7.150 |
| -1 | 300 | D MAP plug-in IPW | -2.548 | 6.162 | 0.2 | 7.490 |
| -1 | 300 | E uncertainty-integrated IPW | -2.306 | 6.089 | 0.3 | 7.946 |
| -1 | 500 | A catalogue-only | -15.220 | 15.326 | 0.0 | 6.548 |
| -1 | 500 | B parent oracle | +2.097 | 4.390 | 0.7 | 7.592 |
| -1 | 500 | C known-probability IPW | +0.607 | 6.174 | 0.4 | 7.150 |
| -1 | 500 | D MAP plug-in IPW | -2.213 | 6.138 | 0.3 | 7.439 |
| -1 | 500 | E uncertainty-integrated IPW | -2.084 | 6.058 | 0.3 | 7.743 |
| 0 | 100 | A catalogue-only | -14.706 | 14.818 | 0.0 | 7.194 |
| 0 | 100 | B parent oracle | +2.097 | 4.390 | 0.7 | 7.592 |
| 0 | 100 | C known-probability IPW | +1.851 | 5.065 | 0.6 | 7.649 |
| 0 | 100 | D MAP plug-in IPW | +1.448 | 6.116 | 0.4 | 7.681 |
| 0 | 100 | E uncertainty-integrated IPW | +1.440 | 5.467 | 0.7 | 9.035 |
| 0 | 300 | A catalogue-only | -14.706 | 14.818 | 0.0 | 7.194 |
| 0 | 300 | B parent oracle | +2.097 | 4.390 | 0.7 | 7.592 |
| 0 | 300 | C known-probability IPW | +1.851 | 5.065 | 0.6 | 7.649 |
| 0 | 300 | D MAP plug-in IPW | +1.303 | 5.122 | 0.6 | 7.776 |
| 0 | 300 | E uncertainty-integrated IPW | +1.438 | 5.028 | 0.6 | 8.376 |
| 0 | 500 | A catalogue-only | -14.706 | 14.818 | 0.0 | 7.194 |
| 0 | 500 | B parent oracle | +2.097 | 4.390 | 0.7 | 7.592 |
| 0 | 500 | C known-probability IPW | +1.851 | 5.065 | 0.6 | 7.649 |
| 0 | 500 | D MAP plug-in IPW | +1.466 | 5.232 | 0.6 | 7.705 |
| 0 | 500 | E uncertainty-integrated IPW | +1.531 | 5.129 | 0.6 | 8.099 |
| +1 | 100 | A catalogue-only | -13.730 | 13.829 | 0.0 | 8.352 |
| +1 | 100 | B parent oracle | +2.097 | 4.390 | 0.7 | 7.592 |
| +1 | 100 | C known-probability IPW | +1.553 | 4.557 | 0.5 | 7.588 |
| +1 | 100 | D MAP plug-in IPW | +3.591 | 6.311 | 0.6 | 7.428 |
| +1 | 100 | E uncertainty-integrated IPW | +3.364 | 5.856 | 0.7 | 8.528 |
| +1 | 300 | A catalogue-only | -13.730 | 13.829 | 0.0 | 8.352 |
| +1 | 300 | B parent oracle | +2.097 | 4.390 | 0.7 | 7.592 |
| +1 | 300 | C known-probability IPW | +1.553 | 4.557 | 0.5 | 7.588 |
| +1 | 300 | D MAP plug-in IPW | +3.765 | 6.036 | 0.4 | 7.384 |
| +1 | 300 | E uncertainty-integrated IPW | +3.672 | 5.914 | 0.4 | 7.878 |
| +1 | 500 | A catalogue-only | -13.730 | 13.829 | 0.0 | 8.352 |
| +1 | 500 | B parent oracle | +2.097 | 4.390 | 0.7 | 7.592 |
| +1 | 500 | C known-probability IPW | +1.553 | 4.557 | 0.5 | 7.588 |
| +1 | 500 | D MAP plug-in IPW | +3.840 | 6.138 | 0.5 | 7.354 |
| +1 | 500 | E uncertainty-integrated IPW | +3.756 | 6.026 | 0.5 | 7.684 |

## Interpretation

1. Catalogue-only A has a very large negative bias (-13.7 to -15.2) and zero observed coverage in all nine cells. This is a strong failure signal for this synthetic setup, not a precise estimate of population coverage.
2. Parent-oracle B is the best numerical reference in this pilot (RMSE 4.39; 7/10 coverage). It is an oracle diagnostic, not an operational catalogue estimator.
3. Known-probability IPW C materially improves over catalogue-only A, but remains less accurate than the parent oracle.
4. Estimated arms D/E are curvature-sensitive. For gamma=+1, MAP and integrated IPW are positively biased by about +3.6 to +3.8, with RMSE about 5.9–6.3.
5. E has lower RMSE than D in the nine cells and equal-or-higher empirical coverage in this pilot, often with wider intervals. This is hypothesis-generating only; it does not show universal improvement.
6. Increasing validation n from 100 to 500 does not produce a clear monotonic improvement in these summaries. The validation samples are nested, so comparisons must use paired replicate-level uncertainty.
7. A/B/C are invariant across validation sizes because those arms do not use the fitted validation inclusion model. Their repeated rows are not additional independent observations.

With only 10 independent replicates, coverage moves in 10-percentage-point increments and uncertainty is large. No estimator superiority, calibrated posterior, or observational claim follows. The fixed Scott multiplier 0.7 and 281-point H0 grid are engineering choices, not a completed training/holdout bandwidth selection or grid-convergence result.

## Required next steps

1. Run at least 100 independently seeded parent/event replicates per gamma × validation-size cell, retaining the complete lifecycle denominator.
2. Add paired replicate-level contrast intervals with parent-population clustering; do not treat nested validation sizes as independent.
3. Reconcile the H0 posterior outputs independently from raw seeds for a subset, beyond child-audit attestations and aggregate arithmetic.
4. Report replicate-level ESS, maximum normalized weight, top-1%/5% mass, support gaps, and inclusion-fit diagnostics jointly with H0 errors; test whether gamma=+1 bias tracks weight instability.
5. Complete H0-grid and KDE bandwidth holdout/convergence checks before any confirmatory preregistration.

No LVK/GWTC data, real galaxy survey, or observational H0 measurement is used.


## Additional raw-posterior boundary diagnostic (10-replicate artifact)

I downloaded and inspected the completed workflow artifact (23 MB) rather than relying only on the printed aggregate table. For each replicate/curvature/arm, I reconstructed the normalized discrete H0 posterior from the archived 281-point log-likelihood array using the same uniform [50,90] prior and trapezoidal normalization.

For the catalogue-only arm A at validation size 500:
- mean posterior mass in H0 ∈ [50,52] was 18.5% across 30 replicate-curvature combinations (median 16.9%; maximum 49.8%);
- mean posterior mass in H0 ∈ [50,55] was 48.3% (median 48.3%; maximum 89.4%);
- by gamma, mean mass in [50,55] was 55.2% at gamma=-1, 49.5% at gamma=0, and 40.1% at gamma=+1;
- the average 2.5th posterior percentile for A was about 50.4, very close to the lower prior boundary.

The other four arms had essentially no mass in [50,52] on average. This indicates that A's large negative median error is accompanied by substantial posterior pile-up near the lower H0 prior boundary, rather than merely a modest shift of an otherwise interior posterior.

This is an important diagnostic, but it is not a numerical-grid validation. The 281-point H0 grid and [50,90] prior truncate the posterior; the next analysis should repeat the A likelihood on an expanded H0 support and finer grids (for example [30,110] with 561/1121 points) to determine whether the catalogue-only posterior mode/median would move farther below 50. The current A bias and coverage results should therefore be interpreted as prior-boundary-truncated outcomes, not an unconstrained estimate of the full bias.


## Independent expanded-prior audit of catalogue-only A

The original [50,90] H0 prior causes substantial boundary truncation for catalogue-only A. To quantify it, I downloaded the complete 10-replicate artifact, reconstructed A's selected-event likelihood directly from the saved KDE density and GW event inputs, and verified it against every archived 281-point A likelihood curve. Maximum absolute log-likelihood discrepancy across all 30 replicate × curvature combinations was **1.71 × 10^-13**.

I then recomputed the A posterior using a uniform expanded H0 prior [30,110], at 561 and 1121 grid points. The maximum 561-versus-1121 differences across all 30 cases were:
- posterior median: **0.000119 km s^-1 Mpc^-1**;
- 2.5th percentile: **0.001963 km s^-1 Mpc^-1**;
- posterior mass in [30,50]: **0.0000265**.

Thus the expanded-grid result is numerically stable at these resolutions. Under the expanded [30,110] prior, the catalogue-only arm has:

| gamma | Mean posterior median | Mean q025 | Mean posterior mass in [30,50] | Mean posterior mass in [30,55] |
|---:|---:|---:|---:|---:|
| -1 | 53.302 | 43.699 | 26.4% | 64.5% |
| 0 | 53.657 | 43.223 | 25.7% | 60.0% |
| +1 | 54.679 | 43.285 | 21.3% | 51.6% |

The expanded-support median is about 1.5–1.6 km s^-1 Mpc^-1 lower than the original [50,90]-truncated median on average. The old q025 values near 50 were therefore boundary artifacts: with a wider support, the 2.5th percentile moves to about 43.2–43.7 and roughly one-fifth to one-quarter of posterior mass lies below H0=50.

This is an independent posterior-boundary diagnostic, not a new inference method. It materially strengthens the interpretation of catalogue-only A: its negative bias is not just a small shift, and the original prior support masked a substantial lower-H0 tail. The five-arm production design should either justify [50,90] as scientifically intended or use a support broad enough to contain the A posterior; any comparison of arms must use the same prior support.

Reproducible implementation: `scripts/audit_catalogue_only_prior_boundary.py`. The dedicated Actions workflow downloads the original 10-replicate artifact and reruns this audit. The dedicated Actions workflow [44464712323](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44464712323) subsequently reran the audit successfully: PASS, 30/30 replicate-curvature rows, with the same reconstruction and grid-refinement errors.
