# Two-stream 200-replicate five-arm readout — 2026-11-19

## Status and scope

Two independently seeded 100-replicate runs completed successfully:

- Run [44460222957](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44460222957), master seed 2026111929.
- Run [44462487127](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44462487127), master seed 2026111930.

Each run used 100 replicates, three curvature settings (-1, 0, +1), validation sizes 100/300/500, N_parent=1200, and 30 detected events. Each run therefore planned 4,500 arm-cells and 300 child replicate-curvature runs.

Both artifacts passed an independent local audit:
- 4,500/4,500 lifecycle cells SUCCESS; 4,500/4,500 audit attestations PASS.
- 300/300 child audit reports PASS.
- 5,708 aggregate manifest entries verified per artifact.
- 4,500 posterior-summary rows, 45 aggregate metric rows, and 7,200 paired contrast rows per artifact.
- Planned/ledger identity, deterministic seed schedule, child/top-level posterior summaries, aggregate point estimates, and archived manifests reconciled.
- Both runs share the same frozen scientific design and have distinct master seeds; all 200 replicate seed pairs are unique.

The 200-replicate pooled readout is reproducible with `scripts/audit_and_analyze_pooled_five_arm.py`. It retains replicate-level pairing and uses Wilson 95% intervals for empirical coverage and Student-t 95% intervals for paired differences.

**Interpretation boundary:** this is a synthetic low-z Euclidean engineering benchmark, not an LVK/GWTC analysis, not an observational H0 measurement, and not a validation against a real galaxy survey. The runs use the original H0 prior support [50,90]. Do not combine these results with the earlier 10-replicate pilot: the seed-2026111929 stream overlaps that pilot.

## Main result at validation size n=500

Pooled results across 200 replicates per curvature:

| Curvature gamma | Arm | Bias | RMSE | 68% coverage | Wilson 95% interval |
|---:|---|---:|---:|---:|---:|
| -1 | A catalogue only | -15.319 | 15.506 | 0.000 | [0.000, 0.019] |
| -1 | B parent oracle | +0.516 | 3.793 | 0.690 | [0.623, 0.750] |
| -1 | C known-probability IPW | +0.052 | 5.109 | 0.530 | [0.461, 0.598] |
| -1 | D MAP plug-in IPW | -2.931 | 5.537 | 0.465 | [0.397, 0.534] |
| -1 | E uncertainty-integrated IPW | -2.720 | 5.370 | 0.490 | [0.422, 0.559] |
| 0 | A catalogue only | -14.639 | 14.866 | 0.000 | [0.000, 0.019] |
| 0 | B parent oracle | +0.516 | 3.793 | 0.690 | [0.623, 0.750] |
| 0 | C known-probability IPW | +0.237 | 4.532 | 0.600 | [0.531, 0.665] |
| 0 | D MAP plug-in IPW | +0.035 | 4.610 | 0.605 | [0.536, 0.670] |
| 0 | E uncertainty-integrated IPW | +0.148 | 4.504 | 0.630 | [0.561, 0.694] |
| +1 | A catalogue only | -13.480 | 13.795 | 0.015 | [0.005, 0.043] |
| +1 | B parent oracle | +0.516 | 3.793 | 0.690 | [0.623, 0.750] |
| +1 | C known-probability IPW | +0.558 | 4.354 | 0.610 | [0.541, 0.675] |
| +1 | D MAP plug-in IPW | +3.122 | 5.459 | 0.490 | [0.422, 0.559] |
| +1 | E uncertainty-integrated IPW | +3.053 | 5.342 | 0.535 | [0.466, 0.603] |

The oracle arm B is the best-performing reference in this design, as expected from its access to the parent population. It is a diagnostic upper-reference arm, not an implementable catalogue estimator. Known-probability IPW (C) does not match the oracle, and estimated-inclusion arms D/E are not uniformly better than C. In particular, at gamma=+1 both D and E have higher RMSE than C. At gamma=0, E is close to C on RMSE and has somewhat higher observed coverage, but the paired evidence should be read alongside the uncertainty intervals and other curvature settings.

### Paired contrasts at n=500

Differences are left arm minus right arm. Selected results:

| gamma | Contrast | Metric | Mean paired difference | 95% t interval |
|---:|---|---|---:|---:|
| -1 | E-D | absolute median error | -0.125 | [-0.193, -0.056] |
| -1 | E-D | 68% coverage | +0.025 | [-0.008, +0.058] |
| 0 | E-D | absolute median error | -0.061 | [-0.126, +0.004] |
| 0 | E-D | 68% coverage | +0.025 | [-0.001, +0.051] |
| +1 | E-D | absolute median error | -0.084 | [-0.154, -0.014] |
| +1 | E-D | 68% coverage | +0.045 | [+0.005, +0.085] |
| +1 | E-C | absolute median error | +0.798 | [+0.465, +1.132] |
| +1 | E-C | 68% coverage | -0.075 | [-0.160, +0.010] |

The E-vs-D absolute-error differences are modest and vary with curvature. The E-vs-C result at positive curvature is adverse. This is not evidence for universal superiority of uncertainty integration.

## Catalogue-only prior-boundary audit

For all 200 replicate populations and all three curvature values (600 replicate-curvature cases), the catalogue-only arm A likelihood was independently reconstructed from archived host-density arrays and observed events on the original [50,90] grid. The maximum absolute log-likelihood discrepancy was below 2.3e-13.

The same selected-event likelihood was then recomputed on an expanded [30,110] H0 grid with 801 points. Results:

| gamma | n | Mean posterior median | Mean q025 | Mean posterior mass [30,50] | Mean mass [30,55] | Mean mass [70,90] |
|---:|---:|---:|---:|---:|---:|---:|
| -1 | 200 | 52.712 | 42.998 | 0.321 | 0.659 | 0.002 |
| 0 | 200 | 53.274 | 42.548 | 0.295 | 0.603 | 0.004 |
| +1 | 200 | 54.352 | 42.377 | 0.257 | 0.520 | 0.010 |

On the expanded prior, the central 68% interval contains the true H0=70 in 0/200 cases for gamma=-1, 0/200 for gamma=0, and 3/200 for gamma=+1. Thus the original lower-bound truncation was a real reporting problem, but widening the prior does not rescue the catalogue-only inference: the posterior remains concentrated far below the generating value in this synthetic design.

The 801-point expanded-grid analysis is a reproducible diagnostic, but this 200-case pass has not yet had a full two-resolution grid-refinement audit. Do not describe the expanded-grid posterior quantiles as fully resolution-certified until that is completed.

## Limitations and next work

1. All 100-replicate runs used H0 prior [50,90]. The expanded-prior analysis was applied only to arm A after the fact. Do not compare expanded-prior A directly with the other arms' [50,90] posteriors as if they shared a prior.
2. The 200 replicates span two master-seed streams. Keep the pilot-overlap caveat: seed 2026111929 overlaps the earlier 10-replicate pilot and must not be added to that pilot as new independent evidence.
3. Coverage is empirical Monte Carlo coverage under one fixed H0 truth and this particular synthetic generator. It is not general SBC certification.
4. Next, complete the expanded-grid refinement audit; run a fully matched wider-prior five-arm experiment; inspect weight-tail and ESS diagnostics by curvature/validation size; and only then decide whether further coefficient-posterior integration work is justified.
5. No estimator ranking or real-data claim is warranted from this experiment.



## Weight concentration and effective sample size diagnostics (n=500)

The paired weight diagnostics reinforce the instability concern. C and D Kish ESS values are bounded by the observed catalogue size (566–660 across the 600 child runs; mean about 611), whereas B's oracle ESS is bounded by the full 1,200-parent population. E uses node-specific diagnostics over the observed catalogue, so its minimum node Kish ESS is not directly interchangeable with a single global Kish ESS.

| gamma | Arm | Mean Kish ESS (C/D) or mean minimum-node Kish ESS (E) | 10th percentile | Mean top-5% mass (or max node top-5% mass for E) |
|---:|---|---:|---:|---:|
| -1 | C known-probability IPW | 106.9 | 51.5 | 0.378 |
| -1 | D MAP plug-in IPW | 182.1 | 110.4 | 0.310 |
| -1 | E uncertainty-integrated IPW | 84.3 | 40.2 | 0.455 |
| 0 | C known-probability IPW | 136.8 | 102.9 | 0.359 |
| 0 | D MAP plug-in IPW | 144.7 | 85.2 | 0.355 |
| 0 | E uncertainty-integrated IPW | 63.3 | 32.0 | 0.518 |
| +1 | C known-probability IPW | 178.1 | 158.3 | 0.332 |
| +1 | D MAP plug-in IPW | 117.4 | 71.6 | 0.405 |
| +1 | E uncertainty-integrated IPW | 51.4 | 28.5 | 0.580 |

For E, the 25-node quadrature ESS averaged about 6.85/25 at n=500. Its node-wise maximum top-5%-weight mass averaged 45.5%, 51.8%, and 58.0% at gamma=-1, 0, +1, respectively. This is a material warning: integrating over inclusion-model uncertainty does not automatically make the weighted redshift representation stable. It also explains why E's posterior can be broader without being better calibrated. These diagnostics are descriptive; they are not by themselves proof that any single cell's posterior is invalid.


## Completed 200-case expanded-prior grid-refinement check

The catalogue-only expanded-prior calculation was repeated for all 200 replicates at both 561- and 1121-point grids over [30,110], for all three curvature settings (600 replicate-curvature cases). The maximum absolute differences across all cases were:
- posterior median: 0.000480;
- q025: 0.002119;
- posterior mass [30,50]: 0.000041.

The archived [50,90] likelihood reconstruction remained within 2.3e-13 maximum absolute log-likelihood error. The 561/1121-grid group means also agree to better than 0.00005 in posterior median. The expanded-prior catalogue-only diagnostic is therefore numerically stable at these two resolutions for the reported summaries. A wrapper, `scripts/audit_pooled_catalogue_boundary.py`, runs the committed single-artifact boundary auditor on both immutable artifacts, verifies design/seed separation, and pools the results.

For completeness, at gamma=+1 and n=500, the paired E-C squared-error difference is +9.576 with 95% paired t interval [+5.927, +13.226], consistent with E's higher RMSE in that setting.


## Exploratory weight-diagnostic association (unadjusted)

Across the 200 pooled replicates, Spearman correlations between absolute posterior-median error and weight concentration change sign with curvature for D and E. At n=500:
- D: rho(top-5% mass, absolute error) = -0.440 at gamma=-1, +0.021 at gamma=0, and +0.203 at gamma=+1.
- E: rho(max node top-5% mass, absolute error) = -0.417 at gamma=-1, +0.033 at gamma=0, and +0.188 at gamma=+1.
- C's corresponding correlations are close to zero (-0.089, -0.030, +0.012).

This sign change is consistent with the broader result that weight concentration is not a universal monotone predictor of posterior error across curvature regimes. These are unadjusted exploratory correlations over multiple arm/size/curvature combinations, not a validated diagnostic rule; no causal interpretation or threshold is claimed.


## Matched wider-prior all-arm follow-up (running)

The five-arm runner and both levels of auditors now accept explicit H0 prior bounds; the default remains [50,90] for backwards compatibility. The failure-inclusive RMSE and interval-width bounds are parameterized by the declared prior support. A regression test exercises the [30,110] path and its audit roundtrip, and another checks the wider failure-bound width.

- Five-replicate wide-prior smoke: [Actions run 44483782186](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44483782186).
- 100-replicate wide-prior stream, master seed 2026111940: [Actions run 44483923518](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44483923518). Its first five replicate seeds overlap the five-replicate smoke; do not pool the smoke separately.
- 100-replicate wide-prior stream, master seed 2026111929: [Actions run 44484212311](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44484212311). This uses the same replicate seed schedule as the original [50,90] run 44460222957 and is intended for paired prior-sensitivity comparisons.

At the time of writing, the regression tests in the first two workflows have passed, while the full wide-prior simulation and audit stages are still running. No results from these new runs are included in the 200-replicate tables above.


## Validation-size pattern

Across n=100, 300, and 500, the known-probability C arm is invariant as expected because its weights use the true inclusion probabilities rather than the validation fit. The estimated arms change with validation size, but more validation does not remove curvature misspecification:
- At gamma=0, D RMSE decreases 5.092 → 4.676 → 4.610, and E remains near C (4.572 → 4.488 → 4.504).
- At gamma=-1, D improves from 6.204 to about 5.5, while E ranges 5.496 → 5.271 → 5.370.
- At gamma=+1, D/E RMSE remain worse than C at all three validation sizes; E moves 5.061 → 5.306 → 5.342.

This is consistent with validation-sample uncertainty and functional-form misspecification being distinct problems: increasing n can reduce estimation noise but cannot make a linear fitted inclusion model represent the quadratic DGP exactly.
