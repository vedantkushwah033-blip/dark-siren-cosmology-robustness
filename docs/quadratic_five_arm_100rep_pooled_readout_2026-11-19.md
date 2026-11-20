# Pooled 200-replicate five-arm 30-event readout

**Status:** completed paired synthetic engineering pilot; descriptive only, not confirmatory and not observational.

## Source runs and independent validation

Two 100-replicate Actions runs completed successfully:

- [Run 44460222957](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44460222957), master seed 2026111929, artifact 14504568415.
- [Run 44462487127](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44462487127), master seed 2026111930, artifact 14506526731.

For each batch, the independent orchestration auditor reported PASS:
- 4,500/4,500 planned lifecycle cells SUCCESS and child-audit PASS;
- 300/300 replicate-by-curvature child runs independently audited;
- 5,708 aggregate manifest entries verified;
- 45 aggregate metric rows reconstructed;
- 7,200 existing paired-contrast rows reconstructed;
- no failed or non-evaluable planned cells.

The two frozen configurations match on all scientific settings checked. Their 100 analysis/validation seed pairs are disjoint (200 unique pairs; zero overlap). The earlier 10-replicate smoke is not counted again: its first ten seed pairs overlap run 44460222957.

The pooled analysis reconstructs the four prespecified contrasts E-D, E-C, D-C, and B-C at replicate level, including median error, absolute error, squared error, 68% coverage, and 68% interval width. It yields 45 cell rows and 180 paired-contrast metric rows. Coverage intervals are Wilson 95%; paired intervals are two-sided Student-t 95% intervals. Intervals are exploratory and are not multiplicity-adjusted.

## Frozen design

Both runs use the same synthetic design:
- 1,200 parent galaxies and 30 detected events; generating H0 = 70;
- curvature gamma = -1, 0, +1;
- validation sizes 100, 300, 500 from a 500-object validation frame;
- five matched arms A-E, sharing the parent/event realization within replicate;
- H0 prior [50,90] on a 281-point grid;
- sigma = 140 Mpc, detection threshold = 250 Mpc;
- Scott bandwidth multiplier 0.7;
- coefficient quadrature order 5 (25 tensor-product nodes).

The independent replicate unit is the parent/event replicate. Gamma conditions share parent/event seeds, and validation sizes are nested. Do not treat arm rows, gamma rows, or validation-size rows as independent replicates.

## Pooled results at validation size 500

All entries use 200 independent parent/event replicates per gamma. Coverage intervals are 95% Wilson intervals.

| gamma | Arm | Bias | RMSE | 68% coverage (95% Wilson CI) | Mean 68% width |
|---:|---|---:|---:|---:|---:|
| -1 | A catalogue-only | -15.319 | 15.506 | 0.000 (0.000, 0.019) | 6.063 |
| -1 | B parent oracle | 0.516 | 3.793 | 0.690 (0.623, 0.750) | 7.246 |
| -1 | C known-pi IPW | 0.052 | 5.109 | 0.530 (0.461, 0.598) | 7.046 |
| -1 | D MAP IPW | -2.931 | 5.537 | 0.465 (0.397, 0.534) | 7.303 |
| -1 | E integrated IPW | -2.720 | 5.370 | 0.490 (0.422, 0.559) | 7.647 |
| 0 | A catalogue-only | -14.639 | 14.866 | 0.000 (0.000, 0.019) | 6.728 |
| 0 | B parent oracle | 0.516 | 3.793 | 0.690 (0.623, 0.750) | 7.246 |
| 0 | C known-pi IPW | 0.237 | 4.532 | 0.600 (0.531, 0.665) | 7.172 |
| 0 | D MAP IPW | 0.035 | 4.610 | 0.605 (0.536, 0.670) | 7.209 |
| 0 | E integrated IPW | 0.148 | 4.504 | 0.630 (0.561, 0.694) | 7.538 |
| +1 | A catalogue-only | -13.480 | 13.795 | 0.015 (0.005, 0.043) | 7.728 |
| +1 | B parent oracle | 0.516 | 3.793 | 0.690 (0.623, 0.750) | 7.246 |
| +1 | C known-pi IPW | 0.558 | 4.354 | 0.610 (0.541, 0.675) | 7.237 |
| +1 | D MAP IPW | 3.122 | 5.459 | 0.490 (0.422, 0.559) | 6.999 |
| +1 | E integrated IPW | 3.053 | 5.342 | 0.535 (0.466, 0.603) | 7.287 |

B is an oracle/reference arm, not an implementable estimator. Its relatively good performance is a diagnostic benchmark only.

## Paired squared-error contrasts at validation size 500

The entries below are mean paired differences in squared posterior-median error, with two-sided 95% t intervals. A negative value favours the left-hand arm on squared error; a positive value favours the reference arm. These are exploratory intervals, without family-wise multiplicity adjustment.

| gamma | Contrast | Mean squared-error difference | 95% paired t interval |
|---:|---|---:|---:|
| -1 | E - D | -1.82 | (-2.53, -1.10) |
| -1 | E - C | +2.74 | (-1.51, +6.99) |
| -1 | D - C | +4.56 | (+0.03, +9.08) |
| -1 | B - C | -11.71 | (-15.45, -7.98) |
| 0 | E - D | -0.96 | (-1.50, -0.43) |
| 0 | E - C | -0.26 | (-2.04, +1.53) |
| 0 | D - C | +0.71 | (-1.09, +2.51) |
| 0 | B - C | -6.16 | (-8.81, -3.50) |
| +1 | E - D | -1.26 | (-1.98, -0.54) |
| +1 | E - C | +9.58 | (+5.93, +13.23) |
| +1 | D - C | +10.84 | (+7.01, +14.67) |
| +1 | B - C | -4.57 | (-6.85, -2.29) |

Interpretation:
- Integrating coefficient uncertainty (E) modestly reduces squared error relative to MAP plug-in (D) in these three curvature cells, but the effect is small compared with the curvature-dependent changes in performance.
- Relative to known-probability IPW (C), E is not uniformly better. At gamma = +1, E has a materially higher squared error; the paired 95% interval is entirely positive.
- The estimated-probability arms D/E retain a curvature-dependent bias pattern. Increasing validation size from 100 to 500 does not remove that pattern.
- No universal estimator ranking follows from this one synthetic design.

## Weight-concentration diagnostics

At gamma = +1 and validation size 500, averaged over the 200 replicate rows:
- C known-pi IPW: mean Kish ESS 178.1; top 5% weight mass 0.332.
- D MAP IPW: mean Kish ESS 117.4; top 5% weight mass 0.405.
- E integrated IPW: mean minimum node-wise Kish ESS 51.4; mean maximum node-wise top-5% mass 0.580; 25 quadrature nodes.

For E, the standard single-weight Kish columns are not applicable and are blank; the reported node-wise diagnostics are the relevant ones. These values indicate appreciable weight concentration, especially across E's coefficient nodes. They should be interpreted jointly with the posterior errors, not as standalone proof of failure.

## Prior-boundary limitation

The H0 prior [50,90] and 281-point grid were held fixed in this experiment. For catalogue-only A at n=500, the mean q025 is 50.47, 50.57, and 50.79 for gamma -1, 0, +1; respectively 73.5%, 68.0%, and 55.5% of replicate q025 values are below 50.5. This is direct evidence that A's lower posterior tail is pressed against the lower prior boundary in this experiment.

A separate expanded-support reconstruction of the earlier 10-replicate smoke found substantial posterior mass below 50 under [30,110]. That result is a boundary diagnostic on the earlier small sample, not a replacement for rerunning the present 200-replicate comparison on a wider common prior. The current table therefore describes the frozen [50,90]-truncated experiment only. Do not interpret A's bias as an unconstrained posterior result, and do not compare arms under different prior supports.

## Limitations and next gates

1. Repeat all five arms on a scientifically justified wider common H0 prior and finer grid; preserve identical support across arms.
2. Perform full-sample H0-grid and KDE-bandwidth convergence audits. The current 281-point grid and Scott multiplier 0.7 are engineering settings, not validated optima.
3. Report replicate-level diagnostics and paired uncertainty in any later confirmatory analysis; account for multiplicity if making formal comparisons.
4. Keep the synthetic low-z Euclidean scope. These results use no LVK/GWTC public data, no real galaxy survey, and are not an observational H0 measurement.

## Reproducibility

The source artifacts are Actions artifacts 14504568415 and 14506526731. The pooling utility is `scripts/analyze_pooled_quadratic_five_arm_precision.py`; it checks the frozen configuration, source audit status, planned denominators, and disjoint seed pairs before producing pooled cell and paired-contrast CSVs.


## Exploratory diagnostic-to-error association

I also joined the archived replicate-level weight diagnostics to the pooled posterior summaries and computed within-curvature Spearman correlations between top-5% normalized weight mass and absolute posterior-median error. These are exploratory, unadjusted associations at the replicate level, not a validated prediction model.

| gamma | Validation n | D MAP IPW: rho(top-5% mass, |error|) | E integrated IPW: rho(top-5% mass, |error|) |
|---:|---:|---:|---:|
| -1 | 100 | -0.46 | -0.37 |
| -1 | 300 | -0.44 | -0.40 |
| -1 | 500 | -0.44 | -0.42 |
| 0 | 100 | -0.16 | -0.10 |
| 0 | 300 | -0.01 | +0.03 |
| 0 | 500 | +0.02 | +0.03 |
| +1 | 100 | +0.24 | +0.15 |
| +1 | 300 | +0.27 | +0.22 |
| +1 | 500 | +0.20 | +0.19 |

The sign reverses between gamma = -1 and gamma = +1: stronger top-weight concentration is associated with *smaller* absolute error in the negative-curvature cells, but with *larger* absolute error in the positive-curvature cells. At gamma = 0 the association is weak or absent. This is a warning against treating ESS or top-weight mass as a universal monotone failure score. Curvature-dependent bias and the direction of misspecification matter.

The p-values were not multiplicity-adjusted and are not evidence of causal effects. The next useful test is to predefine diagnostic thresholds on a training set and evaluate whether they identify posterior error/coverage failures on independent held-out replicate seeds. Until then, these diagnostics are descriptive only.
