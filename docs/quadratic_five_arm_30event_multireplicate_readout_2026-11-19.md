# Paired 10-replicate, 30-event, five-arm H0 pilot

**Workflow:** [44458318640](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44458318640)  
**Status:** COMPLETE; child and orchestration audits PASS.  
**Scope:** synthetic low-redshift Euclidean benchmark only; not confirmatory and not observational.

## Frozen design

- 10 replicate seeds, paired across curvature values `gamma = -1, 0, +1`.
- 1,200 parent galaxies and exactly 30 detected events per replicate.
- Fixed generating H0 = 70 km s^-1 Mpc^-1.
- Quadratic-logit catalogue inclusion with beta = -3.9; degree-one fitted inclusion model.
- Validation frame 500 with nested SRS sizes 100, 300, and 500.
- H0 prior [50,90], 281 grid points; redshift grid 501 points.
- Five arms: A catalogue-only; B full-parent oracle; C known-probability IPW; D MAP plug-in IPW; E uncertainty-integrated IPW.
- Five-node-per-dimension coefficient quadrature (25 nodes); Scott multiplier 0.7.
- Same parent population, GW event realization, validation frame, and sampled validation IDs are paired across curvature within replicate.

The Scott multiplier is a fixed engineering setting, not selected by a completed training/holdout bandwidth experiment. This pilot does not close that gate.

## Execution and integrity

- 30 replicate-by-curvature child runs completed.
- 450/450 planned replicate × curvature × validation-size × arm cells are SUCCESS.
- 450/450 cell audit attestations are PASS.
- 450 posterior-summary rows retained.
- 578 aggregate manifest entries verified.
- 30 cross-curvature paired-validation blocks checked.
- 720 paired-contrast rows and 45 aggregate metric rows independently reconstructed.
- The independent orchestration auditor returned PASS.

A first execution had one child auditor failure caused by the independent damped-Newton implementation stalling on a particular validation seed. The stored mode, nodes, and weights for that seed were independently reconstructed to numerical agreement. The auditor was then given a separate score-equation root-solver fallback, and a regression test was added for that exact seed. This final run includes all 10 replicates and passes the regression and artifact audits; the earlier audit-only failure is not counted as a model failure.

## Posterior-median error, RMSE and 68% coverage at n_validation = 500

Errors are posterior median minus the generating H0, in km s^-1 Mpc^-1. Coverage is the empirical fraction of the 10 replicates whose 68% equal-tail interval contains H0=70.

| gamma | Arm | Mean median error | RMSE | 68% coverage |
|---:|---|---:|---:|---:|
| -1 | A catalogue-only | -15.220 | 15.326 | 0.0 |
| -1 | B parent oracle | +2.097 | 4.390 | 0.7 |
| -1 | C known-probability IPW | +0.607 | 6.174 | 0.4 |
| -1 | D MAP plug-in IPW | -2.213 | 6.138 | 0.3 |
| -1 | E uncertainty-integrated IPW | -2.084 | 6.058 | 0.3 |
| 0 | A catalogue-only | -14.706 | 14.818 | 0.0 |
| 0 | B parent oracle | +2.097 | 4.390 | 0.7 |
| 0 | C known-probability IPW | +1.851 | 5.065 | 0.6 |
| 0 | D MAP plug-in IPW | +1.466 | 5.232 | 0.6 |
| 0 | E uncertainty-integrated IPW | +1.531 | 5.129 | 0.6 |
| +1 | A catalogue-only | -13.730 | 13.829 | 0.0 |
| +1 | B parent oracle | +2.097 | 4.390 | 0.7 |
| +1 | C known-probability IPW | +1.553 | 4.557 | 0.5 |
| +1 | D MAP plug-in IPW | +3.840 | 6.138 | 0.5 |
| +1 | E uncertainty-integrated IPW | +3.756 | 6.026 | 0.5 |

## Validation-size behaviour

A, B, and C are invariant across validation sizes by design: they do not use the fitted validation inclusion model. D and E change with validation size.

The integrated arm E has lower success-only RMSE than MAP plug-in D in all nine gamma-by-validation-size cells in this small pilot. It also has equal or higher empirical 68% coverage than D in all nine cells, generally with wider intervals. This is a useful signal that coefficient uncertainty can soften plug-in instability, but 10 replicates are far too few to establish a general improvement.

E does not consistently beat the known-probability arm C. At gamma=+1, for example, E has larger RMSE than C at all three validation sizes. The parent-oracle arm B remains the reference benchmark with the lowest RMSE in all cells. No estimator ranking or superiority claim is justified.

## Weight-concentration diagnostics

Median diagnostics across the 10 replicates:

- Known-probability IPW C: Kish ESS about 136–185; maximum normalized weight about 0.023–0.051; top-5% weight mass about 0.327–0.347.
- MAP plug-in IPW D: Kish ESS about 130–226; maximum normalized weight about 0.025–0.044; top-5% mass about 0.279–0.374.
- Integrated arm E: coefficient-quadrature ESS about 6.85–6.98 of 25 nodes. Minimum node-wise Kish ESS ranges roughly 22–97, and maximum node-wise top-5% mass roughly 0.44–0.69.

These concentration measures show that integrating over coefficient uncertainty does not eliminate unstable host weights. The coefficient quadrature ESS is an importance-weight concentration diagnostic, not a count of independent posterior draws.

## Interpretation and limits

1. Catalogue-only A is strongly biased low in every cell, showing the effect of catalogue selection in this synthetic setup.
2. Full-parent oracle B is the best reference arm, as expected; it is not an operational estimator because it uses parent-population information unavailable in a real incomplete catalogue.
3. Known-probability IPW C substantially reduces the catalogue-only bias, but weight concentration and RMSE remain material.
4. Estimated-probability arms D/E are sensitive to curvature and validation size. E is modestly more stable than D in this pilot, but neither uniformly improves on C.
5. The sample size is only 10 independent parent/event replicates per curvature. The empirical coverages have large Monte Carlo uncertainty; do not treat them as calibrated coverage estimates.

The experiment is not LVK/GWTC analysis, a real galaxy-survey validation, joint SBC, an observational H0 measurement, or evidence of universal estimator superiority. The next step is a larger repeated-seed study with a bandwidth frozen from independent training/holdout density recovery, explicit uncertainty intervals for paired contrasts, and a prespecified failure threshold.
