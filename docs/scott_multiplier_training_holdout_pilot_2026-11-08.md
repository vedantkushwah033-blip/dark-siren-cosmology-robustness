# Training/held-out Scott-multiplier bandwidth pilot

**Date:** 2026-11-08  
**Status:** isolated local pilot; not a repository benchmark, not CI-verified, not confirmatory  
**Scope:** synthetic parent-host density recovery only; no H0 likelihood or real survey.

## Objective

Test the new `scott_multiplier` parameter without selecting a bandwidth on H0 posterior coverage. The bandwidth candidate is chosen on training simulation seeds using integrated absolute error (IAE) against the known analytic parent-host density, then evaluated on disjoint held-out seeds.

The pilot uses a common bandwidth multiplier across the oracle-parent, catalogue-only, and known-inclusion IPW arms. It does not tune each arm separately.

## Design

- Parent size: 4,000 galaxies.
- Training seeds: 11,000–11,009 (10 independent parent frames).
- Held-out seeds: 21,000–21,009 (10 disjoint parent frames).
- Redshift support: [0.01, 0.40].
- Density evaluation grid: 201 points; parent sampling CDF evaluated on a 20,001-point grid.
- Parent intensity:
  `1 + 2.2 exp(-0.5((z-0.105)/0.025)^2) + 1.6 exp(-0.5((z-0.275)/0.040)^2)`.
- Merger weight: `1 + 1.5 (z-0.01)/0.39`.
- Catalogue inclusion: `0.10 + 0.90/[1 + exp((z-0.22)/0.05)]`.
- Target: normalized analytic parent intensity multiplied by the declared merger weight.
- Arms: parent-sample KDE (oracle diagnostic), unweighted retained-catalogue KDE, and known-inclusion IPW KDE.
- Candidate Scott multipliers: 0.15, 0.40, 1.00.
- Score: integrated absolute error between estimated and analytic host density.

The implementation was an isolated local mirror of the branch KDE algorithm and DGP, not an invocation of the repository's benchmark runner. It used the local environment's NumPy/SciPy versions and is not a substitute for repository CI.

## Training selection

Equal-weight mean IAE pooled across the three arms selected Scott multiplier **0.15**:

| Scott multiplier | Pooled training mean IAE |
|---:|---:|
| 0.15 | 0.2197 |
| 0.40 | 0.2572 |
| 1.00 | 0.3486 |

The oracle-only training objective also selected 0.15. The known-IPW arm selected 0.15. The catalogue-only arm instead preferred 1.00, illustrating why arm-specific bandwidth tuning would make the comparison unfair and could let smoothing partially mask selection bias.

## Held-out results

Mean IAE over the 10 held-out seeds:

| Arm | Scott × 0.15 (selected) | Scott × 0.40 | Scott × 1.00 |
|---|---:|---:|---:|
| Oracle parent KDE | 0.0387 | 0.0942 | 0.2674 |
| Known-inclusion IPW KDE | 0.0622 | 0.1279 | 0.2908 |
| Catalogue-only KDE | 0.5342 | 0.5364 | 0.4846 |

At the training-selected multiplier 0.15, the paired held-out known-IPW minus catalogue-only IAE difference was **-0.4720** (10 paired seeds; percentile-bootstrap 95% interval [-0.4834, -0.4592]). Known-IPW minus oracle IAE was +0.0235 (95% interval [0.0200, 0.0274]).

At Scott multiplier 1.00, known-IPW minus catalogue-only was -0.1937 (95% interval [-0.2010, -0.1863]). In this one controlled DGP, the common narrow multiplier substantially improved recovery for oracle and known-IPW densities, while catalogue-only remained biased by the redshift-dependent thinning.

## Interpretation

1. The explicit Scott-multiplier API behaves coherently in this isolated pilot: 0.15 is a multiplier of the weighted Scott factor, unlike the legacy absolute `bandwidth_factor`.
2. A training/held-out density-recovery protocol can select a bandwidth without using H0 coverage.
3. Under this declared DGP, known-inclusion IPW on held-out seeds recovers the analytic parent-host density much better than catalogue-only KDE at the same selected bandwidth.
4. Catalogue-only prefers a different smoothing level because its redshift distribution is selection-biased. Choosing its bandwidth independently would confound the arm comparison; the common-bandwidth protocol deliberately does not optimize away that bias.
5. This pilot does **not** include finite-validation estimated probabilities. It does not establish the project's central result.

## Limitations and next gate

Only 10 training and 10 held-out seeds, one parent law, one catalogue-inclusion law, one parent size, and three candidate multipliers were tested. The 201-point density grid is a pilot numerical resolution. The intervals quantify only seed variation within this small pilot and should not be treated as confirmatory uncertainty.

Before adopting a default bandwidth rule:
- increase training and held-out seed counts;
- include alternative parent intensities and catalogue-inclusion curves;
- test multiple parent sizes and retention levels;
- refine the multiplier grid and density integration grid;
- retain a single predeclared selection rule and freeze it before the confirmatory comparison;
- run the same protocol on estimated-inclusion arms only after the probability-stage model is properly validated.

**Decision:** the result supports proceeding with a larger training/held-out density-recovery study. It does not justify changing the production benchmark's bandwidth default yet.
