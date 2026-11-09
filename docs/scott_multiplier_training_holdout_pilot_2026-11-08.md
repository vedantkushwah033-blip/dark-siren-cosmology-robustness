# Training/held-out Scott-multiplier bandwidth pilot

> **Methodology status update (2026-11-09): SUPERSEDED PENDING RERUN.** The Scott-multiplier training/holdout pilot used the earlier reflected-sample covariance convention. Its bandwidth selection and held-out IAE values therefore do not validate the documented original-sample physical Scott bandwidth. Draft PR #39 corrects that convention; rerun the training/holdout selection with the corrected estimator and independent reconstruction before reusing the selected multiplier or any reported performance comparison. Results below are retained for provenance only.


**Date:** 2026-11-08  
**Status:** isolated local pilot; not a repository benchmark, not CI-verified, not confirmatory  
**Scope:** synthetic parent-host density recovery only; no H0 likelihood or real survey.

## Objective

Test the new `scott_multiplier` parameter without selecting a bandwidth on H0 posterior coverage. The bandwidth candidate is chosen on training simulation seeds using integrated absolute error (IAE) against the known analytic parent-host density, then evaluated on disjoint held-out seeds.

The pilot uses a common bandwidth multiplier across the oracle-parent, catalogue-only, and known-inclusion IPW arms. It does not tune each arm separately.

## Design

- Parent size: 4,000 galaxies.
- Training seeds: 11,000–11,049 (50 independent parent frames).
- Held-out seeds: 21,000–21,049 (50 disjoint parent frames).
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

The implementation was an isolated local mirror of the branch KDE algorithm and DGP, not an invocation of the repository's benchmark runner. It used the local environment's NumPy/SciPy versions and is not a substitute for repository CI. The expanded 50/50 run wrote separate training and held-out metric CSVs in the local working environment.

## Training selection

Equal-weight mean IAE pooled across the three arms selected Scott multiplier **0.15**:

| Scott multiplier | Pooled training mean IAE |
|---:|---:|
| 0.15 | 0.2179 |
| 0.40 | 0.2575 |
| 1.00 | 0.3487 |

Training mean IAE by arm:

| Arm | 0.15 | 0.40 | 1.00 |
|---|---:|---:|---:|
| Oracle parent KDE | 0.0438 | 0.0975 | 0.2674 |
| Known-inclusion IPW KDE | 0.0698 | 0.1333 | 0.2915 |
| Catalogue-only KDE | 0.5402 | 0.5416 | 0.4872 |

The oracle-only training objective also selected 0.15. The known-IPW arm selected 0.15. The catalogue-only arm instead preferred 1.00, illustrating why arm-specific bandwidth tuning would make the comparison unfair and could let smoothing partially mask selection bias.

## Held-out results

Mean IAE over the 50 held-out seeds:

| Arm | Scott × 0.15 (selected) | Scott × 0.40 | Scott × 1.00 |
|---|---:|---:|---:|
| Oracle parent KDE | 0.0435 | 0.0990 | 0.2681 |
| Known-inclusion IPW KDE | 0.0652 | 0.1310 | 0.2912 |
| Catalogue-only KDE | 0.5311 | 0.5324 | 0.4813 |

At the training-selected multiplier 0.15, the paired held-out known-IPW minus catalogue-only IAE difference was **-0.4659** (50 paired seeds; percentile-bootstrap 95% interval [-0.4712, -0.4605]). Known-IPW minus oracle IAE was +0.0217 (95% interval [0.0179, 0.0258]).

At Scott multiplier 1.00, known-IPW minus catalogue-only was -0.1900 (95% interval [-0.1940, -0.1861]). In this one controlled DGP, the common narrow multiplier substantially improved recovery for oracle and known-IPW densities, while catalogue-only remained biased by the redshift-dependent thinning.

## Interpretation

1. The explicit Scott-multiplier API behaves coherently in this isolated pilot: 0.15 is a multiplier of the weighted Scott factor, unlike the legacy absolute `bandwidth_factor`.
2. A training/held-out density-recovery protocol can select a bandwidth without using H0 coverage.
3. Under this declared DGP, known-inclusion IPW on held-out seeds recovers the analytic parent-host density much better than catalogue-only KDE at the same selected bandwidth.
4. Catalogue-only prefers a different smoothing level because its redshift distribution is selection-biased. Choosing its bandwidth independently would confound the arm comparison; the common-bandwidth protocol deliberately does not optimize away that bias.
5. This pilot does **not** include finite-validation estimated probabilities. It does not establish the project's central result.

## Limitations and next gate

Only 50 training and 50 held-out seeds, one parent law, one catalogue-inclusion law, one parent size, and three candidate multipliers were tested. The 201-point density grid is a pilot numerical resolution. The intervals quantify only seed variation within this small pilot and should not be treated as confirmatory uncertainty.

Before adopting a default bandwidth rule:
- increase training and held-out seed counts;
- include alternative parent intensities and catalogue-inclusion curves;
- test multiple parent sizes and retention levels;
- refine the multiplier grid and density integration grid;
- retain a single predeclared selection rule and freeze it before the confirmatory comparison;
- run the same protocol on estimated-inclusion arms only after the probability-stage model is properly validated.

**Decision:** the result supports proceeding with a larger training/held-out density-recovery study. It does not justify changing the production benchmark's bandwidth default yet.


## 7. Additional held-out parent-law stress

After selecting 0.15 on the original two-overdensity training DGP, I evaluated the frozen candidate on 50 new seeds (31,000–31,049) under a different parent intensity:
`0.8 + 3.0 exp(-0.5((z-0.18)/0.018)^2) + 0.5 exp(-0.5((z-0.32)/0.06)^2)`.
The inclusion curve, merger weight, n_parent=4,000, redshift support, grid and candidate multipliers were unchanged. This is a stress check, not another tuning set.

Mean IAE on the alternative parent law:

| Arm | Scott × 0.15 | Scott × 0.40 | Scott × 1.00 |
|---|---:|---:|---:|
| Oracle parent KDE | 0.0580 | 0.1415 | 0.2664 |
| Known-inclusion IPW KDE | 0.0867 | 0.1775 | 0.2854 |
| Catalogue-only KDE | 0.5061 | 0.5030 | 0.5128 |

At multiplier 0.15, the paired known-IPW minus catalogue-only IAE difference was -0.4194 (50 seeds; percentile-bootstrap 95% interval [-0.4258, -0.4127]); known-IPW minus oracle was +0.0287 (95% interval [0.0254, 0.0321]).

The frozen 0.15 candidate remained best among these three tested values for oracle and known-IPW density recovery under this alternative parent law. Catalogue-only again preferred a broader setting. This is encouraging but not a universal bandwidth result: only two parent laws, one inclusion curve, one parent size and three candidate multipliers have been examined.


## Repository-native corrected-bandwidth rerun protocol (pre-registration)

This section supersedes the isolated local numerical results above as the operational plan. It is a protocol, not a report of completed runs. The corrected reflected-KDE implementation is available on draft PR #39 at commit `a657333f47048bc9ace3bf273cb78beabe798c5c`; its exact-head Python 3.10/3.11/3.12 CI matrix passed. The rerun must use that implementation and its independent auditor.

### Phase A: frozen training and held-out partitions

- Invoke the repository-native `run_density_recovery_benchmark`; do not use the old local mirror.
- Parent size: 4,000; redshift support: [0.01, 0.40]; density grid: 801 points.
- DGP: the runner's declared two-overdensity parent intensity, merger weight, calibrated inclusion intercept, and gamma values -1.2, 0, +1.2 at retention 0.25.
- Nested validation sizes: 100, 300, 1,000; fitted inclusion degrees: 1 and 2.
- Candidate physical Scott multipliers: 0.15, 0.40, 1.00. Treat all as unselected candidates before training.
- Training root seeds: 11,000–11,049 (50 independent replicates).
- Held-out root seeds: 21,000–21,049 (50 disjoint independent replicates).
- Preserve the runner's SeedSequence stream protocol. The independent auditor must reconstruct the raw parent frames and all candidate densities from each partition.

### Phase B: common-bandwidth selection rule

Use a single common multiplier across all arms; never optimize a separate bandwidth for each arm. The predeclared arms are parent-sample KDE (finite-sample oracle diagnostic), catalogue-only KDE, known-inclusion IPW, estimated degree-1 IPW, and estimated degree-2 IPW. The parent-sample KDE is not an operational comparator.

For each candidate, first compute replicate-mean IAE within each prespecified (gamma, validation-size, arm) cell. Then average the cell means with equal cell weights. This prevents large or more frequently successful cells from silently dominating the selection objective. Report the same objective both including and excluding the parent-sample diagnostic, but designate the equal-cell mean over the four non-oracle arms (catalogue-only, known-inclusion, estimated degree 1, estimated degree 2) as the primary selection objective. The oracle-inclusive objective is sensitivity-only.

A candidate is eligible for selection only if every planned training cell has all 50 paired replicate attempts represented and evaluable for that candidate. If any attempt is FAILED, NON_EVALUABLE, or NOT_RUN, do not silently drop it or calculate a successful-only selection score: mark selection blocked, report the complete failure ledger, diagnose the cause, and rerun under a new documented protocol/version if necessary. The selected multiplier is the eligible candidate with the smallest primary equal-cell objective; ties within numerical tolerance are resolved in favour of the broader multiplier. Freeze the selected value and the full candidate table before opening held-out summaries.

### Phase C: held-out evaluation

Evaluate all three candidate multipliers on the untouched held-out seeds, while clearly labelling the training-selected multiplier. Do not reselect, retune, or change the candidate set using held-out outcomes. Report per-cell mean IAE, paired replicate contrasts, Monte Carlo standard errors, and paired intervals for:
- estimated degree-1 IPW minus catalogue-only;
- estimated degree-2 IPW minus catalogue-only;
- known-inclusion IPW minus catalogue-only;
- each estimated-IPW arm minus known-inclusion IPW.

Report all attempt states and planned denominators alongside conditional paired accuracy. A conditional paired contrast is not a replacement for the attempt-level failure rate. Preserve the disjoint-validation interpretation: increasing validation size also reduces the analysis-catalogue size, so this design combines validation information with fewer available analysis hosts. Do not describe its size trend as a pure validation-information effect.

### Phase D: independent audit and decision gate

The auditor must independently verify:
1. non-overlapping train/held-out seed ranges and the declared RNG streams;
2. complete planned design keys, raw parent partitions, nested validation ranks, and catalogue/validation disjointness;
3. truth probabilities, calibrated intercepts, fitted inclusion probabilities, and weight diagnostics;
4. corrected original-sample physical Scott bandwidth and every candidate density reconstruction;
5. attempt statuses, complete denominators, equal-cell training objective, selected multiplier, and frozen held-out evaluation;
6. artifact byte counts and SHA-256 hashes.

The auditor must not import the production runner, fitter, or KDE estimator. A PASS establishes artifact and arithmetic integrity only. It does not establish H0 posterior calibration, a causal effect of validation size, or superiority of estimated IPW.

**Stop rule:** do not launch a large confirmatory run or interpret the superseded positive-curvature failure regime until the corrected training/held-out density study passes its independent audit. Only then proceed to a separately designed H0-stage experiment with grid convergence and fixed-truth calibration/coverage audits.
