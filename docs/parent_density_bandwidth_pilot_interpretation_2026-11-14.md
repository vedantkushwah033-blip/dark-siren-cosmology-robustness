# Parent-density bandwidth pilot: audited result and interpretation

**Date:** 2026-11-14  
**Status:** source and selection audits pass; bandwidth optimum unresolved  
**Scope:** synthetic parent-host density recovery only. This report makes no claim about H0 calibration, estimator ranking, or observational data.

## 1. Primary 200-replicate experiment

Workflow run: [43548695692](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/43548695692)  
Artifact ID: `14095293801`  
Artifact SHA-256: `d35d7081882ccb0b987230eb6115acec5449ba6fa1a0318f96189e87056153fe`  
Source commit: `9f077d151a8eb0a875e44e3505bf31b464aed7bd`

Design:
- 200 independent parent-population replicates, 4,000 parent units each;
- 801-point redshift grid over [0.01, 0.40];
- analytic target proportional to the declared two-component parent intensity times the common merger weight;
- nine weighted-Scott multipliers: 0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30, 0.40;
- deterministic 140/60 training/holdout split;
- selection objective: mean integrated absolute error (IAE) against the analytic target.

### 1.1 Independent source audit

The source auditor independently reconstructed all 1,800 replicate × multiplier attempts and verified all 800,000 parent units. The maximum absolute difference between recorded and independently reconstructed IAE was:

`4.163336342344337e-16`

Source audit status: **PASS**.

### 1.2 Train/holdout results

| Scott multiplier | Training mean IAE (n=140) | Holdout mean IAE (n=60) |
|---:|---:|---:|
| 0.05 | 0.1305253 | 0.1307963 |
| 0.08 | 0.1031205 | 0.1032066 |
| 0.10 | 0.0921284 | 0.0922072 |
| 0.12 | 0.0839649 | 0.0840997 |
| 0.15 | 0.0749159 | 0.0750892 |
| 0.20 | 0.0646492 | 0.0647632 |
| 0.25 | 0.0577219 | 0.0576946 |
| 0.30 | 0.0528231 | 0.0526713 |
| 0.40 | 0.0467434 | 0.0466049 |

The independently audited training selector returned `0.40`, and the selection auditor returned **PASS**. The training and holdout ID sets are disjoint (140 and 60 replicates).

### 1.3 Critical interpretation: boundary minimum

IAE decreases monotonically over the entire candidate grid in both training and holdout. The selected value, 0.40, is the largest tested multiplier.

Therefore, the experiment **does not locate an interior bandwidth optimum**. It establishes only that 0.40 is the best candidate among the nine values tested under this one synthetic parent-density DGP and this particular sample size/grid. It does not justify freezing 0.40 as an optimal, universal, or confirmatory bandwidth.

Recommended status label:

`PASS_AUDIT_BOUNDARY_MINIMUM_UNRESOLVED`

The audit passes; the scientific selection question remains unresolved.

## 2. Overlapping 100-replicate pilot in PR #54

The exact-head PR #54 workflow also completed a 100-replicate pilot and selected 0.40 from a 12-value grid (0.03 through 0.40). Its selection audit passed.

This is **not an independent replication** of the 200-replicate run:
- both use master seed 20261114 and the same parent DGP;
- the 100-replicate pilot uses replicate IDs 0–99;
- the 200-replicate experiment uses replicate IDs 0–199;
- hence the first 100 replicate IDs overlap exactly, while candidate grids and split definitions differ.

The two results are consistent as overlapping sensitivity analyses, but must not be pooled or described as independent confirmation. Both point to the upper edge of their tested range.

## 3. Batched likelihood engineering benchmark in PR #54

Exact-head workflow run: [43550061324](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/43550061324). All Python 3.10, 3.11, and 3.12 jobs passed.

At 321 H0 points × 801 redshift points × 30 events × 23 density nodes (three timing repeats):
- scalar median: 6.07149 s;
- batched median: 0.31343 s;
- observed ratio: 19.37×;
- maximum absolute log-likelihood difference: 1.42e-13;
- maximum absolute log-detection-efficiency difference: 2.66e-15.

This is evidence of numerical equivalence and a useful speed improvement for the isolated likelihood-kernel calculation on that runner. It excludes KDE construction, is hardware/BLAS/load dependent, and is not a universal performance claim. It does not establish scientific calibration or estimator superiority.

## 4. Next required work before freezing a bandwidth

1. Keep the completed 200-replicate source artifact immutable.
2. Predeclare an expanded multiplier grid extending beyond 0.40, with a new master seed and replicate IDs disjoint from both existing pilots.
3. Run the expanded range as an engineering bandwidth-only pilot; retain all candidates and all failures.
4. Check redshift-grid refinement and boundary behavior so apparent improvement is not a finite-grid artifact.
5. If the expanded curve develops an interior minimum, conduct a fresh train/holdout selection on disjoint replicate IDs and independently audit it.
6. Only then freeze a density-only bandwidth rule for the later five-arm H0 study.

The bandwidth must not be selected using H0 posterior error, coverage, PIT, or any event-level outcome. The later H0 study must separately verify H0-grid and coefficient-quadrature convergence.

## 5. Claim boundary

This pilot concerns one synthetic parent population, one parent sample size, and a parent-density IAE objective. It does not evaluate catalogue-inclusion corrections A–E, does not measure H0, does not test real galaxy catalogues, and does not support a claim that any inclusion estimator is superior.


## 6. Expanded-grid sensitivity artifact from PR #54 (same-day update)

**Workflow:** [43554721287](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/43554721287)  
**Artifact:** `five-arm-kde-bandwidth-pilot`, artifact ID `14098956423`  
**Exact source head:** `e8a55f9126d286b2a764d709fc3379f504930fc4`  
**Exact-head CI:** Python 3.10, 3.11, and 3.12 all passed.

The expanded 100-replicate pilot evaluates 17 multipliers from 0.03 to 1.20. Its independent raw-data auditor reconstructed 1,700 replicate-by-multiplier IAE values across 400,000 parent units; maximum absolute reconstruction discrepancy was `7.494005416219807e-16`. The raw and train/holdout audits both returned **PASS**.

The training selector chose multiplier **0.50** (training mean IAE 0.0443212, 70 replicates). The untouched-within-this-artifact holdout mean at 0.50 was 0.0439975 (30 replicates; MCSE 0.0016767). The neighboring 0.60 holdout mean was 0.0441885. The selected value is interior to the expanded candidate grid, unlike the earlier 0.40 boundary result.

**Critical limitation:** this is not a fresh independent sample. The expanded pilot uses master seed `20261114` and replicate IDs 0–99, which overlap the first 100 IDs of the 200-replicate primary experiment and the earlier 100-replicate PR #54 pilot. Its changed grid and split make it a useful expanded-grid sensitivity analysis, but do not turn it into independent confirmation. The 0.50 value is therefore an exploratory density-only candidate for this DGP, N=4,000, grid, and candidate set—not a confirmatory frozen bandwidth.

The same exact-head artifact reports a batched likelihood-kernel numerical-equivalence pass (maximum absolute log-likelihood discrepancy 1.42e-13; maximum log-efficiency discrepancy 2.66e-15). The observed scalar/batched ratio was 19.10× for the isolated kernel benchmark on that runner; KDE construction is excluded and the timing is hardware/BLAS dependent.

### Required gate

Before using 0.50 in a confirmatory five-arm H0 experiment, run a new density-only training/holdout pilot with a master seed and replicate IDs disjoint from every prior bandwidth pilot. Then reserve a separate, disjoint replicate set for H0 evaluation. Preserve the present artifacts unchanged. Do not tune the bandwidth using H0 error, coverage, PIT, or event-level outcomes.

## 7. Independent fresh-seed expanded-grid pilot (PR #44, 2026-11-14)

**Workflow:** [43565461921](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/43565461921)  
**Exact source commit:** `e10c030ef3dbc7373909a26cb7d15729fb11f3ac`  
**Artifact ID:** `14102456385`  
**Artifact ZIP SHA-256:** `c4d758338d05f8755bbeea8722dea3c691102f164fbb469fecc8d03e377bb12f`

The corrected workflow completed successfully. All 13 job steps passed, including the fresh-seed simulation, independent reconstruction, training-only selection, independent split/score audit, provenance gate, and artifact upload. The artifact contains 200 replicates × 4,000 parent units (800,000 total), 2,200 replicate-by-multiplier attempts over the 11 frozen candidates, and the raw parent samples.

The independent raw auditor verified all 200 replicates and 800,000 parent units. It reconstructed all 2,200 IAE attempts with maximum absolute discrepancy **1.18e-16**. The selection auditor independently verified the 140/60 training/holdout split and selected multiplier **0.55**. Both audit reports return **PASS**. The source manifest records SHA-256 hashes for all four source files; the ZIP hash above identifies the downloaded artifact.

### 7.1 Selection and holdout interpretation

Training mean IAE was 0.0447583 at 0.55 and 0.0447867 at 0.50. The 0.55 training minimum is therefore extremely shallow. On the untouched holdout, mean IAE was 0.0434294 at 0.55 and 0.0433576 at 0.50; the neighboring 0.50 candidate is marginally lower on holdout. This is descriptive and does not authorize re-selection.

Because both candidates are evaluated on the same replicate IDs, a paired comparison is more informative than comparing their separate MCSEs. Define the paired difference as IAE(0.55) − IAE(0.50):

| Split | n | Mean paired difference | Paired SE | Approx. 95% t interval |
|---|---:|---:|---:|---:|
| Training | 140 | -2.84e-5 | 9.84e-5 | [-2.23e-4, 1.66e-4] |
| Holdout | 60 | 7.19e-5 | 1.59e-4 | [-2.47e-4, 3.91e-4] |

Both intervals include zero by a wide margin. The paired data do not distinguish 0.50 from 0.55 at useful precision. The result supports a broad low-error region around 0.50–0.55 for this specific synthetic density-recovery setup, not a sharply identified optimum. Keep **0.55 as the preregistered training-selected candidate**, but do not describe it as materially better than 0.50 or freeze it as a universal bandwidth.

### 7.2 What this closes—and what it does not

This is a genuinely fresh parent-sample seed schedule (master seed 20261116; split seed 20261117), disjoint from the earlier 20261114 pilots. It closes the independent density-only expanded-grid audit gate for the declared DGP, N=4,000, grid, and candidate set. It does not establish redshift-grid invariance, H0 posterior calibration, or the best arm-specific physical bandwidth. The holdout remains untouched for selection.

Before confirmatory five-arm H0 evaluation, retain a separate H0 replicate schedule and run a redshift-grid refinement sensitivity using density-only criteria. Then report the selected multiplier as a fixed protocol choice with its selection uncertainty, and evaluate all arms on the same paired H0 replicates. No H0, coverage, PIT, or event-level result may be used to retune this density bandwidth.

