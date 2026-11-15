# Curvature benchmark: proposed training/holdout bandwidth protocol

**Date:** 2026-11-15  
**Status:** proposal only; not yet frozen, run, or audited  
**Depends on:** completion and independent audit of the seven-gamma development pilot  
**Scope:** synthetic low-redshift density recovery only

## 1. Why a separate bandwidth study is required

The existing three-floor experiment and the curvature pilot evaluate several fixed physical bandwidths. Their empirical minima are useful for sensitivity analysis, but choosing a bandwidth and reporting performance on the same replicates would reuse evaluation data. The expanded floor-family grid also placed some catalogue-only minima at the upper boundary, so its minima are not uniformly bracketed.

The next study should separate bandwidth tuning from performance evaluation using disjoint replicate seed sets. It must not use the development-pilot or earlier floor-family replicate outcomes for either training or holdout scoring.

## 2. Scientific estimands

The primary comparison remains estimated-inclusion MAP IPW (Arm D) against catalogue-only density estimation (Arm A). The oracle arms remain diagnostic:

- Arm B: parent-density oracle;
- Arm C: known-inclusion-probability IPW oracle;
- Arm D: degree-one logistic MAP estimated IPW;
- Arm E: posterior-mean density, with posterior-mean density loss primary and posterior expected node-specific loss secondary.

Two bandwidth estimands should be reported.

### Primary: estimator-specific tuning

For each ((\gamma,n,arm)), select a bandwidth using training-replicate mean IAE for that arm, then evaluate that frozen bandwidth on independent holdout replicates. The primary held-out contrast is:
[
\Delta_{D-A}=\operatorname{IAE}_{D}(h_D^*)-\operatorname{IAE}_{A}(h_A^*).
]
This compares the attainable density-recovery performance of the two complete estimator-plus-smoothing procedures under simulation-oracle tuning. It is not an operational bandwidth selector for real galaxy catalogues.

### Secondary: common-bandwidth comparison

For each ((\gamma,n)), select one common bandwidth by minimizing the equally weighted training mean IAE of Arms A and D. Evaluate both arms at that same frozen bandwidth on holdout replicates. This isolates the estimator contrast at matched smoothing and helps reveal whether a result is driven by estimator differences or by the bandwidth selected for each estimator.

Arms B and C should also have estimator-specific training-selected bandwidths for descriptive held-out reference curves. The oracle status of B and C must remain explicit. For Arm E, bandwidth selection is based only on its **primary** estimand, loss of the posterior-mean density. The secondary posterior expected node-specific loss is evaluated at that same frozen Arm E bandwidth; it does not receive a separate bandwidth optimization. This keeps the two Arm E summaries tied to one estimator procedure and prevents the secondary risk from silently becoming a second tuning target.

## 3. Proposed bandwidth grid and boundary rule

Initial physical-bandwidth grid:
[
\{0.004,0.006,0.008,0.010,0.012,0.016,0.024,0.032,
0.048,0.064,0.080,0.100,0.120,0.160\}.
]

The grid is shared across gamma values, validation sizes, and arms. If any primary estimator's training minimum is at either grid boundary, expand that boundary and rerun the training stage only before freezing the selection. A proposed expansion is to halve the lower endpoint or double the upper endpoint. Permit at most two expansion rounds. If a minimum remains at a boundary after two rounds, label that stratum boundary-limited and do not describe its selected bandwidth as an interior optimum.

Tie rule: if candidate mean IAEs differ from the minimum by no more than (10^{-6}) in absolute IAE, choose the largest bandwidth among tied candidates. Record the complete training loss curve, not only the selected value.

## 4. Proposed replicate and seed separation

A resource-aware initial design is:

| Split | Replicates per gamma | Master seed | Purpose |
|---|---:|---:|---|
| Training | 12 | 2026111600 | Select and freeze bandwidths |
| Holdout | 24 | 2026112600 | Estimate held-out performance |

Use the existing seed stride of 10. Within each replicate, all gamma cells share the population and validation seed pair to preserve common-random-number pairing. Replicate seed pairs must be disjoint. Training and holdout seed ranges must be disjoint from each other and from the development pilot seed range.

Validation sizes remain nested at (n=100,300,1000), drawn from the same independent validation frame within a replicate. The largest validation sample is 1,000. Use the seven planned curvature values:
[
\gamma\in\{-2,-1,-0.5,0,0.5,1,2\}.
]

The 12/24 split is a first computational benchmark, not a universal guarantee of adequate Monte Carlo precision. Report the holdout MCSE and interval for every primary contrast. If precision is inadequate, any larger holdout replication must be specified before inspecting the additional outcomes.

### Reproducible seed-partition audit

Before freezing either split, create a machine-readable seed inventory from every prior development and floor-family artifact that could be reused accidentally. Extract the actual population and validation seeds from each replicate/cell ledger and child configuration; where an artifact stores only a master seed and deterministic schedule, independently expand that schedule and reconcile it against the child-level values. Do not infer seed use from directory names or from a single master-seed label.

For the proposed schedule, the population/validation seed pairs are
\[
(s_r,s_r+1),\qquad s_r=M+10r,
\]
with \(M=2026111600\) for training and \(M=2026112600\) for holdout. Enumerate every integer seed actually consumed by each split, including any separately seeded validation frame, bootstrap, or posterior simulation if those are not derived deterministically from the pair. Require exact set intersections to be empty for training versus holdout, training versus every prior development artifact, and holdout versus every prior development artifact. Also check that each replicate's two seeds are distinct and that no replicate pair is reused within a split.

The audit record must contain the source artifact path and hash, the extracted seed pairs, the expanded schedule where applicable, the three pairwise disjointness results, and the code/version used to extract them. If any prior artifact cannot be enumerated reliably, the seed partition is **not yet verified**: do not freeze the proposed seeds or launch training. This is an inventory check, not a claim that distinct seeds make all random streams statistically independent.

## 5. Training-stage algorithm

For each gamma, validation size, and arm:

1. Compute replicate-level IAE at every candidate bandwidth for all successful training replicates.
2. Average IAE across training replicates at each bandwidth.
3. Apply the boundary rule before selection is frozen.
4. Select the minimum-mean-IAE bandwidth, applying the deterministic tie rule.
5. Save all candidate means, replicate-level values, selected bandwidth, boundary flags, training seed manifest, and a hash of the frozen selection file.

For the common-bandwidth estimand, average Arms A and D equally at each candidate bandwidth before applying the same minimization and tie rule.

Do not select on ISE, boundary IAE, ESS, or holdout outcomes. These remain secondary diagnostics. No bandwidth may be changed after holdout scoring begins.

### Training-stage failure and denominator rule

Bandwidth selection requires a complete, fixed training denominator. For each \((\gamma,n,\mathrm{arm})\) selection stratum, every one of the 12 planned training replicates must have a finite IAE at every candidate bandwidth. For the common-bandwidth A/D estimand, both A and D must be present for all 12 replicates at every candidate bandwidth. Do not average over whichever replicates happen to survive, do not allow the effective training denominator to vary by bandwidth, and do not impute missing metrics.

If any required replicate, arm, or candidate-width result is failed, non-evaluable, missing, duplicated, or non-finite, do not freeze a bandwidth for that affected stratum. Preserve the original failure and artifact. A corrected implementation may rerun the complete affected training design using the same prespecified seeds when the failure is demonstrably a software or infrastructure defect; record both attempts and audit the replacement artifact. Do not replace a difficult replicate with a newly chosen seed after examining its outcome. If a complete denominator cannot be recovered, mark the stratum non-evaluable and do not score its corresponding holdout contrast as though a bandwidth had been selected.

The training auditor must independently verify the full replicate × gamma × validation-size × arm × bandwidth grid, one-to-one replicate keys, finite metric values, and the common 12-replicate denominator before accepting any selection file. The selection file must retain the per-replicate training values and the full candidate-width curves, not just the selected bandwidth.

## 6. Holdout-stage algorithm

1. Load the frozen training selection file and verify its hash.
2. Reject any selection file whose seed manifest overlaps the holdout seed range.
3. Run the full seven-gamma holdout design with the frozen physical bandwidths.
4. Independently audit every child raw input, manifest, aggregate, lifecycle record, and paired contrast.
5. Report held-out replicate-level paired differences, mean differences, standard deviations, MCSEs, and uncertainty intervals.
6. Report all planned strata, including failed or non-evaluable cells. Do not silently drop failed pairs.
7. If any planned cell required for a paired contrast fails, suppress the affected paired contrast and mark it non-evaluable; preserve all successful raw records and the failure ledger.

The independent auditor must reconstruct selected-bandwidth assignments from the frozen training file and independently verify that the holdout runner used those exact assignments.

## 7. Reporting

Primary outcomes:
- held-out (\Delta_{D-A}) for estimator-specific tuning;
- held-out common-bandwidth D-minus-A contrast;
- IAE and MCSE by gamma and validation size.

Secondary outcomes:
- ISE and boundary IAE;
- D-minus-C and C-minus-B oracle gaps;
- E-minus-D difference;
- Kish ESS, maximum normalized weight, and top 1%/5% weight mass;
- selected bandwidth stability across training resamples and boundary flags.

Report gamma zero separately as the correctly specified degree-one logistic reference. Nonzero gamma values represent structural misspecification under the degree-one analysis model, while parent-average inclusion rate is held fixed by the calibrated intercept.

Do not interpret a larger Kish ESS as evidence that functional-form bias is absent. Do not combine ESS and density error into an unvalidated scalar score.

### Primary contrast family and uncertainty reporting

The confirmatory family, if this proposal is later frozen as confirmatory, is the 21 estimator-specific held-out D-minus-A contrasts across seven gamma values and three validation sizes. For each stratum, report the 24 planned replicate-level paired differences, their mean, sample standard deviation, MCSE, and a two-sided 95% paired-mean t interval with 23 degrees of freedom. Report the paired t-test p-values with Holm family-wise adjustment across the complete 21-contrast family at 0.05. The unadjusted intervals remain descriptive; do not interpret an isolated unadjusted interval excluding zero as confirmatory evidence after examining all 21 strata.

The common-bandwidth D-minus-A contrasts are secondary and should be presented with the same replicate-level summaries but labelled secondary. Oracle-arm comparisons, Arm E comparisons, ISE, boundary IAE, and weight diagnostics are descriptive secondary analyses; do not promote whichever secondary contrast looks most favourable into the primary claim. If the 24-replicate holdout has a failed planned cell, preserve its failure record and mark the affected primary stratum non-evaluable rather than calculating an interval from a reduced denominator. If the planned primary family is incomplete, report the family as incomplete and do not apply Holm only to the surviving subset as though the original family had been completed.

## 8. Failure and decision gates

The bandwidth protocol is not ready to execute until all of the following hold:

- the seven-gamma development pilot completes;
- its outer audit passes;
- exact-head Python 3.10, 3.11, and 3.12 focused tests pass;
- the independent DGP reconstruction and failure-injection tests pass;
- the training/holdout seed ranges and output schemas are reviewed and frozen.

The training stage must pass its outer audit before the selection file is frozen. The holdout stage must pass its outer audit before any performance conclusion is reported.

Even a successful holdout density study does not authorize an observational (H_0) claim. Density-to-(H_0) propagation remains a separate later stage requiring selected-likelihood normalization, detection-selection treatment, numerical convergence, and calibration coverage audits.

## 9. Current decision

**Do not launch the training/holdout study yet.** First finish and audit the current seven-gamma development pilot and obtain exact-head test results. This document is a proposal for review, not a preregistered or executed experiment.
