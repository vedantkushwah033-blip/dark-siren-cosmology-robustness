# Preregistered Common-Physical-Bandwidth Density Selection: Training and Holdout

**Status:** design specification; not executed.  
**Branch:** `research/five-arm-kde-kernel-reuse-v2-20261114`  
**Scope:** synthetic density recovery only. This is a bandwidth-selection gate before any downstream H0 comparison. It is not an estimator-superiority, posterior-calibration, or observational study.

## 1. Purpose

The 20-replicate absolute-bandwidth pilot is exploratory. It estimates replicate variability for one declared catalogue-inclusion mechanism and checks the full artifact/audit path. It must not be reused to choose a bandwidth.

This protocol defines a fresh, disjoint training/holdout experiment to select one common physical redshift bandwidth for the five density-recovery arms. Selection is based only on density error, never on H0 posterior error, coverage, PIT, or GW event outcomes.

## 2. Target population and density

The parent redshift distribution is generated from the declared latent intensity
`g(z) = 1 + 2.2 exp[-(z-0.105)^2/(2(0.025)^2)] + 1.6 exp[-(z-0.275)^2/(2(0.040)^2)]`
on `[0.01, 0.40]`. The declared common base-population/merger weight is
`b(z) = 1 + 1.5 (z-0.01)/0.39`.

The density target for every arm is the normalized
`f_0(z) proportional to g(z)b(z)`.
The factor `b(z)` is applied exactly once: in raw sample weights A/B as `b(z)`, and C/D/E as `b(z)/pi(z)`. The KDE's post-estimation `merger_weight` is the identity. Applying `b(z)` again after the weighted KDE would change the target to `g(z)b(z)^2` and is prohibited.

The target is an analytic synthetic population target. The finite parent sample in arm B is a parent-population oracle diagnostic, not the exact analytic target.

## 3. Inclusion DGP cells and explicit scope

Use catalogue-inclusion floors `0.00, 0.10, 0.25`, with midpoint 0.22 and width 0.05:
`pi(z) = floor + (1-floor)/(1+exp((z-0.22)/0.05))`.

- Floor 0.00: degree-one logistic truth in the standardized redshift coordinate; correctly specified for the declared fitted model.
- Floors 0.10 and 0.25: the additive floor makes the logit nonlinear; degree-one logistic fitting is structurally misspecified to different degrees.

This factorial studies floor-induced misspecification. It does **not** test the separate quadratic-logit curvature family `expit(alpha + beta*x + gamma*(x^2-c))`. No claim about gamma-curvature robustness may be inferred from this experiment. A separate predeclared DGP extension is required before making such claims.

## 4. Five arms and primary/secondary estimands

For each replicate, generate one parent population, one thinned catalogue, and one independent validation frame. Reuse those exact inputs across all arms, validation sizes, and bandwidths.

- A: catalogue-only, raw weights `b(z)`.
- B: parent-population oracle, raw weights `b(z)`.
- C: known-probability IPW, raw weights `b(z)/pi_true(z)`.
- D: validation-MAP plug-in IPW, raw weights `b(z)/pi_MAP(z)`.
- E: validation-posterior uncertainty-integrated IPW, with one node-specific weighted KDE per shared coefficient node.

**Primary E estimand:** form the posterior-mean density `fbar(z)=sum_k q_k f_k(z)`, then compute its IAE. This is the density estimator delivered by posterior averaging.

**Secondary E estimand:** posterior expected node-specific IAE `sum_k q_k IAE(f_k,f_0)`. This measures node uncertainty/risk and is not interchangeable with primary E. Keep all node-specific metrics for diagnosis. Do not assign a Kish ESS to the posterior-mean mixture density.

Primary selection compares A, B, C, D and E-posterior-mean density. E-node expected loss and node diagnostics are secondary and cannot change the selected bandwidth.

## 5. Fixed design

| Dimension | Frozen value |
|---|---:|
| Parent population per replicate | 4,000 |
| Independent validation frame | 1,500 |
| Nested validation sizes | 100, 300, 1,000 |
| Catalogue floors | 0.00, 0.10, 0.25 |
| Quadrature order | 5 |
| Redshift grid | 801 equally spaced points on [0.01, 0.40] |
| Candidate physical bandwidth SDs | 0.012, 0.016, 0.020, 0.024, 0.028, 0.032, 0.040 |
| Training replicates per floor | 100 |
| Holdout replicates per floor | 100 |
| Primary arms in selection objective | A, B, C, D, E-posterior-mean |
| Primary loss | Integrated absolute error |

Within a replicate, the validation samples are nested prefixes from one fixed random permutation. All candidate widths use the same samples and raw weights. No arm-specific resampling is allowed.

## 6. Disjoint seed schedule

Reserve the following seed ranges before looking at any outcomes:

- Exploratory pilot: master seed 20261201, stride 10 (already reserved; never use for selection).
- Training: master seed 20270101, stride 10, replicate IDs 0–99.
- Holdout: master seed 20280101, stride 10, replicate IDs 0–99.

For replicate `r`, use population seed `master + 10*r` and validation seed `population seed + 1`. The two seed streams within each replicate are distinct. Training and holdout master-seed ranges are disjoint from one another and from the exploratory pilot schedule. Do not change seeds, replicate count, or candidate grid after inspecting any result.

Each floor is run as a separate frozen invocation using the same floor-specific seed schedule. Record the full configuration, seeds, software revision, and artifact hashes.

## 7. Primary selection objective

For candidate width `h`, compute IAE for each replicate, floor, validation size, and primary arm. First average over replicates within each floor × validation-size × arm cell. Then average the cell means with **equal cell weights** across:
`3 floors × 3 validation sizes × 5 primary arms = 45 cells`.

This prevents a larger catalogue or a numerically easier arm from dominating the objective. The objective is:
`R(h) = (1/45) sum_cell mean_replicate IAE(cell, h)`.

Tie-break: choose the smaller physical bandwidth. Report all 7 candidate objective values and all 45 cell-level contributions.

### Failure accounting

The planned denominator is fixed at 100 replicates in every floor × validation-size × primary-arm cell. A failed or non-evaluable primary cell is not silently dropped. Since both estimated and target densities are normalized, IAE is bounded above by 2; assign IAE = 2 to a failed cell for the conservative failure-inclusive training objective, while retaining its exact lifecycle status and error in the ledger. Also report the unpenalized successful-only mean as a clearly secondary descriptive statistic.

A candidate is **selection-eligible** only if every primary cell has at least 95 successful/evaluable replicates out of 100 and its failure-inclusive objective is finite. If no candidate meets this rule, the selection outcome is `NO_SELECTION`; do not relax the threshold after seeing results. The selected candidate must also not be at either endpoint of the candidate grid. If the training minimum is at an endpoint, report `BOUNDARY_LIMITED_EXPAND_GRID` and do not freeze a bandwidth.

The holdout must preserve all 100 planned replicate rows per cell. Apply the same IAE=2 failure penalty and eligibility rule. Holdout results never change the training-selected candidate.

## 8. Holdout reporting

For every candidate width, report:
- failure-inclusive equal-cell objective;
- all 45 cell-level means and successful/evaluable counts;
- replicate-level paired contrasts against the training-selected width;
- paired mean difference, SD, MCSE, and a 95% interval;
- successful-only summaries, clearly labelled secondary;
- arm-specific and validation-size-specific results;
- boundary IAE, ISE, bias/variance decomposition, and weight diagnostics as secondary outcomes;
- the full failure ledger and all raw-data audit outcomes.

The holdout is an independent evaluation of the training selection, not another tuning set. Do not replace or reselect the training choice based on holdout results.

## 9. Independent audit requirements

The independent auditor must not import the production KDE implementation or production selection function. It must independently:

1. Verify every file hash and byte count from each child manifest and the outer manifest.
2. Reconstruct the seed schedule, floor, parent sample, catalogue, validation frame, nested validation IDs, and planned replicate ledger.
3. Recompute the analytic target `g(z)b(z)` and verify that the base measure is applied exactly once.
4. Reconstruct all A–E raw weights from persisted truth/observable inputs; refit the validation logistic posterior independently or explicitly label posterior nodes as unverified. A bandwidth-selection gate requires independent posterior reconstruction for D/E, not merely trusting saved nodes.
5. Reconstruct reflected weighted KDE curves and requested/realized physical widths from raw samples and weights using an independent implementation.
6. Recompute all primary/secondary density metrics and cell-level aggregates.
7. Verify failure penalties, fixed denominators, eligibility, endpoint rule, tie-break, and exact selected candidate.
8. Reject missing, duplicated, extra, or silently dropped replicate × floor × validation-size × arm × bandwidth cells.

Hash verification alone is not an independent statistical audit. Refreshing a manifest after corrupting an output must not make the altered result pass the mathematical reconstruction.

## 10. Interpretation constraints

A selected width is conditional on this synthetic latent intensity, base measure, catalogue inclusion family, validation design, KDE estimator, redshift interval, grid, candidate set, and loss. It is not a universal bandwidth and must not be described as optimal for real galaxy catalogues.

The density-only stage does not establish H0 posterior calibration, coverage, SBC, or estimator superiority. After this gate passes, freeze the selected width and run a separately seeded H0 experiment with a complete selection-conditioned likelihood, shared coefficient node across all events and selection efficiencies, quadrature/grid convergence, and repeated-sampling calibration.

The floor family is not a substitute for the project's quadratic-logit curvature question. Before the final paper, implement and independently audit that separate curvature family (including a principled centering constant under the actual parent measure) or narrow the paper's claims to the floor-based misspecification experiment.

## 11. Execution gate

Do not launch the 600-replicate training/holdout workload until:
- the corrected 20-replicate exploratory pilot completes;
- all 20 child raw-input audits and both aggregate reconstructions pass;
- the pilot's descriptive report and failure analysis are reviewed;
- the runner supports three floor cells without changing any frozen seed or metric definitions;
- the independent posterior-posterior reconstruction requirement for D/E is implemented and tested.

No H0 experiment is authorized by this document alone.
