# Multi-replicate common-physical-bandwidth density pilot

**Status:** orchestration implemented; CI smoke required; no selection or confirmatory inference.  
**Scope:** synthetic low-redshift density recovery only.  
**Parent branch:** research/five-arm-kde-kernel-reuse-v2-20261114.

## Purpose

The single-replicate absolute-bandwidth smoke established that the estimator and an independent raw-input reconstruction can be exercised end to end. It cannot estimate replicate-to-replicate variability or justify a bandwidth choice. This pilot adds repeated independent parent/catalogue/validation realizations while preserving the complete child raw artifacts and audit trail for each replicate.

## Implementation

- Runner: scripts/run_absolute_bandwidth_density_replicates.py.
- Child raw-input auditor: scripts/audit_absolute_bandwidth_density_smoke.py.
- Orchestration auditor: scripts/audit_absolute_bandwidth_density_replicates.py.
- Tests: tests/test_absolute_bandwidth_density_replicates.py.

Each replicate runs the existing one-replicate five-arm density experiment with one parent population, one catalogue realization, and one independently generated validation frame. All arms and candidate widths share those inputs. Arm E node-level densities and posterior-weighted node metrics are retained. The child auditor independently reconstructs weights, reflected KDEs, physical bandwidths, and density metrics from raw source records.

The outer auditor independently checks:
1. every outer manifest hash and byte count;
2. the complete planned replicate ledger and ordered replicate IDs;
3. the deterministic population/validation seed schedule and non-overlap;
4. SUCCESS/FAILED lifecycle accounting, including failure records;
5. a separate raw-input audit for every successful child replicate;
6. exact aggregate metric-key uniqueness and equality to the child records.

## Fresh-seed rule

The runner's default master seed is **20261201**, with stride 10. For replicate r=0,...,R-1:
- population seed = 20261201 + 10r;
- validation seed = population seed + 1.

The two-replicate CI smoke deliberately overrides this with master seed 20261120 and is an engineering test only. Do not include its replicate IDs in a later training or holdout sample. The default seed schedule is reserved for the subsequent fresh density pilot and must not be changed after inspecting outcomes.

## Initial pilot design

Recommended first research run:
- 20 replicates;
- parent population N=4,000;
- validation frame 1,500;
- nested validation sizes 100, 300, 1,000;
- quadrature order 5;
- 801-point redshift grid for the full density-resolution run;
- common physical bandwidth candidates 0.012, 0.016, 0.020, 0.024, 0.028, 0.032, 0.040.

The current CI smoke uses two smaller replicates, N=800, two validation sizes, two widths, and a 201-point grid. It tests orchestration and audit mechanics only. It is not the research pilot.

## Interpretation and exclusions

This pilot is deliberately **not** a bandwidth selector. It has no train/holdout split, no pooled selection objective, and no H0 stage. It does not independently refit the logistic inclusion posterior; the child auditor checks downstream raw weights conditional on the persisted posterior nodes. Posterior quadrature validity remains a separate audit gate.

The pilot can estimate descriptive replicate variability for IAE, ISE, boundary IAE, and the persisted Kish ESS diagnostics. It cannot establish that one arm is superior, that IPW is calibrated, or that a chosen bandwidth generalizes. Failures remain in the planned denominator and must not be silently dropped.

Before a formal bandwidth-selection study, add preregistered training and holdout replicate IDs that are disjoint from all previously inspected seed schedules. Freeze the candidate grid and equal-cell/arm objective before running that study. Only a separately seeded downstream experiment may propagate a density-selected bandwidth into H0 inference.


## Posterior-mean density versus posterior expected node loss

Arm E now reports two distinct density-recovery estimands, both reconstructed by the independent child auditor:

1. **Loss of the posterior-mean density:** for coefficient nodes \(\\theta_b\) and quadrature weights \(q_b\), first form \(\\bar f(z)=\\sum_b q_b f_b(z)\), then evaluate IAE, ISE and boundary IAE of \(\\bar f\) against the analytic target.
2. **Posterior expected node loss:** calculate each node-specific loss and then average it, \(\\sum_b q_b L(f_b,f_0)\). This is an uncertainty/risk summary, not the loss of the posterior-mean density.

These quantities are not interchangeable. Since IAE, ISE and boundary IAE are convex losses, the loss of the posterior-mean density cannot exceed the posterior expected node loss (up to numerical tolerance). The node-specific rows remain available. Kish ESS is not defined for a mixture density itself, so the posterior-mean-density row records it as missing rather than borrowing an arbitrary node's ESS.

This density diagnostic does not replace the Arm E cosmological likelihood: H0 inference must continue to evaluate the complete multi-event selected likelihood at each shared coefficient node and marginalize that joint likelihood once.
