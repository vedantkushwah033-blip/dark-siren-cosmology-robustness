# Absolute Physical Bandwidth Sensitivity Protocol for the Five-Arm Study

**Status:** design specification; not yet implemented or executed.  
**Branch:** `research/five-arm-kde-kernel-reuse-v2-20261114`  
**Purpose:** determine whether conclusions from a shared weighted-Scott multiplier are artifacts of arm-specific realized kernel widths.

## 1. Motivation and estimand

The five-arm KDE implementation applies one Scott multiplier to every arm. This is a common *rule*, not a common physical smoothing scale. The realized redshift-space Gaussian kernel standard deviation is

`h_a = m n_{eff,a}^{-1/5} s_{w,a}`,

where `m` is the shared multiplier, `n_eff,a` is the arm-specific Kish effective sample size, and `s_w,a` is the arm-specific weighted standard deviation. Thus the catalogue-only, oracle, known-IPW, MAP-IPW, and coefficient-node densities can receive different physical bandwidths. A difference between arms may therefore combine a weighting effect with a smoothing effect.

This protocol estimates the sensitivity of **host-redshift density recovery** to that distinction. It does not tune on H0 posterior error, coverage, PIT, or any event-level outcome.

## 2. Questions

1. Does a common physical redshift bandwidth change the relative density-recovery performance of arms A–E compared with the common-Scott-multiplier analysis?
2. Are apparent gains or losses from IPW robust over a predeclared physically common bandwidth range?
3. Does Arm E's shared-coefficient uncertainty integration remain distinguishable from the MAP plug-in arm when every coefficient node uses the same physical kernel width?
4. Are conclusions sensitive to the bandwidth grid endpoints or redshift-grid resolution?

## 3. Fixed data-generating design

Use the already declared synthetic parent-population and catalogue-inclusion DGP from the five-arm preregistration. Preserve common random numbers across all arms and bandwidth values within each replicate:

- one parent population;
- one catalogue-inclusion realization;
- one independently sampled validation frame and its nested validation subsets;
- one detected-event realization, retained only for the later H0 stage;
- one replicate seed schedule shared across all compared bandwidths.

The primary density-recovery stage uses the same parent truth and observed catalogue/validation inputs for A–E. It must not resample the catalogue or validation sample separately for different bandwidths.

The primary grid of validation sizes remains `100, 300, 1000`. The curvature and catalogue-retention factors remain those already declared in the main preregistration. Do not narrow the DGP after seeing the results.

## 4. Bandwidth parameterization

Introduce an explicit `absolute_bandwidth_sd` option in redshift units, mutually exclusive with both `bandwidth_factor` and `scott_multiplier`. It denotes the standard deviation of each Gaussian kernel in the original, unreflected redshift coordinate.

For a reflected weighted KDE, the SciPy covariance factor must be computed from the augmented reflected sample covariance so that the resulting physical kernel SD is exactly the requested `absolute_bandwidth_sd`. The independent auditor must reconstruct this factor from raw samples and raw weights without importing the production estimator.

Do not implement an absolute width by converting it to a Scott multiplier separately for each arm inside the estimator and then reporting only the multiplier. Persist both the requested common physical width and each independently reconstructed realized width; require their agreement within a tight numerical tolerance.

## 5. Candidate bandwidth grid and split

The bandwidth grid must be frozen before density-recovery outcomes are examined. Use a small, log-spaced grid spanning clear undersmoothing through moderate oversmoothing. The exact endpoints should be selected from a bandwidth-only engineering pilot using the distribution of *realized physical widths* under the current Scott rule, not using H0 results. Record the pilot seed range and the rule used to set the endpoints in the configuration.

Use disjoint replicate IDs for:
- **training:** select one common absolute width by the predeclared pooled mean IAE objective, giving equal weight to each arm × validation-size × DGP cell;
- **holdout:** evaluate the frozen selected width and every candidate width without re-selection.

All five arms must have equal eligibility in the training objective. If an arm or planned cell fails, preserve the failure in the planned denominator; do not silently average over whichever arms succeeded. A confirmatory selection cannot proceed unless the predeclared minimum evaluable fraction is met for every arm/cell, or a conservative failure-inclusive objective is explicitly used.

The holdout report must show all candidate widths, not only the selected width. The selected width is chosen from training data only. No holdout re-selection is allowed.

## 6. Primary density metrics

For every replicate, DGP cell, validation size, arm, and bandwidth, compute against the analytic merger-weighted parent-host density:

- integrated absolute error (IAE);
- integrated squared error (ISE);
- integrated bias function and integrated variance across replicates;
- boundary-region IAE, with boundary regions fixed before the run;
- support diagnostics and the arm-specific effective sample size / weight concentration.

The primary metric for bandwidth selection is pooled IAE, with equal cell and arm weighting. ISE and boundary-region IAE are secondary descriptive outcomes. Report per-cell values and paired arm contrasts; do not allow a large/easy cell to dominate by having more successful rows.

For Arm E, compute density recovery for each coefficient-posterior node and the posterior-weighted node-level summary using the same absolute width. Preserve node identities and quadrature weights. Do not average coefficient nodes before constructing the KDE unless the statistical definition explicitly calls for that operation.

## 7. Independent audit requirements

The independent auditor must:

1. verify the source manifest, all raw sample and validation inputs, candidate grid, split, and checksums;
2. independently reconstruct the analytic parent target;
3. independently reconstruct each arm's raw galaxy weights, including the base-population factor and inverse-inclusion factor where applicable;
4. reconstruct each reflected weighted KDE using the requested common physical bandwidth;
5. independently recalculate IAE, ISE, boundary IAE, and all paired contrasts from raw density arrays or raw source samples;
6. verify that every planned replicate × DGP × validation-size × arm × bandwidth key is present exactly once;
7. reconcile every attempt to `SUCCESS`, `NON_EVALUABLE`, `FAILED`, or `NOT_RUN`, retaining all planned attempts in the denominator;
8. reject missing rows, duplicate keys, altered seed schedules, altered candidate grids, manifest mismatch, and any holdout-informed selection.

The audit implementation must not import the production KDE estimator, arm-weight builder, bandwidth-selection function, or aggregate-metric function. A shared low-level constant or documented DGP formula is acceptable only if the auditor independently checks it against the preregistered configuration.

## 8. Numerical gates

Before the density pilot is accepted:

- direct one-dimensional Gaussian-kernel calculations must match the estimator on deterministic weighted fixtures;
- the reflected density must integrate to one on the declared finite grid;
- the reported physical bandwidth must equal the requested absolute width;
- grid refinement must show IAE stability at a predeclared tolerance;
- a scalar reference and batched implementation must agree for all fixed arms and all Arm E nodes;
- zero-weight/underflow, extreme-weight, boundary-concentrated, and small-effective-sample cases must fail explicitly rather than silently change the target bandwidth;
- all Python matrix jobs must pass at the exact PR head.

## 9. H0 propagation gate

Do not run the five-arm H0 comparison merely because a bandwidth is selected. First publish the complete density-recovery training and holdout report and verify that the density stage passes its independent audit.

For the downstream H0 stage:

- freeze the bandwidth using density-only training data;
- use the same absolute physical width across A–E and every Arm E coefficient node;
- preserve the common parent/catalogue/events/validation inputs;
- retain the GW detection-selection normalization separately from catalogue inclusion;
- for Arm E, use one shared coefficient node across all events and marginalize the complete selected joint likelihood exactly once;
- evaluate H0-grid and coefficient-quadrature convergence independently of bandwidth selection;
- keep density-stage selection and H0-stage evaluation on disjoint replicate seeds.

Report both the common-Scott and common-absolute-width analyses as distinct sensitivity analyses. Neither is universally preferred by construction.

## 10. Interpretation limits

A common physical bandwidth is a controlled sensitivity analysis, not automatically the uniquely correct estimator. Conversely, a common Scott multiplier is a coherent adaptive rule but does not isolate weighting from smoothing. The two analyses answer different questions and should be reported side by side.

No result from this protocol alone establishes H0 calibration, estimator superiority, observational validity, or performance on LVK/GWTC data. Any claim must be restricted to the simulated DGPs, sample sizes, validation designs, bandwidth range, and numerical resolution actually audited.
