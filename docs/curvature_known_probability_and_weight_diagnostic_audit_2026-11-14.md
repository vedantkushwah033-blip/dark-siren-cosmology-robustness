# Curvature DGP truth-probability and weight-diagnostic audit

Date: 2026-11-14

## Purpose

This change closes two implementation/audit gaps before the quadratic-logit curvature density experiment is run. It does not launch that experiment and does not authorize bandwidth selection or H0 inference.

## 1. Arm C now accepts the exact known inclusion probabilities

`build_five_arm_weights` now accepts the optional argument
`true_catalogue_inclusion_probability`.

- If omitted, the existing floor/midpoint/width catalogue-selection model is used, preserving all previous callers.
- If supplied, it must be a finite one-dimensional vector with one value per observed catalogue row and every value in `(0, 1]`.
- Arm C uses the supplied truth vector in catalogue row order.
- Arms D and E continue to fit the inclusion model using only the observable validation redshifts and binary retention labels. The supplied truth vector is not passed to the posterior fitter.

This is required for the quadratic-logit stress DGP: the known-probability comparator must use the calibrated quadratic-logit probabilities, while the analysis model remains degree-one logistic under nonzero curvature.

Regression tests cover the explicit-truth path, unchanged validation-only D/E construction, and invalid vector shape/range/finiteness.

## 2. Entropy effective sample size is now reported

`weight_concentration_diagnostics` now reports Shannon entropy ESS:

[
ESS_H = expleft(-sum_i p_ilog p_iight),
qquad p_i = w_i/sum_j w_j.
]

The existing Kish ESS, normalized maximum weight, and top 1%/5% mass diagnostics are retained. For Arm E's aggregate diagnostic row, the minimum across coefficient-quadrature nodes is reported for both ESS measures, while the maximum is reported for concentration-tail measures. Node-level density metrics remain separately available.

## 3. The child auditor independently checks the diagnostic file

The raw-input density auditor now reconstructs all six weight diagnostics directly from saved raw weight arrays and compares them with `results/weight_diagnostics.csv`. It does not import the production diagnostic function.

A regression test modifies the diagnostic file, refreshes its manifest hash, and verifies that the independent mathematical reconstruction still rejects the tampered values. This distinguishes a genuine audit from checksum-only verification.

The multi-replicate runner's configuration text was also corrected: its child auditor independently refits the logistic posterior from observable validation records and reconstructs quadrature nodes/weights. The prior limitation string incorrectly implied those nodes were merely consumed.

## 4. Current validation status

- The pre-existing focused curvature-DGP suite passed on Python 3.10, 3.11, and 3.12 before these changes.
- The new code and regression tests are committed on `research/curvature-known-probability-fix-20261114`.
- Full CI for this exact branch head is still required. No passing status is claimed here.
- The expanded 11-width pilot and three-floor pilot workflows remain separate and are not modified by this branch.
- No seven-gamma density experiment has been launched.

## Scientific guardrails

The quadratic-logit experiment remains a synthetic low-redshift methodology stress test. It is not real-galaxy-survey validation, public LVK/GWTC analysis, or an observational H0 measurement. Oracle probabilities are simulation diagnostics only. A bandwidth must not be selected from the same data used for confirmatory evaluation, and H0 propagation remains blocked until the density-stage audit and preregistered train/holdout design are complete.
