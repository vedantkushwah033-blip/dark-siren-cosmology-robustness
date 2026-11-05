# LP-based logistic separation diagnostic — design and limits

**Date:** 5 November 2026  
**Branch:** `research/logistic-inclusion-v2-newton-2026-11-05`  
**Scope:** diagnostic only; no estimator, objective, coefficient bound, or pilot result is changed.

## Motivation

The V2 logistic fitter reports a box-constrained KKT residual, Hessian rank and coefficient-boundary status. These are useful optimization diagnostics, but they do not formally classify complete, quasi-complete and overlap geometries. A coefficient reaching a finite box face is not itself a formal separation test.

This module adds a separate linear-programming diagnostic:
`src/dark_siren_robustness/logistic_separation_diagnostics.py`.

## Mathematical construction

For design row (x_i^T), binary response (y_i), and sign (s_i=2y_i-1), define the signed margin
[
m_i(b)=s_i x_i^T b.
]

The diagnostic uses (b=b^+-b^-), with (b^+,b^-geq0) and
[
sum_j(b_j^++b_j^-)leq1.
]
This L1 constraint fixes the arbitrary scale of a separating direction; it is not a constraint imposed on the logistic estimator. Nonzero design columns are scaled to unit maximum absolute value before solving. Nonzero column rescaling does not change whether a separating direction exists.

### Complete separation

Solve
[
max_{b^+,b^-,t} t
]
subject to
[
-m_i(b),+tleq0quad	ext{for every }i,qquad tgeq0,
]
and the L1 constraint above. A maximum margin above the declared tolerance indicates complete separation.

### Quasi-complete separation

If the complete-separation margin is not above tolerance, solve
[
max_{b^+,b^-}sum_i m_i(b)
]
subject to
[
m_i(b)geq0quad	ext{for every }i
]
and the same L1 constraint. A positive optimum above tolerance indicates quasi-complete separation: all margins can be nonnegative and at least one is positive, while complete separation was not found.

If both optimized criteria are within tolerance, the diagnostic returns overlap. If an LP fails, it returns indeterminate rather than treating solver failure as evidence of overlap.

Design rank is computed and returned separately. Rank deficiency may coexist with separation or overlap and is not itself a separation category.

## Output

`diagnose_logistic_separation(design, included, tolerance=1e-9)` returns:
- `separation_status`: complete, quasi-complete, overlap, or indeterminate;
- LP status codes/messages and the complete margin/quasi score;
- design rank, parameter dimension and rank-deficiency flag;
- tolerance and column-scaling factors.

## Regression fixtures

`tests/test_logistic_separation_diagnostics.py` checks:
1. complete separation;
2. quasi-complete separation with both labels at the separating boundary;
3. overlap with interleaved outcomes;
4. rank deficiency reported separately from overlap;
5. invariance to nonzero column rescaling;
6. rejection of malformed and single-class inputs.

The test is added to the existing GitHub Actions matrix. Passing CI will establish implementation behavior on these fixtures, not a universal proof of numerical robustness for every pathological design.

## Explicit limitations

- This is a diagnostic of the supplied design/outcome matrix, not a logistic fit and not an automatic estimator switch.
- The LP tolerance is a numerical classification threshold; near-boundary cases should be interpreted with the returned margins and LP status, not forced into a confident binary label.
- Column scaling improves conditioning but does not guarantee every LP is numerically well-conditioned.
- The diagnostic does not establish positivity, validation representativeness, fitted-probability calibration, stable IPW, host-density accuracy or H0 posterior calibration.
- A finite box-constrained logistic solution remains a constrained estimator even when the unconstrained MLE is nonexistent.
- D-01 remains historically unidentifiable without original records.

## Verification status

The implementation and tests have been committed to the research branch. Exact-head GitHub Actions verification is pending until the workflow run for the resulting branch head completes. No local test execution or confirmatory simulation is claimed in this document.
