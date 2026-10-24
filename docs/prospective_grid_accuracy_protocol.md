# Prospective numerical-accuracy protocol for posterior-grid refinement

## Purpose and scope

This protocol defines acceptance criteria for a **new, independent grid-refinement audit** of the synthetic two-host dark-siren model. It is not a claim that the existing 1,001/2,001/4,001-point experiment passed these criteria. The existing experiment was already inspected and includes a known upper-95% endpoint outlier; it cannot honestly be described as a blinded or preregistered test against criteria selected afterward.

This protocol evaluates numerical discretization only. It does not validate the astrophysical assumptions, catalogue selection model, likelihood family, or the interpretation of a synthetic result as an observational H0 measurement.

## Comparison design

- Hold the simulated observation, measurement uncertainty, redshifts, prior, and inference weights fixed across resolutions.
- Use the existing 36 design cells (three generating H0 values, three second-host redshifts, two noise fractions, and two inference-weighting configurations).
- Use a new fixed seed, 20261101, and at least 200 replicates per design cell.
- Compare 1,001, 2,001, and 4,001 grid points for the broad audit.
- For every case whose upper 95% endpoint discrepancy fails the criterion, run a targeted refinement ladder through 8,001, 16,001, 32,001, and 64,001 points.
- Preserve per-observation rows and all configuration metadata. Do not discard outliers.

## Prospective acceptance criteria

For each posterior summary separately (median, lower/upper 68% endpoints, lower/upper 95% endpoints), compare matched observations at 1,001 and 4,001 points.

1. **Cell-level typical-tail criterion:** in every design cell, the empirical 95th percentile of absolute differences must be no greater than 0.01 km s^-1 Mpc^-1.
2. **Individual-case guardrail:** no case may have an absolute difference greater than 0.05 km s^-1 Mpc^-1 without being flagged for targeted refinement. A flagged case is not silently counted as a pass.
3. **Targeted refinement criterion:** for every flagged case, the difference between the 32,001- and 64,001-point summaries must be no greater than 0.001 km s^-1 Mpc^-1, and the full refinement sequence must be reported. If the sequence is non-monotone or the final difference exceeds tolerance, mark the case unresolved and investigate rather than automatically accepting it.
4. **Normalization criterion:** the posterior density must integrate to one on its actual prior support to absolute tolerance 1e-8 using an independent recomputation of the trapezoidal integral.
5. **Credible-interval ordering:** lower_95 <= lower_68 <= median <= upper_68 <= upper_95 for every case, subject only to explicit floating-point tolerance.

The numerical thresholds are operational tolerances for this synthetic analysis, not universal standards in cosmology. They are intended to be small compared with the model's distance/noise-induced uncertainty and the precision of any interpretation. Report results in the native H0 units. If the study's intended scientific precision changes, revise the protocol transparently before running the independent audit.

## Decision rule

Report a summary as numerically adequate only if every design cell passes the 95th-percentile criterion and all cases exceeding the individual-case guardrail are either resolved by targeted refinement or explicitly reported as unresolved. A failed cell cannot be hidden by a pooled median across cells. If any cell fails, do not claim uniform convergence at the 4,001-point resolution.

## Required outputs

Retain:
- complete per-observation CSV,
- per-cell summary CSV,
- targeted-refinement CSV for every flagged case,
- seed, grid sizes, prior bounds, model settings, software commit and Python/NumPy/SciPy versions,
- machine-readable pass/fail status and a plain-language report.

## Limitations

These criteria assess numerical stability under grid refinement only. They do not test repeated-sampling coverage, prior robustness, misspecification robustness, selection effects, or realism relative to LVK data. Those questions require separate experiments and must be reported separately.
