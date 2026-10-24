# Validation status — 2026-10-24

## Scope

This note records the verified state of the synthetic dark-siren robustness repository. It distinguishes software test status, numerical checks, and scientific validation. The project is a controlled synthetic-model study; it is not an LVK observational measurement.

## Verified CI state

The latest fully completed factorial-audit workflow observed during the 2026-10-24 continuation passed on commit `9541cd7a` for Python 3.10, 3.11, and 3.12:

- Factorial audit: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40095985757
- Its uploaded artifacts are listed by that run and expire 2027-01-22.

The two subsequent commits strengthen the artifact validator and add regression tests for invalid host indices and observations outside the declared detection interval:

- `eaa4bcf4`: validator checks observation host index and strict distance bounds.
- `88f8282e`: adds regression tests for both checks.

At the latest status check, CI for both commits was still in progress; do not treat those newest changes as verified until the corresponding runs finish successfully:
- https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40096177016
- https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40096182977

Earlier successful CI runs, including the suite at commit `8634e8aa`, remain evidence only for the code versions tested by those runs. Software test success does not establish astrophysical validity.

## Grid-refinement findings

The existing grid-refinement audit covered 21,600 posterior fits across 36 design cells. Its largest observed absolute difference in the upper 95% credible endpoint between 1,001- and 4,001-point grids was 0.13045331 in the project's H0 units. This is a warning against assuming the coarse grid is uniformly adequate.

A targeted synthetic case was subsequently examined at increasingly fine grids, up to 64,001 points. For that specific case, the upper 95% endpoint stabilized as the grid was refined. This single case is diagnostic, not proof of convergence for all possible observations or design cells.

## Integration cross-check

The targeted convergence output records posterior area from both trapezoidal and Simpson integration and their absolute difference. The test verifies that the quantities are finite/positive where applicable and that the reported difference equals the absolute difference of the recorded areas.

The test intentionally does not force the two integration methods to agree within an arbitrary universal tolerance. Instead, the numerical discrepancy should be inspected and interpreted at each relevant grid resolution. Agreement between two quadrature rules on the same sampled density is useful but is not, by itself, proof that the posterior or its quantiles are correct.

## Factorial artifact validator

The factorial validator checks:

- Expected row counts and exact configuration keys for raw, summary, and paired-contrast tables.
- Complete replicate indices and unique configurations.
- A single simulated observation per truth/replicate shared across inference configurations.
- Finite numeric values and correctly ordered credible intervals.
- Recomputed bias, MAE, RMSE, interval coverage, Wilson coverage intervals, and interval widths.
- Paired mean contrasts and percentile-bootstrap intervals reproduced by a separate code path.

The bootstrap recomputation uses the same resampling design and NumPy random-number-generator family as the experiment implementation. It is an independent code path, not an independent statistical methodology.

The latest test suite exercises the validator on a controlled fixture. **The full 300-replicate factorial artifact has not been verified as completed in this status note.** Do not report its planned row counts as generated results until the manual workflow has completed and its uploaded CSV artifacts have been inspected.

## Next required milestones

1. Run the manual workflow `.github/workflows/full-prior-threshold-factorial.yml` with 300 replicates, 4,001 grid points, and seed 20261028.
2. Confirm the workflow conclusion and download the three CSV artifacts.
3. Inspect the results, including Monte Carlo uncertainty in coverage and paired differences.
4. Perform a final scientific review of assumptions, selection modelling, prior truncation, and interpretation.
5. Keep synthetic conclusions separate from any future real-data analysis.

## Scientific limits

The model uses a simplified Gaussian distance-threshold setup. It does not include the full gravitational-wave detector selection function, SNR-dependent efficiency, duty cycles, waveform effects, multiple-event population inference, redshift uncertainty, or calibrated galaxy-catalogue incompleteness. Passing tests establishes software behaviour against specified checks; it does not establish that these simplifying assumptions are realistic.
