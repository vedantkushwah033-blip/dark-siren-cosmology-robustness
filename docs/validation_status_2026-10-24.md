# Validation status — 2026-10-24

## Scope

This note records the verified state of the synthetic dark-siren robustness repository. It distinguishes software test status, numerical checks, and scientific validation. The project is a controlled synthetic-model study; it is not an LVK observational measurement.

## Verified CI state

The latest completed verification during this continuation includes the focused factorial-audit workflow and focused tests on commit `6ccdd1ad`, with all three Python versions (3.10, 3.11, and 3.12) passing:

- Factorial audit: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40097313503
- Focused tests: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40097313497
- The audit uploaded compact experiment artifacts for all three Python versions.

The README/protocol updates at `91fe3997` and `8c16ac1d` subsequently passed both the focused `tests` workflow and the broader `Test suite` workflow:
- https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40097368888
- https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40097368886
- https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40097372663
- https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40097372541

Later commits `05cbee3f`, `555b9a40`, and `87843974` add an end-to-end compact artifact-validation fixture, a check that summary replicate counts match raw data, and regression tests for both behaviours. Their CI runs were still in progress or queued at the last check; do not treat those newest tests as passed until the corresponding Actions runs finish successfully.

The two subsequent commits strengthen the artifact validator and add regression tests for invalid host indices and observations outside the declared detection interval:

- `eaa4bcf4`: validator checks observation host index and strict distance bounds.
- `88f8282e`: adds regression tests for both checks.

Follow-up CI after those commits:
- The focused `tests` workflow passed on `eaa4bcf4` (run 40096177001) and `88f8282e` (run 40096182919).
- The broader `Test suite` runs and the latest numerical-diagnostics workflow were still running during a later check; inspect the live Actions page before treating the newest commits as verified:
  - https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40096177016
  - https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40096182977
  - https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40096243517

The factorial runner was then optimized to calculate each uncorrected posterior once per observation/prior, reusing that result across the four assumed-threshold cells. This is an algebraically exact reuse because the uncorrected likelihood does not use the assumed threshold; the existing threshold-invariance regression test checks that property. Each completed run now also writes `experiment_manifest.json` with seed, grid, design constants, artifact row counts, software versions, and the CI commit SHA when available. The dispatch workflow validates the manifest against its requested inputs and uploads it alongside the CSV tables.

A workflow-file edit briefly introduced literal escaped newline characters in the artifact path, causing the workflow parser to fail before any job ran. That syntax was corrected in commit `4948be35`; the workflow must be rechecked on the corrected revision before use.

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
