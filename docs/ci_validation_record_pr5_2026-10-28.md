# Exact-commit CI validation record — PR #5 — 2026-10-28

## Provenance

- Repository: `vedantkushwah033-blip/dark-siren-cosmology-robustness`
- Pull request: [#5 — Fix dataset posterior normalization for off-grid prior boundaries](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/pull/5)
- Workflow: `tests`
- Run: [40699584953 (run #1717)](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40699584953)
- Exact tested head SHA: `7e42b7ee38842cc95071800aa4ad7578b17f8e28`
- Head branch: `fix/prior-boundary-grid-normalization`
- Conclusion: **success**
- Matrix: Python 3.10, 3.11, 3.12
- All three `numerical-diagnostics` jobs completed successfully; each reported 45 completed steps and no failed steps.

## Targeted checks relevant to this revision

The passing revision includes the regression test that compares the standalone selection-conditioning helper plus dataset-level normalization against the existing `posterior_h0_selected` implementation for the same finite-interval synthetic model. It also includes the toy Gaussian-threshold efficiency checks and prior-boundary normalization tests.

## What this establishes

The configured GitHub Actions workflow passed on the exact SHA above across the three supported Python versions. This supports software regression confidence for the tested code paths.

## What this does not establish

- It does not establish astrophysical validity or validate a realistic gravitational-wave detector selection function.
- It does not implement event-rate inference or a full selection-aware hierarchical population model.
- It does not validate the proposed selection-aware catalogue-incompleteness benchmark, which remains a protocol/design sketch.
- It does not establish scientific novelty, external CHIMERA v2.2 validation, or a real-data H0 measurement.
- The workflow result applies to SHA `7e42b7ee38842cc95071800aa4ad7578b17f8e28`. Any later commit requires its own CI result before being described as validated.

## Next correctness gate

Before scaling a selection-aware incompleteness experiment, independently verify that the detection-efficiency denominator is integrated over the same population mixture, weights, noise convention, and detection rule as the likelihood numerator. Distinguish an event-specific host-candidate catalogue from the population distribution governing detection. Then test empirical simulated detection frequencies against the analytic efficiency in a small, fixed-seed toy experiment.
