# Validation plan

## Stage 0 — Reproducibility

- Record Python and dependency versions.
- Use explicit random-number generators and record seeds for simulations.
- Store experiment configurations with outputs.
- Keep exploratory runs separate from confirmatory runs.

## Stage 1 — Unit and known-answer tests

- Check cosmological distance functions against limiting cases and an independent numerical reference.
- Check probability functions for valid bounds and invalid-input handling.
- Build analytically tractable likelihood examples before implementing the full model.
- Test numerical stability and normalization.

## Stage 2 — Inference calibration

- Define the exact repeated-sampling experiment and the meaning of each interval.
- Estimate empirical coverage with Monte Carlo uncertainty.
- Run simulation-based calibration under data generated from the inference model.
- Do not interpret SBC as proof of scientific realism or robustness under misspecification.

## Stage 3 — Controlled robustness

- Vary one mechanism at a time before studying interactions.
- Include random incompleteness, distance-dependent incompleteness, and host-weight misspecification.
- Compare methods on matched simulated datasets where feasible.
- Report bias, RMSE, interval width, coverage, failure rate, and computational cost.

## Stage 4 — Independent checks

- Compare key calculations with an independent implementation or analytic limiting case.
- Repeat selected experiments with new random seeds.
- Verify that plots and tables can be regenerated from saved configurations and raw outputs.

## Current status

This plan is a specification, not evidence that the tests or experiments have already passed.
