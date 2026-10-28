# Selection-function misspecification pilot — implementation record

**Status:** implementation added; numerical results are not yet validated in CI or from a separately recorded pilot run. This is a synthetic methodological experiment, not an H0 measurement.

## Implemented design

Module: `src/dark_siren_robustness/selection_misspecification.py`

Runner: `scripts/run_selection_misspecification.py`

Independent artifact auditor: `scripts/audit_selection_misspecification.py`

The generator uses the existing true catalogue inclusion law
`p_true(z) = 0.10 + 0.90 / (1 + exp((z - 0.22) / 0.05))`.
The assumed law affects inverse-inclusion weights only; it does not generate the catalogue. Each replicate shares one parent population and one GW-like detected-distance sample across five arms:

1. `oracle`: parent redshifts, as an idealized reference;
2. `catalogue_only`: observed incomplete catalogue, unweighted;
3. `correct_weighting`: inverse inclusion probabilities under the true declared law;
4. `shifted_transition`: assumed midpoint 0.27 instead of 0.22;
5. `shallow_slope`: assumed width 0.075 instead of 0.05.

The host-density target for integrated squared error is the declared analytic latent galaxy intensity multiplied by the same linear merger-weight function and normalized on the redshift grid. This target is a toy truth, not an empirical galaxy distribution.

## Outputs and reproduction

Run a quick smoke/pilot:

```bash
python scripts/run_selection_misspecification.py \
  --replicates 2 --n-parent 1200 --n-detected 8 \
  --seed 55129 --h0-points 51 --z-points 201 \
  --output-dir artifacts/selection_misspecification_smoke
python scripts/audit_selection_misspecification.py artifacts/selection_misspecification_smoke
```

A 30-replicate pilot at the default grids is the configured initial research pilot, not a precise coverage study:

```bash
python scripts/run_selection_misspecification.py --replicates 30 \
  --output-dir artifacts/selection_misspecification_pilot
python scripts/audit_selection_misspecification.py artifacts/selection_misspecification_pilot
```

The runner records raw replicate/arm metrics, aggregate metrics, runtime, software versions, commit when available, and SHA-256 hashes in a manifest. The independent auditor verifies hashes, paired arm structure, posterior normalization and summaries reconstructed from saved posterior density vectors, coverage flags, and arm aggregates. The tamper test rehashes a corrupted artifact to ensure that a matching manifest alone cannot conceal an internally inconsistent posterior.

## Metrics

Replicate records include density integrated squared error, H0 posterior mean and equal-tail interval, coverage indicator, posterior-grid values, weight ESS and ESS fraction, maximum and 95th-percentile normalized weight, assumed selection probability minimum, and minimum GW-like detection efficiency. Aggregates include bias, RMSE, empirical coverage numerator/denominator, interval width, density ISE, and weight diagnostics.

## Important limitations and gates

- The toy uses low-redshift Euclidean distance scaling, Gaussian distance noise, a fixed detection threshold, and a fixed number of detected events.
- It does not include realistic detector response, event-rate marginalization, peculiar velocities, redshift uncertainty, calibration uncertainty, or real catalogue systematics.
- The oracle is not an available-data method.
- Correct inverse weighting assumes known probabilities and positivity; it is not a learned selection model.
- Severe weight instability is reported rather than silently clipped; no stabilization arm is included.
- The current runner fails closed if a replicate raises an exception. It does not silently drop failures, but it does not yet produce a partial-run failure ledger; a failed run must be documented and debugged before a successful run is claimed.
- Do not claim that any method wins until the pilot has run, its artifacts pass the independent auditor, Monte Carlo uncertainty is reported, and bandwidth/grid sensitivities are examined.
