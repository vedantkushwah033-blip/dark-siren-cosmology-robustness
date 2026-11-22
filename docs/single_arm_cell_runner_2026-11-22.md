# Single-cell execution from frozen shared inputs

**Status:** execution primitive; no coordinator, independent cell auditor, or lifecycle transition is included.

The script `scripts/run_arm_cell_from_bundle.py` consumes a bundle frozen by `shared_input_manifest.json` and evaluates exactly one (replicate, gamma, validation_size, arm) key using the independently tested arm-cell evaluator. It writes a cell-scoped NPZ, JSON summary, and hash manifest through a temporary directory followed by an atomic directory rename. Existing output directories are never overwritten.

The output status is deliberately `COMPUTED_NOT_AUDITED`. The script does not write the four-state checkpoint and cannot declare a cell SUCCESS or NON_EVALUABLE. A future coordinator must independently audit the cell artifact first.

## Invocation

```bash
python scripts/run_arm_cell_from_bundle.py \
  --bundle-dir artifacts/replicate_0001_gamma_0p5 \
  --output-dir artifacts/cells/replicate_0001_gamma_0p5_n300_E \
  --replicate 1 --gamma 0.5 --validation-size 300 \
  --arm E_uncertainty_integrated_ipw
```

The bundle must already contain config, observable inputs, truth inputs, and the verified shared-input manifest emitted by the smoke runner. Requested gamma, validation size, seed, and config are checked against frozen provenance. For D/E, validation labels are loaded; A/B/C do not pass validation labels to the evaluator.

## Explicit non-claims

- This does not yet generate shared inputs without running the legacy five-arm smoke.
- It does not independently audit weights, KDE, likelihood, or posterior quantiles.
- It does not update lifecycle state, schedule sibling cells, or implement restart.
- It does not change the synthetic low-redshift scope or authorize confirmatory calibration.
