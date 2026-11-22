# Preparation-only shared-input bundles

**Status:** engineering workflow documentation. This is not a confirmatory protocol.

The preparation command creates one synthetic replicate/curvature bundle containing the common parent population, observed catalogue, GW events, validation-frame truth, and nested validation samples. It does not fit an inclusion model, construct KDEs, evaluate a likelihood, or calculate an H0 posterior.

## Required independent random-number streams

Supply four distinct, nonnegative seeds:

- `--seed`: parent population and catalogue-thinning stream.
- `--event-seed`: GW event generation stream.
- `--validation-seed`: independent validation-frame generation stream.
- `--validation-sample-seed`: sampling the nested validation prefix.

The command rejects missing or reused stream IDs before creating the bundle. Across paired curvature bundles in one replicate, preserve the exogenous streams and H0 truth. Across independent replicates, use distinct stream IDs and an independently specified H0 truth.

## Prepare a bundle

Example:

```bash
python scripts/prepare_shared_input_bundle.py \
  --seed 20261122 \
  --event-seed 20261123 \
  --validation-seed 20261124 \
  --validation-sample-seed 20261125 \
  --n-parent 1200 \
  --n-detected 30 \
  --validation-frame 500 \
  --validation-sizes 100,300,500 \
  --curvature-gamma 1.0 \
  --inclusion-degree 1 \
  --h0-true 70 \
  --h0-points 1121 \
  --output-dir artifacts/replicate_001_gamma_p1
```

The output directory is immutable: an existing destination is never overwritten. Inputs are first written into a temporary sibling directory, hash-indexed in `shared_input_manifest.json`, then published with an outer `manifest.json`. The config status is `shared_inputs_only_no_inference`; the outer manifest also records `preparation_only=true` and `inference_executed=false`.

## Audit the raw bundle

```bash
python scripts/audit_prepared_shared_input_bundle.py artifacts/replicate_001_gamma_p1
```

The raw-input auditor independently verifies manifest hashes and exact file sets, parent/catalogue retention joins, quadratic-logit inclusion probabilities, validation-frame joins and nested prefixes, GW host/event truth joins, and the selected-event generation trace. It also rejects inference outputs in the preparation-only directory.

## Execute and audit one arm cell

```bash
python scripts/run_arm_cell_from_bundle.py \
  --bundle-dir artifacts/replicate_001_gamma_p1 \
  --output-dir artifacts/cells/r1_g1_n100_A \
  --replicate 1 --gamma 1.0 --validation-size 100 \
  --arm A_catalogue_only

python scripts/audit_arm_cell_from_bundle.py \
  --bundle-dir artifacts/replicate_001_gamma_p1 \
  --cell-dir artifacts/cells/r1_g1_n100_A
```

For a preparation-only bundle, the per-cell audit independently rebuilds the target arm from raw inputs using the separate numerical reference implementation. It checks weights, coefficient posterior quadrature when relevant, reflected KDEs, bandwidths, selected multi-event likelihoods, detection efficiencies, posterior summaries, and weight diagnostics. It does not import the production arm evaluator, KDE, likelihood, or posterior implementation.

The same audit route is used by the serial audit-gated cell coordinator. A successful lifecycle state requires the independent audit to pass.

## Scientific limits

These commands create and validate synthetic low-redshift Euclidean inputs and their numerical execution artifacts only. They do not establish H0 calibration, SBC coverage, estimator superiority, real-survey validity, an observational H0 measurement, or LVK/GWTC public-data analysis. No confirmatory 100/300-truth-per-cell experiment is authorized by this workflow alone. Freeze the scientific design, numerical convergence criteria, lifecycle accounting, and preregistered policy before confirmatory calibration.
