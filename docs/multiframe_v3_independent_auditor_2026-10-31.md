# Independent audit protocol for the v3 multi-frame pilot

**Status:** software-validation protocol and implementation, 2026-10-31  
**Evidence class:** pilot-only; no confirmatory statistical claims

## Purpose

The v3 runner writes finite frames, paired outer samples, three interval-method outputs, frame summaries, and equal-frame-weighted descriptive summaries. This auditor independently reconstructs those outputs rather than trusting the runner's derived fields. It is intended to detect accidental corruption and deliberate artifact tampering, including cases where an attacker refreshes the CSV checksum after editing a value.

Run from the repository root:

```bash
python scripts/run_multiframe_v3_pilot.py \
  --manifest configs/multiframe_v3_smoke_manifest.json \
  --output-dir artifacts/multiframe_v3_smoke
python scripts/audit_multiframe_v3_pilot.py \
  --output-dir artifacts/multiframe_v3_smoke \
  --report artifacts/multiframe_v3_smoke/audit_report.json
pytest -q tests/test_multiframe_v3_runner.py tests/test_multiframe_v3_auditor.py
```

## Independent checks

1. **Artifact contract:** schema, pilot-only status, embedded manifest, exact supported method set, and four output checksums.
2. **Finite frames:** regenerate every scenario/frame from the declared root seed and frame coordinate; verify all unit IDs, values, seed labels, and canonical float64 hashes.
3. **Outer sampling:** reconstruct each SRSWOR sample from its independent sample seed; verify exact ordered indices, no replacement, digest, sample mean, finite-frame target, and frame provenance.
4. **Pairing:** require every method at an outer-sample key to use the same ordered sample indices.
5. **Method streams:** independently derive each method coordinate from SHA-256 and reconstruct the method seed.
6. **Intervals:** independently recompute Wald and ordinary percentile endpoints; separately replay the pseudo-population multinomial draws, SRSWOR resamples, studentized statistics, quantiles, and failure categories.
7. **Completeness:** reject duplicate, omitted, and unexpected outer keys, method rows, frame units, and summary cells.
8. **Summaries:** reconstruct per-frame conditional coverage, interval success, unconditional successful-and-covering probability, mean width, and equal-frame-weighted descriptive aggregates.

The auditor does not import the pilot runner's interval-construction code. It does reuse the manifest validator, finite-frame generator, and public seed-schedule functions; the replay of sampling and interval calculations is independently implemented in the auditor. This is independent artifact reconstruction, not an independent scientific replication or proof of statistical validity.

## Adversarial regression tests

The tests alter a frame value, sample index, interval endpoint, or replicate-row completeness and then refresh the corresponding checksum. The auditor must still reject each altered artifact through its semantic reconstruction checks. A normal checksum-only validator would not detect those edits after checksum refresh.

## Interpretation and limitations

- The smoke manifest is deliberately tiny: three distribution scenarios, two frames per scenario, population size 80, two sample sizes, three outer repetitions per cell, and 25 inner bootstrap attempts.
- This is a software integration test, not an estimate of coverage, a method ranking, or evidence of finite-sample superiority.
- Equal-frame weighting is descriptive across the realized finite frames. With only two frames per scenario, it is not a reliable population-of-frames uncertainty estimate.
- The confirmatory design, number of frames, outer repetitions, inner repetitions, primary contrast, multiplicity strategy, stopping rules, and precision targets remain to be justified and frozen before confirmatory execution.
- The benchmark targets the arithmetic mean of a realized synthetic finite population. It does not establish validity for nonlinear dark-siren posterior functionals, real gravitational-wave events, galaxy catalogues, or an estimate of H0.
