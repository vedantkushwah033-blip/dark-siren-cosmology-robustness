# Fixed-analysis-catalogue pilot RNG-stream incident

**Date:** 2026-11-08  
**Affected workflow run:** 42475754461  
**Affected artifact:** `fixed-analysis-catalogue-density-pilot-py3.12`  
**Affected commit:** `758ba999443133a9703b2bc1d8206592612202a9`  
**Disposition:** technically audited, scientifically invalidated; do not use its numerical results.

## Finding

The first five-replicate FAC pilot completed and its independent auditor returned PASS. A post-run inspection of raw validation-frame records then revealed deterministic dependence between the generated redshifts, validation ordering, and retention uniforms for replicate zero.

The cause was the RNG initialization pattern:

```python
default_rng(seed + stream_offset * replicate)
```

At `replicate == 0`, every stream received exactly the same seed. Consequently, validation redshift draws, validation-order permutation, and retention uniforms were not independent. The first validation subset showed a strong spurious redshift/retention association; the logistic fits failed for all three curvature settings and all validation sizes in replicate zero. The old auditor reconstructed the DGP labels correctly from the stored arrays, but did not detect that the arrays themselves had been generated from aliased random streams.

This is a scientific design defect, not merely a cosmetic RNG issue. The affected artifact is therefore not suitable for inference, estimator comparisons, or project claims, despite its internal artifact-audit PASS.

## Correction

The runner now uses independent NumPy SeedSequence streams:

- stream 1: analysis-frame redshifts;
- stream 2: validation-frame redshifts;
- stream 3: validation-frame permutation;
- stream 4: analysis-frame retention uniforms;
- stream 5: validation-frame retention uniforms.

Each stream is initialized with `SeedSequence([seed, replicate, stream_id])`. The protocol is recorded in the run summary and checked by the independent auditor.

A regression test checks that, for a deterministic 1,000-unit validation sample, validation redshifts are not spuriously correlated with either retention uniforms or validation ranks; the analogous analysis-frame redshift/uniform check is also included. These correlation checks are regression guards, not proofs of independence.

## Re-run requirement

The FAC pilot must be rerun from the corrected commit. Its new artifact must pass:
1. focused and adversarial unit tests;
2. the complete factorial runner;
3. independent raw-frame, DGP, fit, KDE, attempt-ledger and paired-summary audit.

The old artifact remains useful only as a reproducible record of the discovered defect. Its reported IAE contrasts and 80% estimated-fit success rates must not be interpreted scientifically.
