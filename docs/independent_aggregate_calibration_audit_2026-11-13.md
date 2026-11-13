# Independent aggregate calibration audit

**Status:** implementation added on a separate draft branch; no CI or artifact execution is claimed yet.

The inclusion-posterior calibration pilot already stores a replicate ledger and
an aggregate JSON summary. The existing aggregate auditor verifies lifecycle
counts, selected posterior summaries, and per-replicate audit attestations, but
it does not independently reconstruct the paired validation-size contrast
table or the failure-inclusive bounds from the raw ledger.

Run the additional audit after the original per-replicate and aggregate
orchestration audit:

```bash
python scripts/audit_inclusion_posterior_aggregate.py PATH/TO/CALIBRATION_ARTIFACT
```

The script independently rebuilds:

1. Attempt-level lower and upper bounds for each 68%, 90%, and 95% coverage
   (or conditional prior-predictive interval-inclusion) rate. SUCCESS rows
   contribute their observed binary result; FAILED and NON_EVALUABLE attempts
   remain unknown; NOT_RUN rows are excluded from the attempted denominator.
2. Every paired validation-size contrast in the summary, matching on replicate
   identity and retaining only pairs with SUCCESS at both sizes. It recomputes
   the larger-minus-smaller differences, paired Monte Carlo standard errors,
   and normal-approximation intervals from the ledger.

The audit is intentionally separate from the calibration runner and imports no
runner functions. It checks aggregate arithmetic and pairing, not the
scientific validity of the DGP, posterior, or stored per-replicate PASS
attestations. It does not replace raw-artifact reconstruction. A PASS is not
evidence of calibration, estimator superiority, or observational validity.

The tests include a deliberately failed paired attempt and tampering checks to
ensure that changing a stored paired mean or an all-attempt bound is detected.
