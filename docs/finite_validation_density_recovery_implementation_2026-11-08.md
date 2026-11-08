# Finite-validation density-recovery implementation status

**Date:** 2026-11-08  
**Branch:** `research/independent-posterior-audit-claim-fix`  
**Status:** implementation committed; CI execution on the final branch head is not yet verified

## Implemented components

- `src/dark_siren_robustness/finite_validation_density_recovery.py`
  - Declared two-overdensity synthetic parent intensity on (redshift between 0.01 and 0.40).
  - Parent-law curvature centering (the parent-intensity-weighted expectation of squared normalized redshift), not a uniform-law shortcut.
  - Numerically calibrated logistic inclusion intercept for each declared retention target.
  - Common random numbers for parent frames and retention uniforms.
  - Nested simple-random validation samples, independent of retention labels.
  - Validation objects excluded from the analysis catalogue.
  - Catalogue-only, known-inclusion IPW, estimated linear IPW, estimated quadratic IPW, and parent-sample KDE arms.
  - Frozen Scott multipliers 0.15, 0.40 and 1.00.
  - IAE, ISE, maximum absolute density error, Kish ESS, entropy ESS, maximum normalized weight, and top-1%/5% weight mass.
  - Replicate-paired estimated-IPW minus catalogue-only contrasts, percentile bootstrap intervals and paired t-based intervals.
  - Raw parent-unit records, attempt ledger, analytic target, immutable configuration and SHA-256 manifest.

- `src/dark_siren_robustness/finite_validation_density_recovery_audit.py`
  - Does not import the production benchmark runner, inclusion fitter, validation splitter or production KDE estimator.
  - Reconstructs DGP probabilities, retention labels, parent partitions, nested validation and common-random-number invariants.
  - Reconstructs logistic fitted covariance, reflected weighted KDEs, density errors, weight diagnostics and paired contrasts.
  - Verifies manifest checksums, complete factorial attempt counts and summary consistency.
  - Explicitly limits its claim to numerical artifact reconstruction.

- CLI entry points:
  - `scripts/run_finite_validation_density_recovery.py`
  - `scripts/audit_finite_validation_density_recovery.py`

- Tests:
  - `tests/test_finite_validation_density_recovery.py`
  - Includes smoke reconstruction and a checksum-tampering rejection test.

## Reproduction commands

```bash
python -m pytest tests/test_finite_validation_density_recovery.py -q

python scripts/run_finite_validation_density_recovery.py \
  --replicates 1 --n-parent 1200 --grid-points 201 \
  --output-dir artifacts/finite_validation_density_recovery_smoke

python scripts/audit_finite_validation_density_recovery.py \
  --artifact-dir artifacts/finite_validation_density_recovery_smoke \
  --output artifacts/finite_validation_density_recovery_smoke/audit_report.json
```

The default runner uses the 100-replicate primary design and all three predeclared bandwidth multipliers. It is computationally heavier and writes a large raw parent-unit artifact.

## Scientific interpretation limits

This stage evaluates synthetic merger-weighted host-density reconstruction only. It does not implement or validate the estimated-inclusion (H_0) likelihood, posterior calibration, repeated-sampling coverage, real galaxy-survey completeness, or LVK/GWTC event analysis. The inverse observed-information covariance remains a frequentist local covariance estimate, not a Bayesian posterior over inclusion coefficients.

The 0.15 Scott multiplier is a frozen candidate from a separate training/held-out density pilot. It is not a universal optimum; the 0.40 and 1.00 sensitivity arms remain mandatory.

## CI status

The workflow was extended to run the focused tests, a small density-recovery smoke experiment, and the independent artifact auditor, and to upload the resulting smoke artifact. A workflow run associated with an intermediate branch head was later reported as cancelled with no job records. The final branch head has not yet returned a completed CI result. Therefore the new implementation is **not reported as tested or passed** until a completed Actions run and audit artifact are inspected.


## Independent-audit completeness and confirmatory execution gate (2026-11-08)

The density-recovery auditor was strengthened before launching the full design:

- It now reconstructs the complete expected factorial attempt-key set from the immutable summary design and rejects missing, duplicate, or unexpected attempts.
- It accepts only the declared attempt lifecycle states: SUCCESS, NON_EVALUABLE, FAILED, and NOT_RUN.
- Successful estimated-inclusion attempts now store the logistic negative log-likelihood objective, design rank, and design condition number. The auditor independently recomputes these from raw validation labels and stored coefficients, alongside the score, observed information, covariance, and fitted-probability diagnostics.
- Parent-unit partitions are indexed once by experiment key before iterating through attempt rows. This avoids repeatedly scanning the entire multi-million-row raw table for each KDE attempt and makes the full audit computationally feasible in principle.
- Adversarial tests now remove an expected attempt while refreshing the manifest and alter a fitted objective while refreshing the manifest; both must be rejected by the auditor.

A separate workflow, `.github/workflows/finite-validation-density-confirmatory.yml`, has been added for the predeclared 100-replicate primary density-recovery design (N=4,000; gamma -1.2, 0, +1.2; retention 0.25; validation sizes 100, 300, 1,000; linear/quadratic fits; Scott multipliers 0.15, 0.40, 1.00; 801-point density grid). It runs the focused tests, generates the full raw artifact, independently audits it, and uploads the result. The workflow is isolated from the general test matrix and can also be manually dispatched.

**Status at documentation time:** the workflow has been added but the full 100-replicate experiment has not yet completed. The additional completeness/objective tests and the confirmatory workflow must pass before the density stage is declared complete. This remains a density-recovery experiment only; it does not implement H0 inference or posterior calibration.


## Critical RNG-stream defect discovered and corrected (2026-11-08)

Inspection of the Python 3.12 CI smoke artifact exposed a design-breaking defect in the original replicate seeding:

`default_rng(seed + stream_offset * replicate)` gives the **same seed for every stream at replicate 0**, because every offset is multiplied by zero. In the smoke artifact, parent redshifts and retention uniforms consequently had Pearson correlation 0.9969 and Spearman correlation 1.0. All 54 estimated-IPW attempts failed with enormous fitted coefficients; there were zero paired estimated-IPW/catalogue-only pairs. The independent artifact auditor had previously passed because it checked DGP label consistency and common-random-number pairing, but did not independently reconstruct the seed streams.

This is not a benign optimizer failure. The validation/retention design for replicate 0 was not independent, so the smoke artifact cannot support density-recovery conclusions. The same defect affects replicate 0 in the pre-fix 100-replicate workflow run #42474967403 (head `0b637c3cf93c673aee9dcb6b0d046ffd1e6f9b6d`); that run must be treated as **invalid for confirmatory interpretation**, even if its old auditor later reports PASS.

The runner now uses disjoint deterministic streams:
- `SeedSequence([seed, replicate, 1])`: parent redshifts;
- `SeedSequence([seed, replicate, 2])`: validation permutation;
- `SeedSequence([seed, replicate, 3])`: retention uniforms.

The independent auditor separately reconstructs all three streams from the declared seed protocol and verifies every parent redshift, validation rank, and retention uniform against the raw artifact. The stream protocol is recorded in both the immutable config and summary. An adversarial test modifies a raw uniform while preserving its Bernoulli label and refreshing the manifest; the auditor must still reject it.

The objective/coefficient tampering fixtures have been returned to a single replicate because the corrected stream protocol removes the replicate-zero collision; they no longer rely on five replicate attempts to find a successful fit.

**Next gate:** run the focused tests on the corrected code, then rerun the full 100-replicate experiment and inspect its independent RNG-stream audit report. No pre-fix density-recovery results should be pooled with the corrected experiment.
