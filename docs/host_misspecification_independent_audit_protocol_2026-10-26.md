# Independent host-misspecification audit protocol

**Added:** 2026-10-26  
**Scope:** audit of the existing controlled synthetic experiment, not validation against real gravitational-wave or galaxy-catalogue data.

## Why add a second implementation?

The original workflow checked output dimensions, finite summary metrics, and coverage bounds. Those checks could pass even if the random-number stream, posterior calculations, or aggregation were wrong. The independent auditor strengthens this by not importing the production runner or posterior implementation.

## Checks performed

`scripts/audit_host_misspecification_audit.py`:

1. Regenerates the declared NumPy random-number stream from the command-line seed and checks scenario order, replicate index, latent host index/redshift, observed distance, and noise scale against every raw CSV row.
2. Independently reconstructs the four scenario likelihoods from Gaussian mixture densities on the declared H0 grid, including the 30% hidden-host component in the modelled-hidden-host arm.
3. Recomputes posterior medians, central 68% and 95% intervals, interval widths, and coverage for every replicate.
4. Reconstructs all seven reported summary metrics for each scenario.
5. Fails closed on row-count, design, provenance, posterior, or summary discrepancies.

The regression test creates a fresh 200-replicate-per-scenario dataset and verifies the full independent audit. The workflow also runs the auditor on its 1,000-replicate-per-scenario outputs before uploading the artifact.

## Interpretation boundary

This is independent computational reconstruction of the same declared toy model; it is not an independent scientific replication of the model assumptions. Both implementations still share the same scientific design (low-redshift relation, Gaussian distance noise, fixed H0 truth, and a finite H0 grid). A PASS therefore supports provenance and implementation consistency only. It does not establish realistic host-population modelling, calibrated uncertainty under misspecification, or astrophysical validity.

## Reproduction

```bash
python scripts/run_host_misspecification_audit.py --replicates 1000 --seed 20261027 --output-dir results/host-misspecification-audit
python scripts/audit_host_misspecification_audit.py --replicates 1000 --seed 20261027 --input-dir results/host-misspecification-audit
```
