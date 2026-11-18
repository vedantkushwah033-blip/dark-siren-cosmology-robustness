# Cross-gamma audit completion record

**Date:** 2026-11-19  
**Status:** PASS  
**Scope:** paired synthetic low-redshift Euclidean catalogue-inclusion density stress pilot only.

## Provenance

- Workflow run: [44369630799](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44369630799)
- Workflow head: `17ad47e3c4211abb7f0d3317c535d7fd4cbe6186`
- Cross-gamma job: `135042736572`
- Cross-gamma artifact: `14462504250` (`corrected-quadratic-logit-curvature-cross-gamma-summary`)
- Artifact size: 375,146 bytes
- ZIP SHA-256: `c4e86b59dfa37f93b0fd0682b00c440615544fb0f4e10a4989120a7f257b88fc`

The artifact was downloaded after job completion. Its local SHA-256 was independently computed and matched the GitHub artifact digest.

## Audit results

- Gamma cells: 7/7 PASS for gamma in `{-2, -1, -0.5, 0, 0.5, 1, 2}`.
- Replicates: 20/20 successful in each cell; 140/140 total.
- Failed replicates: 0.
- Child raw-input audits: 140/140 passed.
- Manifest entries: 244 per cell; 1,708 total verified.
- Aggregate density metric rows reconstructed: 91,140.
- Realized-Scott-width rows reconstructed: 12,180.
- Cross-gamma summary rows: 4,557.
- Paired arm-contrast rows: 2,205.
- Paired gamma-versus-zero contrast rows: 756.
- Paired seed schedule verified: true.
- Paired parent and validation redshifts verified: true.
- Overall audit JSON status: `PASS`.

The cross-gamma audit independently rechecked configuration/calibration consistency, paired seeds and samples, and per-gamma raw-input reconstruction. It therefore closes the workflow-level computational audit that was pending in the initial density-grid readout.

## Scientific interpretation boundary

This is a synthetic density-recovery stress pilot. It does not perform bandwidth selection, H0 inference, posterior calibration, simulation-based calibration, or fixed-truth interval-coverage validation. The oracle inclusion-probability arm is a simulation diagnostic, not an operational estimator. The completed audit is evidence that the declared computational design and paired artifacts were reconstructed consistently; it is not evidence of estimator superiority or observational validity.

## Next gate: H0 propagation

The repository already contains a shared-node, one-sided selected multi-event likelihood and an independent-validation H0 smoke runner. However, that runner's current generator uses the older floor-logistic catalogue-inclusion law. The completed curvature grid uses a calibrated quadratic-logit inclusion law. Reusing the old H0 smoke as if it tested the curvature grid would silently change the data-generating process.

Before a curvature-specific H0 pilot:
1. Add a separately seeded independent validation-frame generator using the exact same calibrated quadratic-logit inclusion function and parent-redshift law as the density grid.
2. Keep the validation fitter restricted to observed validation redshifts and retention labels; truth probabilities remain scoring-only.
3. Reuse the existing proper logistic posterior quadrature and shared-node likelihood only after matching the generator's host-rate weighting and one-sided GW selection.
4. Add deterministic tests for the quadratic DGP, truth firewall, event-selection arithmetic, and likelihood marginalization ordering.
5. Run a small H0 smoke with explicit numerical convergence checks, then independent raw-artifact audit.
6. Only after that, design repeated fixed-truth coverage and SBC as distinct experiments.

No H0 result is claimed by this audit record.
