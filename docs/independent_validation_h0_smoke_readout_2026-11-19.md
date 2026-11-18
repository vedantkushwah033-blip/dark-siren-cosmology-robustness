# Independent-validation selected-H0 smoke readout

**Date:** 2026-11-19  
**Status:** workflow PASS; independent raw-artifact audit PASS; coefficient-quadrature convergence PASS.  
**Scope:** one low-redshift Euclidean synthetic realization. This is a numerical/integration smoke, not SBC, fixed-truth coverage, or an H0 measurement.

## Exact execution record

- Workflow run: 44404971527
- Commit: 9faf7ec62ecb38bf3fb9afc08a589304b3c42dd8
- Artifact: 14474930145 (`independent-validation-h0-smoke-py3.12`)
- Artifact ZIP size: 242,750 bytes
- ZIP SHA-256: `2383b6ea1ce97a4b586a29f2a3db7234330916825f5104f778180179c0e76e27`
- GitHub artifact digest: `sha256:2383b6ea1ce97a4b586a29f2a3db7234330916825f5104f778180179c0e76e27`
- Independently downloaded ZIP digest matched the GitHub digest.
- All 18 files listed in the internal manifest were independently verified by byte count and SHA-256.

## Design actually run

- Analysis parent: 1,200; observed catalogue: 707.
- Independent validation frame: 1,500; SRS validation sample: 200.
- Separate seeds: 2026111901 (analysis), 2026111902 (validation).
- Four detected events; 1,000 generation attempts; 936 detected candidates; 64 rejected candidates.
- Catalogue inclusion DGP: floor-logistic with floor 0.10, midpoint 0.22, width 0.05.
- Inference: degree-one logistic inclusion model with independent Normal coefficient prior, SD (5.0, 2.5).
- H0 prior/grid: uniform [50, 90], 41 points.
- Shared inclusion-coefficient nodes across events; node-specific detection efficiency and selection correction evaluated before marginalization.
- Prior slope-scale sensitivity multipliers: 0.5, 1, 2.

**Important design boundary:** this smoke uses the earlier floor-logistic catalogue-inclusion DGP. It does not yet propagate the quadratic-logit gamma design from the seven-curvature density grid into H0. The separate quadratic-logit H0 work is tracked in PRs #69 and #70.

## Numerical and independent audit results

- Independent auditor: `status=PASS`.
- 18 manifest entries verified.
- 3 H0 posterior normalization checks passed.
- Independent selected-likelihood reconstruction maximum absolute posterior-density differences: exactly 0 at quadrature orders 5, 7, and 9.
- Prior-scale sensitivity posterior-density reconstruction differences: exactly 0 for multipliers 0.5, 1, and 2.
- Quadrature convergence: `CONVERGED` at order 9, after two consecutive passes.
- Order 7 → 9 maximum absolute density difference: (2.42\times10^{-8}).
- Order 7 → 9 maximum absolute quantile shift: (6.82\times10^{-6}) km s\(^{-1}\) Mpc\(^{-1}\).

At order 9:
- Posterior median: 73.535 km s\(^{-1}\) Mpc\(^{-1}\).
- Central 68% interval: [61.703, 83.429].
- Central 95% interval: [53.129, 88.744].
- Generating value: 70 km s\(^{-1}\) Mpc\(^{-1}\).

Prior slope-scale sensitivity, holding the data and H0 grid fixed:
- Multiplier 0.5: median 72.948 (shift -0.587 from reference).
- Multiplier 1.0: median 73.535.
- Multiplier 2.0: median 73.708 (shift +0.173).

## Interpretation and limitations

This run establishes that the existing independent-validation coefficient-posterior → shared selected-H0 likelihood → marginalization → raw-artifact auditor pipeline executes coherently for this small floor-logistic synthetic design. It is not a repeated-sampling calibration result. One posterior containing the generating H0 does not establish coverage, and quadrature convergence does not establish H0-grid convergence or model adequacy.

The initial H0 grid had 41 points over [50, 90]. A same-data 81-point refinement has now passed; see Section 5 below. This establishes grid stability for this one realization, not across the full parameter or DGP space. The quadratic-logit matched-DGP smoke and repeated fixed-truth calibration pilot remain separate gates.

## Next gate

1. Repeat the H0-grid calculation on a finer grid while holding all data and coefficient quadrature settings fixed.
2. Complete the curvature-matched quadratic-logit H0 smoke and independent audit.
3. Finish the paired fixed-truth repeated-sampling pilot; report success, non-evaluable, failed, and not-run counts alongside coverage.
4. Keep conditional H0 PIT distinct from joint SBC; the latter requires generating inclusion coefficients from the declared coefficient prior.


## 5. Same-realization H0-grid refinement (completed)

**Workflow:** 44407318606  
**Artifact:** 14476521866 (`independent-validation-h0-smoke-py3.12`)  
**ZIP SHA-256:** `31dd6c95e58392cac93f1ed72e3d36cc7a11570b7d56e09ccf5a62ee4d7b4877` (matches GitHub artifact digest).

The exact same analysis and validation seeds, parent catalogue, observed GW events, validation sample, coefficient prior, and quadrature orders were run on 41- and 81-point H0 grids. Both independent raw-artifact audits passed; both internal manifests (18 entries each) were independently checked and all entries matched.

At quadrature order 9:
- 41-point posterior median: 73.53528667.
- 81-point posterior median: 73.53898610.
- Absolute median shift: 0.00369943 km s\(^{-1}\) Mpc\(^{-1}\).
- Maximum absolute shift across the seven reported quantiles (2.5%, 5%, 16%, 50%, 84%, 95%, 97.5%): 0.01035126 km s\(^{-1}\) Mpc\(^{-1}\).
- L1 difference between posterior densities after interpolation to the coarse grid: 0.00017062.
- Comparison status: PASS.

The workflow also reran focused coefficient-posterior, shared selected-likelihood, H0 normalization, quadrature-convergence, and prior-sensitivity tests; all passed. Coefficient quadrature converged at order 9 on both grids (orders 5, 7, 9).

This is a numerical grid-resolution result for one floor-logistic synthetic realization. It does not establish calibration, SBC, fixed-truth coverage, or robustness under the quadratic-logit curvature DGP.
