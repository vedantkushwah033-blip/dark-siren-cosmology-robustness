# CI validation record — 2026-10-31

## Scope

This record documents the GitHub Actions validation of the draft linear finite-population bootstrap benchmark changes at commit `b29f80faf934ac9a3f736a12aa4abdc1e58e6574` on the PR #17 head branch.

## Exact-head workflow results

Both exact-head workflow runs completed successfully:

- Pull-request run: [41159462185](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/41159462185) — success.
- Push run: [41159456933](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/41159456933) — success.

The PR run completed the numerical-diagnostics jobs on Python 3.10, 3.11, and 3.12. The Python 3.12 job also completed the v2 linear finite-population bootstrap pilot and independent audit.

## Independently audited pilot artifact

The Python 3.12 workflow log reports:

- Audit status: `PASS`.
- Finite frames independently verified: 3.
- Frame-seed reconstruction: `PASS`.
- Inner-attempt accounting: `PASS`.
- Outer method attempts independently reconstructed: 1,800.
- Summary cells verified: 18.
- Reconstructed estimands: conditional coverage, interval success rate, unconditional successful-and-covering probability, and conditional/unconditional Monte Carlo standard errors.
- Scientific scope explicitly limited to the linear finite-population mean pilot; no transfer to nonlinear dark-siren posterior inference.

The uploaded artifact is `numerical-diagnostics-py3.12` (artifact ID `13026512154`) in the PR workflow run. Artifact availability is subject to GitHub retention.

## Interpretation

This is an exploratory pilot validation, not a confirmatory multi-frame result. Three independently generated frames are insufficient to establish stable between-frame variability or broad method superiority. The successful audit establishes reproducibility of the prespecified pilot design and its recorded outputs; it does not establish general finite-sample guarantees.

The precision planner and its sensitivity analysis remain planning tools. A confirmatory campaign requires a separately justified frame count, frozen scenarios and estimands, frame-aware uncertainty, and a new independently audited run.

## Status

- Exact-head CI: passed on Python 3.10, 3.11, and 3.12.
- v2 pilot artifact audit: passed.
- Multi-frame confirmatory experiment: not yet executed.
- PR #17: remains draft and unmerged.
