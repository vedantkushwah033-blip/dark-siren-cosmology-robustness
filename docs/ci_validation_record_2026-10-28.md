# Exact-commit CI validation record — 2026-10-28

## Scope and provenance

This record captures the GitHub Actions run observed for research-readiness branch commit `f88c7b09793b8232efefa3f71c845a369126e78f`.

- Workflow: `tests`
- Run: [40640604729](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40640604729)
- Run conclusion: **success**
- Run-associated head branch: `audit/research-readiness-2026-10-28`
- Run-associated head SHA, as reported in the uploaded artifact metadata: `f88c7b09793b8232efefa3f71c845a369126e78f`
- Python matrix: 3.10, 3.11, 3.12; all three numerical-diagnostics jobs concluded successfully.
- Artifact set: `numerical-diagnostics-py3.10`, `numerical-diagnostics-py3.11`, `numerical-diagnostics-py3.12`.

## What passed

All three jobs completed the common test and diagnostic steps, including:
- Explicit distance-uncertainty likelihood tests and paired sensitivity smoke experiment.
- Independent reconstruction of distance-uncertainty artifacts.
- Host-weight normalization and joint-selection/missing-host tests.
- Selected-likelihood PIT calibration and independent audit.
- Single- and multi-seed grid-sensitivity checks.
- Host-misspecification and prior-support factorial checks.
- CHIMERA benchmark-data preparation tests.
- Extreme posterior-grid stress audit and fixed-truth reference-interval audit tests.
- Posterior-grid stability and targeted upper-tail convergence diagnostics.
- Artifact upload.

The Python 3.12 job also completed the full 300-replicate distance-uncertainty experiment, independently audited its 32,400-row output, ran the fixed-truth reference-interval comparison, summarized paired contrasts, and recomputed posteriors on a finer grid. Those five full-run steps were correctly skipped in Python 3.10 and 3.11 by workflow conditions.

## Interpretation limits

This is genuine CI evidence for the exact SHA listed above. It is not evidence that the draft PR has been merged, that the simplified model is astrophysically valid, or that the project has established scientific novelty. Passing tests cannot replace independent model validation or a realistic detector-selection comparison.

The full experiment was run in the Python 3.12 job only. The 3.10/3.11 jobs passed the common suite and smoke diagnostics but did not rerun the 300-replicate full experiment.

## Artifact digests reported by GitHub

- Python 3.10 artifact: `sha256:17177dc006f42c81658c473ffd0288d669b43c4fc81869684b9a57d1c6d95cee`
- Python 3.11 artifact: `sha256:6ccf9bd76699f35cf1addf3cd87a41ee00d78386d535bdfdfa9f86a91135606f`
- Python 3.12 artifact: `sha256:ae7417f3ab58e5172e02d502bd35229b83529edae197133d0ec7045b9ba2b049`
