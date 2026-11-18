# Execution audit addendum — 2026-11-18

## Purpose

This addendum records a fresh read-only inspection of the live curvature-density implementation and the historical CI failure. It does not report a new simulation run, and it does not alter production code.

## Current implementation inspected

On branch `research/curvature-known-probability-fix-20261114`, I fetched and inspected:

- `.github/workflows/quadratic-logit-curvature-density-pilot.yml`
- `tests/test_absolute_bandwidth_density_replicates.py`
- `scripts/run_absolute_bandwidth_density_replicates.py`
- `scripts/audit_absolute_bandwidth_density_replicates.py`
- `scripts/audit_quadratic_logit_curvature_pilot.py`
- `src/dark_siren_robustness/five_arm_comparison.py`

The fetched file blob SHAs were respectively:

- workflow: `867e4977b99bf09e3ca479f5bd83404469c23717`
- replicate tests: `d1dfa8a1e6aef178e49e208073faf38c16de8cd8`
- replicate runner: `9d3dfbc28dcc6fca147ffa0b26978996b33bff18`
- five-arm comparison: `b495a0370c692b22966c1c7f2996df5fd00147d4`

The current manual workflow is configured for seven gamma values, 20 replicates per gamma, N_parent=4000, validation frame=1500, validation sizes 100/300/1000, 11 absolute bandwidths, and a 401-point redshift grid. It has only a `workflow_dispatch` trigger. The workflow file is configuration, not evidence that the experiment ran. No dispatch operation is exposed by the connected GitHub actions tools inspected in this session; consequently no new workflow was launched.

The current replicate test module is `tests/test_absolute_bandwidth_density_replicates.py`. Its source was fetched successfully and contains valid Python syntax. Its tests cover:
- a two-replicate baseline density pilot plus independent outer audit;
- a two-gamma quadratic-logit smoke plus cross-gamma paired audit;
- rejection of duplicate aggregate keys;
- preservation of failed replicates in the planned denominator;
- rejection of corrupted aggregate Scott widths even when the outer manifest is refreshed.

These are test definitions, not test execution. No new CI run was produced for this branch.

## Historical syntax failure: exact provenance clarified

I fetched job logs for workflow run `43634450957`, job `132568058719` (Python 3.12). The job checked out merge ref `refs/pull/57/merge` at SHA `b45a540b54f49a5fd50a4c9e7597f6fc5d458d24`.

The failing collection step was running a different, older module:
`tests/test_quadratic_logit_curvature_density_replicates.py`.

The log identifies a SyntaxError at line 148 in:
`child_manifest.write_text("{\\\\n  \\\\"tampered\\\\": true\\\\n}\\\\n")`
(the source had over-escaped string delimiters/backslashes). The same collection step failed on Python 3.10 and 3.11. Earlier smoke and single-replicate test groups passed; subsequent steps were skipped.

I fetched that exact old test file from the run's merge SHA and confirmed the malformed source at the reported line. This is a genuine historical failure on the PR #57 merge revision.

**Important correction in interpretation:** the failed old test module is not the current PR #55 replicate test module. The current PR #55 branch uses the newer `test_absolute_bandwidth_density_replicates.py` orchestration and its corresponding runner/auditors. The historical failure therefore must not be described as proof that the current PR #55 source has the same syntax error. Conversely, the current source being syntactically valid does not prove that its tests pass. Exact-head CI remains unverified.

## Code-level observations

1. The current multi-replicate runner labels a child `SUCCESS` when the single-replicate runner returns without raising. It records the child output, then the independent outer auditor re-runs the child raw-input audit. The outer audit is therefore the critical integrity gate; the runner's SUCCESS label alone is not a child audit pass.
2. The runner preserves failed replicate rows in `replicate_ledger.csv` and aggregates metrics only for successful children. The outer audit checks planned counts and reconstruction; incomplete cells must not be promoted to confirmatory results.
3. The cross-gamma auditor explicitly requires all expected gamma artifacts and complete replicate cells, verifies common configuration and calibrated intercepts, and checks paired parent/validation redshifts. This is a valuable safeguard against comparing unmatched Monte Carlo samples.
4. The current workflow is density-only. Its own configuration explicitly excludes H0 likelihood, posterior, calibration, PIT and coverage. No density-only output can support a cosmological calibration claim.
5. In `five_arm_comparison.py`, Arm C can receive the simulation-truth inclusion probabilities, while Arms D/E fit from validation redshifts and binary retention labels. The docstring explicitly warns that the true-probability vector is simulation truth and must not enter D/E fitting. This separation is appropriate and should remain protected by tests.

## Execution status after this inspection

- Current workflow/test/runner/auditor source inspection: **DONE**.
- Historical CI failure reproduced from its exact log and source revision: **CONFIRMED**.
- Current PR #55 exact-head CI: **NOT VERIFIED**.
- Current seven-gamma workflow dispatch: **NOT RUN** (no dispatch tool available in the connected GitHub tool set).
- Seven-gamma density experiment: **NOT RUN / NO NEW RESULTS**.
- Cross-gamma raw artifact audit for the current design: **NOT RUN**.
- H0 inference/calibration: **NOT RUN**.
- Production-code changes: **NONE**.

## Next execution gate

Do not start H0 work. First obtain a passing exact-head test run for the current PR #55 branch (including the new replicate module), then dispatch the seven-gamma density workflow, preserve every raw artifact, run the independent per-cell and cross-gamma auditors, and report failures/non-evaluable attempts with the planned denominator intact. Only after the density endpoint, bandwidth policy, and Monte Carlo precision gate are satisfied should an H0 calibration study be designed.
