# Research completion status

**Updated:** 5 November 2026  
**Branch:** `research/bootstrap-failure-path-repairs-2026-11-02`  
**Recorded source head before this status-file update:** `9c16b93dded047a7eb5adf760535c6ab7da3cc90`  
**Default branch:** `main` (not modified)

This is a remote GitHub status record. It does not inspect or certify the user's local working tree.

## Session objective

Resolve the exact CI failure in the configuration-driven estimated-inclusion curvature runner test, document the current evidence boundary, and define the next executable release gate.

## Changes made in this session

1. Updated `tests/test_estimated_inclusion_config_runner.py`.
   - Added a deterministic monkeypatch for `configured.platform.platform()`.
   - Rationale: the production provenance collector calls `platform.platform()`, which may itself invoke subprocess (for example `uname -p`). The test's subprocess mock is intended to validate the runner invocation and Git revision lookup, not to intercept Python's platform introspection.
   - Commit: `97e1b6ee6ec1e195409c3ca5951161d6ac4785d9` — `Isolate platform provenance probe in config-runner test`.

2. Added `docs/research_completion_audit.md` and updated it with the smoke failure and repair.
   - Initial audit commit: `221312794531525d02300702b6b078296c8b708d`.
   - Follow-up audit commit: `9c16b93dded047a7eb5adf760535c6ab7da3cc90`.

3. Fixed the direct-CLI import path in `scripts/run_estimated_inclusion_curvature_screen.py` and added a subprocess `--help` regression test in `tests/test_estimated_inclusion_curvature_screen.py`.
   - Commits: `82732e32e484b17d2233731e51cd7b903f239554` and `2c356da161e701cade8d1f8c0df4a56e1b5388bd`.

No scientific estimator, likelihood, simulation DGP, or scientific estimand was changed. The CI smoke attempted a reduced experiment but failed during script import before simulation. No full experiment was launched. No PR was merged and `main` was not modified.

## Tests and CI

### Exact CI evidence before this session's test fix

Workflow run `41976309029` on prior head `eccaeb8513508cc07ba8ebb8a1b83c4c0e873f47` failed.

Python 3.11 focused steps:
- `tests/test_estimated_inclusion_study.py`: 13 passed.
- `tests/test_audit_estimated_inclusion_pilot.py`: 30 passed.
- `tests/test_estimated_inclusion_curvature_screen.py`: 5 passed.
- `tests/test_estimated_inclusion_config_runner.py`: 8 passed, 1 failed due to the platform subprocess mock interaction described above.

Python 3.10 and 3.12 jobs were cancelled by fail-fast; later steps were skipped. This is not a passing matrix.

Historical curvature regression workflow `41976283260` passed on `fa05e2fc488831166266f0fe1c01a16994344445`; this is not evidence for the new head.

### After the fix

- Test-only fix committed as `97e1b6ee6ec1e195409c3ca5951161d6ac4785d9`; its focused test passed (9 passed on Python 3.12), but the overall workflow failed during configured smoke.
- The smoke exposed a direct-CLI import defect; code/test repairs are in `82732e32e484b17d2233731e51cd7b903f239554` and `2c356da161e701cade8d1f8c0df4a56e1b5388bd`.
- The exact-head run for the final code/doc commits must be refreshed from GitHub Actions before claiming success.
- No local pytest execution was performed in this session.

## Experiment status

| Work item | Status |
|---|---|
| Config schema and canonical hash implementation | Implemented/source-inspected |
| Curvature runner guards | Implemented/source-inspected |
| Config-runner mock repair | Implemented/committed; focused test passed, full workflow then failed in smoke |
| Direct-CLI import repair and regression test | Implemented/committed; exact-head CI pending |
| Reduced configuration-wrapper smoke | Not run |
| Full five-curvature × 20-replicate screen | Not run/verified |
| Full screen artifact audit | Not run |
| Confirmatory pilot and frozen protocol | Not run/not frozen |
| Publication figures from audited screen data | Not available |
| Completed manuscript results | Pending |

The versioned manifest is `configs/estimated_inclusion_curvature_screen_v1.json`. Its declared design is 5 curvature values × 20 replicates, parent population 4,000, 30 detected events, validation sizes 100/300/1,000, H0 grid 41, and redshift grid 201. It remains an exploratory screen.

## Completion-level distinctions

- **Implemented:** the code and configuration wrapper exist.
- **Locally tested:** not performed in this session.
- **CI tested:** earlier commit CI is mixed/failing; exact post-fix head requires refresh.
- **Experiment executed:** no curvature screen executed in this session.
- **Independently audited:** no full-screen artifacts available to audit in this session.
- **Scientifically validated:** not established by this session.
- **Publication ready:** no.

## Known blockers

1. Exact-head Python 3.10/3.11/3.12 CI must pass after the test fix.
2. The configuration-wrapper smoke must be run in a new, empty artifact directory and audited.
3. The full 20-replicate screen must then be run at the exact manifest settings; reduced smoke is not a substitute.
4. The independent auditor must verify all condition outputs, raw-to-summary reconstruction, hashes, weight diagnostics, posterior normalization, and failure denominators.
5. Confirmatory design requires a prespecified primary contrast, uncertainty/power assumptions, and multiplicity policy.
6. Uncertainty propagation for estimated inclusion probabilities through the nonlinear H0 posterior is not established.
7. Finite-population bootstrap validity for the nonlinear posterior functional remains an open methodological question.
8. All results remain synthetic and toy-model-specific; no LVK or survey validation is claimed.

## Next highest-priority task

Refresh GitHub Actions for the current branch head. If the focused test fix passes, execute the reduced config-driven smoke, inspect all emitted artifacts, and independently audit it before starting the full exploratory screen.
