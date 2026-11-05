# Research completion status

**Updated:** 5 November 2026  
**Branch:** `research/bootstrap-failure-path-repairs-2026-11-02`  
**Recorded source head before this status-file update:** `0361a862a9fe2bc774b864a019f09ae93d7fb181`  
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

3. Fixed the direct-CLI import path in `scripts/run_estimated_inclusion_curvature_screen.py` and added a subprocess `--help` regression test.
   - Commits: `82732e32e484b17d2233731e51cd7b903f239554` and `2c356da161e701cade8d1f8c0df4a56e1b5388bd`.
4. Fixed the configuration wrapper's negative-curvature CLI argument by passing it atomically as `--curvature-levels=-2,0,2`; updated the config-runner assertion.
   - Commits: `cb32c98b711d9cab19ef5ee961891ee6db830b67` and `55099af8665c8d1b8b92bbcebf5a1231e33e35f3`.

No scientific estimator, likelihood, simulation DGP, or scientific estimand was changed. The CI smoke was attempted twice on earlier commits: first failed before simulation due to direct-CLI import resolution, then after that repair failed because argparse treated a negative curvature string as an option. Both defects have been repaired in later commits. The smoke has not yet been shown to pass on the current head. No full experiment was launched. No PR was merged and `main` was not modified.

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
- A second smoke exposed negative curvature argument parsing; wrapper/test repairs are in `cb32c98b711d9cab19ef5ee961891ee6db830b67` and `55099af8665c8d1b8b92bbcebf5a1231e33e35f3`.
- Dedicated curvature regression workflow passed on `2c356da161e701cade8d1f8c0df4a56e1b5388bd`; full CI on that head failed in the smoke step before the negative-argument fix. Exact-head CI for the latest fix is required.
- No local pytest execution was performed in this session.

## Experiment status

| Work item | Status |
|---|---|
| Config schema and canonical hash implementation | Implemented/source-inspected |
| Curvature runner guards | Implemented/source-inspected |
| Config-runner mock repair | Implemented/committed; focused test passed, full workflow then failed in smoke |
| Direct-CLI import repair and regression test | Implemented/committed; dedicated regression workflow passed on 2c356da161e701cade8d1f8c0df4a56e1b5388bd |
| Negative-curvature argparse fix | Implemented/committed; exact-head CI pending |
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

## Verified update — 5 November 2026, exact-head CI and configured smoke

**Verified source revision:** `d7dcc79f966fc35a6eec6db60620eff4a3122927`  
**Workflow run:** [41983103045](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/41983103045)  
**Workflow result:** `success`  
**Matrix:** Python 3.10, 3.11, and 3.12 all completed successfully.

The run completed all three matrix jobs (105 steps in each job). The Python 3.12 job also ran the configured curvature-screen smoke and its condition-level independent audits.

### Configured smoke artifact evidence

The Python 3.12 workflow artifact `numerical-diagnostics-py3.12` was downloaded and inspected. It contains:

- Configuration: `configured_inclusion_screen_smoke`.
- Source revision recorded in provenance: `d7dcc79f966fc35a6eec6db60620eff4a3122927`.
- Configuration SHA-256: `d81b6fe0365f215348abce0bbb4618ac748a0e66cf682fb9452b7227703242b2`.
- Curvature levels: `[-2, 0, 2]`.
- Replicates per curvature: 2.
- Parent population: 1,200.
- Detected events: 8.
- Validation sizes: 100 and 300.
- H0 grid: 41 points; redshift grid: 201 points.
- Aggregate records: 36.
- Failed records: 0.
- All three condition-level independent audits: `PASS`.

The persisted aggregate CSV was inspected: it contains the expected 3 × 2 × 6 = 36 condition/replicate/arm records, with all 36 marked successful. Each independent audit JSON exactly matched the corresponding audit embedded in the screen report.

The independent audit's `comparison_diagnostics` also records, for every comparison arm and condition, the number of attempted common replicates, jointly successful pairs, baseline-success/comparison-failure cases, baseline-failure/comparison-success cases, and jointly failed cases. The paired estimator contrast remains explicitly conditional on joint fit success. This differential-failure accounting is present in the independent audit; it should be retained and surfaced clearly in the final full-screen report.

### Evidence boundary and next gate

This is a **reduced smoke test**, not the full exploratory experiment. Its 2 replicates and 3 curvature levels do not establish estimator performance, calibration, or scientific conclusions. It does establish that the configured runner completed on the verified source revision and that its emitted smoke artifacts passed the existing independent condition-level audits.

The full manifest `configs/estimated_inclusion_curvature_screen_v1.json` remains unexecuted in this evidence record. Its five-curvature × 20-replicate screen must run into a fresh output directory and pass independent audits before confirmatory design work is based on it.

No scientific estimator, likelihood, DGP, or estimand was changed by this verification. No PR was merged and `main` was not modified.
