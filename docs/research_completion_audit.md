# Research completion audit

**Audit date:** 5 November 2026  
**Branch audited:** `research/bootstrap-failure-path-repairs-2026-11-02`  
**Scope:** remote GitHub source, documentation, workflow and run inspection. This is not a local checkout audit; no claim is made about local uncommitted files.

## Executive finding

The repository already contains a broad synthetic dark-siren software and simulation codebase. The immediate release bottleneck is not lack of additional estimator families: it is completion and verification of the externally validated inclusion-probability curvature experiment, a coherent uncertainty treatment for its nonlinear downstream posterior functional, and exact-head CI/reproducibility evidence.

The configuration-driven curvature runner and independent artifact auditor exist. The prescribed five-curvature, 20-replicate experiment has **not** been established as executed by this audit. A recent full test workflow failed in a test double, not in the scientific simulation: `platform.platform()` internally invoked a subprocess that the test had replaced with a narrowly scoped mock. A targeted test-only repair was committed during this session; its resulting exact-head CI must be checked separately.

## 1. Repository state and evidence boundary

- Default branch: `main`.
- Active working branch: `research/bootstrap-failure-path-repairs-2026-11-02`.
- The branch history includes the estimated-inclusion runner, its configuration wrapper, the independent pilot auditor, pre-posterior diagnostics, and the recent failure-path/configuration test work.
- Main was not modified. No PR was merged.
- Open PR #15 remains a draft and is not treated as a release.
- GitHub's remote API does not expose the user's local working-tree state; no statement about local uncommitted changes is possible.
- The earlier `docs/research_release_baseline_2026-11-05.md` is useful background but its recorded inspected SHA and CI snapshot predate the present audit and must not be reused as the current head record.

## 2. Implemented functionality found by source inspection

### Estimated inclusion and host reconstruction

`src/dark_siren_robustness/estimated_inclusion_study.py` implements a controlled synthetic study with:
- A parent population and a separate catalogue-retention draw.
- A logistic inclusion law with linear and quadratic terms.
- Quadratic centering calculated under the non-uniform parent-redshift intensity, rather than assuming the uniform[-1,1] moment 1/3.
- A simple random validation subset of the finite parent frame, sampled without replacement and independent of catalogue retention.
- Catalogue-only, oracle-diagnostic, and estimated-linear/quadratic arms at each validation size.
- Probability MAE scoring against known synthetic truth, kept separate from estimator fitting.
- Pre-posterior warning diagnostics, weight ESS and maximum inverse weight, density ISE, posterior summaries, and explicit per-arm failure rows.
- Persisted numerical payloads intended to support independent downstream reconstruction.

The oracle arm is a simulation diagnostic because it uses true inclusion probabilities; it is not an implementable estimator and does not imply an upper bound for every target metric.

### Configuration and experiment runner

- `configs/estimated_inclusion_curvature_screen_v1.json` declares 20 replicates, parent size 4,000, 30 detected events, validation sizes 100/300/1,000, seed 20261031, 41 H0 grid points, 201 redshift grid points, and curvature levels -2,-1,0,1,2.
- `scripts/run_estimated_inclusion_curvature_screen_from_config.py` validates a strict JSON schema, canonicalizes and hashes the configuration, invokes the single underlying screen runner, records package/platform/source provenance, and hashes output files.
- `scripts/run_estimated_inclusion_curvature_screen.py` rejects invalid or colliding curvature labels, refuses to overwrite a nonempty output directory, writes condition-level and aggregate raw outputs, and runs the independent pilot auditor per condition.
- `scripts/audit_estimated_inclusion_pilot.py` reads persisted artifacts, independently refits validation logistic models, reconstructs weighted KDE and posterior calculations, and recomputes summary and paired absolute-error contrasts.

### Other research tracks

The source tree also contains separate likelihood/selection, host-mixture, calibration, grid-stability, finite-parent sampling, finite-population bootstrap, and nested validation-bootstrap modules and workflows. These tracks have distinct estimands and assumptions. Their existence must not be presented as a single validated end-to-end astrophysical analysis. In particular, the linear finite-population mean benchmark does not establish validity of a resampling method for the nonlinear posterior functional.

## 3. Scientific assumptions and interpretation limits

1. The inclusion model is `P(R=1 | z)` under a specified synthetic parent-frame law. The validation sample is SRSWOR from that same finite parent frame. This does not establish transportability to a different survey, covariate distribution, or validation design.
2. Inverse-probability weighting requires positivity and adequate estimation of small inclusion probabilities. Large weights can dominate the KDE even when a fitted optimizer reports convergence.
3. Probability MAE is a model-scoring diagnostic; it is not itself evidence of lower H0 error or calibrated credible intervals.
4. Kish/entropy ESS and top-weight mass describe weight concentration, not spatial independence of galaxies.
5. Posterior-mean absolute error, credible-interval coverage, and posterior width are different outcomes. A positive paired reduction in absolute error does not establish calibration.
6. The toy selection-aware likelihood is conditional on its declared generative model. It is not a realistic LVK detector-selection or event-rate likelihood.
7. Synthetic results cannot be described as a real H0 measurement, LVK result, or observational validation.
8. The 20-replicate curvature design is an exploratory pipeline/overlap screen, not a confirmatory study.

## 4. CI evidence inspected

GitHub Actions run `41976309029` (`tests`, branch `research/bootstrap-failure-path-repairs-2026-11-02`, commit `eccaeb8513508cc07ba8ebb8a1b83c4c0e873f47`) completed with failure.

On Python 3.11:
- The estimated-inclusion study tests passed (13).
- The estimated-inclusion artifact auditor tests passed (30).
- The curvature-screen tests passed (5).
- The configuration-runner tests had 8 passes and 1 failure.

The failing test was `test_configured_runner_records_manifest_and_hash`. Its subprocess mock was intended to intercept the runner invocation and `git rev-parse`, but the production provenance collector calls `platform.platform()`; Python's platform implementation may call `subprocess.run(["uname", "-p"], ...)`. The mock incorrectly treated that unrelated platform probe as the runner invocation. This is a test-isolation defect. The repair patches `configured.platform.platform` to a deterministic test value before installing the subprocess mock.

The Python 3.10 and 3.12 jobs were cancelled by the matrix's fail-fast behavior after the Python 3.11 failure; later workflow steps were skipped. Therefore the run is a **failed CI run**, not a passing matrix.

Separate historical evidence:
- Curvature-screen regression workflow `41976283260` passed on commit `fa05e2fc488831166266f0fe1c01a16994344445`.
- Full test workflow `41976283105` failed on that earlier commit.
- These do not certify the current branch head.

## 5. Defects and gaps

### Immediate software/CI defect
- The config-runner test mock was not isolated from `platform.platform()`. A targeted test-only fix has been committed; exact-head CI after that fix remains to be verified.

### Scientific and experiment blockers
- The full prescribed 5 × 20 curvature screen has not been verified as run.
- No raw screen artifact has been independently audited in this session.
- The experiment's SRSWOR validation design estimates inclusion under the same parent-frame mechanism; it is not a transportability experiment.
- The primary contrast and confirmatory estimand need a separately frozen protocol before any large run.
- The proposed pilot sample-size rule requires checking its precision target, pilot uncertainty, paired-contrast distribution, and desired power; it must not be applied mechanically.
- Inclusion-model uncertainty is not yet established as fully propagated into the final H0 posterior uncertainty.
- Paired-record and pseudo-population resampling methods need validation against the complete nonlinear posterior target; finite-population-mean results do not transfer automatically.
- The toy likelihood and synthetic population are not realistic detector/galaxy-survey validation.

### Numerical and reporting limitations
- Grid stability, prior-boundary behavior, and selection normalization must be checked for the exact curvature-screen settings, not inferred from tests of adjacent modules.
- Failed, non-evaluable, and not-run attempts must remain distinguishable in the experiment ledger. Existing pilot rows use a success/failed convention in some paths; the wider four-status taxonomy should be checked before confirmatory use.
- The artifact auditor has substantial independent reconstruction, but an auditor PASS is a software/artifact-integrity result, not scientific validation of the DGP or estimator.

## 6. Verified versus unverified checklist

| Item | Status |
|---|---|
| Source files and configuration inspected remotely | Implemented/source-inspected |
| Curvature CLI collision and overwrite guards | Source-inspected; prior focused regression workflow passed |
| Config schema validation and canonical hash | Source-inspected; focused tests partly passed in failed CI |
| Config-runner test mock isolation repair | Committed; exact-head CI pending at audit time |
| Full supported-version test matrix | Not verified |
| Five-curvature × 20-replicate screen | Not run/verified |
| Screen raw artifacts | Not available in this audit |
| Independent audit of full screen | Not run |
| Confirmatory protocol and frozen contrast | Not complete |
| Full inclusion-uncertainty propagation | Not established |
| Nonlinear posterior bootstrap validity | Not established |
| Real LVK/survey validation | Out of scope/not performed |
| Publication-ready manuscript with completed results | Not ready |

## 7. Next release sequence

1. Verify exact-head CI for the test-isolation repair; inspect all Python-version jobs.
2. Run focused config-runner and curvature tests on the exact head.
3. Run a reduced smoke through the versioned configuration wrapper in a fresh artifact directory; verify manifest, hashes, all conditions/arms, deterministic replay, and independent audit.
4. Only after smoke PASS, execute the prescribed five-curvature, 20-replicate exploratory screen at the declared settings.
5. Independently audit all condition and aggregate artifacts; report every failure denominator and weight/overlap diagnostic.
6. Freeze a separate confirmatory protocol only after examining exploratory feasibility and scientific estimands.
7. Resolve nonlinear inclusion uncertainty and finite-validation uncertainty before making claims about uncertainty calibration.
8. Draft the manuscript around the narrow synthetic benchmark and its limitations; leave numerical results pending until independently audited artifacts exist.

**Audit conclusion:** the repository has substantial implementation, but the central estimated-inclusion experiment remains unexecuted/unverified and the research release is not publication-ready. The immediate next gate is exact-head CI, then a reproducible smoke and independently audited full exploratory screen.
