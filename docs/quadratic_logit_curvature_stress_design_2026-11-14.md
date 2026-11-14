# Parent-calibrated quadratic-logit inclusion stress design (2026-11-14)

## Research purpose

The completed floor-family density pilot changes the shape of the inclusion curve through an additive floor. A separate quadratic-logit family is needed to isolate curvature misspecification while holding parent-averaged catalogue completeness fixed.

The experiment asks whether finite-sample estimation of a degree-one logistic inclusion model deteriorates as the true inclusion mechanism acquires quadratic curvature, even when the mean inclusion rate is held constant.

This is a synthetic statistical stress test, not a real survey model.

## Parent measure and coordinate

The parent redshift distribution is proportional to the existing latent intensity:

[
g(z)=1+2.2expleft[-rac{(z-0.105)^2}{2(0.025)^2}ight]
+1.6expleft[-rac{(z-0.275)^2}{2(0.040)^2}ight],
quad zin[0.01,0.40].
]

Use:

[
x=(z-0.205)/0.195.
]

Numerical integration on the same 20,001-point grid used by the implementation gives:

[
c=E_g[x^2]=0.2845005579.
]

The centering is with respect to the parent sampling measure (g(z)), not the host/merger-weighted target density (g(z)b(z)). The host/merger weight must not enter catalogue inclusion calibration.

## Inclusion mechanism

[
pi_gamma(z)=
operatorname{expit}left[
alpha_gamma-3.9x+gamma(x^2-c)
ight].
]

The reference model is gamma=0 with alpha_0=0. Its parent-average inclusion rate is:

[
r_0=E_g[operatorname{expit}(-3.9x)]
=0.5101919962.
]

**Correction to an earlier preliminary note:** 0.5487448435 was an incorrect preliminary rate and is superseded. It must not be used in configs or analysis.

For every gamma, alpha_gamma is deterministically calibrated so:

[
E_g[pi_gamma(z)]=r_0.
]

The calibrated intercepts on the implementation's 20,001-point quadrature grid are:

| gamma | alpha_gamma | achieved parent-average rate |
|---:|---:|---:|
| -2.0 | -0.2404784703 | 0.5101919962 |
| -1.0 | -0.1247216935 | 0.5101919962 |
| -0.5 | -0.0629707145 | 0.5101919962 |
| 0.0 | 0.0000000000 | 0.5101919962 |
| 0.5 | 0.0631226953 | 0.5101919962 |
| 1.0 | 0.1253055057 | 0.5101919962 |
| 2.0 | 0.2424413834 | 0.5101919962 |

Gamma zero is exactly degree-one logistic truth. Nonzero gamma values are structurally misspecified by a degree-one logistic analysis model. Equal parent-average inclusion prevents average completeness from being a confounder of the curvature comparison.

## Implemented DGP component

Added to `src/dark_siren_robustness/generative_selection.py`:

- `quadratic_logit_parent_measure`: computes the parent-measure centering constant.
- `calibrate_quadratic_logit_intercept`: calibrates alpha for a specified curvature.
- `quadratic_logit_inclusion_probability`: evaluates the calibrated probability.
- `generate_quadratic_logit_inclusion_samples`: creates the parent catalogue and independent validation frame.

The independent validation sample is generated with a separate seed. Its observable table contains only validation ID, redshift, and retention outcome. Inclusion probabilities and the full validation frame remain truth-only. The generated configuration records the model, calibration constants, seeds, and sample sizes.

Regression tests check:
- parent-measure centering;
- equal parent-average rates across all seven gamma values;
- the gamma-zero degree-one logistic identity;
- nonzero curvature changing the inclusion curve;
- separation of observable and truth-only fields;
- deterministic reproduction under fixed seeds.

## Density smoke implementation (added 2026-11-15)

Two separate scripts now exercise the five-arm density calculation under the quadratic-logit DGP:

- `scripts/run_quadratic_logit_curvature_density_smoke.py` creates one gamma-specific synthetic replicate, stores observable and truth-only records separately, passes true quadratic inclusion probabilities to Arm C as an oracle diagnostic, and fits the degree-one logistic analysis model for Arms D/E.
- `scripts/audit_quadratic_logit_curvature_density_smoke.py` verifies the manifest, regenerates the population and independent validation frame from recorded seeds, checks the observable/truth artifacts, then independently reconstructs density metrics, weight diagnostics, and Scott widths from raw inputs.
- `tests/test_quadratic_logit_curvature_density_smoke.py` checks the gamma-zero and curved cases, nested validation prefixes, observable-field separation, oracle-weight construction, and output structure.

This is a one-replicate engineering smoke path, not the seven-gamma experiment. Its CI result is pending at the time of this design update; no curvature-family performance conclusion is made from it.

## Paired multi-replicate orchestration (implemented 2026-11-15; CI pending)

Added:

- `scripts/run_quadratic_logit_curvature_density_replicates.py`: enumerates the replicate × gamma design, runs the audited single-cell density runner, preserves every cell attempt in a lifecycle ledger, and aggregates density, weight-diagnostic, and Scott-width records.
- `scripts/audit_quadratic_logit_curvature_density_replicates.py`: verifies the outer manifest and complete planned cell grid, checks paired seed schedules and failure preservation, invokes the independent child raw-input auditor for every successful cell, and reconstructs all aggregate tables from child outputs.
- `tests/test_quadratic_logit_curvature_density_replicates.py`: checks a compact paired 2-replicate × 3-gamma run and rejects aggregate corruption even when the outer manifest is refreshed.
- The focused curvature workflow now includes both single-cell and multi-replicate tests plus a compact Python 3.12 end-to-end run.

The runner uses common random numbers across gamma values within each replicate: all gamma cells share that replicate's population and independent-validation seeds. Replicate seed pairs remain disjoint. This improves precision of within-replicate gamma contrasts, but does not make the cells statistically independent. All gamma-specific DGPs still recalibrate the intercept to the same parent-average inclusion rate.

The new orchestration has not yet completed CI at this exact head. The full seven-gamma × 20-replicate experiment has not been run, and no curvature-family result is reported.

## Proposed full density experiment (not yet run)

The production-scale design should use:

- gamma values: -2, -1, -0.5, 0, 0.5, 1, 2;
- validation sizes: 100, 300, 1,000;
- independent validation frame: 1,500;
- parent population: 4,000;
- five arms A-E, preserving the existing definitions;
- primary Arm E estimand: loss of posterior-mean density;
- secondary Arm E estimand: posterior expected node-specific loss;
- common physical bandwidth grid, initially inherited from the expanded-grid review;
- paired replicate seeds across gamma values and validation sizes;
- explicit lifecycle ledger and preservation of all failures;
- independent child reconstruction and outer aggregate audit.

The analysis should report paired contrasts against catalogue-only, known-probability IPW, and parent oracle, along with ESS and weight-tail diagnostics. Because mean inclusion is fixed across gamma, differences can be interpreted as curvature stress rather than a shift in overall catalogue retention.

Do not use gamma values or validation sizes to tune the model after inspecting holdout results. Freeze the density estimand and bandwidth-training/holdout split before confirmatory runs.

## Scope boundaries

This experiment is density-only. It does not evaluate the GW detection-selection function, H0 likelihood, H0 posterior calibration, real galaxy surveys, or public LVK/GWTC data. The separate GW host-weighting model and GW detection-selection normalization remain distinct from catalogue inclusion.

No curvature-family performance claim is made until the multi-replicate runner, independent auditor, and artifacts are complete.

## Paired multi-gamma pilot runner and outer auditor (2026-11-15)

The GitHub Actions experiment workflow uses the paired-replicate implementation:

- `scripts/run_quadratic_logit_curvature_density_replicates.py` enumerates the replicate × gamma design, runs each single-cell density experiment, preserves every attempt in `cell_ledger.csv`, and reconstructs aggregate metric, diagnostic, Scott-width, and paired-contrast tables.
- `scripts/audit_quadratic_logit_curvature_density_replicates.py` verifies the outer manifest, complete planned cell grid, seed schedule, child raw-input audits, and exact aggregate/contrast reconstruction.
- `tests/test_quadratic_logit_curvature_density_replicates.py` exercises paired runs, aggregate tampering, one-cell failure, and all-cell failure.
- `tests/test_quadratic_logit_curvature_incomplete_reference.py` additionally injects failure of every gamma-zero reference cell and checks that the remaining experiment remains auditable.
- `.github/workflows/quadratic-logit-curvature-density-pilot.yml` runs the configured seven-gamma, three-replicate engineering pilot and uploads the artifact even when an earlier step fails.
- `.github/workflows/quadratic-logit-curvature-density-pilot-tests.yml` runs focused regression tests across Python 3.10, 3.11, and 3.12.

The runner uses common population and validation seeds across gamma values within each replicate. The validation sizes remain nested prefixes of one independent validation frame. Arm C consumes truth-only quadratic probabilities and remains an oracle diagnostic; Arms D/E fit the degree-one logistic analysis model.

Lifecycle semantics are explicit: successful cells are marked `PENDING_OUTER_AUDIT` in the cell ledger until the separate outer auditor independently reconstructs their child outputs; failed cells remain `NOT_RUN` for child audit. If any cell fails, the outer manifest and audit report are `INCOMPLETE`, and the auditor exits non-zero. If any planned replicate × gamma cell fails—including every gamma-zero reference cell—the common-random-number grid is incomplete. The runner therefore suppresses all paired contrast rows and summaries, including within-gamma arm-pair contrasts, rather than reporting a selectively surviving subset. Successful child metrics remain available descriptively, and the complete failure ledger and child directories remain auditable. The experiment must never be presented as a complete PASS in any incomplete case.

The workflow's configured three replicates per gamma are an engineering pilot only. At the latest source review, the updated runner, auditor, and regression-test files were fetched from GitHub, but the exact-head Actions results and uploaded artifact had not yet been retrieved and independently verified. Therefore no seven-gamma performance result is claimed.

The full seven-gamma × 20-replicate experiment has not been run. No curvature-family performance conclusion, leakage-free bandwidth selection, held-out density validation, or H0 propagation is authorized by this pilot.


## Outer-manifest inventory hardening (2026-11-15)

A follow-up audit review identified that the outer auditor verified every file named in the manifest, including byte counts and SHA-256 hashes, but did not reject additional files that were present in the artifact directory and absent from the manifest. Such extra files could not silently change the reconstructed CSV metrics, but their presence meant the manifest was not an exact inventory of the artifact.

The outer auditor now compares the complete set of files under the artifact root (excluding the manifest itself) with the manifest's listed paths. It rejects both missing listed files and unmanifested extra files before proceeding to child audits and numerical reconstruction. A regression test creates a valid paired pilot, adds an extra file without updating the manifest, and requires the auditor to reject the artifact.

Implementation commits on the research branch:
- `f434904391247c9a43a9c75d80dff7d9c7daec8f`: exact artifact-inventory check.
- `fc713b0aa106d47c1e9134510346c5eba52925e9`: unmanifested-file regression test.

The exact branch head was re-fetched after both writes. The GitHub combined-status endpoint returned no status records for that head at inspection time; therefore this change is implemented but its CI result is **not yet verified**. This audit hardening does not change the statistical estimands or provide any new experimental result.


## Nested child-manifest coverage hardening (2026-11-15)

A second inventory review found that the outer manifest builder and auditor excluded every file whose basename was `manifest.json`, rather than excluding only the top-level outer manifest. Each successful child cell has its own `manifest.json`; consequently, those child manifests were not themselves covered by the outer SHA-256 inventory.

The runner now excludes only `<output-root>/manifest.json` when constructing the outer manifest. Nested child manifests are included in the outer file inventory. The outer auditor applies the same exact-path rule, so changes to a child manifest are detected before the child audit is trusted.

A regression test now tampers with a nested child manifest and requires the outer auditor to reject it with an outer-manifest verification failure. The existing extra-file test continues to check exact directory inventory.

Commits:
- `40543c2f2cd76165662e5df4739ca01a2798c75e`: include nested child manifests in the outer manifest.
- `e6f73948364ab6c2290dba9605705412547cad18`: verify nested child manifests in the outer auditor.
- `18889829ddacb05f4fbd5a219fb60071b0d089da`: regression test for nested child-manifest tampering.

The latest branch head was re-fetched after these changes. CI is still pending/cancelled on earlier intermediate heads because subsequent commits triggered new runs; no passing CI claim is made for the latest head yet. This change is artifact-integrity hardening only and does not alter the simulation or its estimands.


## Independent posterior auditor: objective-resolution follow-up (2026-11-15)

The seven-gamma engineering smoke reached the independent raw-input audit for its positive-curvature cells and exposed a numerical issue in the auditor's separately coded Newton polish. The reference implementation accepted a Newton step only when the objective satisfied an Armijo decrease. Near floating-point objective resolution, that decrease can round away even when the independently evaluated posterior score is still strictly improving. The production fitter already permits either Armijo decrease or strict score-norm improvement during its damped Newton refinement.

The independent reference now uses the same mathematical acceptance alternatives while retaining its own BFGS initialization, independently evaluated gradient/Hessian, damped Newton solve, and final score/Hessian checks. This is a convergence-control alignment, not a change to the posterior target, prior, quadrature rule, or analysis model. The final independent score threshold remains 1e-6; no threshold was loosened.

A focused regression test now regenerates the paired-seed gamma=0.5 and gamma=1.0 validation samples at N=4,000 and n_validation=1,000, then exercises the independent posterior reconstruction and checks finite mode, quadrature nodes, normalized weights, and positive output structure.

The paired multi-replicate runner also had a missing NumPy import in its contrast-construction path; that import was added separately. The prior exact-head CI logs confirmed both issues: the runner failed with NameError before contrast generation, while the seven-gamma smoke failed because the independent Newton reference did not converge. Both code paths have now been patched. The newest branch commits have not yet been confirmed by a completed exact-head Actions run, so these fixes remain **pending CI verification**. The earlier engineering artifact remains failed and is not retroactively converted into a passing result.


## Failure-inclusive paired-contrast gate (2026-11-15)

A further review of the multi-replicate runner found that it could calculate paired contrasts from the intersection of successful gamma cells when one planned cell had failed. Although the lifecycle ledger retained that failure and the auditor reconstructed the same intersection, the resulting paired sample could be selectively incomplete. That is not an acceptable basis for the predeclared common-random-number comparison.

The runner now suppresses both per-cell paired contrasts and their summary whenever any planned replicate × gamma cell fails. Successful child metrics, diagnostics, and Scott-width records remain available for descriptive inspection; the complete ledger and all failed child directories remain preserved. The outer auditor independently enforces the same complete-grid rule and rejects a nonempty paired-contrast table for an incomplete design.

The existing injected-failure regression now additionally requires both paired output tables to be empty. The all-failed regression continues to require a fully auditable incomplete artifact with zero paired rows. This change deliberately favors a clear “not estimable under the preregistered paired design” outcome over a misleading partial estimate.


## Inclusion-model specification audit (2026-11-15)

The paired outer configuration now records the fixed degree-one logistic coordinate transform (center 0.205, scale 0.195) and independent Gaussian coefficient prior (mean [0, 0], standard deviations [5, 2.5]). The independent outer auditor checks these fixed specifications and verifies that every successful child artifact agrees on the coordinate transform, coefficient prior, analysis-model label, and estimator-arm labels. A regression test changes the outer prior while refreshing its manifest hash and requires the auditor to reject the artifact. This closes a configuration-drift route in which a validly hashed outer manifest could otherwise describe a different inferential prior from the one used in its child runs.

This remains an engineering-integrity improvement, not a scientific result. The pilot still requires execution, independent artifact audit, a frozen training-only bandwidth-selection protocol, and held-out density validation before confirmatory interpretation or any H0 propagation.


## Explicit inclusion-prior handoff (2026-11-15)

A reproducibility review found that the single-cell runner recorded the fixed Gaussian coefficient prior in its configuration, while relying on the fitter's matching defaults when constructing Arms D/E. The numerical target was the same at this revision, but the implementation had an avoidable implicit dependency: a future change to fitter defaults could make the recorded prior differ from the prior actually used.

The single-cell runner now defines the coordinate transform and Gaussian prior once as named constants, records those constants in each child configuration, and passes the prior mean and standard deviations explicitly to the posterior fitter. The paired runner imports the same constants when writing its outer configuration. A regression test intercepts the five-arm builder and verifies that the declared center, scale, prior mean, and prior standard deviations are passed explicitly.

This is a specification-to-implementation consistency fix. It does not change the current numerical prior (mean [0, 0], standard deviations [5, 2.5]), the posterior target, or any reported experimental result.


## Raw-input common-random-number audit and artifact naming guard (2026-11-15)

The outer auditor now verifies the common-random-number design from persisted raw inputs rather than relying only on repeated seed labels. Within each replicate it compares parent object IDs/redshifts, validation-frame IDs/redshifts, validation-frame selection indicators, and selected validation IDs/redshifts across successful gamma cells. Retention outcomes and catalogue membership are deliberately allowed to differ because gamma changes the inclusion probabilities. The compact paired regression checks these raw-input invariants directly.

The paired runner also rejects gamma values that collide under its two-decimal child-directory naming convention. This prevents two distinct requested curvature values from attempting to write to the same artifact directory.

These changes strengthen reproducibility and artifact integrity only. They do not change the DGP, inferential target, estimator arms, or scientific estimands. Exact-head CI is still required before treating the updated code as validated.


## All-failed and missing-reference lifecycle hardening (2026-11-15)

A static edge-case review found that the all-cells-failed path constructed an empty metric table with object-typed columns and then passed its gamma column to numpy.isclose. The runner and independent outer auditor now explicitly coerce empty-or-populated gamma and node columns to numeric types before contrast reconstruction. This allows an all-failed artifact to remain auditable and correctly report zero paired contrasts rather than failing inside the audit machinery.

The dedicated tests/test_quadratic_logit_curvature_incomplete_reference.py test injects failure of every gamma-zero reference cell while allowing the other gamma cells to complete. It requires the artifact to remain auditable as INCOMPLETE, with no gamma-versus-zero or within-gamma paired summaries. The focused curvature CI workflow now explicitly runs this test and includes it in the workflow path filters.


## Secondary Arm E estimand in paired curvature summaries (2026-11-15)

The paired runner previously retained the posterior expected node-specific loss in the raw child and aggregate metric tables but omitted its gamma-versus-zero paired summary. The primary Arm E estimand remains loss of the posterior-mean density. The secondary estimand is posterior expected node-specific loss; it is now included as its own gamma-versus-zero contrast series, with no within-gamma arm-pair comparison that could blur the distinction between the posterior-mean density and the posterior distribution over node-specific densities. The outer auditor independently reconstructs this additional contrast series. The compact complete-grid regression now expects 576 paired rows and 288 summary rows, including 48 secondary-estimand gamma contrasts. This changes reporting completeness, not the underlying density calculations or estimands.


## Child lifecycle row-count reconciliation (2026-11-15)

The outer auditor now reconciles all three per-cell output counts against the lifecycle ledger: density metric rows, weight-diagnostic rows, and realized Scott-width rows. Previously only the density metric count was checked explicitly, even though the diagnostic and width tables were reconstructed at the aggregate level. A regression test changes the ledger's diagnostic count, refreshes the outer manifest hash, and requires the auditor to reject the inconsistency. This closes a ledger-consistency gap without changing any density result.


## Arm E convex-loss ordering audit (2026-11-15)

The child-level independent density auditor now checks the mathematical ordering between the two distinct Arm E estimands. For each validation size and common physical bandwidth, the primary quantity is loss of the posterior-mean density, while the secondary quantity is posterior expected node-specific loss. Since IAE, ISE, and boundary IAE are convex losses and the posterior quadrature weights are nonnegative and normalized, Jensen's inequality requires the secondary quantity to be at least the primary quantity (up to numerical tolerance). The auditor now rejects violations separately for all three losses, and a regression test injects a violation to verify rejection. This is a consistency invariant, not new empirical evidence or a change to either estimand.
