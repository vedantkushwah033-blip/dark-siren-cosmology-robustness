# Repository baseline and scientific implementation audit

**Audit date:** 5 November 2026  
**Inspected branch:** `research/bootstrap-failure-path-repairs-2026-11-02`  
**Inspected branch HEAD:** `4227ed8aaf3007359b3a52acff1e6077e32753f0`  
**Scope:** source/documentation/CI inspection through the GitHub repository API. This is not a local working-tree inspection.

## Executive status

The repository contains several distinct synthetic research tracks: the original low-redshift host-mixture calibration work; toy event-detection selection; catalogue-thinning and host-density estimation; externally estimated catalogue-inclusion probabilities; and a separate finite-population bootstrap methods track. They must not be collapsed into one claim of a single validated end-to-end dark-siren analysis.

The estimated-inclusion branch contains linear and quadratic logistic inclusion fits, paired estimator arms, pre-posterior diagnostics, a curvature-screen runner, and an independent artifact auditor. The five-curvature 20-replicate screen is documented but has not been established as executed by this audit.

## 1. Repository and workflow state

- The inspected branch ref resolves to the HEAD above.
- The default branch is `main`; it is not the inspected research branch.
- Open draft PR #15 remains open and unmerged. Its historical head/base are not treated as the active branch.
- Open draft PR #17 is a separate linear finite-population bootstrap methods track. Its scope and results do not validate nonlinear dark-siren posterior functionals.
- No merge was performed.
- The GitHub Actions run triggered by inspected HEAD, run `41972443506`, was pending when inspected; its job list was empty. Exact-HEAD CI is therefore **not verified**.
- The GitHub API provides remote repository state, not the user's local uncommitted working tree. Local `git status` cannot be asserted from this audit.

## 2. Software layout and dependencies

The branch includes package modules for cosmological inference, event and dataset likelihoods, toy selection, catalogue completeness, KDE host-density estimation, estimated inclusion, validation/bootstrap sensitivity, and independent audits. It also includes experiment scripts, a substantial pytest suite, and extensive protocol/results documentation.

`pyproject.toml` declares Python `>=3.10`, NumPy, SciPy, pandas, and matplotlib; development extras include pytest and Ruff. The inspected CI workflow tests Python 3.10, 3.11, and 3.12. Ruff is declared but the inspected workflow does not show a Ruff lint step. No static type-checking step was found in the inspected workflow.

## 3. Likelihood architecture findings

### 3.1 Missing-host mixture exists in dedicated likelihood code

The mixture is implemented in `src/dark_siren_robustness/dataset_likelihood.py` and `selection.py`, rather than in `selection_toy.py`.

For event (e), the dataset-level event likelihood uses normalized within-component weights and a pre-detection population mass (f):

[
L_e(H_0)=(1-f)sum_j 	ilde w_{c,j},p(d_emid z_{c,j},H_0)
+fsum_k 	ilde w_{m,k},p(d_emid z_{m,k},H_0).
]

The event-detection-efficiency helper uses the same component masses and within-component weights. The selected fixed-count log likelihood applies (sum_elog L_e(H_0)-Nlogeta(H_0)). The dedicated selected-posterior helper likewise combines catalogue and missing components before applying the selection correction.

This is internally coherent for the declared finite-mixture toy model when (f) is interpreted as the pre-detection population mass and the same population mixture, noise model, and detection rule define both numerator and (eta(H_0)). It is not a generic correction for an event-specific catalogue, and the code documentation correctly warns that the supplied finite mixture must represent the population used in the selection integral.

### 3.2 Separate toy likelihoods must not be conflated

- `selection_toy.py` integrates a supplied normalized continuous host-redshift density and conditions on a fixed number of events above a one-sided distance threshold. It does not expose an explicit missing-host fraction; incompleteness must already be represented in the supplied density.
- `selection.py` implements a finite host mixture with a finite noisy-distance interval and includes a dedicated missing-host wrapper.
- `dataset_likelihood.py` provides event-level finite-mixture composition and a separate one-sided detection-efficiency helper. Its own documentation explicitly says it is a minimal prototype without event-rate inference or shared nuisance-parameter marginalization.
- `inference.py` is an unselected low-redshift known-answer baseline, not a selection-aware hierarchical likelihood.

The likelihoods have related mathematics but different data and selection conventions. Cross-module agreement is established only where explicit tests compare them.

### 3.3 Selection normalization and limitations

For conditionally independent detected events under a fixed-size sample and the declared toy assumptions, the correction (-Nlogeta(H_0)) is appropriate when the event numerator is unconditional. It must not be applied a second time to an event likelihood that is already conditioned on detection. The implementation documents this distinction.

The code does not establish a realistic GW detector-selection model, event-rate likelihood, observing-run response, waveform/SNR selection, sky-position effects, or a fully hierarchical model for uncertain population parameters. Those are out of scope, not implicitly solved by the toy (eta(H_0)).

## 4. Estimated-inclusion implementation

The inspected source contains:
- Linear logistic inclusion fitting.
- Quadratic logistic inclusion fitting.
- A synthetic DGP with a quadratic perturbation.
- Paired linear and quadratic estimated-inclusion arms for each validation size.
- An independently sampled parent-frame validation subset with both redshift and catalogue-retention indicator.
- A separate independent auditor that refits the inclusion model and reconstructs the downstream likelihood/posterior from persisted numerical inputs.
- A five-level curvature-screen runner at (g\in\{-2,-1,0,1,2\}).

The DGP uses (x^2-c_{parent}), where (c_{parent}=E_{p(parent)}[x^2]) is computed under the actual non-uniform parent-redshift law. It intentionally does **not** use (1/3), which would center under a uniform distribution on ([-1,1]), not under this parent law. The fitted quadratic model may use the uncentered (x^2) basis because a fixed centering shift is absorbed by its intercept.

The validation design is simple random sampling without replacement from the complete parent frame. This supports estimating (P(R=1mid z)) under the declared design. It does not test transportability to a different survey, population, covariate distribution, or validation mechanism.

## 5. Diagnostic definitions

For fitted inclusion probabilities (hatpi_i) on retained catalogue objects, raw inverse weights are (w_i=1/hatpi_i), with normalized weights (p_i=w_i/\sum_jw_j).

- **Kish ESS:** ((\sum_iw_i)^2/\sum_iw_i^2). It measures weight concentration, not the number of independent galaxies after spatial clustering.
- **Entropy ESS:** (exp[-\sum_i p_i\log p_i]). It is the effective support size implied by Shannon entropy.
- **ESS fraction:** Kish ESS divided by catalogue size.
- **Entropy ESS fraction:** entropy ESS divided by catalogue size.
- **Maximum normalized weight share:** (max_i p_i).
- **Top-(q) weight mass:** sum of the largest (lceil qn\rceil) normalized weights; the implementation uses (q=0.01) and (0.05).
- **Inverse-weight quantiles:** empirical 95th and 99th percentiles, plus the maximum, of (1/hatpi_i).
- **Validation support-gap fraction:** fraction of retained-catalogue redshifts outside the observed validation redshift minimum/maximum.
- **Validation boundary-proximity fraction:** fraction of catalogue redshifts within the validation range and within 5% of its observed span from either endpoint.
- **Coefficient-boundary flag:** whether the bounded logistic optimizer reached a declared coefficient bound.
- **Fit status:** convergence or boundary-limited status from the bounded fit; failures are retained in the pilot ledger.

The support-gap and boundary-proximity features are observable warning proxies, not formal tests of covariate shift or positivity. The current diagnostics are computed before posterior construction and do not accept posterior outcomes or simulation truth.

## 6. Confirmed gaps / release blockers

1. **Exact-head CI pending:** run `41972443506` had not produced jobs/results at audit time.
2. **No completed 20-replicate curvature screen established:** the runner and reduced regression test exist, but the protocol says the full five-condition screen is still to be run and audited.
3. **No complete uncertainty propagation for estimated inclusion probabilities established:** the current validation-resampling sensitivity is conditional on the fixed catalogue and GW events and is not a full (H_0) posterior uncertainty interval.
4. **No config-first experiment contract established for the curvature runner:** several experiment factors and curvature levels are CLI defaults or module constants; a versioned validated manifest with canonical hashing should precede larger runs.
5. **No end-to-end realistic GW likelihood:** the current selection models are controlled toy models.
6. **No real LVK/GWTC or survey validation:** all conclusions remain synthetic.
7. **CI quality gap:** Python matrix tests are present; the inspected workflow does not run Ruff or a type checker.

## 7. Recommended staged release sequence

1. Obtain exact-HEAD CI success and inspect every job.
2. Freeze a validated, versioned curvature-screen configuration and its hash.
3. Execute the five-curvature, 20-replicate screen; independently audit all raw artifacts and report failure denominators and overlap diagnostics.
4. Resolve the uncertainty-propagation design before any expanded simulation: compare a design-respecting repeated-sampling method with a joint inclusion-parameter model, and validate each against known-answer synthetic cases.
5. Only then preregister and run a larger simulation with one primary contrast, multiplicity control, explicit failure accounting, and grid/prior sensitivity.
6. Generate figures only from independently audited completed artifacts.

## 8. Claim boundary

This audit establishes source-level implementation facts and identifies release gaps. It does not establish estimator superiority, confirmatory simulation outcomes, observational validation, or a cosmological measurement. Passing software tests would establish software properties only, not scientific adequacy of the simplified generative model.
