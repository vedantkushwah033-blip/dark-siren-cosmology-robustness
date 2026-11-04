# Independent Reconciliation of the Master Computational Research Dossier

**Review date:** 4 November 2026  
**Repository:** `vedantkushwah033-blip/dark-siren-cosmology-robustness`  
**Reviewed branch:** `research/bootstrap-failure-path-repairs-2026-11-02`  
**Reviewed branch HEAD at review:** `335c0dad543016a65e95d042b8e2c98e71faea11`  
**Review status:** source-level reconciliation; not an experimental result or a release approval.

## Executive conclusion

The 2 November dossier is useful as a planning document, but it must not be treated as the current verified baseline or as authorization to launch the proposed confirmatory campaign. Several statements conflict with the live repository or overstate what the code and CI currently establish.

The current repair branch contains a parent-law quadratic-centering correction and a non-confirmatory estimated-inclusion pilot. Its latest GitHub Actions run at the reviewed HEAD was pending with no jobs reported. The confirmatory experiments have not been independently verified from persisted artifacts in this review.

## 1. Live repository facts checked

- Default branch: `main`, at `d461740e9106356f7fdc1c4ec7ed96fb938645a0` at the time of inspection.
- Repair branch: `research/bootstrap-failure-path-repairs-2026-11-02`, at `335c0dad543016a65e95d042b8e2c98e71faea11`.
- PR #15 remains open and draft: `research/estimated-inclusion-audit` → `research/audit-posterior-grid-integrity`; its head is `4f700fcf120ce547de598e6199d3829bc88c8498`.
- PR #16 and #17 also remain draft and unmerged. PR #17 is the separate linear finite-population-mean bootstrap benchmark; it is not validation of the nonlinear dark-siren posterior bootstrap.
- Workflow run #41775011181 targets the repair-branch HEAD above. At inspection it was pending and its jobs endpoint returned an empty list. No exact-HEAD green CI result is established here.
- The current branch has a dedicated selection toy likelihood in `src/dark_siren_robustness/selection_toy.py`, not the `src/dark_siren_robustness/likelihoods.py` path named in the supplied dossier.
- The current estimated-inclusion implementation is in `src/dark_siren_robustness/estimated_inclusion.py`; its fitted model is linear logistic inclusion, `logit pi(z) = intercept + slope * standardized_z`.

## 2. Material corrections to the supplied dossier

| Dossier statement | Source-level reconciliation |
|---|---|
| Quadratic DGP uses `x^2 - 1/3` | The current repair-branch pilot now uses `x^2 - c_parent`, with `c_parent = E_parent[x^2]` computed from the actual latent parent-redshift intensity. The recorded value for the current domain/intensity is approximately 0.284500557. The old 1/3 specification must not be mixed with new artifacts. |
| Quadratic logistic estimator is already available | The inspected inclusion fitter currently fits only an intercept and linear standardized-redshift slope. The quadratic-fitted estimator arm is not yet implemented in this module. |
| A complete two-component catalogue/missing-host mixture likelihood is implemented in the current toy likelihood | The inspected `selection_aware_toy_likelihood` accepts one normalized host density and conditions on a fixed number of detected events. It does not itself implement a separate catalogue/missing-host mixture parameter (f). |
| Validation-bootstrap likelihood integration is an accepted compromise | Not established. A nonparametric/design-aware bootstrap over validation samples approximates a repeated-sampling distribution under its resampling assumptions. It is not automatically a posterior distribution for inclusion-model parameters. Averaging unnormalized likelihood curves over bootstrap fits is not, by itself, Bayesian marginalization. |
| 3,201 H0 grid points eliminate quantile jitter | No fixed grid size proves convergence. A convergence claim requires nested-grid comparisons of posterior summaries and interval endpoints against prespecified tolerances, over representative difficult cases. |
| Prior fixed-truth numerical metrics are verified | The numerical values in the supplied dossier were not independently tied to raw artifacts and an independent reconstruction during this review. Treat them as unverified historical claims until the exact artifact, configuration, seed lineage, audit output, and reconstruction are available. |
| CI PASS across Python 3.10–3.12 | Historical run claims do not establish current exact-HEAD status. The latest inspected repair-branch run was pending with no jobs reported. |

## 3. Statistical reconciliation

### 3.1 Inclusion-model uncertainty

Three targets must remain distinct:

1. **Frequentist repeated-sampling uncertainty:** variability of the fitted inclusion model and downstream estimator across validation samples generated under a specified sampling design.
2. **Bayesian inclusion-parameter uncertainty:** a posterior distribution such as (p(\boldsymbol\beta\mid D_{val})), obtained from an explicit likelihood and prior.
3. **Final H0 posterior uncertainty:** uncertainty conditional on the adopted generative model and all modeled latent quantities.

A validation bootstrap can address (1), subject to a correct resampling design and demonstrated performance for the nonlinear functional. It does not automatically provide (2). The current code explicitly labels its ordinary bootstrap as approximate IID validation-sample bootstrap with the catalogue fixed and no finite-population correction; the design-aware pseudo-population SRSWOR routine is separately labelled exploratory and does not claim validity for the nonlinear logistic→IPW→KDE functional.

**Required decision:** do not implement a class named `BootstrapInclusionPropagator` that averages unnormalized likelihoods until its target estimand and probability interpretation are derived. Preferred candidates for evaluation are (a) a joint Bayesian model for inclusion parameters and H0, or (b) an explicitly frequentist bootstrap distribution of downstream posterior summaries, reported as repeated-sampling uncertainty rather than folded into a Bayesian posterior. A parametric bootstrap may be compared as a third method. The methods must not be conflated.

### 3.2 Selection normalization

For a fixed number of detected events under the inspected toy model,

[
\ell_{sel}(H_0)=\sum_{k=1}^{N_{det}}\log
\left[\int h(z)\,p(d_k\mid z,H_0)\,dz\right]
-N_{det}\log\beta(H_0),
]

where

[
\beta(H_0)=\int h(z)P(\mathrm{det}\mid z,H_0)\,dz.
]

The inspected implementation computes the event numerator from the supplied normalized host density, computes the detection efficiency from the same density and the same distance-noise/detection-threshold model, and subtracts (N_{det}\log\beta). This is internally aligned for the single-host-density toy model, subject to numerical and generative checks.

For a two-component host model, define normalized component densities (h_{cat}(z)) and (h_{miss}(z)), and explicitly define what (f) means. If it is the pre-detection mixture fraction, then the event density and selection efficiency must be derived from the same mixture:
[
q(z)=(1-f)h_{cat}(z)+f h_{miss}(z),
quad
\beta(H_0,f)=\int q(z)P(\mathrm{det}\mid z,H_0)dz.
]
The per-event numerator must likewise integrate against (q(z)). Do not normalize a component twice, apply a catalogue inclusion factor both in (h_{cat}) and again in the mixture, or treat (f) as an observed post-detection fraction without deriving the corresponding transformation. This derivation is a prerequisite to implementing the dossier's mixture formula.

### 3.3 Normalized IPW weights

Normalizing inverse-inclusion weights to sum to one is appropriate for constructing a density estimate, but it does not remove finite-sample instability or uncertainty in fitted inclusion probabilities. The estimand represented by the weighted KDE, its bandwidth selection, boundary behavior, and the role of the merger-rate factor must be stated explicitly. ESS and maximum normalized weight are diagnostics, not guarantees of calibration.

### 3.4 Primary estimand and multiplicity

The supplied dossier alternates between interior coverage and absolute-error reduction as primary outcomes. They are different estimands and cannot both be primary without an explicit hierarchy. The current agreed research protocol selects paired absolute posterior-mean error reduction at the prespecified curvature/validation-size condition as the primary contrast; coverage, interval width, density ISE, and weight diagnostics are secondary. Any change to this hierarchy must be documented before confirmatory runs and accompanied by a multiplicity plan.

## 4. Current experiment scope

The current pilot implementation uses a parent population of 4,000 and 30 detected events by default. Its current arms are catalogue-only, oracle, and validation-estimated linear logistic inclusion at the requested validation sizes. The validation records are sampled without replacement from the complete parent frame and include both retention outcomes. The pilot's own thinning law is separate from the generator's default catalogue thinning; the generated event distances are retained from the parent-population simulation.

The current exploratory nonlinear protocol specifies curvature levels (0,\pm1,\pm2), validation sizes 100, 300, and 1,000, and a staged smoke/200-replicate plan. It does **not** establish the later eight-arm confirmatory Study I or the independent-cohort Study II. Those designs require a distinct, reconciled implementation and preregistration.

The supplied dossier's (N_{det}=300), 2,000-replicate full factorial, and interior-coverage primary outcome are not silently adopted. They conflict with the current pilot scale and the latest agreed staged protocol. In particular, increasing event count by a factor of ten changes the experiment and its runtime/precision properties; it requires an explicit scientific rationale and new planning pilot.

## 5. Release gates before confirmatory execution

1. **DGP gate:** parent-law centering is tested; true inclusion probabilities and the validation sampling design are independently reconstructed; old and new DGP artifacts are separated.
2. **Estimator gate:** implement and test linear and quadratic fitted inclusion models, including boundary/failure states, without truth leakage.
3. **Uncertainty gate:** specify the target of inclusion-uncertainty propagation; compare candidate procedures on controlled data; do not label bootstrap sampling distributions as Bayesian posteriors.
4. **Likelihood gate:** derive the one-component and any two-component selection normalization from the same generative process; test analytic/simple limiting cases and simulation agreement.
5. **Numerical gate:** nested H0/redshift grids with prespecified tolerances; prior-domain sensitivity; independent posterior reconstruction.
6. **Accounting gate:** preserve SUCCESS, NON_EVALUABLE, FAILED, NOT_RUN; reconcile planned counts and wrapper-level failures without fabricated inner-resample counts.
7. **Artifact gate:** immutable output directory, full seed/config provenance, hashes, and an independent auditor that reconstructs outputs rather than merely replaying producer code.
8. **CI gate:** exact branch HEAD must pass all required supported-Python jobs. A run with no jobs or pending status is not a pass.
9. **Planning gate:** use a pilot to estimate Monte Carlo precision and runtime before authorizing large runs. Do not assume (R=2000) or (N_{det}=300) without this evidence.

## 6. Recommended immediate work order

1. Wait for and inspect exact-HEAD CI for `335c0dad543016a65e95d042b8e2c98e71faea11`; fix any failures.
2. Add a quadratic logistic inclusion fitter as a separately named model, with stable optimization, explicit coefficient/boundary diagnostics, and adversarial tests.
3. Extend the independent pilot auditor to reconstruct both linear and quadratic fits and scoring-only DGP probabilities without using truth in estimation.
4. Write and test the selection-likelihood derivation before implementing any missing-host mixture.
5. Only then run a small paired planning pilot; do not launch the 2,000-replicate campaign yet.

## 7. Claims boundary

This review verifies source structure and repository metadata only. It does not verify any scientific performance result, confirmatory simulation, external CHIMERA comparison, LVK data analysis, observational H0 measurement, or publication readiness. The CHIMERA v2.2 benchmark remains a separate external benchmark workstream and is not an arm of the synthetic factorial study.
