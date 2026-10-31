# Independent Inclusion-Audit Cohort: Protocol Draft v1

**Status:** design proposal only; not implemented, preregistered, or executed.  
**Scope:** controlled synthetic methodology; no real survey, LVK, or observational $H_0$ claim.

## 1. Motivation

The existing estimated-inclusion pilot samples audit records from the same finite parent frame that generates the retained catalogue. The audit selection is random, but the audit labels are the same realized retention indicators that define the catalogue. Audit records and catalogue objects can overlap, and the two data sources are dependent.

This protocol defines a separate experiment to test whether the estimated-inclusion correction behaves differently when its inclusion model is trained on an independently generated audit cohort. It does not assume that independence must improve inference.

## 2. Target question

Under a correctly shared redshift population and catalogue-retention mechanism, how does training the inclusion model on an independent audit cohort, rather than an internal subsample of the catalogue's parent frame, change host-density reconstruction and downstream fixed-truth $H_0$ inference?

The experiment isolates sampling dependence. It does **not** test population shift, different survey footprints, different targeting policies, or misspecified transport between surveys.

## 3. Generative design

For each replicate:

1. Generate the **source frame** and GW-event realization using the existing estimated-inclusion pilot's source-population and event-generation procedure.
2. Apply the prespecified logistic catalogue-retention law to the source frame. This creates the source catalogue and its observed GW-event distances.
3. Construct the **shared-frame audit** exactly as in the existing pilot: randomly permute the source frame and take nested prefixes of sizes 100, 300, and 1,000, retaining each sampled object's realized retention label.
4. Independently generate an **audit frame** of redshifts from the same normalized latent galaxy-intensity distribution and the same redshift bounds as the source frame. Draw a fresh Bernoulli retention indicator for each audit-frame object using the *same fixed* logistic retention law. The audit frame must use an independent random stream and contain no source-frame object IDs.
5. Take nested prefixes of the independent audit frame at the same sizes (100, 300, and 1,000). Fit the same linear-logistic and quadratic-logistic inclusion models using only audit redshift and its observed retention indicator.
6. Apply each fitted inclusion model to the unchanged source catalogue to estimate its host density and infer $H_0$ using the same GW distances and likelihood settings as the other arms.

The independent audit cohort is independent of the source catalogue **conditional on the fixed shared data-generating law**. It is not a model of survey-to-survey transport under changing inclusion mechanisms.

## 4. Arms and pairing

Retain the existing source realization and compare:

- Catalogue-only source baseline.
- Oracle source correction using the source population's true inclusion probabilities (simulation benchmark only).
- Shared-frame estimated-linear and estimated-quadratic arms at each audit size.
- Independent-cohort estimated-linear and estimated-quadratic arms at each audit size.

Within a replicate, all arms must share the exact same source parent population, source catalogue-retention realization, detected GW-event distances, redshift grid, $H_0$ grid, and likelihood configuration. The independent audit cohort is generated once per replicate and its nested prefixes are reused across model families and sizes.

Use separate deterministic child seeds for source generation, source catalogue retention, shared-frame audit ordering, independent audit-frame redshifts, and independent audit retention labels. Record all seeds and the audit-cohort provenance. The independent cohort must never reuse source object IDs or source retention labels.

## 5. Estimands

For method $a$ and audit size $m$, retain the existing catalogue-only versus method absolute-posterior-mean-error contrast.

The new primary dependence contrast is:

$$
I_{r,m,a} =
|\widehat H_{0,r,m,a}^{shared}-H_0^{true}|
-
|\widehat H_{0,r,m,a}^{independent}-H_0^{true}|.
$$

Positive values favor the independent-cohort arm **only in absolute posterior-mean error**. They do not establish better interval calibration or universal superiority.

Report this paired contrast separately for linear and quadratic model families, each audit size, and each DGP scenario. Also report paired differences in fixed-truth coverage, interval width, host-density ISE, inclusion-probability MAE, failure rate, boundary-limited fits, ESS, maximum inverse weight, and weight concentration.

Keep attempted, successful, and failed denominators explicit. Successful-pair estimates condition on both fits succeeding; they must not conceal differential failures.

## 6. Data-generating scenarios

Start with the existing two prespecified inclusion laws:

1. Correctly specified linear-logistic DGP (quadratic coefficient 0).
2. Nonlinear quadratic-logistic stress DGP (quadratic coefficient 1.2 in the existing centered basis).

Report scenarios separately. The independent cohort uses the same true inclusion function as the source population in each scenario. Do not pool them.

## 7. Scoring and leakage controls

- Inclusion-model fitting may use only the relevant audit cohort's redshifts and observed retention indicators.
- True inclusion probabilities, source parent truth, GW host truth, and $H_0^{true}$ are scoring-only or oracle-only inputs.
- The independent audit frame's labels must not be generated by conditioning on source catalogue retention.
- The source catalogue and GW-event realization must not change between shared-frame and independent-cohort arms.
- Record source/audit ID disjointness and fail the audit if any IDs or source retention outcomes leak into the independent cohort.
- Audit provenance must independently reconstruct the independent-cohort redshifts, retention labels, and nested-prefix membership from primitive seeds and the declared DGP.
- Report the independent cohort's empirical redshift-distribution discrepancy relative to the source-frame generating distribution. This is a Monte Carlo diagnostic, not permission to reweight after seeing inference outcomes.

## 8. Interpretation limits

This design isolates finite-sample sampling dependence under a common, correctly shared covariate and inclusion law. It does not establish performance under covariate shift, different survey selection functions, measurement error in redshift, clustering dependence, or detector-realistic GW selection.

A difference between shared and independent arms may reflect both dependence structure and ordinary finite-sample variation in the audit redshift/retention draws. Pairing on the source realization reduces unrelated variation but does not make the two audit datasets identical. Interpret uncertainty at the replicate level.

Do not infer that independence necessarily improves estimates. The independent cohort could perform similarly or worse due to its own finite-sample realization.

## 9. Stop/go gates before implementation

1. Review and freeze this protocol before inspecting comparative outcomes.
2. Unit-test independent RNG stream separation and source/audit ID disjointness.
3. Test exact nested-prefix membership and reproducibility.
4. Add adversarial tests proving that source labels/truth cannot enter the independent audit fit.
5. Add independent artifact/provenance reconstruction for both audit designs.
6. Run tiny smoke tests first; label them implementation checks only.
7. Run a full-configuration planning pilot and inspect runtime, failures, positivity, ESS, weight concentration, and paired-difference variability.
8. Set confirmatory replicate count from prespecified Monte Carlo precision targets, then freeze seeds and code commit.

**Execution status (updated):** a standalone independent-cohort generator now exists in `src/dark_siren_robustness/independent_audit.py`, with unit tests in `tests/test_independent_audit.py`. The generator has not yet been integrated into the main pilot runner, and the tests have not yet been executed in a verified CI run. No smoke run, planning run, or scientific result is claimed.
