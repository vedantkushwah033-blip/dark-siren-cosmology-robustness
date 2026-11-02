# Follow-up audit: nested validation-bootstrap replay independence and failure completeness

**Reviewed source reference:** PR #17 head `7c612efb1b075ea2df0e665b93091fa78c77fe84`  
**Review branch:** `research/hostile-referee-audit-response-2026-11-02`  
**Scope:** static source audit of the nested runner and replay auditor. No simulation was executed and no scientific result was generated.

## Executive finding

The independent replay architecture is materially stronger than a simple rerun of the production bootstrap. For successful inner attempts, the replay auditor reconstructs resampling indices from seeds and independently refits the logistic inclusion model, rebuilds the reflected weighted KDE, recalculates the selection-conditioned likelihood, and reconstructs posterior summaries. It does not call the production estimator or likelihood for those successful attempts.

However, the audit is not completely independent for wrapper-level failures: `scripts/audit_nested_validation_bootstrap_replay.py` imports and reruns `bootstrap_validation_h0_posterior` when a whole bootstrap wrapper is recorded as failed. That is a deterministic failure-reproduction check, not an independent reconstruction. This limited circular dependency does not invalidate successful-attempt replay, but the report must not describe wrapper-failure verification as independent.

## 1. What is genuinely independent

### 1.1 Independent scenario reconstruction

`scripts/audit_nested_validation_bootstrap_replay.py::_reconstruct_generator` independently reconstructs:
- the parent redshift frame by inverse-CDF sampling;
- the generator's catalogue-thinning RNG consumption;
- host probabilities and GW distance rejection sampling.

It then verifies parent, catalogue, and GW-distance fingerprints. This is a substantive independent replay, subject to matching the generator's exact algorithm and numerical-library behavior.

### 1.2 Independent outer validation samples and posterior

The replay auditor regenerates the outer random permutation from the scenario seed and checks each saved validation-index prefix. It reconstructs the finite-parent target and each successful outer posterior using `posterior_validation_bootstrap_audit._independent_posterior`.

### 1.3 Independent successful inner bootstrap attempts

`src/dark_siren_robustness/posterior_validation_bootstrap_audit.py::_independent_posterior` separately implements:
- bounded logistic maximum-likelihood fitting using `scipy.optimize.minimize`;
- inverse-inclusion weighting and Kish ESS;
- reflected weighted KDE;
- selection-conditioned likelihood;
- normalized posterior and posterior summaries.

The audit module explicitly avoids importing the production estimator or likelihood functions. It independently regenerates paired-bootstrap or pseudo-population/SRSWOR indices and checks every successful attempt's posterior density, posterior mean, quantiles, diagnostics, and reference mean.

This supports the description **independent numerical replay of successful attempts**, not independent validation of the scientific appropriateness of the estimator or resampling method.

## 2. Limited circularity in wrapper-level failures

At the top of `scripts/audit_nested_validation_bootstrap_replay.py`, the replay script imports both:
- the independent audit functions from `posterior_validation_bootstrap_audit.py`; and
- the production `bootstrap_validation_h0_posterior` function.

The production function is used in the branch for rows whose entire bootstrap wrapper is marked failed. The auditor reruns the same production wrapper with the same inputs and seed and compares exception class/message.

**Interpretation:** this verifies deterministic reproducibility of the wrapper failure, but it is not independent verification that the failure was scientifically or computationally appropriate. A bug shared by producer and replay call can reproduce itself.

**Required correction before a confirmatory artifact audit:** independently reconstruct wrapper-level failure conditions from the saved inputs and attempt trace, or explicitly downgrade the claim to “deterministic wrapper-failure replay.” Do not claim fully independent wrapper-failure verification.

## 3. Failure-accounting gap in the runner

In `scripts/run_nested_validation_bootstrap_comparison.py::run_comparison`, several per-draw failures are preserved, but a scenario-level failure while constructing the finite-parent target is raised as a `RuntimeError` before rows and attempt records for that scenario are written. The scenario gates include insufficient catalogue size/unique redshifts and a non-finite full-frame target.

This means the runner does not always produce a complete artifact with an explicit failure record for every planned scenario/draw. The current protocol requires failures and denominators to be retained. A failed scenario may instead abort the run and leave partial outputs absent or incomplete.

**Required correction:** add scenario-level status records and ensure the manifest records planned versus completed scenarios, draws, and method attempts. Do not silently treat an aborted run as a completed benchmark. The independent auditor should reject incomplete manifests unless an explicitly supported partial-run schema is used.

A second edge case exists inside `posterior_validation_bootstrap.py`: the paired-bootstrap branch appends an empty index placeholder if an exception occurs before index persistence, but the design-aware branch does not have the analogous placeholder logic. The normal validated path should generate indices before model fitting, so this is an exceptional-path integrity gap rather than evidence of a routine failure. Add adversarial tests that force an exception at the index-generation boundary and verify attempt-aligned arrays remain complete for both methods.

## 4. Important statistical interpretation checks

1. The outer finite-parent target is the posterior mean obtained by fitting the declared inclusion pipeline to all actual parent-frame records. It is not the injected (H_0).
2. The outer SRSWOR estimates measure validation-sample variability conditional on the fixed scenario (parent frame, catalogue, and GW distances).
3. The paired nonparametric bootstrap and multinomial pseudo-population/SRSWOR comparator are not proven valid uncertainty estimators for this nonlinear posterior functional.
4. The FPC-normal interval is a diagnostic comparator only. Its use of (sqrt{1-n/N}) does not establish correct variance scaling for the fitted logistic + weighted KDE + selection-normalized posterior-mean functional.
5. The independent replay validates arithmetic and artifact reconstruction under the declared algorithm. It does not establish frequentist coverage, model adequacy, real-survey validity, or real detector realism.

## 5. Prioritized next changes

### P0
- Remove or explicitly scope the production-function dependency for wrapper-level failure replay.
- Preserve scenario-level failures and planned/completed denominators in partial-run artifacts.
- Add adversarial tests for design-aware index-generation failure and attempt alignment.
- Ensure the replay auditor rejects truncated/incomplete artifact sets.

### P1
- Independently validate the generator reconstruction against a separately implemented statistical distribution check, not only exact seeded replay.
- Add scenario-level summaries of failure rates and verify denominators against the manifest.
- Add tests that deliberately tamper with summary CSVs, manifest counts, missing rows, duplicated keys, and scenario omissions.

### P2
- Keep all method labels explicitly exploratory until repeated-sampling calibration supports stronger claims.
- Continue to report empirical outer SRSWOR SD beside inner SE and target-coverage summaries.
- Keep scenario as the top-level independent unit in uncertainty intervals.

## 6. Execution status

This was a static source inspection at PR #17 head `7c612efb1b075ea2df0e665b93091fa78c77fe84`. No local Python environment was available through the GitHub source-inspection workflow. No tests were run, no code was modified, and no confirmatory simulation was launched.
