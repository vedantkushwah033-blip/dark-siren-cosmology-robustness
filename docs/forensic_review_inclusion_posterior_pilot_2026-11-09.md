# Forensic review of the finite-validation inclusion-posterior pilot

**Review date:** 2026-11-09  
**Reviewed reference:** PR #38, head `1b3e3769ecca93f6ebb944b1e631cb622277938b`  
**Scope:** synthetic low-redshift Euclidean dark-siren methodology only.  
**Purpose:** record what the accessible source, workflow metadata, and CI logs establish; define the next evidence gates. This document does not claim a new scientific result.

## Executive finding

PR #38 is a substantial implementation milestone: it adds a proper Gaussian-prior logistic inclusion posterior approximated by Laplace-centred, importance-corrected tensor Gauss-Hermite quadrature; propagates a shared inclusion-coefficient node through the joint multi-event likelihood and the matching GW detection-selection normalization; and provides a 360-cell calibration orchestration pilot with a planned-attempt ledger.

However, the pilot is engineering-scale and is not confirmatory evidence. Its own PR description records 360 planned/attempted cells, 359 SUCCESS and 1 FAILED. The aggregate auditor reports PASS, but its source explicitly says it checks orchestration and arithmetic integrity and trusts saved per-replicate PASS reports. The PR description also acknowledges that the aggregate audit does not independently reconstruct all paired contrasts or all-attempt bounds. The H0 grid has 81 points, and the pilot uses one Bayesian estimated-inclusion arm rather than a matched comparison of the five scientifically relevant procedures.

**Decision:** retain PR #38 as a draft research prototype. Do not merge it as a completed scientific result and do not describe the pilot as confirmatory posterior calibration.

## Evidence inspected

1. Repository metadata and branch/PR inventory on 2026-11-09.
2. PR #38 metadata and full PR description.
3. PR #38 source files:
   - `src/dark_siren_robustness/inclusion_posterior.py`
   - `src/dark_siren_robustness/shared_inclusion_likelihood.py`
   - `src/dark_siren_robustness/h0_posterior.py`
   - `scripts/run_inclusion_posterior_calibration.py`
   - `scripts/audit_inclusion_posterior_calibration.py`
   - `docs/finite_validation_h0_marginalization_design_contract_2026-11-08.md`
4. GitHub Actions run 42505936724: completed successfully on the exact PR #38 head; the PR description records Python 3.10, 3.11, and 3.12 test-matrix success.
5. GitHub Actions run 42505936771: completed successfully for the dedicated calibration pilot. Its job log reports the final planned cell, the aggregate 360 attempts, 359 SUCCESS and 1 FAILED, followed by an aggregate audit status PASS.

No local execution of the repository was performed during this review. No raw artifact ZIP was downloaded or independently recomputed here. The accessible workflow log does not expose the failed cell's underlying exception or cause.

## What is implemented in the reviewed branch

### Inclusion-coefficient uncertainty

The posterior target is the Bernoulli validation likelihood multiplied by a declared proper Gaussian prior. The Laplace approximation supplies an importance proposal; normalized importance-corrected quadrature weights target the declared posterior. This is materially different from treating an inverse-Hessian covariance as if it generated exact posterior draws.

The implementation supports degree-one and degree-two logistic models, checks design rank and Hessian positive-definiteness, and records quadrature diagnostics. Quadrature-order convergence is evaluated separately from downstream H0-grid convergence.

### Shared-coefficient event and selection likelihood

The shared likelihood computes node-specific inverse-inclusion weights once, uses those weights for every event, computes the detection efficiency from the same node-specific population mixture, forms the joint selected log likelihood, and only then marginalizes over the validation-posterior quadrature weights.

This ordering is essential. It avoids independently marginalizing the catalogue-inclusion coefficient for each event and avoids applying the selection correction after prematurely averaging detection efficiencies.

The model is explicitly modular/cut: validation data determine the inclusion-coefficient posterior, and GW data do not feed back into those coefficients. It is not presented as a fully joint posterior.

### Calibration orchestration

The runner creates paired validation-size cells with a fixed analysis catalogue and event dataset across nested validation samples. It preserves the four-state attempt vocabulary SUCCESS, NON_EVALUABLE, FAILED, and NOT_RUN. It reports conditional H0 prior-predictive PIT with fixed inclusion truth and fixed-truth frequentist coverage as separate estimands.

## Findings and unresolved risks

### F1 — One failed cell remains unresolved

The pilot log identifies:

`fixed_truth_floor_0.25_n_1000_rep_0012: FAILED (RUN_SMOKE)`

The overall pilot still accounts for all 360 planned attempts: 359 SUCCESS and 1 FAILED. The accessible log does not reveal the underlying exception. The failed cell must remain in the denominator and must not be silently rerun or replaced. Obtain its raw replicate directory and failure record from the workflow artifact, inspect the exception, and classify the cause before deciding whether a code repair is justified.

### F2 — Aggregate PASS is not an independent raw-data audit

The aggregate auditor checks ledger/summary arithmetic, per-replicate audit attestations, selected posterior summary fields, PIT/coverage flags, and nested validation structure. It does not re-run the per-replicate raw likelihood reconstruction itself. Its own module documentation explicitly limits the claim.

The next auditor should independently reconstruct each successful replicate from its immutable observable inputs and saved model/configuration, rather than trusting only the stored `audit_report.json` PASS field. It should independently recompute:
- coefficient posterior/quadrature outputs or a separately coded reference on selected adversarial cases;
- node-specific catalogue weights and their normalization;
- event-mixture log likelihoods;
- node-specific GW detection efficiencies;
- shared-node selected log likelihoods and log-sum-exp marginalization;
- H0 posterior normalization and quantiles from raw grid-level likelihood/density arrays;
- PIT and fixed-truth interval indicators;
- aggregate summaries and paired contrasts.

The independent implementation should avoid importing the production function being audited for the critical calculation.

### F3 — The pilot does not establish H0-grid convergence

The pilot uses an 81-point H0 grid. Coefficient quadrature convergence is not a substitute for H0-grid convergence. Before interpreting posterior coverage or PIT, repeat selected identical datasets on a predeclared sequence such as 81, 161, 321, and a finer reference grid if required. Compare posterior quantiles, density shape, normalization, and boundary mass. Report non-convergence rather than silently accepting a coarse grid.

### F4 — One inference arm cannot answer the comparative research question

The pilot uses one Bayesian estimated-inclusion arm. The central research question concerns when correction helps or harms. A matched comparison should preserve the same generated parent frame, catalogue, GW events, and validation sample across:
1. catalogue-only unweighted inference;
2. oracle-inclusion IPW (simulation-only diagnostic);
3. estimated-inclusion plug-in IPW;
4. estimated-inclusion modular posterior marginalization;
5. a clearly specified parent/oracle population reference where appropriate.

Each arm must use the same event and detection-selection convention. The oracle arm must never be described as an operational estimator. The inferential target and comparison estimand must be stated before running the experiment.

### F5 — The current PIT is conditional, not joint SBC

The pilot's H0 prior-predictive PIT randomizes H0 while holding inclusion-truth coefficients fixed within a design cell. It is a conditional diagnostic under that fixed inclusion data-generating process; it is not joint simulation-based calibration for the inclusion coefficients and H0. Do not label it ordinary joint SBC. A genuine joint SBC experiment would need to draw all relevant parameters from the declared joint prior and generate the validation, catalogue, and GW data under the same generative model, while preserving the distinction between modular and fully joint inference.

### F6 — The modular/cut target must remain explicit

The validation posterior for inclusion coefficients is deliberately insulated from GW feedback. This is a defensible modular target for a controlled benchmark, but it is not interchangeable with a full joint posterior. Manuscript language, artifact metadata, and comparisons must preserve that distinction.

### F7 — Model support and failure behaviour need adversarial review

The quadrature implementation should be stress-tested for:
- near-separation and strongly prior-dominated validation samples;
- validation designs with poor covariate spread;
- quadrature importance-weight concentration and tail-node influence;
- probabilities close to zero and extreme inverse weights;
- normalized host-weight underflow;
- coefficient-prior scale sensitivity;
- H0 likelihood concentration near prior boundaries;
- non-finite or effectively zero detection efficiencies.

A successful numerical return is not automatically a reliable approximation. The convergence ladder and diagnostics must be interpreted together.

## Required next gates

| Gate | Required evidence | Release condition |
|---|---|---|
| G1: failure forensics | Raw artifact and exception for the single failed cell | Cause documented; failure retained in original ledger |
| G2: raw independent audit | Recompute replicate likelihood/posterior outputs from raw observable inputs using independent code | Zero unexplained discrepancies; adversarial tests pass |
| G3: H0-grid convergence | Matched 81/161/321 (or finer) calculations | Predeclared quantile/density tolerances met |
| G4: matched estimator arms | Same datasets/validation samples across five arms | Complete paired attempt ledger and auditable raw outputs |
| G5: uncertainty sensitivity | Prior-scale and quadrature-order sensitivity, including difficult cells | Stable or explicitly non-evaluable; no hidden fallback |
| G6: powered confirmatory design | Replication selected from target Monte Carlo precision, with all-attempt bounds | Protocol frozen before confirmatory execution |
| G7: final synthesis | Independent audit, negative findings, uncertainty intervals, limitations | Claims limited to the synthetic estimand and tested conditions |

## Interpretation boundary

Even after these gates pass, the result will remain a synthetic low-redshift methodological benchmark. It will not establish performance on LVK/GWTC events, validate a real galaxy survey, model a realistic detector pipeline, or produce an observational H0 measurement.

The most defensible possible contribution is an independently audited comparison of how finite validation, logistic inclusion-model misspecification, and IPW weight instability affect host-mixture reconstruction and selection-normalized H0 posterior behaviour under explicitly stated synthetic assumptions.

## Immediate next action

Retrieve and inspect the raw artifact for run 42505936771, resolve the failed replicate's cause without deleting or replacing its attempt, and upgrade the independent auditor so it reconstructs the likelihood and posterior from raw inputs. Only after those steps should the matched five-arm H0-grid-converged pilot be expanded.
