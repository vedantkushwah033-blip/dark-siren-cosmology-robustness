# v3 multi-frame pilot: CI-verified smoke and independent audit report

**Report date:** 2026-10-31  
**PR:** #17 (open, draft, unmerged)  
**Source head:** `f869524c3372c1285e6daa212cbfbfa6af4498ea`  
**Workflow run:** [41176574199](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/41176574199)  
**Scope:** software validation of the linear finite-population arithmetic-mean benchmark. This is not confirmatory scientific evidence.

## 1. Verification evidence

The GitHub Actions run completed successfully on the PR merge commit containing source head `f869524c3372c1285e6daa212cbfbfa6af4498ea`.

- Python 3.10 numerical-diagnostics job: success.
- Python 3.11 numerical-diagnostics job: success.
- Python 3.12 numerical-diagnostics job: success.
- The dedicated v3 runner/auditor tests ran in the Python 3.12 job: **7 passed**.
- The v3 smoke runner and independent auditor both executed in the Python 3.12 job and completed successfully.

The Python 3.10 and 3.11 jobs ran the v3 unit tests but skipped the artifact-generating smoke workflow, which is configured to run only on Python 3.12. Do not describe the v3 artifact itself as having been generated independently under all three Python versions.

## 2. Smoke artifact dimensions

The frozen smoke manifest produced:

| Quantity | Value |
|---|---:|
| Distribution regimes | 3 |
| Independent frames per regime | 2 |
| Total finite frames | 6 |
| Population size per frame | 80 |
| Outer repetitions per frame/sample-size cell | 3 |
| Sample sizes | 10 and 30 |
| Interval methods | 3 |
| Outer sample keys | 36 |
| Method-level replicate rows | 108 |
| Frame-summary cells verified | 36 |
| Equal-frame summary cells verified | 18 |
| Inner bootstrap replicates | 25 |
| Nominal interval level | 95% |

The three regimes are Normal, Gamma(shape 0.7), and Lognormal. The manifest seed is 20261031. The artifact status is `PILOT_NOT_CONFIRMATORY`.

## 3. Independent audit outcome

The auditor reported:

- artifact checksums: PASS;
- six finite frames verified;
- 36 outer sample keys verified;
- 108 replicate rows independently reconstructed;
- paired method inputs: PASS;
- seed and sampling replay: PASS;
- inner-attempt accounting: PASS;
- interval endpoints: PASS;
- frame-summary reconstruction: PASS;
- equal-frame summary reconstruction: PASS;
- overall audit status: PASS.

The workflow log records the interpretation as **pilot software-validation evidence only; not confirmatory scientific evidence**.

## 4. What this establishes

This run provides evidence that, for the checked smoke manifest and tested software environment, the v3 workflow can:

1. generate deterministic independent finite frames under the declared seed schedule;
2. preserve common outer samples for paired method comparisons;
3. serialize the expected frame and replicate structure;
4. account for inner bootstrap attempts; and
5. independently reconstruct the sampled indices, interval endpoints, and summaries from the persisted artifact.

This is an important reproducibility milestone. It is not evidence that any interval method attains nominal coverage in a scientifically meaningful population regime. With only three outer repetitions per frame/sample-size cell and two frames per regime, the smoke artifact is far too small for coverage comparison, method ranking, or generalization across finite populations.

## 5. Next scientific gates

Before a confirmatory run:

1. resolve every OPEN item in `confirmatory_design_lock_register_2026-10-31.md`;
2. obtain a defensible estimate of between-frame variability from a planning study with enough independent frames, keeping it separate from confirmatory inference;
3. prespecify a primary contrast, estimand, practical relevance threshold, precision target, multiplicity plan, and stopping rule;
4. select outer repetitions from Monte Carlo precision requirements and independent frame count from frame-level uncertainty;
5. freeze the complete design manifest, source commit, dependency environment, and seed schedule;
6. execute the confirmatory study and independently reproduce its full artifact.

No confirmatory result, preferred method, or publication claim is inferred from this smoke test.

## 6. Scientific boundary

The benchmark concerns interval estimation for the arithmetic mean of a realized finite population under SRSWOR. It does not validate the nonlinear dark-siren posterior, galaxy-catalogue completeness model, host-galaxy weighting, gravitational-wave selection function, or a measurement of the Hubble constant. Those are separate inferential targets and require their own validation evidence.
