# Paired nonlinear inclusion-curvature screen runner

**Status:** implementation added; the five-condition screen has not yet been executed at the 20-replicate protocol setting.  
**Scope:** synthetic host-catalogue retention only. This is not LVK/real-survey validation or an observational $H_0$ result.

## What the runner does

`scripts/run_estimated_inclusion_curvature_screen.py` executes the prespecified curvature levels
$g\in\{-2,-1,0,+1,+2\}$ from
`docs/estimated_inclusion_nonlinear_stress_protocol_2026-10-30.md`.

For each condition it:

1. Calls the existing estimated-inclusion pilot without changing the fitted linear-logistic estimator.
2. Uses the same root seed across curvature conditions. This creates common random numbers for parent populations, GW events, catalogue-uniform draws, and validation-sample ordering; the retention probabilities themselves change with $g$.
3. Preserves the catalogue-only, oracle diagnostic, and all requested validation-size arms.
4. Saves the raw replicate table and the existing summary separately for each curvature.
5. Calls the independent persisted-artifact auditor for each condition and stops immediately if any audit fails.
6. Produces a combined raw table and a screen report with failure counts, fit-boundary counts, true-probability ranges, low-probability mass, and maximum true inverse probability.

No probabilities are clipped and no inverse weights are capped. Oracle weighting remains a diagnostic upper-bound arm, not an implementable method.

## Reproducible command

After installing the project and test dependencies, run from the repository root:

```bash
python scripts/run_estimated_inclusion_curvature_screen.py \
  --replicates 20 \
  --n-parent 4000 \
  --n-detected 30 \
  --validation-sizes 100,300,1000 \
  --seed 20261031 \
  --h0-points 41 \
  --z-points 201 \
  --output-dir artifacts/estimated_inclusion_curvature_screen
```

The 20-replicate run is a **screening/smoke gate**, not a basis for method ranking or confirmatory inference. Review every condition's overlap diagnostics and failure accounting before considering the 200-replicate stage. A cell with severe positivity problems must remain a separately labelled stress condition; do not pool it with moderate-overlap cells.

## Outputs

Each condition directory contains:

- `replicates.csv`
- `summary.json`
- `independent_audit.json`

The root output directory contains:

- `all_conditions_replicates.csv`
- `screen_report.json`

The test suite exercises a reduced one-replicate version of all five conditions and requires all five independent audits to pass. That regression test validates orchestration and artifact handling only; it does not substitute for the 20-replicate screen.

## Interpretation and stop/go

Do not interpret the output until all five independent audits pass. Then review:

- true inclusion-probability minimum/maximum and fractions below 0.05 and 0.10;
- true inverse-probability tails;
- retained catalogue sizes;
- validation-model fit failures and coefficient-boundary hits;
- arm-specific failure rates and paired-success denominators;
- estimated-arm probability error, ESS, and weight tails;
- paired absolute-error contrasts with uncertainty intervals.

The pilot summary's normal intervals are not a confirmatory uncertainty procedure. The five curvature levels and three validation sizes form a factorial family; report the complete pattern and label cell-wise contrasts exploratory. Do not proceed to 200 replicates until the 20-replicate screen is audited and the overlap/failure review is documented.

**Claim boundary:** this screen tests a deliberately simplified one-dimensional synthetic retention model with an independent parent-frame validation sample. It does not establish performance under covariate shift, validation nonresponse, redshift error, spatial clustering, unobserved selection covariates, realistic survey masks, or uncertain GW detection efficiency.
