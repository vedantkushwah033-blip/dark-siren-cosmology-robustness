# Joint catalogue incompleteness × event-detection selection: prospective protocol

**Status: design proposal only; not implemented or executed.** This is not a result.

## Research question

When GW-event detection selection and galaxy-catalogue incompleteness operate together, how do their separate and interacting effects change posterior bias, RMSE, and credible-interval coverage? Under what conditions does explicitly modelling each mechanism improve calibration when the model is correctly specified?

These are different mechanisms. Event selection determines which GW events enter the sample; catalogue incompleteness determines which potential host galaxies are available to inference. Host-population/weight mismatch is a third factor.

The current inverse-redshift catalogue mask is a synthetic stress mechanism, not a calibrated survey selection function. The current hard observed-distance threshold is a toy event-detection model, not a detector-network sensitivity model.

## Stage 0: compatibility gate (required)

1. Pin source revisions and write a joint generative graph: latent source/host population → event properties → GW detection decision → observed GW data → catalogue inclusion/missingness → inference inputs.
2. State whether catalogue inclusion depends on GW detection or on shared latent properties; implement that dependence explicitly if present.
3. Derive the selection normalization for the actual detected-data likelihood. Do not reuse a correction formula from another likelihood without deriving why it applies.
4. In a complete-catalogue, correctly weighted, selection-aware control, compare numerical calculations with independent quadrature/brute-force checks.
5. Demonstrate prior-predictive SBC for a genuinely matched generator and inference model. Keep fixed-truth coverage as a separate estimand.
6. Stop if the simulator and inference likelihood do not describe the same joint model. Do not run the full factorial until this gate passes.

## Staged design

First validate each factor alone, then one combined matched control:
- event selection on/off with a complete catalogue and matched host weights;
- catalogue completeness on/off without event-detection truncation;
- matched versus deliberately reversed host weights on shared observations;
- both mechanisms active, with both corrections correctly specified.

Only after these controls pass, consider the confirmatory factors:

| Factor | Initial levels |
|---|---|
| True H0 | 55, 70, 85 km s⁻¹ Mpc⁻¹ (interior values; boundary stress tests separate) |
| Event detection | no truncation; the specified toy detected-event rule |
| Catalogue completeness | 100%; 70% random missingness; 70% inverse-redshift-weighted synthetic missingness |
| Selection treatment | matched normalization; omitted normalization when detection is active |
| Catalogue treatment | matched model; explicitly defined misspecified/ignored-incompleteness control |
| Host weights | matched; reversed negative control |

This is a design map, not permission to cross incompatible likelihoods. The catalogue correction must be derived before the factor is run.

## Estimands

Report truth/design cells before pooling:
- mean error and median error, labelled separately;
- RMSE and mean absolute error;
- 68% and 95% interval coverage with binomial uncertainty intervals;
- interval widths and prior-boundary accumulation;
- failed/non-finite fits;
- paired replicate-level method contrasts and a prespecified interaction contrast between selection treatment and catalogue treatment.

SBC ranks are prior-predictive diagnostics; fixed-truth coverage is repeated-sampling behavior at a fixed truth. Do not combine them into a single score. State the pooling weights explicitly.

## Audit and reproducibility gates

Freeze protocol, seeds, grid, contrasts, and stopping rules before confirmatory runs. Use disjoint seed blocks across true-H0 strata and common random numbers within meaningful paired contrasts. Save raw replicate-level outputs, metadata, environment versions, and hashes. An independent audit must reconstruct posterior summaries and contrasts without importing primary summary functions; reject missing/duplicate cells, invalid masks, inconsistent detection flags, non-finite values, tampered inputs, and summaries inconsistent with raw rows.

Report Monte Carlo uncertainty for coverage and continuous metrics. Run focused and full CI on Python 3.10, 3.11, and 3.12. Keep large raw artifacts in workflow artifacts rather than Git.

## Stop and interpretation rules

- If the matched selected-data control fails SBC, debug the joint model before interpreting misspecification arms.
- If SBC passes but fixed-truth coverage varies near prior boundaries, report the heterogeneity.
- If the joint likelihood cannot be derived and independently checked, keep the two existing experiment strands separate and state that their interaction remains untested.
- Do not select a correction after seeing results and then describe that result as confirmatory.
- Results remain conditional on the synthetic model; they do not establish that the assumed selection mechanisms match real LVK data.

**Go/no-go:** no full factorial until the joint likelihood derivation and matched-control tests pass. A trustworthy joint model is more valuable than a larger run count.
