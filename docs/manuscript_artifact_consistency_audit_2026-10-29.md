# Manuscript-to-artifact consistency audit — 29 October 2026

**Scope:** read-only review of the manuscript draft, Stage B report, Study C confirmatory report, independent reproduction note, literature-positioning audit, README and research-readiness report on `research/audit-posterior-grid-integrity`.

**Disposition: CONDITIONAL GO for continued internal review; NO-GO for submission until the checks below are resolved.** This is a document-level consistency review. It is not a fresh execution of simulations or tests.

## Checks that agree

1. The manuscript describes three separate studies and explicitly avoids pooling their estimands.
2. Study A's primary paired contrast is consistently defined as catalogue-only absolute error minus correct-weighting absolute error; positive values favor correct weighting.
3. The Stage B primary result is consistently reported near 5.9519 km s⁻¹ Mpc⁻¹ with bootstrap 95% interval [5.5410, 6.3573], for the declared fixed-H0 synthetic setup.
4. Study C's three confirmatory condition estimates are consistent to rounding across the confirmatory report and independent-reproduction note. Small bootstrap-interval differences are expected because the independent reproduction uses a different bootstrap seed.
5. The manuscript explicitly distinguishes checksum/posterior-summary reconstruction from regenerating Study C's simulation from its random seeds. This limitation must remain in the abstract and methods.
6. The literature-positioning audit correctly states that incomplete-catalogue reconstruction, Gaussian-process methods, clustering-based completion and practical completeness tests predate this project. A first-method or superiority claim is not supported.
7. The external CHIMERA v2.2 benchmark remains blocked; no real-data inference is claimed.

## Corrections and checks required before submission

### 1. Make the Study C audit wording uniform

The manuscript says Study C audits checked artifact integrity and reconstructed posterior summaries and aggregates, but did not independently regenerate the simulations from seeds. The confirmatory report says the author independently recomputed aggregate metrics from downloaded CSVs. These statements can coexist, but the manuscript should distinguish three operations consistently:

- artifact hash/byte-count verification;
- independent reconstruction of posterior grids and raw-table aggregates from saved rows;
- independent regeneration of the original simulated data from the declared seeds.

The available documentation supports the first two for Study C, not the third. Avoid the unqualified phrase “independently reproduced the simulation” for Study C.

### 2. Do not call the three Study C conditions independent replications

The condition seed ranges overlap. The manuscript and reports disclose this correctly. Abstracts, cover letters, figures and conclusions must call them three tested conditions in one related simulation campaign, not three independent replications. The within-condition paired contrasts remain meaningful under their stated design.

### 3. Keep coverage claims bounded by Monte Carlo precision

The Study C coverage estimates use 250 replicates per condition. Report covered counts and Wilson intervals alongside percentages where space permits. Do not imply that 65.6%, 69.6% or 70.8% proves exact nominal calibration. Study A's 40-replicate bandwidth screen is exploratory and is too small for strong coverage comparisons.

### 4. Separate known-law correction from practical estimation

The main weighting arm is supplied the true inclusion probability by construction. The paper should repeatedly make clear that this is a diagnostic benchmark, not evidence that real-survey completeness can be known accurately or that inverse weighting is robust when inclusion probabilities are small or misspecified.

### 5. Recheck every table against immutable artifacts

The documents identify artifact IDs, hashes and workflow runs, but this audit did not download and rerun every artifact. Before external submission, create a release manifest mapping each manuscript table/figure to:
- exact workflow run and source commit;
- artifact identifier and SHA-256;
- script and configuration;
- replicate count, seed range and estimand;
- independent auditor result and its scope.

Do not cite a workflow's green status as evidence of scientific validity.

### 6. Clarify the manuscript's intended contribution

The defensible contribution is a controlled, reproducible, low-redshift synthetic robustness benchmark plus software/reproducibility diagnostics. The current evidence does not establish a new astrophysical result, realistic detector performance, superiority to GP/clustering-based completion, or an observational H0 measurement.

### 7. Pull-request governance

PR #7 remains open and unmerged, has a large review surface, and targets an experimental branch. No submitted reviews or inline review threads were returned in the checked review endpoints. Do not merge automatically. Have a maintainer confirm the intended base branch and review the changed-file set; consider splitting the PR before merge.

## Release gates

- [ ] Exact current code-bearing commit passes supported CI.
- [ ] All manuscript tables and figures are mapped to immutable artifacts and rechecked.
- [ ] Study C language clearly separates artifact/posterior reconstruction from seed-level simulation replay.
- [ ] Overlapping condition seeds are disclosed everywhere conditions are summarized.
- [ ] Coverage uncertainty and exploratory-screen limitations are visible.
- [ ] A domain expert reviews likelihood conditioning, selection normalization, host-population assumptions and interpretation.
- [ ] CHIMERA v2.2 remains marked NO-GO until compatible, version-pinned inputs and an auditable end-to-end execution are available.

**Bottom line:** the documentation supports a credible internal synthetic-methodology draft. It does not yet support calling the manuscript publication-ready or observationally validated.
