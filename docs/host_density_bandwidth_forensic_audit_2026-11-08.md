# Forensic audit: KDE bandwidth semantics and H0 baseline sensitivity

**Date:** 2026-11-08  
**Status:** local diagnostic; requires repository-level implementation and CI follow-up  
**Branch:** `research/independent-posterior-audit-claim-fix`

## Executive finding

The existing host-density benchmark passes `bandwidth_factor=0.8` to `scipy.stats.gaussian_kde`. SciPy interprets a scalar `bw_method` as the **absolute covariance factor**, not as a multiplier of Scott's rule. Thus 0.8 means a kernel standard deviation of approximately (0.8\,s_z), not (0.8\times n_{\rm eff}^{-1/5}s_z).

For (n=4{,}000), Scott's factor is approximately 0.19. The existing 0.8 setting is therefore about four times Scott's factor. A controlled local diagnostic found substantial density-recovery and H0-posterior sensitivity to this setting. This is a major confounder for the existing oracle/known-completeness H0 baseline and must be resolved before interpreting any estimated-inclusion H0 contrast.

No repository source code was changed by this diagnostic. The design document was updated to remove 0.8 as a recommended default.

## Code-path evidence

The current `estimate_host_density` implementation constructs the reflected sample and then calls:

```python
gaussian_kde(augmented_sample, weights=augmented_weights, bw_method=bandwidth_factor)
```

When `bandwidth_factor` is a scalar, SciPy uses that scalar directly as its covariance factor. The no-argument path separately computes the weighted original-sample effective size and uses its Scott factor. These are different semantics.

Consequently:
- `bandwidth_factor=None`: use the computed weighted Scott factor.
- `bandwidth_factor=0.8`: use the absolute SciPy factor 0.8.
- It does **not** currently mean “80% of Scott bandwidth.”

The existing host-density benchmark passes 0.8 explicitly. The current unit tests check normalization and a relative mean-shift property, but do not assert the intended scalar-bandwidth semantics or assess recovery of the declared parent density.

## Local H0 diagnostic design

A targeted local experiment used the existing estimator, validation-split, and selection-conditioned likelihood code paths copied from the current branch, with a small standalone driver:

- 20 paired simulation replicates.
- 4,000 parent galaxies per replicate.
- Parent redshifts drawn from the repository's declared two-overdensity intensity on (0.01\le z\le0.40).
- Merger weight (m(z)=1+1.5(z-0.01)/0.39).
- 15 detected events per replicate.
- (H_{0,\rm true}=70\;\mathrm{km\,s^{-1}\,Mpc^{-1}}), Gaussian distance noise (sigma=140\) Mpc, threshold 250 Mpc.
- H0 grid: 201 points on [55, 85]; redshift grid: 201 points on [0.01, 0.40].
- Identical parent/event realizations across the bandwidth comparisons.
- The analytic reference density was the normalized declared intensity times the declared merger weight. The KDE arms used the same parent samples and fixed bandwidth choices.

This is a diagnostic pilot, not a full repository benchmark, not an independent artifact audit, and not a confirmatory coverage study.

## H0 posterior sensitivity

Across the 20 paired replicates:

| Host-density construction | Mean posterior-median bias | Empirical 68% coverage | Mean 68% width |
|---|---:|---:|---:|
| Analytic parent density | +0.691 | 14/20 = 0.70 | 9.660 |
| Parent KDE, factor 0.03 | +0.753 | 14/20 = 0.70 | 9.626 |
| Parent KDE, factor 0.05 | +0.827 | 14/20 = 0.70 | 9.663 |
| Parent KDE, factor 0.08 | +1.004 | 14/20 = 0.70 | 9.722 |
| Parent KDE, factor 0.10 | +1.161 | 14/20 = 0.70 | 9.750 |
| Parent KDE, factor 0.15 | +1.618 | 12/20 = 0.60 | 9.764 |
| Parent KDE, factor 0.20 | +2.019 | 12/20 = 0.60 | 9.721 |
| Parent KDE, Scott default | +1.951 | 12/20 = 0.60 | 9.732 |
| Parent KDE, factor 0.80 | +2.775 | 10/20 = 0.50 | 9.490 |

The 20-replicate uncertainty is substantial; these are not precise coverage estimates. The important finding is the large systematic movement across bandwidth choices, including between the analytic density and the current factor-0.8 KDE.

In the same local plug-in smoke, the known-inclusion IPW and estimated-inclusion IPW arms at factor 0.8 were close to the parent-KDE oracle arm. That does **not** demonstrate that inclusion correction is ineffective. It shows that the current smoothing choice can dominate or obscure the inclusion-model contrast.

## Parent-density recovery

For each of the same 20 parent frames, the KDE host density was compared with the normalized analytic parent intensity times merger weight using integrated absolute error (IAE):

| SciPy covariance factor | Mean IAE |
|---:|---:|
| 0.03 | 0.04044 |
| 0.05 | 0.05349 |
| 0.08 | 0.10461 |
| 0.10 | 0.14311 |
| 0.15 | 0.22528 |
| 0.20 | 0.27478 |
| Scott default | 0.26769 |
| 0.80 | 0.31305 |

These results explain the H0 sensitivity: broad kernels erase the narrow overdensity structure in the declared parent law. Factor 0.03 performed best in this particular finite sample, but it must **not** be selected as a universal default from these same 20 seeds. It may under-smooth in other parent laws or sample sizes.

## Implications for the project

1. **The existing H0 grid-convergence result is not a density-estimation validation.** It establishes numerical integration stability conditional on the supplied host density. It does not establish that the KDE represents the parent host population.
2. **The current factor-0.8 oracle arm is not an adequate calibration reference until bandwidth sensitivity is resolved.** Oracle input redshifts do not make a heavily smoothed KDE equal to the parent density.
3. **Do not attribute the observed H0 bias to finite-validation inclusion estimation.** The parent-KDE oracle itself exhibits the bandwidth-associated shift.
4. **Do not choose a bandwidth by maximizing H0 coverage on the same seeds.** Use a predeclared density-recovery criterion on training simulation seeds and evaluate the chosen rule on held-out seeds.
5. Compare at least: analytic parent-density reference, parent-sample KDE, catalogue-only KDE, known-inclusion IPW KDE, and estimated-inclusion IPW KDE under matched bandwidth rules.
6. Report bandwidth semantics explicitly in every manifest and figure. Distinguish absolute SciPy covariance factors from multipliers of Scott's factor.

## Required remediation before estimated-inclusion H0 claims

- Clarify the `bandwidth_factor` argument in the code docstring and tests.
- Add an explicitly named Scott-multiplier option (or retain the current scalar API but rename it clearly); avoid silently changing the old API's meaning.
- Add density-recovery tests against analytic truth on training seeds, then verify on held-out seeds.
- Run a bandwidth sensitivity grid for catalogue-only, oracle, known-inclusion, and estimated-inclusion arms.
- Freeze the bandwidth-selection rule before the confirmatory H0 campaign.
- Re-run the H0 baseline and independent artifact audit under the corrected bandwidth protocol.
- Only after the oracle/known-inclusion baseline is adequately calibrated should finite-validation coefficient uncertainty be propagated into H0 inference.

## Decision

**NO-GO for interpreting the current factor-0.8 H0 benchmark as a calibrated reference for the finite-validation project.**  
**GO for a dedicated KDE bandwidth forensic and density-recovery phase.**

The probability-stage finite-validation pilot remains independently audited. This bandwidth issue is a separate downstream host-density problem and does not invalidate the raw-record audit of that probability-stage experiment.


## IPW density-recovery diagnostic under finite validation

A second local experiment used 20 paired parent frames, (n_{parent}=4{,}000), (n_{val}=300), target retention 0.25, and curvature (gamma\in\{-1.2,0,1.2\}). Validation was sampled SRS from the full parent frame and excluded from the retained analysis catalogue. Linear and quadratic inclusion models were fitted only on validation labels. Each estimated/known IPW host density was compared with the analytic parent-host density using integrated absolute error (IAE).

Mean IAE by absolute SciPy bandwidth factor:

| Curvature | Arm | Factor 0.03 | Factor 0.05 | Factor 0.80 |
|---:|---|---:|---:|---:|
| -1.2 | Catalogue-only | 0.3716 | 0.3712 | 0.3104 |
| -1.2 | Known inclusion IPW | 0.0977 | 0.0888 | 0.3132 |
| -1.2 | Estimated degree-1 IPW | 0.1564 | 0.1465 | 0.3135 |
| -1.2 | Estimated degree-2 IPW | 0.1871 | 0.1781 | 0.3138 |
| 0.0 | Catalogue-only | 0.1760 | 0.1723 | 0.3115 |
| 0.0 | Known inclusion IPW | 0.0849 | 0.0779 | 0.3131 |
| 0.0 | Estimated degree-1 IPW | 0.1289 | 0.1200 | 0.3135 |
| 0.0 | Estimated degree-2 IPW | 0.1571 | 0.1485 | 0.3135 |
| +1.2 | Catalogue-only | 0.1030 | 0.0988 | 0.3134 |
| +1.2 | Known inclusion IPW | 0.0790 | 0.0709 | 0.3131 |
| +1.2 | Estimated degree-1 IPW | 0.1182 | 0.1123 | 0.3131 |
| +1.2 | Estimated degree-2 IPW | 0.1304 | 0.1215 | 0.3132 |

Interpretation:

- At factors 0.03 and 0.05, known-inclusion IPW improved density recovery relative to catalogue-only in all three tested curvature scenarios.
- Estimated inclusion was less accurate than known inclusion, as expected with finite validation. In the (gamma=+1.2) scenario, both estimated models were worse than catalogue-only at these settings; this is a useful candidate failure condition, but only a 20-replicate pilot.
- At factor 0.80, every arm had IAE around 0.31, and differences between catalogue-only, known, and estimated IPW were almost completely erased by smoothing.
- The lower IAE at factor 0.03/0.05 does not establish those as universally optimal bandwidths. They must be evaluated on held-out seeds and alternative parent laws.

This directly demonstrates that the current factor-0.8 smoothing choice can mask the very IPW differences the project is designed to measure. A density-recovery bandwidth audit is therefore a prerequisite, not a cosmetic sensitivity check.
