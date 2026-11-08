# Fixed-analysis-catalogue confirmatory density-recovery results

**Date:** 2026-11-08  
**Scope:** synthetic host-density recovery only. These results are not an H0 posterior-calibration or real-data result.

## 1. Artifact provenance and design

The archived 100-replicate fixed-analysis-catalogue artifact was downloaded from GitHub Actions run [#42476989338](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/42476989338), artifact ID `13616760895`. Its manifest records generating commit `c18de7b013cdf583a710dba552daf796fad0e994`.

Frozen design in that artifact:
- 100 replicates;
- 4,000 independent analysis-frame parent units per replicate/scenario;
- independent validation frame of maximum size 1,000, with nested sizes 100, 300 and 1,000;
- curvature \\(\gamma\in\{-1.2,0,+1.2\}\\);
- target retention 0.25;
- degree-1 and degree-2 logistic inclusion fits;
- Scott multiplier 0.15;
- 801-point redshift grid on [0.01, 0.40];
- five arms: catalogue-only, known-inclusion IPW, parent-frame KDE, estimated degree-1 IPW and estimated degree-2 IPW.

Unlike the earlier disjoint-validation design, the validation frame here is independent of the analysis frame. Increasing validation size therefore does not remove units from the analysis catalogue. The two designs answer related but distinct questions and must not be pooled as if they were one experiment.

## 2. Raw-artifact audit

The archived raw artifact's independent auditor reports:
- 4,500/4,500 density attempts reconstructed;
- 300 analysis/validation frame-scenario partitions checked;
- 100 RNG stream sets checked;
- all 4,500 attempts SUCCESS;
- audit status PASS.

The manifest's four raw/configuration files (attempt ledger, frame units, configuration and summary) were independently re-hashed locally; all byte counts and SHA-256 values matched.

A separate direct Gaussian-kernel calculation was also performed locally on 15 selected raw attempts: one replicate for each of the three curvature values and each of the five estimator arms. This calculation used an explicit weighted Gaussian sum and weighted covariance formula, not `scipy.stats.gaussian_kde`. Maximum absolute difference from the stored density vectors was \\(1.69\times10^{-14}\\); maximum absolute IAE discrepancy was below \\(3.4\times10^{-16}\\).

The direct-formula spot check is now being incorporated into both repository auditors so future smoke and confirmatory runs execute it automatically. The full current-head workflow for those new audit changes is a separate verification gate; the local reconstruction above is not a claim that the updated auditor has already passed CI.

## 3. Predeclared primary interaction

The analysis JSON defines

\\[
\Delta(\gamma,n_V,d)=
\operatorname{IAE}_{\mathrm{estimated\ IPW}}-
\operatorname{IAE}_{\mathrm{catalogue-only}},
\\]

and the pilot-informed primary interaction

\\[
\Psi=
\Delta(+1.2,300,1)-\Delta(-1.2,300,1),
\\]

at Scott multiplier 0.15.

Independent reconstruction from the raw attempt CSV gave:

| Quantity | Reconstructed value |
|---|---:|
| Mean interaction \\(\hat\Psi\\) | +0.2778930 |
| Monte Carlo SE | 0.0054182 |
| Paired bootstrap 95% interval | [0.2674495, 0.2886512] |
| Paired-t 95% interval | [0.2671421, 0.2886440] |
| Complete paired replicates | 100/100 |

The independently recomputed primary values exactly match the stored analysis JSON at floating-point precision.

**Selection caveat:** the symmetric-curvature interaction was motivated by earlier exploratory work. It is therefore a pilot-informed confirmatory-style benchmark contrast, not a hypothesis selected independently of all prior results. The interval quantifies Monte Carlo variability under this declared synthetic DGP; it is not a real-universe uncertainty interval.

## 4. Estimated-IPW versus catalogue-only density error

Paired mean IAE differences at Scott multiplier 0.15 are shown below. Negative values favour estimated IPW; positive values mean estimated IPW has larger IAE than catalogue-only.

### Degree 1

| Curvature \\(\gamma\\) | Validation 100 | Validation 300 | Validation 1,000 |
|---:|---:|---:|---:|
| -1.2 | -0.22186 (SE 0.00778) | -0.25927 (0.00427) | -0.27864 (0.00357) |
| 0.0 | -0.02773 (0.00878) | -0.07394 (0.00466) | -0.09350 (0.00377) |
| +1.2 | +0.06619 (0.00850) | +0.01863 (0.00366) | +0.00245 (0.00223) |

### Degree 2

| Curvature \\(\gamma\\) | Validation 100 | Validation 300 | Validation 1,000 |
|---:|---:|---:|---:|
| -1.2 | -0.10225 (SE 0.02199) | -0.21422 (0.00856) | -0.26932 (0.00441) |
| 0.0 | +0.08388 (0.02106) | -0.04093 (0.00788) | -0.08152 (0.00419) |
| +1.2 | +0.13209 (0.01895) | +0.02982 (0.00596) | -0.00158 (0.00285) |

For the positive-curvature degree-1 cell at validation size 1,000, the paired-t 95% interval is approximately [-0.00198, 0.00688]. For degree 2 it is approximately [-0.00723, 0.00407]. These intervals include zero. At positive curvature and validation sizes 100 and 300, both fitted models have positive paired IAE differences in this DGP.

The 18 secondary paired contrasts and their Holm adjustments were independently reconstructed from the raw attempt CSV. Maximum absolute discrepancy across the reported means, Monte Carlo SEs, p-values and Holm-adjusted p-values was \\(2.3\times10^{-16}\\).

## 5. Oracle/reference arms

At validation size 300, paired differences against catalogue-only were independently calculated for the two reference arms:

| Curvature | Known-inclusion IPW minus catalogue-only | Parent-frame KDE minus catalogue-only |
|---:|---:|---:|
| -1.2 | -0.29857 (SE 0.00343) | -0.34058 (SE 0.00315) |
| 0.0 | -0.10540 (0.00329) | -0.13964 (0.00305) |
| +1.2 | -0.02065 (0.00158) | -0.05162 (0.00229) |

The known-inclusion arm is an oracle diagnostic because it uses true inclusion probabilities. The parent-frame KDE is also an oracle diagnostic because it uses the full analysis-frame parent redshifts. Neither is an operational estimator when the corresponding truth is unavailable.

These reference arms show that IPW itself is not universally harmful in this DGP: known-probability IPW improves density IAE relative to catalogue-only at all three curvature values. The estimated-probability arms can nevertheless be worse in some finite-validation/misspecification regimes.

## 6. Interpretation

The defensible synthetic density-level observation is conditional:
- estimated IPW can improve density recovery substantially in some curvature regimes;
- finite validation and model complexity can produce a variance penalty;
- under positive curvature, the fitted-IPW arms are worse than catalogue-only at validation sizes 100 and 300 in this configuration;
- the positive-curvature penalty becomes statistically indistinguishable from zero at validation size 1,000 for both fitted degrees;
- the oracle known-inclusion arm improves density recovery in every tested curvature condition.

This is evidence of regime dependence for the declared synthetic DGP, not a universal ranking of IPW methods. The fixed catalogue across validation sizes makes the validation-size comparisons cleaner than the earlier design in which validation units were removed from the analysis catalogue.

## 7. Analysis-audit provenance and remaining gates

The archived artifact was generated at commit `c18de7b...`, whose workflow version uploaded the raw density audit and confirmatory analysis JSON but did not yet run the separate confirmatory-analysis arithmetic auditor. That JSON was therefore independently recomputed locally as described above. The current workflow definition now includes a separate analysis-audit step, and the current repository tests exercise the arithmetic auditor; a new current-head confirmatory artifact and its audit report must still be verified before claiming that the revised workflow has completed end to end.

The next scientific gates remain:
1. verify the revised full density artifacts and direct-kernel audit crosschecks on current-head workflows;
2. run a repeated-sampling H0 propagation experiment using the same finite-validation design;
3. assess H0-grid/redshift-grid convergence and prior-scale sensitivity;
4. report H0-specific likelihood-contribution ESS and marginalization MCSE;
5. conduct fixed-truth coverage and, separately, prior-drawn SBC;
6. retain the synthetic-only scope and avoid claims about LVK/GWTC data or real galaxy surveys.

**No H0 posterior-calibration, coverage, SBC, detector-response, or real-data conclusion follows from this density-recovery result.**
