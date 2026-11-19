# Paired H0 prior-support sensitivity readout — 2026-11-20

## Scope and frozen comparison

This is a paired synthetic sensitivity analysis of the same 100 replicate realizations under two uniform H0 prior supports:

- Narrow: Uniform(50, 90), source workflow `44460222957`, artifact `14504568415`.
- Expanded: Uniform(30, 110), source workflow `44484212311`, artifact `14520551336`.
- Both use master seed `2026111929`, 30 detected events, and five inference arms.

The paired audit workflow `44600400410` completed successfully. It verified:

- 5,708 manifest entries in each artifact;
- 4,500/4,500 lifecycle cells in each artifact;
- 300 replicate-by-curvature child pairs;
- 2,700 raw parent/event/validation input file hashes;
- 600 weight/bandwidth file hashes;
- 4,500 paired posterior rows;
- maximum common-support log-likelihood difference (2.27\times10^{-13}).

The common-support check establishes that the paired posterior changes are attributable to prior support and its normalization/truncation, not changed simulated inputs or likelihood calculations.

## Paired results at validation size 500

All differences below are expanded-prior minus narrow-prior. The reported intervals are paired 95% t intervals over 100 replicate pairs. Median-error difference is the signed posterior-median error shift; absolute-error difference is the change in absolute median error. Width is the change in the central 68% interval width. Coverage difference is the paired change in the 68% coverage indicator.

| Curvature | Arm | Median-error change | 95% interval | Absolute-error change | 95% interval | 68% width change | 95% interval |
|---:|---|---:|---:|---:|---:|---:|---:|
| -1 | Catalogue-only A | -2.281 | [-2.763, -1.799] | +2.281 | [1.799, 2.763] | +2.958 | [2.511, 3.405] |
| -1 | Parent oracle B | +0.001 | [-0.0003, 0.0030] | +0.001 | [-0.0003, 0.0030] | +0.003 | [-0.0006, 0.0066] |
| -1 | Known-π IPW C | +0.010 | [-0.0087, 0.0280] | +0.011 | [-0.0068, 0.0298] | +0.020 | [-0.0105, 0.0514] |
| -1 | MAP IPW D | -0.001 | [-0.0063, 0.0038] | +0.006 | [0.0009, 0.0107] | +0.013 | [0.0024, 0.0230] |
| -1 | Integrated IPW E | -0.001 | [-0.0071, 0.0051] | +0.007 | [0.0008, 0.0127] | +0.015 | [0.0027, 0.0271] |
| 0 | Catalogue-only A | -2.477 | [-3.026, -1.929] | +2.477 | [1.929, 3.026] | +3.339 | [2.767, 3.911] |
| 0 | Parent oracle B | +0.001 | [-0.0003, 0.0030] | +0.001 | [-0.0003, 0.0030] | +0.003 | [-0.0006, 0.0066] |
| 0 | Known-π IPW C | +0.004 | [-0.0037, 0.0119] | +0.005 | [-0.0029, 0.0127] | +0.010 | [-0.0052, 0.0248] |
| 0 | MAP IPW D | +0.005 | [-0.0059, 0.0156] | +0.007 | [-0.0033, 0.0181] | +0.014 | [-0.0056, 0.0343] |
| 0 | Integrated IPW E | +0.005 | [-0.0064, 0.0174] | +0.008 | [-0.0035, 0.0202] | +0.016 | [-0.0056, 0.0378] |
| +1 | Catalogue-only A | -2.661 | [-3.338, -1.984] | +2.661 | [1.984, 3.338] | +3.595 | [2.906, 4.283] |
| +1 | Parent oracle B | +0.001 | [-0.0003, 0.0030] | +0.001 | [-0.0003, 0.0030] | +0.003 | [-0.0006, 0.0066] |
| +1 | Known-π IPW C | +0.004 | [-0.0001, 0.0087] | +0.005 | [0.0001, 0.0089] | +0.009 | [0.0003, 0.0187] |
| +1 | MAP IPW D | +0.029 | [-0.0014, 0.0596] | +0.029 | [-0.0012, 0.0597] | +0.050 | [-0.0001, 0.1000] |
| +1 | Integrated IPW E | +0.029 | [-0.0026, 0.0596] | +0.029 | [-0.0024, 0.0598] | +0.048 | [-0.0012, 0.0971] |

The paired 68% coverage indicator did not change for any arm/curvature cell in this 100-pair n=500 comparison. This means no replicate crossed the central-interval coverage boundary when changing support; it does not imply that the coverage level itself is nominal.

## Interpretation

1. **Catalogue-only A is strongly prior-support sensitive.** Expanding the support shifts its posterior median error downward by about 2.3–2.7 km/s/Mpc and increases its absolute median error by the same amount. Its central 68% interval also widens by about 3.0–3.6 km/s/Mpc. All paired 95% intervals exclude zero.
2. **Parent-oracle B is essentially invariant** to this support change in the paired sample; all paired intervals include zero and effects are numerically tiny.
3. **Known-π IPW C is comparatively stable.** Changes are small; only the positive-curvature absolute-error and interval-width contrasts narrowly exclude zero.
4. **Estimated-IPW D/E show small-to-moderate changes**, especially at positive curvature, but paired intervals are wide and include zero for signed median error and interval width. There is no evidence here that integration E is more prior-robust than MAP D.
5. The observed catalogue-only sensitivity is not a universal H0-likelihood boundary effect: it is arm-dependent. The same simulated events and common-support likelihood curves were used across the pair.

## Limitations

- Fixed-truth synthetic experiment; not prior-predictive SBC under the expanded prior.
- One matched 100-replicate seed stream, three curvature settings, and 30 detected events per replicate.
- This readout reports n=500; n=100 and n=300 paired outputs remain in the uploaded audit artifact and should be summarized before any pooled paper-level claim.
- Paired t intervals describe Monte Carlo uncertainty for the mean paired contrasts; they do not certify frequentist coverage or estimator optimality.
- No observational H0 or real galaxy-catalogue claim follows.

## Reproducibility

- Audit PR: [#125](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/pull/125).
- Audit workflow: `44600400410`.
- Audit artifact: `14567533099`, SHA-256 `627311b9c186bd4f5dfd4af6d935de3c9ce847b804b4cf08b86fb5afd294a986`.
- Source artifact hashes are verified by the audit script before any paired summary is calculated.

**Status:** paired prior-support attribution gate PASS. Full repository test matrix and scientific review remain separate gates. Draft/unmerged; no main-branch changes.


## Complete validation-size follow-up (n=100, 300, 500)

The updated audit workflow on commit `f61d6affd6fbe53290c0e5561556f815faa5d31e` printed all 180 paired metric rows (3 curvatures × 3 validation sizes × 5 arms × 4 metrics). The audit and its focused regression tests both completed successfully in workflow `44601966016`. The independent full repository test workflow was still pending at the time of this readout.

The main additional result is that the prior-support effect on the estimated-inclusion arms is not a simple monotone function of validation size. The table reports the change in absolute posterior-median error (expanded minus narrow prior), in km/s/Mpc, with paired 95% t intervals over 100 matched replicates:

| Curvature | n_validation | MAP IPW D | Integrated IPW E |
|---:|---:|---:|---:|
| -1 | 100 | +0.0139 [0.0065, 0.0212] | +0.0158 [0.0075, 0.0240] |
| -1 | 300 | +0.0082 [0.0010, 0.0155] | +0.0083 [0.0023, 0.0144] |
| -1 | 500 | +0.0058 [0.0009, 0.0107] | +0.0067 [0.0008, 0.0127] |
| 0 | 100 | +0.0121 [-0.0055, 0.0297] | +0.0130 [-0.0037, 0.0298] |
| 0 | 300 | +0.0062 [-0.0022, 0.0146] | +0.0073 [-0.0028, 0.0174] |
| 0 | 500 | +0.0074 [-0.0033, 0.0181] | +0.0084 [-0.0035, 0.0202] |
| +1 | 100 | +0.0268 [-0.0022, 0.0558] | +0.0252 [-0.0011, 0.0514] |
| +1 | 300 | +0.0225 [0.0012, 0.0438] | +0.0231 [0.0002, 0.0460] |
| +1 | 500 | +0.0293 [-0.0012, 0.0597] | +0.0287 [-0.0024, 0.0598] |

Additional controls:
- Catalogue-only A's absolute-error increase is validation-size invariant: +2.281 [1.799, 2.763] at gamma=-1; +2.477 [1.929, 3.026] at gamma=0; +2.661 [1.984, 3.338] at gamma=+1.
- Known-probability IPW C is also invariant across validation sizes, as expected because its inclusion probabilities are known rather than estimated: +0.0115 [-0.0068, 0.0298], +0.0049 [-0.0029, 0.0127], and +0.0045 [0.0001, 0.0089] for gamma=-1, 0, +1 respectively.
- Parent-oracle B changes are numerically tiny: +0.0014 [-0.0003, 0.0030] across these paired cells.
- Coverage-indicator changes are zero in all but one of the 45 arm-curvature-validation cells: E at gamma=0, n=100 changes by -1 percentage point, with paired interval [-2.98, +0.98] percentage points.

Interpretation: widening the prior materially worsens the catalogue-only arm, while the oracle and known-probability IPW controls are nearly invariant. D/E show much smaller shifts, with some curvature/validation-size cells excluding zero, but the effect is modest and does not show a consistent improvement for posterior integration E over MAP D. The 45 cell-wise intervals are exploratory and unadjusted for multiplicity; do not read isolated exclusions as confirmatory discoveries.

This is still one fixed-truth, one matched 100-replicate seed stream. It is a prior-support sensitivity audit, not an expanded-prior SBC calibration study or evidence about observational H0.


## Expanded validation-size readout (n=100, 300, 500)

Workflow 44602206708 reran the complete paired audit after the fail-closed pairing checks were added. Its unit tests, artifact downloads, full reconstruction, and report upload all completed successfully. The output covers every validation size. For estimated-IPW arms D/E, the table gives the change in absolute median error and central 68% interval width (expanded minus narrow prior), each with paired 95% t intervals over 100 matched replicates. Values are km/s/Mpc.

| Curvature | n | Arm | Absolute-error change [95% CI] | 68% width change [95% CI] |
|---:|---:|---|---:|---:|
| -1 | 100 | D MAP | +0.0139 [0.0065, 0.0212] | +0.0298 [0.0144, 0.0452] |
| -1 | 100 | E integrated | +0.0158 [0.0075, 0.0240] | +0.0374 [0.0199, 0.0550] |
| -1 | 300 | D MAP | +0.0082 [0.0010, 0.0155] | +0.0192 [0.0015, 0.0370] |
| -1 | 300 | E integrated | +0.0083 [0.0023, 0.0144] | +0.0215 [0.0046, 0.0384] |
| -1 | 500 | D MAP | +0.0058 [0.0009, 0.0107] | +0.0127 [0.0024, 0.0230] |
| -1 | 500 | E integrated | +0.0067 [0.0008, 0.0127] | +0.0149 [0.0027, 0.0271] |
| 0 | 100 | D MAP | +0.0121 [-0.0055, 0.0297] | +0.0226 [-0.0075, 0.0526] |
| 0 | 100 | E integrated | +0.0130 [-0.0037, 0.0298] | +0.0255 [-0.0026, 0.0537] |
| 0 | 300 | D MAP | +0.0062 [-0.0022, 0.0146] | +0.0131 [-0.0033, 0.0295] |
| 0 | 300 | E integrated | +0.0073 [-0.0028, 0.0174] | +0.0151 [-0.0038, 0.0341] |
| 0 | 500 | D MAP | +0.0074 [-0.0033, 0.0181] | +0.0144 [-0.0056, 0.0343] |
| 0 | 500 | E integrated | +0.0084 [-0.0035, 0.0202] | +0.0161 [-0.0056, 0.0378] |
| +1 | 100 | D MAP | +0.0268 [-0.0022, 0.0558] | +0.0461 [-0.0014, 0.0937] |
| +1 | 100 | E integrated | +0.0252 [-0.0011, 0.0514] | +0.0426 [0.0010, 0.0841] |
| +1 | 300 | D MAP | +0.0225 [0.0012, 0.0438] | +0.0414 [0.0029, 0.0800] |
| +1 | 300 | E integrated | +0.0231 [0.0002, 0.0460] | +0.0409 [0.0019, 0.0799] |
| +1 | 500 | D MAP | +0.0293 [-0.0012, 0.0597] | +0.0499 [-0.0001, 0.1000] |
| +1 | 500 | E integrated | +0.0287 [-0.0024, 0.0598] | +0.0479 [-0.0012, 0.0971] |

Interpretation remains cautious:

- At negative curvature, widening the prior produces small but consistently positive changes in absolute error and interval width for D/E; the changes shrink as validation size increases.
- At zero curvature, paired intervals for these D/E changes include zero at every validation size.
- At positive curvature, changes are larger, particularly at n=300, but many intervals are close to zero. This is not evidence of a robust D-versus-E difference.
- Catalogue-only A's prior effect is identical across validation sizes because A does not use the validation sample. Its absolute-error increase is +2.281, +2.477, and +2.661 for curvature -1, 0, and +1, respectively; all paired 95% intervals exclude zero.
- No 68% coverage indicator changed between the two priors in any paired cell. This is a statement about crossing the interval boundary, not about nominal coverage.

The paired-artifact integrity gate again passed: 5,708 manifest entries per artifact, 4,500 lifecycle cells per artifact, 300 child pairs, 2,700 raw-input hashes, 600 weight/bandwidth hashes, 4,500 posterior rows, and maximum common-support log-likelihood discrepancy 2.27e-13. The new pairing-key and truth-consistency unit tests passed before the reconstruction.

Workflow: 44602206708 (success). Separate repository-wide test workflow: 44602206997 (status pending at report time; do not treat the audit's embedded unit-test step as a substitute for the repository-wide matrix).
