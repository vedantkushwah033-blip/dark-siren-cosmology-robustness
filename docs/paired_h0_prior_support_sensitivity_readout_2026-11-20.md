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
- One matched 100-replicate seed stream, three curvature settings, four events? **Correction: this experiment uses 30 detected events per replicate.**
- This readout reports n=500; n=100 and n=300 paired outputs remain in the uploaded audit artifact and should be summarized before any pooled paper-level claim.
- Paired t intervals describe Monte Carlo uncertainty for the mean paired contrasts; they do not certify frequentist coverage or estimator optimality.
- No observational H0 or real galaxy-catalogue claim follows.

## Reproducibility

- Audit PR: [#125](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/pull/125).
- Audit workflow: `44600400410`.
- Audit artifact: `14567533099`, SHA-256 `627311b9c186bd4f5dfd4af6d935de3c9ce847b804b4cf08b86fb5afd294a986`.
- Source artifact hashes are verified by the audit script before any paired summary is calculated.

**Status:** paired prior-support attribution gate PASS. Full repository test matrix and scientific review remain separate gates. Draft/unmerged; no main-branch changes.
