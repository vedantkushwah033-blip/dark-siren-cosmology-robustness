# Corrected gamma=0.5 20-replicate density pilot readout

**Date:** 2026-11-18  
**Status:** Independently audited exploratory density pilot; not confirmatory.  
**Scope:** Synthetic low-redshift Euclidean catalogue-inclusion experiment only.

## Provenance and audit

- Workflow: [GitHub Actions run 44366794840](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44366794840).
- Source branch: `research/gamma05-20rep-density-pilot-20261118`.
- Workflow head: `46ad562fbc03cae28c51f8373a3ab93b70431e46`.
- Design: 20 replicates; N_parent=4000; validation frame=1500; validation sizes 100, 300, 1000; gamma=0.5; beta=-3.9; quadrature order 5; redshift grid 401; Scott multiplier 0.7; absolute bandwidth candidates 0.012, 0.016, 0.020, 0.024, 0.028, 0.032, 0.040.
- Artifact ID: `14457626969`; size 70,214,078 bytes; GitHub artifact SHA-256 `90ab6a5cb0ad1e68d631ea6b8945387c882df38b154ad8106222efb850059007`.
- Outer audit: PASS; 20/20 successful; 0 failed; 20/20 child raw-input audits passed; 244 manifest files verified; 13,020 aggregate metric rows and 1,740 realized-Scott-width rows reconstructed.
- Independent post-download verification recomputed all 244 manifest byte counts and SHA-256 digests: 244/244 matched. The downloaded ZIP digest also matched GitHub's artifact digest.
- A/B sample weights were exactly one. C weights matched independently reconstructed inverse true inclusion probabilities to maximum absolute difference 7.1e-14. The common host-rate factor is retained as a downstream grid factor, not duplicated in sample-level A/B/C weights.

## Exploratory density results

For each replicate and validation size, I first averaged IAE over the seven fixed absolute-bandwidth candidates. I then summarized those replicate-level averages across the 20 paired replicates. This is a descriptive grid-average, **not** bandwidth selection and not a claim that the seven candidates are equally scientifically optimal.

| Validation n | A: catalogue only | B: parent oracle | C: known-probability IPW | D: degree-one MAP IPW | E: posterior-mean density |
|---:|---:|---:|---:|---:|---:|
| 100 | 0.7191 | 0.1009 | 0.1281 | 0.2223 | 0.2354 |
| 300 | 0.7191 | 0.1009 | 0.1281 | 0.1903 | 0.1993 |
| 1000 | 0.7191 | 0.1009 | 0.1281 | 0.1715 | 0.1746 |

Paired replicate-level contrasts against C (D minus C):

| Validation n | Mean paired IAE difference | 95% paired t interval | Replicates with D worse than C |
|---:|---:|---:|---:|
| 100 | 0.0943 | [0.0394, 0.1492] | 19/20 |
| 300 | 0.0622 | [0.0313, 0.0932] | 18/20 |
| 1000 | 0.0434 | [0.0268, 0.0601] | 19/20 |

E minus C was positive in 20/20, 19/20, and 19/20 replicates for n=100, 300, and 1000, respectively. These are exploratory paired t intervals over 20 replicates, not confirmatory inferential guarantees.

## Weight-instability diagnostics

Mean Kish ESS across replicates (same catalogue realization across validation sizes for A/B/C):

| Validation n | C: known pi | D: MAP plug-in | E: posterior-integrated node weights |
|---:|---:|---:|---:|
| 100 | 779 | 752 | 128 |
| 300 | 779 | 611 | 183 |
| 1000 | 779 | 603 | 298 |

For D, the mean top-1% normalized weight mass was approximately 12.5%, 13.5%, and 13.2% across n=100, 300, and 1000, compared with 10.3% for C. E's posterior-integrated node weights were more concentrated still (mean top-1% mass about 36.8%, 28.5%, and 20.9%). This supports inspecting weight concentration and ESS alongside density loss rather than assuming coefficient uncertainty integration stabilizes IPW.

## Interpretation and limitations

At gamma=0.5, this pilot does **not** show an advantage for the fitted degree-one estimated-inclusion IPW arms over the known-probability IPW reference. The MAP and posterior-mean-density arms have higher average IAE in this design, although their gap to C narrows as validation size increases. This is consistent with the risk that degree-one logistic misspecification and finite-sample coefficient estimation can offset the intended benefit of estimated inclusion correction. It does not establish a universal ranking.

This is density recovery only. There is no H0 likelihood, H0 posterior, selection-conditioned inference, SBC, coverage, or calibration analysis here. No bandwidth was selected. Oracle/reference probabilities are available only because this is a synthetic simulation. No public LVK/GWTC or observational galaxy-catalogue data are used.

**Metadata caveat:** child-level configs still inherit the single-replicate runner's phrase “one engineering replicate only” in their limitations field, even though the outer config, ledger, 20 child artifacts, and independent orchestration audit describe the 20-replicate pilot. This stale wording should be corrected in the runner before the artifacts are treated as publication-ready provenance; it does not change the recorded numerical outputs or the successful raw-input audit.

## Next step

A separate corrected seven-gamma paired 20-replicate workflow is running from a branch based directly on corrected PR #59 head `390ba806c42376b1a2d4b8d083b92b8d206c8a88`. Its per-cell workflow requires 20/20 successful replicates and independent child audits, preserves partial artifacts on failure, and runs the cross-gamma paired-seed/calibration audit only after all seven cells pass. The older seven-gamma run from the pre-correction branch remains quarantined and must not be used.
