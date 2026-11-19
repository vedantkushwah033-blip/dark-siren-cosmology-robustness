# H0 grid-refinement readout for the 100-replicate fixed-truth pilot

**Source run:** `44442239902`  
**Artifact:** `14493625716` (`h0-100rep-grid-refinement-281-561-1121`)  
**Interpretation:** numerical posterior-grid sensitivity only; not additional calibration replicates.

## Integrity and execution

- The archived baseline artifact contains 14,106 manifest entries; all were verified by the workflow.
- The refinement output contains 1,281 manifest entries; an independent SHA-256 pass found zero mismatches.
- 18 baseline successful cells were selected (three replicates in each of six gamma × validation-size cells).
- 54 refined child runs completed successfully; all 54 child audits passed.
- The same archived observables/truth were used across 281-, 561-, and 1,121-point reruns.

## Comparisons

Predeclared tolerances: maximum posterior-density difference ≤ 0.001 and maximum absolute shift among reported posterior quantiles ≤ 0.05 km s^-1 Mpc^-1.

| Grid comparison | Paired comparisons | Within both tolerances | Maximum density difference | Maximum quantile shift |
|---|---:|---:|---:|---:|
| 281 → 561 | 36 | 36/36 | 0.00029696 | 0.04378 |
| 561 → 1,121 | 18 | 18/18 | 0.00007444 | 0.01179 |
| 281 → 1,121 | 18 | 15/18 | 0.00029694 | 0.05430 |

The 281→561 and 561→1,121 adjacent comparisons pass throughout the tested subset. Three direct 281→1,121 comparisons narrowly exceed the 0.05 quantile threshold (maximum 0.05430). This is consistent with small cumulative discretization shifts even when each adjacent refinement is within tolerance; the direct comparison should not be hidden.

## Interpretation and limitations

For this selected subset, 561 points appears numerically stable relative to 1,121 points under the stated tolerances. This is not a full audit of all 600 cells in the 100-replicate experiment, and the 18 selected cells do not establish uniform grid convergence across the full design. It also does not validate the likelihood, calibrate the posterior, or change fixed-truth coverage. Retain the 1,121-point grid as a reference for further posterior-grid checks, and report the three direct 281→1,121 tolerance exceedances.

The workflow report's `status: PASS` means all planned refinement runs and child audits completed; it does **not** mean every comparison passed. The machine-readable report correctly records 87/90 comparisons within both tolerances and `all_comparisons_within_tolerance: false`.
