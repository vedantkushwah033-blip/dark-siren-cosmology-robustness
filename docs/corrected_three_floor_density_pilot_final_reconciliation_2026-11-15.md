# Corrected three-floor density pilot: final audit reconciliation

Date: 2026-11-15

## Status

The corrected five-width, 20-replicate density-only pilot (workflow `43619594762`, run head `17466326a4be59ca7698c6c1f34529801e6e5507`) completed successfully in all three catalogue-floor cells. Each cell's outer artifact audit passed and independently reconstructed all 20 successful replicate children.

This report records the completed workflow and its archived artifacts. The branch has advanced 22 commits beyond the workflow head; those later commits include documentation and further audit/test changes. The workflow result must therefore be cited with its exact run SHA rather than represented as a fresh run of the current branch head.

## Design and audit

- Parent population: 4,000 objects.
- Independent validation frame: 1,500; nested validation sizes: 100, 300, 1,000.
- 20 paired replicate seeds per floor, master seed 20261215, stride 10; validation seed is population seed + 1.
- Catalogue floors: 0.00, 0.10, 0.25.
- Five physical bandwidths: 0.012, 0.016, 0.024, 0.032, 0.040.
- Degree-one logistic inclusion model; quadrature order 5; 401-point redshift grid.
- Each floor artifact contains 9,300 aggregate metric rows and 1,740 realized-Scott-width rows.

| Floor | Planned | Successful | Failed | Child raw audits | Manifest files verified | Aggregate metric rows reconstructed | Outer audit |
|---:|---:|---:|---:|---:|---:|---:|---|
| 0.00 | 20 | 20 | 0 | 20 | 244 | 9,300 | PASS |
| 0.10 | 20 | 20 | 0 | 20 | 244 | 9,300 | PASS |
| 0.25 | 20 | 20 | 0 | 20 | 244 | 9,300 | PASS |

The earlier floor-0.25 optimizer failure did not recur in this corrected run. All three corrected artifacts have 20/20 successful replicates; do not carry forward the earlier 19/20 count.

## Results at h = 0.016, validation n = 1,000

Mean integrated absolute error (IAE), averaged over the 20 paired replicates:

| Arm | Floor 0.00 | Floor 0.10 | Floor 0.25 |
|---|---:|---:|---:|
| A catalogue-only | 0.53049 | 0.37809 | 0.21349 |
| B parent oracle | 0.18378 | 0.18378 | 0.18378 |
| C known-probability IPW | 0.20391 | 0.18653 | 0.18065 |
| D MAP estimated IPW | 0.20720 | 0.29618 | 0.27331 |
| E posterior-mean density | 0.21032 | 0.30006 | 0.27605 |

Paired replicate-level contrasts (mean difference; Monte Carlo standard error):

| Contrast | Floor 0.00 | Floor 0.10 | Floor 0.25 |
|---|---:|---:|---:|
| D − A | −0.32330 (0.01145) | −0.08192 (0.01470) | +0.05983 (0.00907) |
| D − C | +0.00329 (0.01189) | +0.10964 (0.01359) | +0.09266 (0.01012) |
| C − B | +0.02013 (0.00654) | +0.00275 (0.00401) | −0.00313 (0.00256) |
| E − D | +0.00312 (0.00035) | +0.00388 (0.00012) | +0.00274 (0.00008) |

Positive differences indicate larger IAE (worse density recovery) for the first arm named relative to the second.

## Interpretation limits

1. At floor 0.00, the degree-one logistic analysis model is correctly specified. MAP estimated IPW substantially improves on the catalogue-only density and is close to known-probability IPW at n=1,000.
2. At floors 0.10 and 0.25, the floor-offset inclusion curve is not represented by the fitted degree-one logistic family. Estimated IPW remains materially worse than known-probability IPW; at floor 0.25 it is also worse than the catalogue-only baseline at this fixed bandwidth.
3. The oracle arm is available only because this is a simulation. It is a diagnostic benchmark, not an ordinary-data estimator.
4. Arm E's primary metric is the loss of its posterior-mean density. Its node-averaged posterior expected loss is a separate secondary estimand; do not substitute one for the other.
5. The expanded 11-width workflow re-evaluated the same paired seeds. It is not an independent second set of replicates and must not be pooled with this run as additional sample size.

## What this does not establish

This is a low-redshift synthetic density-recovery pilot. It does not establish a selected bandwidth, confirmatory estimator superiority, H0 accuracy or calibration, performance on real galaxy catalogues, or validity for LVK/GWTC data. No H0 propagation should be treated as validated until the density-stage bandwidth protocol and the separate selection-conditioned likelihood/calibration gates are met.

## Provenance

- Workflow: [absolute-width-distribution-pilot #20](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/43619594762)
- Exact workflow head: `17466326a4be59ca7698c6c1f34529801e6e5507`
- Current branch head when this report was written: `9894307b78c569adcc137c27b530afbdaf2247bb`
- Expanded 11-width companion workflow: [run 43617052462](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/43617052462)
- Full numerical-diagnostics workflow on Python 3.10, 3.11, and 3.12: [run 43619577181](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/43619577181), all three matrix jobs passed.
