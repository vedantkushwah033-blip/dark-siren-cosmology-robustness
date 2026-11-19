# Full-parent paired event-count calibration readout

**Workflow:** [44447250706](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44447250706)  
**Status:** successful; all planned rows and artifact hashes verified.  
**Scope:** synthetic full-parent fixed-truth diagnostic; no catalogue thinning or inclusion fitting.

## Design

- 500 independent parent-population replicates.
- 1,200 parent galaxies per replicate.
- Fixed generating H0 = 70 km s^-1 Mpc^-1.
- Uniform H0 prior [20,300], 561-point H0 grid.
- For each replicate, the first four detected events are compared with the nested first thirty events from the same generated 30-event sequence.
- The same full-parent host mixture and selected-event likelihood are used in both event-count conditions.
- 1,000 raw rows total; 500 paired replicate comparisons.

## Results

| Event count | Replicates | Mean posterior-median bias | Bias MCSE | RMSE | 68% coverage | 90% coverage | 95% coverage | Mean 68% width |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 500 | +5.060 | 0.717 | 16.793 | 69.6% | 91.0% | 95.4% | 23.436 |
| 30 | 500 | +0.866 | 0.182 | 4.160 | 66.2% | 86.0% | 92.0% | 7.355 |

The paired 30-event minus 4-event median-bias difference is **-4.194 km s^-1 Mpc^-1** (paired MCSE 0.699; 95% t interval [-5.567, -2.820]).

The 4-event case has broad posterior intervals and a positive finite-event median bias, while 30 events materially reduce both bias and RMSE. The 68% coverage estimate is near nominal in both conditions given its Monte Carlo uncertainty; the 30-event Wilson interval is [61.9%, 70.2%].

## Integrity

The workflow passed its nested-design and paired-summary tests, completed all 500 paired replicates (500 successes at each event count), verified all four aggregate manifest entries, reconciled 1,000 raw rows, and checked the paired contrast denominator.

## Interpretation

This experiment establishes that finite detected-event count alone can generate a positive posterior-median bias of roughly +5 km s^-1 Mpc^-1 in this full-parent toy setup with four events, while the paired 30-event bias is much smaller. It provides an important control for interpreting the separate four-event catalogue-inclusion calibration.

It does **not** validate any inclusion correction. The full-parent event-count experiment uses the exact parent host mixture; the inclusion-calibration pilot uses a finite catalogue, fitted inclusion probabilities, IPW weights, and a different posterior computation. Differences between those studies cannot be attributed solely to catalogue correction. The paired 30-event five-arm study is the appropriate next test of the inclusion pipeline.
