# Matched host-structure misspecification control protocol

**Status:** implemented; results are not considered available until the workflow artifact has completed and been inspected.  
**Scope:** a synthetic statistical control, not a gravitational-wave detector or real galaxy-catalogue analysis.

## Question

How much does collapsing a two-host redshift mixture into one effective host change fixed-truth inference when both analyses use exactly the same detected distance observations?

## Data-generating process

For each fixed truth (H_0in\{52,55,60,65,70,75,80,85,88\}), sample a host from redshifts (z=(0.02,0.04)) with probabilities (w=(0.7,0.3)). Generate a noisy distance from a Gaussian with mean (cz/H_0) and standard deviation 20 Mpc, retaining only observations strictly within (0,120) Mpc. Use 300 detected observations per truth, a fixed seed, and a common grid of 4,001 points from 25 to 125 with a uniform prior on [50,90].

## Analysis arms

Each observed distance is fitted by four arms:
1. two-host mixture without selection normalization;
2. two-host mixture with the toy selection normalization;
3. one collapsed host at the mixture-weighted redshift (z_*=0.7(0.02)+0.3(0.04)=0.026), without selection normalization;
4. the collapsed host with the same toy selection normalization.

The true and assumed detection intervals are the same in this experiment. Thus this is not a threshold-misspecification experiment. It isolates host-mixture representation under the specified toy likelihood, while retaining the fixed prior and selection-conditioned data-generation process.

## Outcomes and uncertainty

Report bias, mean absolute error, RMSE, 68% and 95% interval coverage with Wilson binomial intervals, and mean interval widths for each truth and method. For each truth and correction status, compute the paired difference in absolute error (collapsed minus two-host) using the same observations and a percentile bootstrap with 5,000 resamples. Positive paired differences mean the collapsed model has greater absolute error on average; negative values favor the collapsed model for that metric.

## Interpretation limits

- A weighted-mean redshift is not generally likelihood-equivalent to a mixture of host redshifts; differences here are a model-compression diagnostic.
- The experiment does not prove that host mixture is the cause of previous edge undercoverage; prior boundaries, noise, and selection conditioning remain relevant.
- No result should be described as a real (H_0) measurement or validation of LVK selection effects.
