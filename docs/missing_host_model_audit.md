# Scientific audit checkpoint: missing-host mixture control

## What this experiment tests

The missing-host control compares an observed-catalogue-only likelihood with a mixture likelihood that adds redshifts for the hidden part of a simulated catalogue. The original simulator is pinned to commit `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`.

The `oracle_true` option uses the simulation's true host-generation weights to construct the hidden-population component. This is deliberately an oracle control: a real analysis cannot ordinarily inspect the unobserved galaxy population or know its true host weights.

## Interpretation constraint

The current control keeps the observed-catalogue inference weights uniform in all six scenarios, while the oracle scenario uses true host-generation weights to calculate the hidden component's prior mass and internal mixture weights. Therefore, the oracle comparison is a hybrid diagnostic: it changes the missing component using privileged simulation information, rather than providing a fully oracle-weighted posterior over both observed and hidden galaxies.

Do not describe this as a direct comparison between a practical estimator and a fully oracle-optimal estimator. It is a controlled sensitivity test of how the hidden component is constructed while the observed component remains misspecified.

## Recommended follow-up before publication

Add a separately named, internally consistent oracle benchmark in which the same true host-generation weights are used for both observed and hidden galaxies, and the total prior mass is normalized across the complete simulated population. Compare at least:

1. observed-only inference with uniform weights;
2. missing-host term using inference-model weights;
3. hidden-component oracle weights with the observed component unchanged (current hybrid control);
4. full-population oracle weighting over observed plus hidden galaxies.

Keep the full-population oracle result labelled as a simulation-only upper-information benchmark. Report coverage, bias, RMSE, interval width, and Monte Carlo uncertainty together. The oracle benchmark should not be used as evidence that the same information is available in real catalogues.

## Workflow status

The automated unit-test workflow passing does not mean this experiment has run. Before interpreting results, run the manual `factorial-sensitivity` workflow, verify all steps succeed, and inspect the artifact files. No scientific conclusion should be drawn from unexecuted scenario definitions.
