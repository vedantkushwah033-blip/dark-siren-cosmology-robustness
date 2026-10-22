# Reproducing the original complete-catalogue coverage result

## What this workflow does

The manually triggered workflow imports the original project's `run_repeated`
and `summarize_runs` functions from a pinned commit and reruns the
complete-catalogue, uniform-host-generation, uniform-inference-weighting
scenario. It saves per-run rows and a summary JSON artifact.

Calibration report's source commit: `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`.

This is intended to test whether the reported high coverage can be reproduced
from the source code, rather than to reimplement the original likelihood in
the new repository.

## Run it

1. Open the repository's Actions tab.
2. Choose `original-coverage-reproduction`.
3. Click `Run workflow`.
4. Use 1,000 replicates and seed 2026 for the first reproduction.
5. After completion, download the `original-coverage-reproduction` artifact.

GitHub's manual workflow controls require the workflow file to exist on the
default branch. See the official [GitHub manual workflow documentation](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/manually-run-a-workflow).

## What to compare

Compare the reproduced `coverage_68`, `coverage_95`, interval widths,
bias, RMSE, and pull diagnostics with `docs/baseline_calibration.md` in the
original repository. Small numerical differences can arise from code changes,
dependency versions, or different commits; the pinned commit and seed make
the comparison more controlled.

Do not call this a reproduction until the workflow completes and the output
artifact has been inspected. A matching coverage value alone does not prove
that the model is scientifically realistic or correctly specified.

## Interpretation

The source-level review indicates that the original posterior is a mixture
over candidate redshifts using a low-redshift distance relation and a uniform
prior over H0. Such a mixture can be broad because many catalogue galaxies
can partially explain the observed distance at different H0 values. However,
that is a possible mechanism, not yet a demonstrated explanation for the
reported 94.35% coverage. The rerun is the next evidence-gathering step.
