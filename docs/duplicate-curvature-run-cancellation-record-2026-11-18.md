# Duplicate curvature run cancellation record

The standalone corrected density workflow `44367512399` is retained as the sole 11-bandwidth execution for this branch. The duplicate `tests.yml` run `44367673265` uses the same seed schedule, gamma grid, replicate count, and bandwidth grid, so it adds no independent simulation information. The branch-specific `tests.yml` concurrency policy was changed to cancel that duplicate and its pending replacement. The separate corrected 7-bandwidth run `44369630799` remains the primary paired-curvature pilot.
