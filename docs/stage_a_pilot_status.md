# Stage A selection-misspecification pilot

This is an exploratory synthetic pilot, not real LVK inference or a measurement of H0.

The 200-replicate artifact audit passed. Repository CI passed on Python 3.10, 3.11, and 3.12. The conditional long-run sensitivity steps were skipped in this run.

The catalogue-only arm showed mean H0 bias -9.468 km s^-1 Mpc^-1, RMSE 9.984, 68% interval coverage 7.5%, and mean density ISE 0.87242. Correct weighting showed bias +0.351, RMSE 3.785, coverage 71.0%, and ISE 0.03829. The analytic-truth control had RMSE 3.754 despite zero density ISE, so finite-event information and other inference approximations remain to be separated.

Do not interpret these numbers as real-world H0 performance. Oracle and analytic-truth arms use truth information and are diagnostic controls. Next steps are event-level posterior reconstruction, finite-parent versus continuous-target checks, posterior-grid refinement, and regeneration of host-density artifacts under the corrected weighted-Scott bandwidth convention.
