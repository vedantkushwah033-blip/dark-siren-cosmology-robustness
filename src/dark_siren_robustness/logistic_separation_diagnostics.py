"""Linear-programming diagnostics for separation in binary logistic designs.

This module diagnoses the geometry of a supplied design/outcome pair. It does
not fit a logistic model, alter an estimator, or decide whether inverse weights
are scientifically usable. Coefficient directions are L1-normalized only
inside the diagnostic LP to make the separation margin finite.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import linprog


def diagnose_logistic_separation(
    design: ArrayLike,
    included: ArrayLike,
    *,
    tolerance: float = 1e-9,
) -> dict[str, object]:
    """Classify overlap, complete separation, or quasi-complete separation.

    For y_i in {0,1}, let s_i = 2*y_i - 1. Complete separation exists when
    some coefficient direction b has s_i * x_i.T @ b > 0 for every row.
    Quasi-complete separation exists when all signed margins can be made
    nonnegative and at least one can be made strictly positive, but no
    direction makes every margin strictly positive.

    Two linear programs test these conditions after harmless per-column
    scaling. The L1 constraint on the positive/negative split of b fixes its
    arbitrary scale. A failed LP is reported as indeterminate rather than
    being mislabeled overlap.

    Rank deficiency is reported separately: it can coexist with any
    separation category and is not itself evidence of separation.
    """
    x = np.asarray(design, dtype=float)
    raw_y = np.asarray(included)
    if x.ndim != 2 or x.shape[0] < 2 or x.shape[1] < 1:
        raise ValueError("design must be a 2D matrix with at least two rows and one column")
    if not np.all(np.isfinite(x)):
        raise ValueError("design must contain only finite values")
    if raw_y.ndim != 1 or raw_y.shape[0] != x.shape[0]:
        raise ValueError("included must be a 1D binary vector matching design rows")
    if not np.all(np.isin(raw_y, [0, 1, False, True])):
        raise ValueError("included must contain only binary 0/1 values")
    y = raw_y.astype(float)
    if np.unique(y).size != 2:
        raise ValueError("included must contain both outcome classes")
    if not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("tolerance must be finite and positive")

    column_scale = np.max(np.abs(x), axis=0)
    safe_scale = np.where(column_scale > 0.0, column_scale, 1.0)
    scaled_x = x / safe_scale
    rank = int(np.linalg.matrix_rank(scaled_x))
    dimension = int(scaled_x.shape[1])
    signed_x = (2.0 * y - 1.0)[:, None] * scaled_x
    n, p = signed_x.shape

    # Variables are (b_plus[p], b_minus[p], margin), with b=b_plus-b_minus.
    # max margin subject to every signed margin >= margin and ||b||_1 <= 1.
    complete_c = np.r_[np.zeros(2 * p), -1.0]
    complete_a_ub = np.vstack([
        np.column_stack([-signed_x, signed_x, np.ones(n)]),
        np.r_[np.ones(2 * p), 0.0][None, :],
    ])
    complete_b_ub = np.r_[np.zeros(n), 1.0]
    complete_bounds = [(0.0, None)] * (2 * p) + [(0.0, None)]
    complete = linprog(
        complete_c, A_ub=complete_a_ub, b_ub=complete_b_ub,
        bounds=complete_bounds, method="highs",
    )

    if not complete.success:
        return {
            "schema": "logistic_separation_diagnostic_v1",
            "separation_status": "indeterminate",
            "complete_lp_status": int(complete.status),
            "complete_lp_message": str(complete.message),
            "quasi_lp_status": None,
            "quasi_lp_message": None,
            "complete_margin": None,
            "quasi_score": None,
            "design_rank": rank,
            "parameter_dimension": dimension,
            "design_rank_deficient": rank < dimension,
            "tolerance": float(tolerance),
            "column_scaling": safe_scale.tolist(),
        }

    complete_margin = float(complete.x[-1])
    if complete_margin > tolerance:
        return {
            "schema": "logistic_separation_diagnostic_v1",
            "separation_status": "complete_separation",
            "complete_lp_status": int(complete.status),
            "complete_lp_message": str(complete.message),
            "quasi_lp_status": None,
            "quasi_lp_message": None,
            "complete_margin": complete_margin,
            "quasi_score": None,
            "design_rank": rank,
            "parameter_dimension": dimension,
            "design_rank_deficient": rank < dimension,
            "tolerance": float(tolerance),
            "column_scaling": safe_scale.tolist(),
        }

    # max sum_i signed_margin_i subject to all signed_margin_i >= 0
    # and ||b||_1 <= 1. Since complete separation was ruled out above, a
    # positive optimum means at least one strict and at least one zero margin.
    quasi_c = np.r_[-np.sum(signed_x, axis=0), np.sum(signed_x, axis=0)]
    quasi_a_ub = np.vstack([
        np.column_stack([-signed_x, signed_x]),
        np.ones((1, 2 * p)),
    ])
    quasi_b_ub = np.r_[np.zeros(n), 1.0]
    quasi = linprog(
        quasi_c, A_ub=quasi_a_ub, b_ub=quasi_b_ub,
        bounds=[(0.0, None)] * (2 * p), method="highs",
    )
    if not quasi.success:
        status = "indeterminate"
        quasi_score = None
    else:
        quasi_score = float(-quasi.fun)
        status = (
            "quasi_complete_separation"
            if quasi_score > tolerance else "overlap"
        )

    return {
        "schema": "logistic_separation_diagnostic_v1",
        "separation_status": status,
        "complete_lp_status": int(complete.status),
        "complete_lp_message": str(complete.message),
        "quasi_lp_status": int(quasi.status),
        "quasi_lp_message": str(quasi.message),
        "complete_margin": complete_margin,
        "quasi_score": quasi_score,
        "design_rank": rank,
        "parameter_dimension": dimension,
        "design_rank_deficient": rank < dimension,
        "tolerance": float(tolerance),
        "column_scaling": safe_scale.tolist(),
    }
