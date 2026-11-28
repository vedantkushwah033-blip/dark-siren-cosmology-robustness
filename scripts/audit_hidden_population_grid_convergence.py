"""Paired grid and quadrature convergence audit for the continuous hidden model."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from dark_siren_robustness.hidden_population import posterior_grid, posterior_marginals
from scripts.run_continuous_hidden_population_recovery import (
    SIGMA_U, TRUTH, generate_dataset, LISTED_Z, LISTED_W, Z_BOUNDS,
)


def _fit(u, h_points, f_points, mu_points, order):
    posterior = posterior_grid(
        u, SIGMA_U,
        np.linspace(50, 90, h_points),
        np.linspace(0, 1, f_points),
        np.linspace(0.03, 0.095, mu_points),
        [TRUTH["sigma_h"]],
        listed_redshifts=LISTED_Z, listed_weights=LISTED_W,
        z_bounds=Z_BOUNDS, quadrature_order=order,
    )
    return posterior_marginals(posterior)


def run(n_events=30, seed=20261129, coarse=(11, 7, 7), fine=(21, 13, 13),
        quadrature_orders=(24, 48)):
    if n_events < 1:
        raise ValueError("n_events must be positive")
    if any(a < 5 for a in coarse) or any(f <= c for c, f in zip(coarse, fine)):
        raise ValueError("fine grids must be strictly finer and each grid dimension >= 5")
    if len(quadrature_orders) != 2 or quadrature_orders[0] < 8 or quadrature_orders[1] <= quadrature_orders[0]:
        raise ValueError("quadrature_orders must be two increasing values >= 8")
    u, _, _ = generate_dataset(n_events, np.random.default_rng(seed))
    coarse_fit = _fit(u, *coarse, quadrature_orders[0])
    fine_same_order = _fit(u, *fine, quadrature_orders[0])
    fine_higher_order = _fit(u, *fine, quadrature_orders[1])
    fields = ("median", "lower_68", "upper_68", "lower_95", "upper_95")
    comparisons = {}
    for parameter in ("H0", "fraction", "mu_h", "sigma_h"):
        comparisons[parameter] = {}
        for label, other in (
            ("grid_refinement", fine_same_order),
            ("quadrature_refinement", fine_higher_order),
        ):
            differences = {
                field: float(other[parameter][field]-coarse_fit[parameter][field])
                for field in fields
            }
            comparisons[parameter][label] = {
                "absolute_max_quantile_change": max(abs(v) for v in differences.values()),
                "signed_changes": differences,
            }
    return {
        "seed": seed, "n_events": n_events,
        "coarse_grid_points": list(coarse), "fine_grid_points": list(fine),
        "quadrature_orders": list(quadrature_orders),
        "coarse_fit": {
            p: {f: float(coarse_fit[p][f]) for f in fields}
            for p in ("H0", "fraction", "mu_h", "sigma_h")
        },
        "comparisons": comparisons,
        "interpretation": (
            "Paired numerical sensitivity for one fixed synthetic dataset only. "
            "This is not a population-level convergence guarantee."
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-events", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20261129)
    parser.add_argument("--coarse", default="11,7,7")
    parser.add_argument("--fine", default="21,13,13")
    parser.add_argument("--quadrature-orders", default="24,48")
    parser.add_argument("--output", default="results/continuous-hidden-population-recovery/grid_convergence.json")
    args = parser.parse_args()
    result = run(
        args.n_events, args.seed,
        tuple(int(x) for x in args.coarse.split(",")),
        tuple(int(x) for x in args.fine.split(",")),
        tuple(int(x) for x in args.quadrature_orders.split(",")),
    )
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result["comparisons"], indent=2))


if __name__ == "__main__":
    main()
