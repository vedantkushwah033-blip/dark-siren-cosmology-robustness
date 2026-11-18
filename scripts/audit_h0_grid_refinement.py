"""Compare H0 posterior discretization on identical synthetic observations."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def summarize(root: Path, points: int) -> dict:
    conv=json.loads((root/"quadrature_convergence.json").read_text())
    if conv.get("status")!="CONVERGED":
        raise RuntimeError(f"{points}-point run did not converge in coefficient quadrature")
    orders=conv["orders_evaluated"]
    record=next(x for x in conv["posterior_by_order"] if x["order"]==orders[-1])
    frame=pd.read_csv(root/"h0_posterior_density_by_order.csv")
    h0=frame["H0"].to_numpy(float)
    density=frame[f"density_order_{orders[-1]}"].to_numpy(float)
    if len(h0) != points or not np.isfinite(h0).all() or not np.isfinite(density).all():
        raise AssertionError(f"{points}-point H0 grid has invalid length or non-finite values")
    if np.any(np.diff(h0) <= 0) or np.any(density < 0):
        raise AssertionError(f"{points}-point H0 grid must be increasing with non-negative density")
    config=json.loads((root/"config.json").read_text())
    prior_bounds=config["inference"]["uniform_h0_prior"]
    if not (np.isclose(h0[0],prior_bounds[0]) and np.isclose(h0[-1],prior_bounds[1])):
        raise AssertionError(f"{points}-point H0 grid does not span its configured prior bounds")
    # Grid resolution is the only allowed configuration difference.
    normalized_config=json.loads(json.dumps(config))
    normalized_config["inference"].pop("h0_grid_points",None)
    config_hash=hashlib.sha256(json.dumps(
        normalized_config,sort_keys=True,separators=(",",":")
    ).encode()).hexdigest()
    return {"points":points,"step":280.0/(points-1),"root":str(root),
        "model_config_hash":config_hash,"prior_bounds":prior_bounds,
        "observable_hashes":{name:sha(root/name) for name in [
            "observable/analysis_catalogue.csv","observable/gw_events_observed.csv",
            "observable/validation_observed.csv"]},
        "quadrature_order":orders[-1],"quadrature_orders":orders,
        "median":record["h0_median"],"q16":record["h0_q16"],
        "q84":record["h0_q84"],"q025":record["h0_q025"],
        "q975":record["h0_q975"],"h0":h0.tolist(),"density":density.tolist()}


def compare_observable_inputs(left_root: Path, right_root: Path, *, atol: float = 1e-12) -> dict:
    """Compare the same generated observations, allowing only machine-roundoff noise."""
    key_columns = {
        "observable/analysis_catalogue.csv": "object_id",
        "observable/gw_events_observed.csv": "event_id",
        "observable/validation_observed.csv": "validation_id",
    }
    diagnostics = {}
    for relative, key in key_columns.items():
        left = pd.read_csv(left_root / relative)
        right = pd.read_csv(right_root / relative)
        if list(left.columns) != list(right.columns) or len(left) != len(right):
            raise AssertionError(f"observable schema/row-count differs for {relative}")
        if key not in left.columns:
            raise AssertionError(f"observable key {key} missing from {relative}")
        left = left.sort_values(key).reset_index(drop=True)
        right = right.sort_values(key).reset_index(drop=True)
        if not np.array_equal(left[key].to_numpy(), right[key].to_numpy()):
            raise AssertionError(f"observable identifiers differ for {relative}")
        max_abs_difference = 0.0
        for column in left.columns:
            if column == key:
                continue
            a = left[column].to_numpy()
            b = right[column].to_numpy()
            if pd.api.types.is_float_dtype(left[column].dtype):
                difference = np.abs(a.astype(float) - b.astype(float))
                max_abs_difference = max(max_abs_difference, float(difference.max(initial=0.0)))
                if not np.allclose(a, b, rtol=0.0, atol=atol, equal_nan=True):
                    raise AssertionError(
                        f"observable float values differ beyond atol={atol:g} for {relative}:{column}"
                    )
            elif not np.array_equal(a, b):
                raise AssertionError(f"observable values differ for {relative}:{column}")
        diagnostics[relative] = {
            "rows": len(left),
            "key_column": key,
            "max_abs_float_difference": max_abs_difference,
            "absolute_tolerance": atol,
            "non_float_columns_exact": True,
        }
    return diagnostics


def compare(records: list[dict], out: Path) -> dict:
    reference=max(records,key=lambda x:x["points"])
    observable_diagnostics = {}
    reference_root = Path(reference["root"])
    for r in records:
        observable_diagnostics[str(r["points"])] = compare_observable_inputs(
            reference_root, Path(r["root"]), atol=1e-12
        )
    if any(r["model_config_hash"]!=reference["model_config_hash"] for r in records):
        raise AssertionError("inference or generator configuration differs across H0 grids")
    if any(r["prior_bounds"]!=reference["prior_bounds"] for r in records):
        raise AssertionError("H0 prior bounds differ across grids")
    common=np.linspace(20.0,300.0,11201)
    ref=np.interp(common,reference["h0"],reference["density"])
    ref/=np.trapezoid(ref,common)
    ref_cdf=np.concatenate([[0.0],np.cumsum((ref[1:]+ref[:-1])*np.diff(common)/2)])
    ref_cdf/=ref_cdf[-1]
    ref_quantiles={q:float(np.interp(q,ref_cdf,common)) for q in (.025,.16,.5,.84,.975)}
    rows=[]
    for r in sorted(records,key=lambda x:x["points"]):
        den=np.interp(common,r["h0"],r["density"])
        den/=np.trapezoid(den,common)
        cdf=np.concatenate([[0.0],np.cumsum((den[1:]+den[:-1])*np.diff(common)/2)])
        cdf/=cdf[-1]
        quantiles={q:float(np.interp(q,cdf,common)) for q in (.025,.16,.5,.84,.975)}
        rows.append({"points":r["points"],"step":r["step"],
            "quadrature_order":r["quadrature_order"],
            "max_abs_density_difference_vs_finest":float(np.max(np.abs(den-ref))),
            "total_variation_vs_finest":float(.5*np.trapezoid(np.abs(den-ref),common)),
            "max_abs_quantile_difference_vs_finest":float(max(abs(quantiles[q]-ref_quantiles[q]) for q in quantiles)),
            "median":quantiles[.5],"q16":quantiles[.16],"q84":quantiles[.84],
            "quantiles":{str(k):v for k,v in quantiles.items()}})
    result={"status":"PASS","reference_grid_points":reference["points"],
        "common_comparison_grid_points":len(common),"same_observable_inputs":True,
        "grid_summaries":rows,
        "reference_quantiles":{str(k):v for k,v in ref_quantiles.items()},
        "interpretation":"Numerical H0-grid sensitivity only for one fixed synthetic realization. This does not establish calibration, SBC, coverage, or estimator performance."}
    out.write_text(json.dumps(result,indent=2)+"\n")
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    records=[summarize(a.root/f"grid_{n}",n) for n in (141,281,561)]
    print(json.dumps(compare(records,a.output),indent=2))


if __name__=="__main__": main()
