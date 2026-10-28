"""Write separated observable and truth artifacts for the synthetic selection toy."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from dark_siren_robustness.generative_selection import generate_synthetic_population

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261028)
    parser.add_argument("--n-parent", type=int, default=12000)
    parser.add_argument("--n-detected", type=int, default=250)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/synthetic_selection"))
    args = parser.parse_args()
    output = generate_synthetic_population(seed=args.seed, n_parent=args.n_parent, n_detected=args.n_detected)
    observable_dir = args.output_dir / "observable"
    truth_dir = args.output_dir / "truth"
    observable_dir.mkdir(parents=True, exist_ok=True)
    truth_dir.mkdir(parents=True, exist_ok=True)
    output["catalogue_observed"].to_csv(observable_dir / "catalogue_observed.csv", index=False)
    output["gw_events_observed"].to_csv(observable_dir / "gw_events_observed.csv", index=False)
    output["parent_truth"].to_csv(truth_dir / "parent_truth.csv", index=False)
    output["gw_event_truth"].to_csv(truth_dir / "gw_event_truth.csv", index=False)
    (args.output_dir / "config.json").write_text(json.dumps(output["config"], indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": output["config"]["status"], "n_parent": len(output["parent_truth"]), "n_catalogue_observed": len(output["catalogue_observed"]), "n_detected_events": len(output["gw_events_observed"]), "output_dir": str(args.output_dir), "truth_is_separate": True}, indent=2))

if __name__ == "__main__":
    main()
