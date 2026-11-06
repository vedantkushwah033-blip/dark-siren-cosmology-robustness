"""Command-line entry point for the independent scenario lifecycle audit."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from dark_siren_robustness.scenario_lifecycle_audit import audit_scenario_lifecycle


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit_scenario_lifecycle(args.artifact_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
