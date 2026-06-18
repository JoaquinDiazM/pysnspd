"""PRE-run placeholder pipeline."""

import argparse

from pysnspd.config import load_config, validate_config
from pysnspd.runs.prerun import run_prerun


def main():
    """Launch the PRE-run scaffold."""
    parser = argparse.ArgumentParser(description="Run PRE-run template.")
    parser.add_argument("--config", required=True, help="Path to project YAML configuration.")
    parser.add_argument("--run-name", required=True, help="Canonical run name.")
    args = parser.parse_args()

    config = load_config(args.config)
    validate_config(config)
    run_prerun(config, args.run_name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
