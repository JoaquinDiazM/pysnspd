"""Project-configuration placeholder pipeline."""

import argparse

from pysnspd.config import load_config, validate_config, summarize_config
from pysnspd.io.manager import initialize_project_storage


def main():
    """Validate configuration and initialize the external big-data folder."""
    parser = argparse.ArgumentParser(description="Configure a pySNSPD project.")
    parser.add_argument("--config", required=True, help="Path to project YAML configuration.")
    args = parser.parse_args()

    config = load_config(args.config)
    validate_config(config)
    summarize_config(config)
    initialize_project_storage(config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
