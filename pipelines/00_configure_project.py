"""
Configure a pySNSPD project.

This pipeline validates a YAML configuration file, creates the external
big-data folder layout, creates the default run layout and writes a project
manifest. It performs no physics calculation.
"""

from __future__ import annotations

import argparse

from pysnspd.config import load_config, validate_config, summarize_config
from pysnspd.io.manager import (
    initialize_project_storage,
    create_run_layout,
    write_manifest,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate pySNSPD configuration and initialize storage."
    )
    parser.add_argument(
        "--config",
        required=True,
        help="Path to the YAML project configuration.",
    )
    parser.add_argument(
        "--run-name",
        default=None,
        help=(
            "Optional run name. If omitted, project.default_run_name "
            "from the config is used."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    cfg = load_config(args.config)
    cfg = validate_config(cfg)

    print(summarize_config(cfg))
    print()

    base_layout = initialize_project_storage(cfg)
    run_layout = create_run_layout(cfg, args.run_name)

    manifest_path = write_manifest(
        cfg,
        args.run_name,
        stage="project",
        extra={
            "pipeline": "00_configure_project.py",
            "purpose": "OE1 project configuration and storage validation",
        },
    )

    print("Base storage layout")
    for key, value in base_layout.items():
        print(f"  {key}: {value}")

    print()
    print("Run storage layout")
    for key, value in run_layout.items():
        print(f"  {key}: {value}")

    print()
    print(f"Project manifest written to: {manifest_path}")
    print("Status: OK")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())