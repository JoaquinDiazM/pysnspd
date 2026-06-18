"""Plotting placeholder pipeline."""

import argparse

from pysnspd.config import load_config, validate_config
from pysnspd.plotting.figures import (
    plot_catalog_diagnostics,
    plot_mesh,
    plot_photon_transient,
    plot_stationary_state,
)


def main():
    """Generate placeholder plots for a run."""
    parser = argparse.ArgumentParser(description="Plot a pySNSPD run template.")
    parser.add_argument("--config", required=True, help="Path to project YAML configuration.")
    parser.add_argument("--run-name", required=True, help="Canonical run name.")
    args = parser.parse_args()

    config = load_config(args.config)
    validate_config(config)
    plot_mesh(config, args.run_name)
    plot_catalog_diagnostics(config, args.run_name)
    plot_stationary_state(config, args.run_name)
    plot_photon_transient(config, args.run_name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
