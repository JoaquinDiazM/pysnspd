"""Centralized file and run-directory manager.

No solver should manually build output paths. This module will define the
canonical layout

    big_data_root/raw/<run_name>/
    big_data_root/plots/<run_name>/
    big_data_root/logs/<run_name>/

and will be responsible for manifests, metadata, cache discovery, and safe
read/write policies.
"""


def initialize_project_storage(config):
    """Initialize the external big-data directory for a project."""
    return 0


def create_run_layout(config, run_name):
    """Create raw, plots, and logs folders for a run using one shared run name."""
    return 0


def write_manifest(config, run_name, stage):
    """Write a manifest for a simulation stage."""
    return 0


def read_manifest(config, run_name, stage):
    """Read a manifest for a simulation stage."""
    return 0


def resolve_stage_path(config, run_name, stage):
    """Resolve the canonical raw-data path for a stage of a run."""
    return 0


def resolve_plot_path(config, run_name):
    """Resolve the canonical plot path for a run."""
    return 0
