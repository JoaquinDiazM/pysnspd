"""
File and folder management for pySNSPD runs.

This module owns the external data layout. The repository should remain light;
large raw data, catalogs, logs and figures live under ``project.big_data_root``.
"""

from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping
import os
import platform
import re
import sys


try:
    import yaml
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "pySNSPD requires PyYAML to write .yaml manifests. "
        "Install it with: python -m pip install pyyaml"
    ) from exc


from pysnspd import __version__
from pysnspd.config import validate_config


class StorageError(RuntimeError):
    """Raised when pySNSPD cannot create or write to its data folders."""


_RUN_NAME_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.+-]*$")

_BASE_FOLDERS = ("raw", "plots", "logs", "catalogs", "tmp")
_RUN_STAGES = ("pre", "ss", "photon")


def initialize_project_storage(config: Mapping[str, Any]) -> dict[str, str]:
    """
    Create and validate the base external data folders.

    Parameters
    ----------
    config:
        Valid pySNSPD configuration dictionary.

    Returns
    -------
    dict
        Dictionary with absolute paths to the base folders.
    """
    cfg = validate_config(config, require_big_data_root_exists=False)
    root = Path(cfg["project"]["big_data_root"]).expanduser().resolve()

    try:
        root.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise StorageError(f"Could not create big_data_root: {root}") from exc

    if not root.is_dir():
        raise StorageError(f"big_data_root is not a directory: {root}")

    layout = {"big_data_root": str(root)}

    for name in _BASE_FOLDERS:
        path = root / name
        try:
            path.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise StorageError(f"Could not create storage folder: {path}") from exc
        layout[name] = str(path)

    _assert_writable(root)

    for name in _BASE_FOLDERS:
        _assert_writable(Path(layout[name]))

    return layout


def create_run_layout(
    config: Mapping[str, Any],
    run_name: str | None = None,
) -> dict[str, str]:
    """
    Create the folder structure for one simulation run.

    The same ``run_name`` is used for raw data, plots and logs. Renaming runs
    after producing data is discouraged because it breaks the implicit link
    between raw outputs and figures.

    Parameters
    ----------
    config:
        Valid pySNSPD configuration dictionary.
    run_name:
        Run identifier. If None, ``project.default_run_name`` is used.

    Returns
    -------
    dict
        Dictionary with absolute paths to all run folders.
    """
    cfg = validate_config(config, require_big_data_root_exists=False)

    if run_name is None:
        run_name = str(cfg["project"]["default_run_name"])

    run_name = validate_run_name(run_name)
    base = initialize_project_storage(cfg)

    raw_run = Path(base["raw"]) / run_name
    plots_run = Path(base["plots"]) / run_name
    logs_run = Path(base["logs"]) / run_name

    folders = {
        "raw_run": raw_run,
        "raw_pre": raw_run / "pre",
        "raw_ss": raw_run / "ss",
        "raw_photon": raw_run / "photon",
        "plots_run": plots_run,
        "plots_figures": plots_run / "figures",
        "plots_mesh": plots_run / "mesh",
        "plots_diagnostics": plots_run / "diagnostics",
        "plots_comparisons": plots_run / "comparisons",
        "logs_run": logs_run,
    }

    for path in folders.values():
        try:
            path.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise StorageError(f"Could not create run folder: {path}") from exc
        _assert_writable(path)

    layout = dict(base)
    layout["run_name"] = run_name
    layout.update({key: str(path) for key, path in folders.items()})

    return layout


def write_manifest(
    config: Mapping[str, Any],
    run_name: str | None = None,
    *,
    stage: str = "project",
    extra: Mapping[str, Any] | None = None,
) -> Path:
    """
    Write a YAML manifest for a run or a specific stage.

    Parameters
    ----------
    config:
        Valid pySNSPD configuration dictionary.
    run_name:
        Run identifier. If None, ``project.default_run_name`` is used.
    stage:
        One of ``project``, ``pre``, ``ss``, ``photon`` or ``plots``.
    extra:
        Optional extra metadata to include in the manifest.

    Returns
    -------
    pathlib.Path
        Path to the manifest written.
    """
    cfg = validate_config(config, require_big_data_root_exists=False)
    layout = create_run_layout(cfg, run_name)
    stage = _normalize_stage(stage)

    manifest = build_manifest(cfg, layout, stage=stage, extra=extra)

    manifest_path = _manifest_path_from_layout(layout, stage)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    with manifest_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(
            manifest,
            f,
            sort_keys=False,
            allow_unicode=True,
            default_flow_style=False,
        )

    if stage == "project":
        plots_manifest = Path(layout["plots_run"]) / "manifest.yaml"
        with plots_manifest.open("w", encoding="utf-8") as f:
            yaml.safe_dump(
                manifest,
                f,
                sort_keys=False,
                allow_unicode=True,
                default_flow_style=False,
            )

    return manifest_path


def read_manifest(
    config: Mapping[str, Any],
    run_name: str | None = None,
    *,
    stage: str = "project",
) -> dict[str, Any]:
    """
    Read a previously written YAML manifest.

    Parameters
    ----------
    config:
        Valid pySNSPD configuration dictionary.
    run_name:
        Run identifier. If None, ``project.default_run_name`` is used.
    stage:
        One of ``project``, ``pre``, ``ss``, ``photon`` or ``plots``.

    Returns
    -------
    dict
        Manifest contents.
    """
    cfg = validate_config(config, require_big_data_root_exists=False)
    layout = create_run_layout(cfg, run_name)
    stage = _normalize_stage(stage)

    path = _manifest_path_from_layout(layout, stage)

    if not path.exists():
        raise StorageError(f"Manifest does not exist: {path}")

    with path.open("r", encoding="utf-8") as f:
        manifest = yaml.safe_load(f)

    if manifest is None:
        raise StorageError(f"Manifest is empty: {path}")

    if not isinstance(manifest, dict):
        raise StorageError(f"Manifest must contain a YAML mapping: {path}")

    return manifest


def build_manifest(
    config: Mapping[str, Any],
    layout: Mapping[str, str],
    *,
    stage: str,
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Build the manifest dictionary without writing it to disk.
    """
    cfg = deepcopy(dict(config))
    stage = _normalize_stage(stage)

    manifest: dict[str, Any] = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "package": {
            "name": "pysnspd",
            "version": __version__,
        },
        "python": {
            "version": sys.version,
            "executable": sys.executable,
            "platform": platform.platform(),
        },
        "project": {
            "name": cfg["project"]["name"],
            "big_data_root": cfg["project"]["big_data_root"],
        },
        "run": {
            "name": layout["run_name"],
            "stage": stage,
        },
        "paths": dict(layout),
        "config": cfg,
    }

    if extra is not None:
        manifest["extra"] = deepcopy(dict(extra))

    return manifest


def resolve_stage_path(
    config: Mapping[str, Any],
    run_name: str | None = None,
    *,
    stage: str = "project",
) -> Path:
    """
    Resolve the raw-data path associated with a run stage.
    """
    cfg = validate_config(config, require_big_data_root_exists=False)
    layout = create_run_layout(cfg, run_name)
    stage = _normalize_stage(stage)

    if stage == "project":
        return Path(layout["raw_run"])
    if stage == "pre":
        return Path(layout["raw_pre"])
    if stage == "ss":
        return Path(layout["raw_ss"])
    if stage == "photon":
        return Path(layout["raw_photon"])
    if stage == "plots":
        return Path(layout["plots_run"])

    raise StorageError(f"Unknown stage: {stage}")


def resolve_plot_path(
    config: Mapping[str, Any],
    run_name: str | None = None,
    *,
    subfolder: str | None = None,
) -> Path:
    """
    Resolve the plot path associated with a run.

    Parameters
    ----------
    subfolder:
        Optional subfolder inside the run plot directory. Common values are
        ``figures``, ``mesh``, ``diagnostics`` and ``comparisons``.
    """
    cfg = validate_config(config, require_big_data_root_exists=False)
    layout = create_run_layout(cfg, run_name)
    base = Path(layout["plots_run"])

    if subfolder is None:
        return base

    safe = validate_run_name(subfolder)
    path = base / safe
    path.mkdir(parents=True, exist_ok=True)
    _assert_writable(path)

    return path


def validate_run_name(run_name: str) -> str:
    """
    Validate a run name.

    Allowed characters are letters, numbers, underscores, hyphens, dots and
    plus signs. Slashes and spaces are intentionally disallowed.
    """
    if not isinstance(run_name, str):
        raise StorageError("run_name must be a string.")

    cleaned = run_name.strip()

    if cleaned == "":
        raise StorageError("run_name must not be empty.")

    if not _RUN_NAME_PATTERN.match(cleaned):
        raise StorageError(
            "Invalid run_name. Use only letters, numbers, underscores, "
            "hyphens, dots and plus signs. The name must start with a letter "
            f"or number. Got: {run_name!r}"
        )

    return cleaned


def _manifest_path_from_layout(layout: Mapping[str, str], stage: str) -> Path:
    stage = _normalize_stage(stage)

    if stage == "project":
        return Path(layout["raw_run"]) / "manifest.yaml"
    if stage == "pre":
        return Path(layout["raw_pre"]) / "manifest.yaml"
    if stage == "ss":
        return Path(layout["raw_ss"]) / "manifest.yaml"
    if stage == "photon":
        return Path(layout["raw_photon"]) / "manifest.yaml"
    if stage == "plots":
        return Path(layout["plots_run"]) / "manifest.yaml"

    raise StorageError(f"Unknown manifest stage: {stage}")


def _normalize_stage(stage: str) -> str:
    if not isinstance(stage, str):
        raise StorageError("stage must be a string.")

    normalized = stage.strip().lower()

    aliases = {
        "project": "project",
        "run": "project",
        "pre": "pre",
        "prerun": "pre",
        "pre-run": "pre",
        "ss": "ss",
        "stationary": "ss",
        "ss-run": "ss",
        "photon": "photon",
        "photonrun": "photon",
        "photon-run": "photon",
        "plot": "plots",
        "plots": "plots",
        "plotting": "plots",
    }

    if normalized not in aliases:
        raise StorageError(
            "stage must be one of: project, pre, ss, photon, plots. "
            f"Got: {stage}"
        )

    return aliases[normalized]


def _assert_writable(path: Path) -> None:
    try:
        path.mkdir(parents=True, exist_ok=True)
        test_file = path / ".pysnspd_write_test"
        test_file.write_text("ok\n", encoding="utf-8")
        test_file.unlink()
    except OSError as exc:
        raise StorageError(f"Path is not writable: {path}") from exc