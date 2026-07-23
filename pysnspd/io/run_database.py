"""Lightweight discovery utilities for pySNSPD raw-run databases.

The Z-series plot pipelines compare multiple runs.  They need a first pass that
can answer a simple question before any physics-level plotting is attempted:
what files, summaries, manifests and NPZ keys are available in
``project.big_data_root/raw``?

This module intentionally reads metadata only.  It does not load full heavy
arrays into memory, except for tiny scalar previews and string metadata fields.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import yaml

from pysnspd.config import load_config, validate_config


STAGE_NAMES: tuple[str, ...] = ("pre", "ss", "photon")
SUMMARY_LIKE_SUFFIXES: tuple[str, ...] = (".yaml", ".yml", ".json")


def load_database_config(config_path: str | Path) -> dict[str, Any]:
    """Load and validate a pySNSPD config for read-only database scanning."""
    return validate_config(load_config(config_path), require_big_data_root_exists=True)


def big_data_root_from_config(config: Mapping[str, Any]) -> Path:
    """Return the configured external data root."""
    cfg = validate_config(config, require_big_data_root_exists=True)
    return Path(cfg["project"]["big_data_root"]).expanduser().resolve()


def raw_root_from_config(config: Mapping[str, Any]) -> Path:
    """Return ``project.big_data_root/raw``."""
    return big_data_root_from_config(config) / "raw"


def discover_raw_run_records(
    config: Mapping[str, Any],
    *,
    run_names: Sequence[str] | None = None,
    run_prefixes: Sequence[str] | None = None,
    stages: Sequence[str] | None = None,
    include_npz_keys: bool = True,
    include_yaml_data: bool = True,
) -> list[dict[str, Any]]:
    """Discover available raw run data under ``big_data_root/raw``.

    Parameters
    ----------
    config:
        Valid pySNSPD config dictionary.
    run_names:
        Optional explicit run names.  If omitted, every folder under ``raw`` is
        scanned, optionally filtered by ``run_prefixes``.
    run_prefixes:
        Optional prefixes used to filter the discovered raw-run directories.
    stages:
        Stage names to inspect.  Defaults to ``("pre", "ss", "photon")``.
    include_npz_keys:
        If true, open every NPZ file and record key, shape and dtype metadata.
    include_yaml_data:
        If true, parse YAML/JSON-like text summaries and manifests.

    Returns
    -------
    list[dict]
        YAML/JSON-serializable run records.
    """
    cfg = validate_config(config, require_big_data_root_exists=True)
    raw_root = raw_root_from_config(cfg)
    selected_stages = _normalize_stages(stages)
    selected_names = _select_run_names(
        raw_root,
        run_names=run_names,
        run_prefixes=run_prefixes,
    )

    records: list[dict[str, Any]] = []
    for run_name in selected_names:
        raw_run = raw_root / run_name
        record: dict[str, Any] = {
            "run_name": str(run_name),
            "raw_run": str(raw_run),
            "raw_run_exists": bool(raw_run.exists()),
            "stages": {},
        }
        for stage in selected_stages:
            record["stages"][stage] = _scan_stage_dir(
                raw_run / stage,
                stage_name=stage,
                raw_run=raw_run,
                include_npz_keys=include_npz_keys,
                include_yaml_data=include_yaml_data,
            )
        records.append(record)
    return records


def write_database_inventory(
    records: Sequence[Mapping[str, Any]],
    output_dir: str | Path,
    *,
    basename: str = "Z1_current_sweep_inventory",
) -> dict[str, Path]:
    """Write inventory records as YAML, JSON and compact CSV files."""
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    serializable = _to_jsonable(list(records))
    yaml_path = out_dir / f"{basename}.yaml"
    json_path = out_dir / f"{basename}.json"
    csv_path = out_dir / f"{basename}_files.csv"

    with yaml_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(serializable, f, sort_keys=False, allow_unicode=True, default_flow_style=False)

    with json_path.open("w", encoding="utf-8") as f:
        json.dump(serializable, f, indent=2, ensure_ascii=False)

    write_database_inventory_csv(records, csv_path)

    return {"yaml": yaml_path, "json": json_path, "csv": csv_path}


def write_database_inventory_csv(
    records: Sequence[Mapping[str, Any]],
    output_path: str | Path,
) -> Path:
    """Write one compact row per discovered NPZ/YAML/JSON file."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = inventory_file_rows(records)
    fieldnames = [
        "run_name",
        "stage",
        "kind",
        "relative_path",
        "exists",
        "size_bytes",
        "n_keys",
        "keys",
    ]
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})
    return path


def inventory_file_rows(records: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    """Flatten database records into file-level rows."""
    rows: list[dict[str, Any]] = []
    for record in records:
        run_name = str(record.get("run_name", ""))
        stages = record.get("stages", {})
        if not isinstance(stages, Mapping):
            continue
        for stage_name, stage in stages.items():
            if not isinstance(stage, Mapping):
                continue
            for npz in stage.get("npz_files", []) or []:
                if not isinstance(npz, Mapping):
                    continue
                keys = npz.get("keys", {})
                key_names = list(keys.keys()) if isinstance(keys, Mapping) else []
                rows.append(
                    {
                        "run_name": run_name,
                        "stage": str(stage_name),
                        "kind": "npz",
                        "relative_path": npz.get("relative_path", ""),
                        "exists": True,
                        "size_bytes": npz.get("size_bytes", ""),
                        "n_keys": len(key_names),
                        "keys": ";".join(key_names),
                    }
                )
            for text_file in stage.get("summary_files", []) or []:
                if not isinstance(text_file, Mapping):
                    continue
                rows.append(
                    {
                        "run_name": run_name,
                        "stage": str(stage_name),
                        "kind": str(text_file.get("kind", "summary")),
                        "relative_path": text_file.get("relative_path", ""),
                        "exists": True,
                        "size_bytes": text_file.get("size_bytes", ""),
                        "n_keys": text_file.get("n_top_level_keys", ""),
                        "keys": ";".join(text_file.get("top_level_keys", []) or []),
                    }
                )
    return rows


def summarize_inventory(records: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    """Return compact counts for terminal reporting and manifests."""
    n_runs = len(records)
    n_stage_dirs = 0
    n_npz = 0
    n_summary = 0
    stage_counts: dict[str, dict[str, int]] = {}
    for record in records:
        stages = record.get("stages", {})
        if not isinstance(stages, Mapping):
            continue
        for stage_name, stage in stages.items():
            if not isinstance(stage, Mapping):
                continue
            counts = stage_counts.setdefault(str(stage_name), {"dirs": 0, "npz": 0, "summary": 0})
            if stage.get("exists"):
                n_stage_dirs += 1
                counts["dirs"] += 1
            n_npz_stage = len(stage.get("npz_files", []) or [])
            n_summary_stage = len(stage.get("summary_files", []) or [])
            n_npz += n_npz_stage
            n_summary += n_summary_stage
            counts["npz"] += n_npz_stage
            counts["summary"] += n_summary_stage
    return {
        "n_runs": int(n_runs),
        "n_stage_dirs": int(n_stage_dirs),
        "n_npz_files": int(n_npz),
        "n_summary_files": int(n_summary),
        "stage_counts": stage_counts,
    }


def _scan_stage_dir(
    stage_dir: Path,
    *,
    stage_name: str,
    raw_run: Path,
    include_npz_keys: bool,
    include_yaml_data: bool,
) -> dict[str, Any]:
    stage: dict[str, Any] = {
        "stage": stage_name,
        "path": str(stage_dir),
        "exists": bool(stage_dir.exists()),
        "npz_files": [],
        "summary_files": [],
    }
    if not stage_dir.exists():
        return stage

    npz_paths = sorted(stage_dir.rglob("*.npz"))
    for path in npz_paths:
        stage["npz_files"].append(
            _summarize_npz_file(
                path,
                raw_run=raw_run,
                include_keys=include_npz_keys,
            )
        )

    summary_paths = [
        path
        for path in sorted(stage_dir.rglob("*"))
        if path.is_file() and path.suffix.lower() in SUMMARY_LIKE_SUFFIXES
    ]
    for path in summary_paths:
        stage["summary_files"].append(
            _summarize_text_summary(
                path,
                raw_run=raw_run,
                include_data=include_yaml_data,
            )
        )
    return stage


def _summarize_npz_file(path: Path, *, raw_run: Path, include_keys: bool) -> dict[str, Any]:
    out: dict[str, Any] = {
        "path": str(path),
        "relative_path": _safe_relative(path, raw_run),
        "size_bytes": _file_size(path),
        "keys": {},
    }
    if not include_keys:
        return out
    try:
        with np.load(path, allow_pickle=False) as data:
            for key in data.files:
                arr = np.asarray(data[key])
                out["keys"][str(key)] = _summarize_array(arr)
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
    return out


def _summarize_array(arr: np.ndarray) -> dict[str, Any]:
    summary: dict[str, Any] = {
        "shape": list(arr.shape),
        "dtype": str(arr.dtype),
        "size": int(arr.size),
    }
    preview = _array_preview(arr)
    if preview is not None:
        summary["preview"] = preview
    return summary


def _array_preview(arr: np.ndarray) -> Any | None:
    if arr.size == 0:
        return None
    if arr.shape == ():
        return _to_jsonable(arr.item())
    if arr.dtype.kind in {"U", "S"} and arr.size <= 4:
        return _to_jsonable(arr.reshape(-1)[:4].tolist())
    if arr.dtype.kind in {"i", "u", "f", "b"} and arr.size <= 6:
        return _to_jsonable(arr.reshape(-1).tolist())
    return None


def _summarize_text_summary(path: Path, *, raw_run: Path, include_data: bool) -> dict[str, Any]:
    out: dict[str, Any] = {
        "path": str(path),
        "relative_path": _safe_relative(path, raw_run),
        "kind": _summary_kind(path),
        "size_bytes": _file_size(path),
    }
    if not include_data:
        return out
    try:
        data = _load_text_data(path)
        out["data"] = _to_jsonable(data)
        if isinstance(data, Mapping):
            keys = [str(key) for key in data.keys()]
            out["top_level_keys"] = keys
            out["n_top_level_keys"] = len(keys)
    except Exception as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
    return out


def _load_text_data(path: Path) -> Any:
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() == ".json":
        return json.loads(text)
    return yaml.safe_load(text)


def _summary_kind(path: Path) -> str:
    name = path.name.lower()
    if "manifest" in name:
        return "manifest"
    if "summary" in name:
        return "summary"
    if "config" in name:
        return "config"
    return path.suffix.lower().lstrip(".") or "text"


def _select_run_names(
    raw_root: Path,
    *,
    run_names: Sequence[str] | None,
    run_prefixes: Sequence[str] | None,
) -> list[str]:
    if run_names:
        names = [str(name) for name in run_names]
    else:
        names = sorted(path.name for path in raw_root.iterdir() if path.is_dir()) if raw_root.exists() else []

    prefixes = [str(prefix) for prefix in (run_prefixes or []) if str(prefix)]
    if prefixes:
        names = [name for name in names if any(name.startswith(prefix) for prefix in prefixes)]
    return names


def _normalize_stages(stages: Sequence[str] | None) -> tuple[str, ...]:
    if stages is None:
        return STAGE_NAMES
    cleaned = tuple(str(stage).strip().lower() for stage in stages if str(stage).strip())
    if not cleaned or "all" in cleaned:
        return STAGE_NAMES
    invalid = [stage for stage in cleaned if stage not in STAGE_NAMES]
    if invalid:
        raise ValueError(f"Invalid stage names: {invalid}. Expected one of {STAGE_NAMES}.")
    return cleaned


def _file_size(path: Path) -> int:
    try:
        return int(path.stat().st_size)
    except OSError:
        return 0


def _safe_relative(path: Path, parent: Path) -> str:
    try:
        return str(path.relative_to(parent))
    except ValueError:
        return str(path)


def _to_jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _to_jsonable(val) for key, val in value.items()}
    if isinstance(value, (list, tuple)):
        return [_to_jsonable(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return _to_jsonable(value.tolist())
    if isinstance(value, np.generic):
        return _to_jsonable(value.item())
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


__all__ = [
    "STAGE_NAMES",
    "big_data_root_from_config",
    "discover_raw_run_records",
    "inventory_file_rows",
    "load_database_config",
    "plots_root_from_config",
    "raw_root_from_config",
    "summarize_inventory",
    "write_database_inventory",
    "write_database_inventory_csv",
]
