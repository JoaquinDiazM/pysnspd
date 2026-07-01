"""Data analysis utilities for stationary SS gTDGL runs.

This module intentionally contains no matplotlib calls.  It reads the raw SS
``.npz`` files, extracts physical arrays, builds masks/profiles, and returns
plain dictionaries that plotting functions can consume.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import yaml

from pysnspd.config import load_config, validate_config
from pysnspd.io.manager import create_run_layout
from pysnspd.mesh.delaunay import load_mesh_npz
from pysnspd.mesh.edges import load_edges_npz

MEV_J = 1.602176634e-22


@dataclass(frozen=True)
class SSRunData:
    """Loaded raw data for one stationary SS run."""

    run_name: str
    pre_run_name: str | None
    raw_ss: Path
    figures_dir: Path
    mesh: Any
    edge_data: Any
    state: dict[str, np.ndarray]
    history: dict[str, np.ndarray]
    summary: dict[str, Any]


def load_ss_run(
    *,
    config_path: str | Path,
    run_name: str,
    pre_run_name: str | None = None,
) -> SSRunData:
    """Load mesh, edge table, stationary state, relaxation history and summary.

    Parameters
    ----------
    config_path:
        YAML project config.
    run_name:
        SS run to analyze.
    pre_run_name:
        PRE run that provided ``mesh.npz`` and ``edges.npz``.  If omitted,
        the function first reads it from ``ss_summary.yaml`` and otherwise
        falls back to ``run_name``.
    """

    cfg = validate_config(load_config(config_path))
    ss_layout = create_run_layout(cfg, run_name)
    raw_ss = Path(ss_layout["raw_ss"])
    summary = _load_yaml(raw_ss / "ss_summary.yaml")
    resolved_pre = pre_run_name or str(summary.get("pre_run_name") or run_name)
    pre_layout = create_run_layout(cfg, resolved_pre)
    raw_pre = Path(pre_layout["raw_pre"])

    state_path = Path(summary.get("outputs", {}).get("stationary_state_npz", raw_ss / "stationary_state.npz"))
    history_path = Path(summary.get("outputs", {}).get("relaxation_history_npz", raw_ss / "relaxation_history.npz"))

    mesh = load_mesh_npz(raw_pre / "mesh.npz")
    edge_data = load_edges_npz(raw_pre / "edges.npz")
    state = _load_npz_dict(state_path)
    history = _load_npz_dict(history_path)
    figures_dir = Path(ss_layout["plots_figures"])
    figures_dir.mkdir(parents=True, exist_ok=True)

    return SSRunData(
        run_name=run_name,
        pre_run_name=resolved_pre,
        raw_ss=raw_ss,
        figures_dir=figures_dir,
        mesh=mesh,
        edge_data=edge_data,
        state=state,
        history=history,
        summary=summary,
    )


def build_ss_plot_dataset(run: SSRunData) -> dict[str, Any]:
    """Build physical plotting arrays from raw SS state/history files."""

    mesh = run.mesh
    nodes = np.asarray(mesh.nodes, dtype=float)
    x_m = nodes[:, 0]
    y_m = nodes[:, 1]
    state = run.state
    history = run.history
    summary = run.summary
    solver_summary = _as_mapping(summary.get("solver", {}))

    psi_J = _complex_from_state(state, "psi_real_J", "psi_imag_J")
    delta_meV = np.abs(psi_J) / MEV_J
    delta0_meV = _scalar_from_sources(
        history,
        solver_summary,
        "delta0_meV",
        default=_metadata_scalar(state, "delta0_meV", default=np.nan),
    )
    delta_over_delta0 = delta_meV / delta0_meV if np.isfinite(delta0_meV) and delta0_meV > 0.0 else delta_meV

    phi_V = _array_or_zeros(state, "phi_V", x_m.size)
    phi_mV = 1.0e3 * phi_V

    jtot_x = _array_or_zeros(state, "node_jtot_x_A_m2", x_m.size)
    jtot_y = _array_or_zeros(state, "node_jtot_y_A_m2", x_m.size)
    js_x = _array_or_zeros(state, "node_js_us_x_A_m2", x_m.size)
    js_y = _array_or_zeros(state, "node_js_us_y_A_m2", x_m.size)
    jn_x = _array_or_zeros(state, "node_jn_x_A_m2", x_m.size)
    jn_y = _array_or_zeros(state, "node_jn_y_A_m2", x_m.size)
    jtot_mag = _magnitude(jtot_x, jtot_y)
    js_mag = _magnitude(js_x, js_y)
    jn_mag = _magnitude(jn_x, jn_y)

    javg = _scalar_from_sources(history, solver_summary, "javg_A_m2", default=np.nan)
    if not np.isfinite(javg) or abs(javg) <= 0.0:
        target_current_A = float(solver_summary.get("target_current_A", np.nan))
        width = float(getattr(mesh, "width_m", np.nan))
        thickness = _infer_thickness_from_summary(summary, default=np.nan)
        if np.isfinite(target_current_A) and np.isfinite(width) and np.isfinite(thickness) and width > 0 and thickness > 0:
            javg = abs(target_current_A / (width * thickness))
    jscale = abs(javg) if np.isfinite(javg) and abs(javg) > 0.0 else max(float(np.nanmax(jtot_mag)), 1.0)

    pairbreaking = _array_or_zeros(state, "node_pairbreaking_ratio", x_m.size)
    div_j = _array_or_zeros(state, "node_div_jtot_A_m3", x_m.size)

    normal_terminal_mask = _history_bool_mask(history, "normal_terminal_node_mask", x_m.size)
    bulk_mask = _bulk_node_mask_from_summary(x_m, solver_summary, fallback=np.logical_not(normal_terminal_mask))

    t_s = _history_time_s(history)
    t_ps = t_s / 1.0e-12
    dt_s = _history_array(history, "dt_s")
    if dt_s.size == 0:
        dt_s = np.diff(np.r_[0.0, t_s]) if t_s.size else np.array([], dtype=float)
    dt_fs = dt_s / 1.0e-15

    dataset: dict[str, Any] = {
        "run_name": run.run_name,
        "pre_run_name": run.pre_run_name,
        "x_m": x_m,
        "y_m": y_m,
        "x_nm": x_m * 1.0e9,
        "y_nm": y_m * 1.0e9,
        "triangles": np.asarray(mesh.triangles, dtype=np.int64),
        "psi_J": psi_J,
        "delta_meV": delta_meV,
        "delta0_meV": float(delta0_meV) if np.isfinite(delta0_meV) else np.nan,
        "delta_over_delta0": delta_over_delta0,
        "phi_V": phi_V,
        "phi_mV": phi_mV,
        "javg_A_m2": float(jscale),
        "jtot_x_A_m2": jtot_x,
        "jtot_y_A_m2": jtot_y,
        "jtot_mag_A_m2": jtot_mag,
        "js_us_x_A_m2": js_x,
        "js_us_y_A_m2": js_y,
        "js_us_mag_A_m2": js_mag,
        "jn_x_A_m2": jn_x,
        "jn_y_A_m2": jn_y,
        "jn_mag_A_m2": jn_mag,
        "pairbreaking_ratio": pairbreaking,
        "div_jtot_A_m3": div_j,
        "normal_terminal_node_mask": normal_terminal_mask,
        "bulk_node_mask": bulk_mask,
        "t_ps": t_ps,
        "dt_fs": dt_fs,
        "eta_R": _resize_to_time(_history_array(history, "eta_R"), t_ps),
        "terminal_voltage_mV": 1.0e3 * _resize_to_time(_history_array(history, "terminal_voltage_V"), t_ps),
        "normal_current_fraction": _normal_current_fraction_history(history, t_ps),
        "pairbreaking_max_history": _resize_to_time(_history_array(history, "pairbreaking_max"), t_ps),
        "adaptive_retries": _resize_to_time(_history_array(history, "adaptive_retries"), t_ps),
        "adaptive_rejected_attempts": _resize_to_time(_history_array(history, "adaptive_rejected_attempts"), t_ps),
        "adaptive_window_mean_d_abs_sq": _resize_to_time(_history_array(history, "adaptive_window_mean_d_abs_sq"), t_ps),
        "dt_attempt_fs": _resize_to_time(_history_array(history, "dt_attempt_s") / 1.0e-15, t_ps),
        "dt_accepted_fs": _resize_to_time(_history_array(history, "dt_accepted_s") / 1.0e-15, t_ps),
        "dt_next_fs": _resize_to_time(_history_array(history, "dt_next_s") / 1.0e-15, t_ps),
        "adaptive_target_dt_fs": _resize_to_time(_history_array(history, "adaptive_target_dt_s") / 1.0e-15, t_ps),
        "summary_scalars": _summary_scalars(summary),
        "npz_keys": summarize_ss_npz_contents(state=state, history=history),
    }

    # Backward compatibility: if the adaptive diagnostics were not saved yet,
    # still provide useful curves from the accepted dt history.
    if dataset["dt_accepted_fs"].size == 0 and dt_fs.size:
        dataset["dt_accepted_fs"] = _resize_to_time(dt_fs, t_ps)
    if dataset["dt_attempt_fs"].size == 0 and dt_fs.size:
        dataset["dt_attempt_fs"] = _resize_to_time(dt_fs, t_ps)
    if dataset["dt_next_fs"].size == 0 and dt_fs.size:
        dataset["dt_next_fs"] = _resize_to_time(dt_fs, t_ps)

    dataset["x_profile_nm"], dataset["profiles"] = compute_x_profiles(dataset)
    return dataset


def compute_x_profiles(dataset: Mapping[str, Any], *, n_bins: int = 51) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    """Compute simple x-binned mean profiles for final SS fields."""

    x_nm = np.asarray(dataset["x_nm"], dtype=float)
    if x_nm.size == 0:
        return np.array([], dtype=float), {}
    bins = np.linspace(float(np.nanmin(x_nm)), float(np.nanmax(x_nm)), int(n_bins) + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])
    idx = np.digitize(x_nm, bins) - 1
    idx = np.clip(idx, 0, centers.size - 1)

    profiles: dict[str, np.ndarray] = {}
    keys = [
        "delta_over_delta0",
        "phi_mV",
        "jtot_mag_A_m2",
        "js_us_mag_A_m2",
        "jn_mag_A_m2",
        "pairbreaking_ratio",
    ]
    scale = float(dataset.get("javg_A_m2", 1.0))
    for key in keys:
        values = np.asarray(dataset.get(key, []), dtype=float)
        if values.size != x_nm.size:
            continue
        out = np.full(centers.shape, np.nan, dtype=float)
        for k in range(centers.size):
            mask = idx == k
            if np.any(mask):
                out[k] = float(np.nanmean(values[mask]))
        if key.endswith("_A_m2") and scale > 0.0:
            profiles[key + "_over_javg"] = out / scale
        else:
            profiles[key] = out
    return centers, profiles


def summarize_ss_npz_contents(
    *,
    state: Mapping[str, np.ndarray] | None = None,
    history: Mapping[str, np.ndarray] | None = None,
) -> dict[str, Any]:
    """Return key/shape/dtype summaries for SS ``.npz`` files."""

    out: dict[str, Any] = {}
    for name, data in (("state", state), ("history", history)):
        if data is None:
            continue
        out[name] = {
            str(key): {"shape": list(np.asarray(value).shape), "dtype": str(np.asarray(value).dtype)}
            for key, value in sorted(data.items())
        }
    return out


def _load_npz_dict(path: str | Path) -> dict[str, np.ndarray]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Required SS NPZ file does not exist: {p}")
    with np.load(p, allow_pickle=True) as data:
        return {key: np.asarray(data[key]) for key in data.files}


def _load_yaml(path: str | Path) -> dict[str, Any]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Required SS summary does not exist: {p}")
    with p.open("r", encoding="utf-8") as f:
        obj = yaml.safe_load(f)
    if not isinstance(obj, dict):
        raise ValueError(f"YAML summary must contain a mapping: {p}")
    return obj


def _as_mapping(obj: Any) -> Mapping[str, Any]:
    return obj if isinstance(obj, Mapping) else {}


def _complex_from_state(state: Mapping[str, np.ndarray], real_key: str, imag_key: str) -> np.ndarray:
    real = np.asarray(state.get(real_key, []), dtype=float)
    imag = np.asarray(state.get(imag_key, np.zeros_like(real)), dtype=float)
    if real.size == 0:
        return np.array([], dtype=np.complex128)
    return real + 1j * imag


def _array_or_zeros(data: Mapping[str, np.ndarray], key: str, n: int) -> np.ndarray:
    arr = np.asarray(data.get(key, np.zeros(int(n), dtype=float)), dtype=float)
    if arr.size != int(n):
        arr = np.resize(arr, int(n))
    return arr


def _history_array(history: Mapping[str, np.ndarray], key: str) -> np.ndarray:
    arr = np.asarray(history.get(key, []), dtype=float)
    return arr.reshape(-1) if arr.size else np.array([], dtype=float)


def _history_time_s(history: Mapping[str, np.ndarray]) -> np.ndarray:
    t = _history_array(history, "t_s")
    if t.size:
        return t
    dt = _history_array(history, "dt_s")
    if dt.size:
        return np.cumsum(dt)
    return np.array([], dtype=float)


def _resize_to_time(values: np.ndarray, t_ps: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float).reshape(-1)
    n = np.asarray(t_ps).size
    if arr.size == 0:
        return np.array([], dtype=float)
    if arr.size == n:
        return arr
    if n == 0:
        return arr
    return np.resize(arr, n)


def _history_bool_mask(history: Mapping[str, np.ndarray], key: str, n: int) -> np.ndarray:
    raw = np.asarray(history.get(key, np.zeros(int(n), dtype=bool)))
    if raw.size != int(n):
        raw = np.resize(raw, int(n))
    return raw.astype(bool)


def _magnitude(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.sqrt(np.asarray(x, dtype=float) ** 2 + np.asarray(y, dtype=float) ** 2)


def _scalar_from_sources(
    history: Mapping[str, np.ndarray],
    summary: Mapping[str, Any],
    key: str,
    *,
    default: float = np.nan,
) -> float:
    if key in summary:
        try:
            return float(summary[key])
        except Exception:
            pass
    arr = np.asarray(history.get(key, []), dtype=float).reshape(-1)
    if arr.size:
        try:
            return float(arr[-1])
        except Exception:
            pass
    return float(default)


def _metadata_scalar(state: Mapping[str, np.ndarray], key: str, *, default: float = np.nan) -> float:
    raw = state.get("metadata_json")
    if raw is None:
        return float(default)
    try:
        text = str(np.asarray(raw).reshape(()).item())
        obj = json.loads(text)
        return float(obj.get(key, default))
    except Exception:
        return float(default)


def _infer_thickness_from_summary(summary: Mapping[str, Any], *, default: float = np.nan) -> float:
    # The SS summary stores most physical scalars under solver, while the
    # manifest/config path is not guaranteed to be present in ad-hoc tests.
    cfg = _as_mapping(summary.get("config", {}))
    mat = _as_mapping(cfg.get("material", {}))
    try:
        return float(mat.get("thickness_m", default))
    except Exception:
        return float(default)


def _normal_current_fraction_history(history: Mapping[str, np.ndarray], t_ps: np.ndarray) -> np.ndarray:
    jn = _history_array(history, "normal_current_max_A_m2")
    jt = _history_array(history, "total_current_max_A_m2")
    jn = _resize_to_time(jn, t_ps)
    jt = _resize_to_time(jt, t_ps)
    if jn.size == 0 or jt.size == 0:
        return np.array([], dtype=float)
    return jn / np.maximum(np.abs(jt), 1.0e-300)


def _bulk_node_mask_from_summary(x_m: np.ndarray, solver_summary: Mapping[str, Any], *, fallback: np.ndarray) -> np.ndarray:
    try:
        stationarity = _as_mapping(solver_summary.get("stationarity", {}))
        exclusion_m = float(stationarity.get("bulk_exclusion_length_m", np.nan))
    except Exception:
        exclusion_m = np.nan
    if not np.isfinite(exclusion_m) or exclusion_m <= 0.0:
        return np.asarray(fallback, dtype=bool)
    xmin = float(np.nanmin(x_m))
    xmax = float(np.nanmax(x_m))
    return (x_m >= xmin + exclusion_m) & (x_m <= xmax - exclusion_m)


def _summary_scalars(summary: Mapping[str, Any]) -> dict[str, Any]:
    solver = _as_mapping(summary.get("solver", {}))
    stationarity = _as_mapping(solver.get("stationarity", {}))
    contact = _as_mapping(solver.get("contact_recovery", {}))
    continuity = _as_mapping(solver.get("continuity", {}))
    keys = [
        "first_magic_ready",
        "accepted_steps",
        "rejected_steps",
        "final_time_ps",
        "target_current_A",
        "terminal_voltage_V",
        "max_pairbreaking_ratio",
        "normal_current_fraction_max",
        "usadel_current_backend",
    ]
    out = {key: solver.get(key) for key in keys if key in solver}
    out["stationarity_passes"] = stationarity.get("passes")
    out["contact_recovery_passes"] = contact.get("passes")
    out["continuity_passes"] = continuity.get("passes")
    return out
