"""Power, energy, and snapshot figures for stationary SS runs.

The SS solver writes two optional post-processing files:

``stationary_snapshots.npz``
    Mesoscopic fields sampled only at requested physical snapshot times.

``snapshot_power_energy_diagnostics.npz``
    Runtime lookup of PRE power/energy/transport catalogues at the same
    snapshot times.

This module deliberately treats both files as diagnostics.  It never changes
solver state and it gracefully returns no figures for older runs that do not
have the new files yet.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import numpy as np

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.colors import LogNorm, SymLogNorm

from pysnspd.analysis.snapshots import compute_snapshot_joule_power_density
from pysnspd.plotting.style import THESIS_WIDTH_IN, apply_thesis_style

apply_thesis_style()

MEV_J = 1.602176634e-22

def _with_recomputed_joule_power(
    power: Mapping[str, np.ndarray],
    *,
    snapshots: Mapping[str, np.ndarray],
    dataset: Mapping[str, Any],
) -> dict[str, np.ndarray]:
    """Return power maps with positive-definite Joule recomputed from snapshots.

    Older SS runs may have ``joule_snapshot_W_m3`` saved as
    ``j_tot*j_n/sigma_n``.  The plotting pipeline has enough frozen snapshot
    data to repair the diagnostic in memory, so figures can be regenerated
    without rerunning PRE or SS.
    """
    out = dict(power)
    if not out or not snapshots:
        return out

    sigma_n = _dataset_sigma_n(dataset)
    if sigma_n is None:
        return out

    p_ep = _snapshot_array(out, ("P_total_snapshot_W_m3",))
    if p_ep.size == 0:
        return out

    joule = compute_snapshot_joule_power_density(
        snapshots,
        sigma_n_S_m=sigma_n,
        n_snap=int(p_ep.shape[0]),
        n_nodes=int(p_ep.shape[1]),
    )
    if joule is not None and joule.shape == p_ep.shape:
        out["joule_snapshot_W_m3"] = np.asarray(joule, dtype=float)
        out["joule_formula"] = np.asarray(["jn_squared_over_sigma_n"], dtype=object)

    return out


def _dataset_sigma_n(dataset: Mapping[str, Any]) -> float | None:
    for key in ("sigma_n_S_m", "sigma_n"):
        if key in dataset:
            try:
                val = float(np.asarray(dataset[key]).reshape(-1)[0])
            except Exception:
                continue
            if np.isfinite(val) and val > 0.0:
                return val
    return None


def _load_npz_if_exists(path: Path) -> dict[str, np.ndarray]:
    if not path.exists():
        return {}
    with np.load(path, allow_pickle=True) as data:
        return {key: np.asarray(data[key]) for key in data.files}


def _prepare_output(output_path: str | Path) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def _mesh_x_nm(mesh: Any, dataset: Mapping[str, Any]) -> np.ndarray:
    if "x_nm" in dataset:
        return np.asarray(dataset["x_nm"], dtype=float)
    nodes = np.asarray(getattr(mesh, "nodes", getattr(mesh, "sites", [])), dtype=float)
    return nodes[:, 0] * 1.0e9


def _mesh_y_nm(mesh: Any, dataset: Mapping[str, Any]) -> np.ndarray:
    if "y_nm" in dataset:
        return np.asarray(dataset["y_nm"], dtype=float)
    nodes = np.asarray(getattr(mesh, "nodes", getattr(mesh, "sites", [])), dtype=float)
    return nodes[:, 1] * 1.0e9


def _mesh_nodes_m(mesh: Any, dataset: Mapping[str, Any]) -> np.ndarray:
    nodes = np.asarray(getattr(mesh, "nodes", getattr(mesh, "sites", [])), dtype=float)
    if nodes.ndim == 2 and nodes.shape[1] >= 2 and nodes.shape[0] > 0:
        return nodes[:, :2]
    if "x_m" in dataset and "y_m" in dataset:
        return np.column_stack((np.asarray(dataset["x_m"], dtype=float), np.asarray(dataset["y_m"], dtype=float)))
    return np.column_stack((_mesh_x_nm(mesh, dataset) * 1.0e-9, _mesh_y_nm(mesh, dataset) * 1.0e-9))


def _mesh_triangles(mesh: Any, dataset: Mapping[str, Any]) -> np.ndarray:
    return np.asarray(
        dataset.get("triangles", getattr(mesh, "triangles", getattr(mesh, "elements", []))),
        dtype=np.int64,
    )


def _triangulation(mesh: Any, dataset: Mapping[str, Any]) -> mtri.Triangulation:
    x = _mesh_x_nm(mesh, dataset)
    y = _mesh_y_nm(mesh, dataset)
    triangles = _mesh_triangles(mesh, dataset)
    return mtri.Triangulation(x, y, triangles)


def _delta0_mev(dataset: Mapping[str, Any], snapshots: Mapping[str, np.ndarray] | None = None) -> float:
    for source in (dataset, snapshots or {}):
        if "delta0_meV" in source:
            arr = np.asarray(source["delta0_meV"], dtype=float).reshape(-1)
            if arr.size and np.isfinite(arr[-1]) and arr[-1] > 0.0:
                return float(arr[-1])
    return 1.0


def _javg(dataset: Mapping[str, Any], snapshots: Mapping[str, np.ndarray] | None = None) -> float:
    for source in (dataset, snapshots or {}):
        if "javg_A_m2" in source:
            arr = np.asarray(source["javg_A_m2"], dtype=float).reshape(-1)
            if arr.size and np.isfinite(arr[-1]) and abs(arr[-1]) > 0.0:
                return abs(float(arr[-1]))
    return 1.0


def _snapshot_array(
    data: Mapping[str, np.ndarray] | None,
    keys: tuple[str, ...],
    *,
    fallback: np.ndarray | None = None,
    shape_like: np.ndarray | None = None,
) -> np.ndarray:
    if data:
        for key in keys:
            if key in data:
                arr = np.asarray(data[key], dtype=float)
                if arr.ndim == 1:
                    arr = arr[None, :]
                return arr
    if fallback is not None:
        arr = np.asarray(fallback, dtype=float)
        if arr.ndim == 1:
            arr = arr[None, :]
        return arr
    if shape_like is not None and np.asarray(shape_like).size:
        return np.zeros_like(np.asarray(shape_like, dtype=float))
    return np.empty((0, 0), dtype=float)


def _resize_snapshot_field(arr: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    a = np.asarray(arr, dtype=float)
    if a.size == 0:
        return np.zeros(shape, dtype=float)
    if a.shape == shape:
        return a
    return np.resize(a, shape)


def _delta_mev_from_psi_snapshots(snapshots: Mapping[str, np.ndarray] | None) -> np.ndarray:
    if not snapshots:
        return np.empty((0, 0), dtype=float)
    real = snapshots.get("psi_snapshot_real_J")
    imag = snapshots.get("psi_snapshot_imag_J")
    if real is None:
        return np.empty((0, 0), dtype=float)
    r = np.asarray(real, dtype=float)
    i = np.asarray(imag if imag is not None else np.zeros_like(r), dtype=float)
    return np.sqrt(r * r + i * i) / MEV_J


def _snapshot_current_mag(snapshots: Mapping[str, np.ndarray], *, family: str) -> np.ndarray:
    if family == "jtot":
        mag_keys = ("jtot_snapshot_mag_A_m2", "jmag_snapshot_A_m2", "current_density_snapshot_A_m2")
        x_keys = (
            "jtot_snapshot_x_A_m2",
            "current_density_snapshot_x_A_m2",
            "node_jtot_x_snapshot_A_m2",
            "jx_snapshot_A_m2",
        )
        y_keys = (
            "jtot_snapshot_y_A_m2",
            "current_density_snapshot_y_A_m2",
            "node_jtot_y_snapshot_A_m2",
            "jy_snapshot_A_m2",
        )
    elif family == "jn":
        mag_keys = ("jn_snapshot_mag_A_m2", "normal_current_density_snapshot_A_m2")
        x_keys = ("jn_snapshot_x_A_m2", "normal_current_density_snapshot_x_A_m2")
        y_keys = ("jn_snapshot_y_A_m2", "normal_current_density_snapshot_y_A_m2")
    else:
        raise ValueError(f"unknown current family: {family}")

    mag = _snapshot_array(snapshots, mag_keys)
    if mag.size:
        return mag
    x = _snapshot_array(snapshots, x_keys)
    y = _snapshot_array(snapshots, y_keys, shape_like=x)
    if x.size:
        return np.sqrt(x * x + y * y)
    return np.empty((0, 0), dtype=float)


def _snapshot_times_ps(data: Mapping[str, np.ndarray], *, preferred: tuple[str, ...], n: int) -> np.ndarray:
    for key in preferred:
        if key in data:
            arr = np.asarray(data[key], dtype=float).reshape(-1)
            if arr.size:
                if arr.size != n:
                    arr = np.resize(arr, n)
                return arr / 1.0e-12
    return np.arange(int(n), dtype=float)


def _representative_snapshot_indices(n: int, *, max_panels: int) -> np.ndarray:
    n = int(n)
    if n <= 0:
        return np.array([], dtype=int)
    if n <= max_panels:
        return np.arange(n, dtype=int)
    return np.unique(np.linspace(0, n - 1, int(max_panels)).round().astype(int))


def _snapshot_grid_figure(
    *,
    nrows: int,
    ncols: int,
    title: str,
    left: float = 0.055,
    right: float = 0.910,
    bottom: float = 0.080,
    top: float = 0.900,
    wspace: float = 0.10,
    hspace: float = 0.24,
):
    height = max(3.0, 1.80 * max(nrows, 1) + 1.1)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(THESIS_WIDTH_IN, height),
        squeeze=False,
        constrained_layout=False,
    )
    fig.subplots_adjust(left=left, right=right, bottom=bottom, top=top, wspace=wspace, hspace=hspace)
    fig.suptitle(title, y=0.975)
    return fig, axes


def _format_map_axis(
    ax,
    *,
    show_xlabel: bool = False,
    show_xlabel_top: bool = False,
    show_ylabel: bool = True,
) -> None:
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)
    if show_xlabel_top:
        ax.set_xlabel("x [nm]")
        ax.xaxis.set_label_position("top")
        ax.xaxis.tick_top()
        ax.tick_params(axis="x", labeltop=True, top=True, labelbottom=False, bottom=False, pad=2)
    elif show_xlabel:
        ax.set_xlabel("x [nm]")
        ax.tick_params(axis="x", labeltop=False, top=False, labelbottom=True, bottom=True)
    else:
        ax.tick_params(axis="x", labeltop=False, top=False, labelbottom=False, bottom=True)
    ax.set_xticklabels([])
    if show_ylabel:
        ax.set_ylabel("y [nm]")
        ax.tick_params(axis="y", labelleft=True, left=True)
    else:
        ax.tick_params(axis="y", labelleft=False, left=True)
    ax.set_yticklabels([])


def _annotate_snapshot_time(ax, t_ps: float) -> None:
    ax.text(
        0.97,
        0.94,
        f"t={t_ps:.3g} ps",
        transform=ax.transAxes,
        ha="right",
        va="top",
        color="white",
        fontsize=8.5,
        bbox={
            "boxstyle": "round,pad=0.18",
            "facecolor": "#c00000",
            "edgecolor": "white",
            "linewidth": 0.9,
            "alpha": 0.98,
        },
        zorder=5,
    )


def _add_row_colorbar(fig, axes_row, mappable, label: str, *, width_fraction: float = 0.80) -> None:
    row_axes = list(np.ravel(np.asarray(axes_row, dtype=object)))
    left = min(ax.get_position().x0 for ax in row_axes)
    right = max(ax.get_position().x1 for ax in row_axes)
    bottom = min(ax.get_position().y0 for ax in row_axes)
    row_width = right - left
    cb_width = row_width * float(width_fraction)
    cb_left = left + 0.5 * (row_width - cb_width)
    cb_height = 0.018
    cb_bottom = max(bottom - 0.022, 0.020)
    cax = fig.add_axes([cb_left, cb_bottom, cb_width, cb_height])
    cb = fig.colorbar(mappable, cax=cax, orientation="horizontal")
    cb.set_label(label)


def _wrap_values_to_range(values: np.ndarray, vmin: float | None, vmax: float | None) -> np.ndarray:
    """Wrap diagnostic values into a finite color range instead of clipping."""
    z = np.asarray(values, dtype=float)
    if vmin is None or vmax is None or not np.isfinite(vmin) or not np.isfinite(vmax):
        return z
    width = float(vmax) - float(vmin)
    if not np.isfinite(width) or width <= 0.0:
        return z
    out = np.array(z, copy=True)
    finite = np.isfinite(out)
    out[finite] = ((out[finite] - float(vmin)) % width) + float(vmin)
    return out


def _plot_values_for_mode(values: np.ndarray, *, mode: str, norm: Any) -> np.ndarray:
    z = np.asarray(values, dtype=float)
    if mode == "positive_log" and isinstance(norm, LogNorm):
        floor = float(norm.vmin) if norm.vmin is not None else 1.0e-300
        return np.where(z > 0.0, z, floor)
    return z


def _node_color_limits(
    values: np.ndarray,
    *,
    symmetric: bool,
    vmin: float | None,
    vmax: float | None,
):
    z = np.asarray(values, dtype=float)
    finite = z[np.isfinite(z)]
    if finite.size == 0:
        finite = np.array([0.0])
    norm = None
    if symmetric:
        vm = float(np.nanpercentile(np.abs(finite), 99.5))
        vm = max(vm, 1.0e-30)
        return norm, -vm, vm
    if vmax is None:
        vmax = float(np.nanpercentile(finite, 99.5))
    if vmin is None:
        vmin = float(np.nanpercentile(finite, 0.5))
    if not np.isfinite(vmax) or not np.isfinite(vmin) or vmax <= vmin:
        vmax = float(vmin) + 1.0
    return norm, vmin, vmax


def _norm_for_mode(values: np.ndarray, mode: str):
    z = np.asarray(values, dtype=float)
    finite = z[np.isfinite(z)]
    if finite.size == 0:
        finite = np.array([0.0])
    if mode == "signed":
        vmax = float(np.nanpercentile(np.abs(finite), 99.2))
        vmax = max(vmax, 1.0e-30)
        linthresh = max(vmax * 1.0e-6, 1.0e8)
        return SymLogNorm(linthresh=linthresh, vmin=-vmax, vmax=vmax), None, None
    if mode == "positive_log":
        pos = finite[finite > 0.0]
        if pos.size == 0:
            return None, 0.0, 1.0
        vmin = float(np.nanpercentile(pos, 1.0))
        vmax = float(np.nanpercentile(pos, 99.2))
        if not np.isfinite(vmin) or vmin <= 0.0:
            vmin = max(float(np.nanmin(pos)), 1.0e-300)
        if not np.isfinite(vmax) or vmax <= vmin:
            vmax = vmin * 10.0
        return LogNorm(vmin=vmin, vmax=vmax), None, None
    return _node_color_limits(z, symmetric=False, vmin=None, vmax=None)


def _snapshot_diffusion_power_density(
    mesh: Any | None,
    *,
    snapshots: Mapping[str, np.ndarray] | None,
    power: Mapping[str, np.ndarray] | None,
    dataset: Mapping[str, Any],
    shape_like: np.ndarray,
) -> np.ndarray:
    """Return a diagnostic electron-diffusion power-density map.

    Preference order:
      1. use a saved diffusion map if a future SS writer provides one;
      2. reconstruct ``div(kappa_s grad T_e)`` from saved snapshot topology;
      3. return NaNs with the requested shape so the row renders as unavailable.

    This helper is intentionally plotting-side only.  It does not modify the SS
    solver, the raw diagnostics NPZ, or the PRE power catalogue.
    """
    shape_arr = np.asarray(shape_like, dtype=float)
    if shape_arr.ndim != 2 or shape_arr.size == 0:
        return np.empty((0, 0), dtype=float)
    shape = shape_arr.shape

    saved = _snapshot_array(
        power,
        (
            "P_diff_snapshot_W_m3",
            "P_diffusion_snapshot_W_m3",
            "diffusion_snapshot_W_m3",
            "thermal_diffusion_snapshot_W_m3",
            "electron_diffusion_snapshot_W_m3",
        ),
    )
    if saved.size:
        return _resize_snapshot_field(saved, shape)

    if not snapshots:
        return np.full(shape, np.nan, dtype=float)

    Te = _snapshot_array(snapshots, ("Te_snapshot_K",), shape_like=shape_arr)
    if Te.size == 0 or not np.any(np.isfinite(Te)):
        return np.full(shape, np.nan, dtype=float)
    Te = _resize_snapshot_field(Te, shape)

    kappa = _snapshot_array(power, ("kappa_s_snapshot_W_m_K",), shape_like=shape_arr)
    if kappa.size == 0 or not np.any(np.isfinite(kappa)):
        return np.full(shape, np.nan, dtype=float)
    kappa = _resize_snapshot_field(kappa, shape)

    edge_i, edge_j = _edge_indices_from_snapshots_or_mesh(snapshots, mesh, dataset)
    if edge_i.size == 0 or edge_i.size != edge_j.size:
        return np.full(shape, np.nan, dtype=float)

    n_snap, n_nodes = shape
    if int(np.max(edge_i, initial=-1)) >= n_nodes or int(np.max(edge_j, initial=-1)) >= n_nodes:
        return np.full(shape, np.nan, dtype=float)

    edge_length = _edge_length_from_snapshots_or_mesh(snapshots, mesh, dataset, edge_i=edge_i, edge_j=edge_j)
    dual_length = _dual_length_from_snapshots(snapshots, edge_length)
    node_area = _node_control_areas_m2(mesh, dataset, n_nodes=n_nodes)
    if node_area.size != n_nodes or not np.all(np.isfinite(node_area)):
        return np.full(shape, np.nan, dtype=float)

    edge_length = np.maximum(np.asarray(edge_length, dtype=float), 1.0e-300)
    dual_length = np.maximum(np.asarray(dual_length, dtype=float), 0.0)
    node_area = np.maximum(np.asarray(node_area, dtype=float), 1.0e-300)

    out = np.zeros(shape, dtype=float)
    geom = dual_length / edge_length
    for s in range(n_snap):
        k_edge = 0.5 * (kappa[s, edge_i] + kappa[s, edge_j])
        dT = Te[s, edge_j] - Te[s, edge_i]
        flux_into_i = k_edge * geom * dT
        acc = np.zeros(n_nodes, dtype=float)
        np.add.at(acc, edge_i, flux_into_i)
        np.add.at(acc, edge_j, -flux_into_i)
        out[s] = acc / node_area
    return out


def _edge_indices_from_snapshots_or_mesh(
    snapshots: Mapping[str, np.ndarray] | None,
    mesh: Any | None,
    dataset: Mapping[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    if snapshots and "edge_i" in snapshots and "edge_j" in snapshots:
        return (
            np.asarray(snapshots["edge_i"], dtype=np.int64).reshape(-1),
            np.asarray(snapshots["edge_j"], dtype=np.int64).reshape(-1),
        )

    triangles = _mesh_triangles(mesh, dataset) if mesh is not None else _mesh_triangles(_NullMesh(), dataset)
    if triangles.size == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    edge_set: set[tuple[int, int]] = set()
    for tri in np.asarray(triangles, dtype=np.int64):
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            ia, ib = int(a), int(b)
            edge_set.add((ia, ib) if ia < ib else (ib, ia))
    if not edge_set:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    edges = np.array(sorted(edge_set), dtype=np.int64)
    return edges[:, 0], edges[:, 1]


class _NullMesh:
    nodes = np.empty((0, 2), dtype=float)
    triangles = np.empty((0, 3), dtype=np.int64)


def _edge_length_from_snapshots_or_mesh(
    snapshots: Mapping[str, np.ndarray] | None,
    mesh: Any | None,
    dataset: Mapping[str, Any],
    *,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
) -> np.ndarray:
    if snapshots and "edge_length_m" in snapshots:
        arr = np.asarray(snapshots["edge_length_m"], dtype=float).reshape(-1)
        if arr.size == edge_i.size:
            return arr
    if mesh is None:
        return np.ones(edge_i.size, dtype=float)
    nodes = _mesh_nodes_m(mesh, dataset)
    if nodes.shape[0] <= max(int(np.max(edge_i, initial=0)), int(np.max(edge_j, initial=0))):
        return np.ones(edge_i.size, dtype=float)
    return np.linalg.norm(nodes[edge_j] - nodes[edge_i], axis=1)


def _dual_length_from_snapshots(snapshots: Mapping[str, np.ndarray] | None, edge_length: np.ndarray) -> np.ndarray:
    if snapshots and "dual_face_length_m" in snapshots:
        arr = np.asarray(snapshots["dual_face_length_m"], dtype=float).reshape(-1)
        if arr.size == np.asarray(edge_length).size:
            return arr
    return np.asarray(edge_length, dtype=float)


def _node_control_areas_m2(mesh: Any | None, dataset: Mapping[str, Any], *, n_nodes: int) -> np.ndarray:
    for key in (
        "node_control_area_m2",
        "control_area_m2",
        "node_area_m2",
        "dual_area_m2",
        "site_areas_m2",
    ):
        if key in dataset:
            arr = np.asarray(dataset[key], dtype=float).reshape(-1)
            if arr.size == n_nodes:
                return arr

    if mesh is None:
        return np.full(n_nodes, np.nan, dtype=float)
    nodes = _mesh_nodes_m(mesh, dataset)
    triangles = _mesh_triangles(mesh, dataset)
    if nodes.shape[0] < n_nodes or triangles.size == 0:
        return np.full(n_nodes, np.nan, dtype=float)

    tri = np.asarray(triangles, dtype=np.int64)
    p0 = nodes[tri[:, 0]]
    p1 = nodes[tri[:, 1]]
    p2 = nodes[tri[:, 2]]
    area = 0.5 * np.abs(
        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
    )
    out = np.zeros(n_nodes, dtype=float)
    for local in range(3):
        valid = (tri[:, local] >= 0) & (tri[:, local] < n_nodes)
        np.add.at(out, tri[valid, local], area[valid] / 3.0)
    return out


def _plot_snapshot_metric(ax, t_ps: np.ndarray, values: np.ndarray, *, reducer: str, label: str) -> None:
    arr = np.asarray(values, dtype=float)
    if arr.ndim != 2 or arr.size == 0:
        return
    finite_by_row = np.any(np.isfinite(arr), axis=1)
    if not np.any(finite_by_row):
        return
    with np.errstate(all="ignore"):
        if reducer == "max":
            y = np.nanmax(arr, axis=1)
        elif reducer == "max_abs":
            y = np.nanmax(np.abs(arr), axis=1)
        elif reducer == "min":
            y = np.nanmin(arr, axis=1)
        elif reducer == "mean":
            y = np.nanmean(arr, axis=1)
        elif reducer == "p99":
            y = np.nanpercentile(arr, 99.0, axis=1)
        else:
            raise ValueError(f"unknown reducer: {reducer}")
    n = min(np.asarray(t_ps).size, y.size)
    if n:
        ax.plot(np.asarray(t_ps)[:n], y[:n], marker="o", linewidth=1.2, label=label)


def _plot_series_if_any(ax, x, y, *, label: str) -> None:
    x_arr = np.asarray(x, dtype=float).reshape(-1)
    y_arr = np.asarray(y if y is not None else [], dtype=float).reshape(-1)
    if x_arr.size == 0 or y_arr.size == 0:
        return
    n = min(x_arr.size, y_arr.size)
    if n <= 0:
        return
    mask = np.isfinite(x_arr[:n]) & np.isfinite(y_arr[:n])
    if not np.any(mask):
        return
    ax.plot(x_arr[:n][mask], y_arr[:n][mask], linewidth=1.2, label=label)


def _legend_if_labels(ax, *, frameon: bool = False, loc: str = "best") -> None:
    handles, labels = ax.get_legend_handles_labels()
    filtered = [
        (handle, label)
        for handle, label in zip(handles, labels)
        if label and not label.startswith("_")
    ]
    if filtered:
        handles, labels = zip(*filtered)
        ax.legend(handles, labels, frameon=frameon, loc=loc)


def _plot_binned_profile(ax, x_nm: np.ndarray, z: np.ndarray, *, label: str, n_bins: int = 80) -> None:
    x = np.asarray(x_nm, dtype=float)
    y = np.asarray(z, dtype=float)
    if x.size == 0 or y.size == 0:
        return
    if y.size != x.size:
        y = np.resize(y, x.size)
    bins = np.linspace(float(np.nanmin(x)), float(np.nanmax(x)), int(n_bins) + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])
    which = np.clip(np.digitize(x, bins) - 1, 0, centers.size - 1)
    prof = np.full(centers.size, np.nan, dtype=float)
    for k in range(centers.size):
        mask = which == k
        if np.any(mask):
            prof[k] = float(np.nanmean(y[mask]))
    ax.plot(centers, prof, marker="o", markersize=2.2, linewidth=1.2, label=label)
