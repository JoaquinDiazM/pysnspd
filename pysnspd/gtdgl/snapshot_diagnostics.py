"""Snapshot-level SS diagnostics for frozen-thermal gTDGL runs.

This module is deliberately diagnostic-only: it does not feed any power or
energy quantity back into the SS solver.  It evaluates the PRE power/energy
catalogue only at the already requested SS snapshot times, so the runtime cost
is negligible compared with the relaxation itself.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

import numpy as np

MEV_J = 1.602176634e-22


def save_ss_snapshot_bundle_npz(history: Mapping[str, Any], output_path: str | Path) -> Path:
    """Save only snapshot and static-FV topology arrays from a relaxation history.

    ``relaxation_history.npz`` remains the complete compact history.  This file is
    a convenience product for plotting/inspection pipelines that should not have
    to know which arrays are final scalars and which arrays are full snapshot
    maps.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    keep: dict[str, np.ndarray] = {}
    for key, value in history.items():
        if _is_snapshot_or_static_topology_key(str(key)):
            keep[str(key)] = np.asarray(value)
    if "snapshot_t_s" not in keep and "delta_snapshot_t_s" in history:
        keep["snapshot_t_s"] = np.asarray(history["delta_snapshot_t_s"])
    np.savez_compressed(output, **keep)
    return output


def write_ss_snapshot_power_diagnostics(
    *,
    history: Mapping[str, Any],
    state: Any,
    power_table_npz: str | Path,
    output_path: str | Path,
    sigma_n_S_m: float | None = None,
) -> Path:
    """Evaluate PRE energy/power tables on all stored SS snapshots and save NPZ.

    The SS solver is still frozen-thermal.  Therefore these maps are diagnostics:
    they answer what local ``P_ep``, ``u_e``, ``C_e``, ``kappa_s`` and phonon
    escape terms would be for the instantaneous ``(|Delta|, q, Te, Tph)`` state.
    """
    diagnostics = compute_ss_snapshot_power_diagnostics(
        history=history,
        state=state,
        power_table_npz=power_table_npz,
        sigma_n_S_m=sigma_n_S_m,
    )
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **diagnostics)
    return output


def compute_ss_snapshot_power_diagnostics(
    *,
    history: Mapping[str, Any],
    state: Any,
    power_table_npz: str | Path,
    sigma_n_S_m: float | None = None,
) -> dict[str, np.ndarray]:
    """Return node maps of PRE power/energy quantities at SS snapshot times."""
    cat = _load_power_table_npz(power_table_npz)
    delta_meV = np.asarray(history["delta_snapshot_meV"], dtype=float)
    if delta_meV.ndim != 2:
        raise ValueError("history['delta_snapshot_meV'] must have shape (n_snap, n_nodes).")
    n_snap, n_nodes = delta_meV.shape
    delta_J = np.maximum(delta_meV * MEV_J, 0.0)

    snapshot_t_s = np.asarray(
        history.get("snapshot_t_s", history.get("delta_snapshot_t_s", np.arange(n_snap, dtype=float))),
        dtype=float,
    )
    if snapshot_t_s.shape != (n_snap,):
        snapshot_t_s = np.asarray(snapshot_t_s[:n_snap], dtype=float)

    Te_node = _broadcast_node_temperature(getattr(state, "Te_K"), n_snap=n_snap, n_nodes=n_nodes, name="Te_K")
    Tph_node = _broadcast_node_temperature(getattr(state, "Tph_K"), n_snap=n_snap, n_nodes=n_nodes, name="Tph_K")

    q_node = _snapshot_node_q_abs(history, n_snap=n_snap, n_nodes=n_nodes)

    iTe = _nearest_indices(cat["Te_values_K"], Te_node)
    iTph = _nearest_indices(cat["Tph_values_K"], Tph_node)
    iDelta = _nearest_indices(cat["delta_values_J"], delta_J)
    iQ = _nearest_indices(cat["q_values_m_inv"], q_node)

    P_S = cat["P_S_W_m3"][iTe, iTph, iDelta, iQ]
    P_R = cat["P_R_W_m3"][iTe, iTph, iDelta, iQ]
    P_total = cat["P_total_W_m3"][iTe, iTph, iDelta, iQ]
    u_e = cat["u_e_J_m3"][iTe, iDelta, iQ]
    C_e = cat["C_e_J_m3_K"][iTe, iDelta, iQ]

    if cat["kappa_s_W_m_K"].size:
        kappa_s = cat["kappa_s_W_m_K"][iTe, iDelta]
    else:
        kappa_s = np.full_like(P_total, np.nan, dtype=float)

    iTph_only = _nearest_indices(cat["Tph_values_K"], Tph_node)
    u_ph = _table_1d_or_nan(cat.get("u_ph_J_m3"), iTph_only, shape=(n_snap, n_nodes))
    C_ph = _table_1d_or_nan(cat.get("C_ph_J_m3_K"), iTph_only, shape=(n_snap, n_nodes))
    P_esc = _table_1d_or_nan(cat.get("P_esc_W_m3"), iTph_only, shape=(n_snap, n_nodes))

    out = {
        "snapshot_t_s": snapshot_t_s,
        "snapshot_t_ps": snapshot_t_s / 1.0e-12,
        "Te_snapshot_K": Te_node,
        "Tph_snapshot_K": Tph_node,
        "delta_snapshot_J": delta_J,
        "delta_snapshot_meV": delta_meV,
        "q_abs_snapshot_m_inv": q_node,
        "power_table_iTe": iTe.astype(np.int64),
        "power_table_iTph": iTph.astype(np.int64),
        "power_table_iDelta": iDelta.astype(np.int64),
        "power_table_iQ": iQ.astype(np.int64),
        "P_S_snapshot_W_m3": np.asarray(P_S, dtype=float),
        "P_R_snapshot_W_m3": np.asarray(P_R, dtype=float),
        "P_total_snapshot_W_m3": np.asarray(P_total, dtype=float),
        "u_e_snapshot_J_m3": np.asarray(u_e, dtype=float),
        "C_e_snapshot_J_m3_K": np.asarray(C_e, dtype=float),
        "kappa_s_snapshot_W_m_K": np.asarray(kappa_s, dtype=float),
        "u_ph_snapshot_J_m3": np.asarray(u_ph, dtype=float),
        "C_ph_snapshot_J_m3_K": np.asarray(C_ph, dtype=float),
        "P_esc_snapshot_W_m3": np.asarray(P_esc, dtype=float),
        "metadata_json": np.asarray([_diagnostic_metadata_json(cat, power_table_npz)], dtype=object),
    }

    joule = _snapshot_joule_power_density(history, sigma_n_S_m=sigma_n_S_m, n_snap=n_snap, n_nodes=n_nodes)
    if joule is not None:
        out["joule_snapshot_W_m3"] = joule
    return out


def _load_power_table_npz(path: str | Path) -> dict[str, np.ndarray]:
    with np.load(Path(path), allow_pickle=True) as data:
        required = [
            "Te_values_K",
            "Tph_values_K",
            "delta_values_J",
            "q_values_m_inv",
            "P_S_W_m3",
            "P_R_W_m3",
            "P_total_W_m3",
            "u_e_J_m3",
            "C_e_J_m3_K",
        ]
        missing = [key for key in required if key not in data.files]
        if missing:
            raise ValueError(f"power_table_catalog.npz is missing keys: {missing}")
        out = {key: np.asarray(data[key]) for key in required}
        out["kappa_s_W_m_K"] = np.asarray(data.get("kappa_s_W_m_K", np.array([], dtype=float)), dtype=float)
        out["u_ph_J_m3"] = np.asarray(data.get("u_ph_J_m3", np.array([], dtype=float)), dtype=float)
        out["C_ph_J_m3_K"] = np.asarray(data.get("C_ph_J_m3_K", np.array([], dtype=float)), dtype=float)
        out["P_esc_W_m3"] = np.asarray(data.get("P_esc_W_m3", np.array([], dtype=float)), dtype=float)
        if "metadata" in data.files:
            try:
                out["metadata"] = data["metadata"].item()
            except Exception:
                out["metadata"] = {}
        else:
            out["metadata"] = {}
        return out


def _broadcast_node_temperature(values: Any, *, n_snap: int, n_nodes: int, name: str) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.shape == (n_nodes,):
        return np.broadcast_to(arr[None, :], (n_snap, n_nodes)).copy()
    if arr.shape == (n_snap, n_nodes):
        return arr.copy()
    if arr.size == 1:
        return np.full((n_snap, n_nodes), float(arr.ravel()[0]), dtype=float)
    raise ValueError(f"state.{name} must be scalar, (n_nodes,), or (n_snap,n_nodes); got {arr.shape}.")


def _snapshot_node_q_abs(history: Mapping[str, Any], *, n_snap: int, n_nodes: int) -> np.ndarray:
    if "edge_Q_snapshot_m_inv" not in history:
        return np.zeros((n_snap, n_nodes), dtype=float)
    edge_q = np.abs(np.asarray(history["edge_Q_snapshot_m_inv"], dtype=float))
    if edge_q.ndim != 2 or edge_q.shape[0] != n_snap:
        raise ValueError("edge_Q_snapshot_m_inv must have shape (n_snap,n_edges).")
    return _edge_to_node_snapshots(edge_q, history=history, n_nodes=n_nodes)


def _snapshot_joule_power_density(
    history: Mapping[str, Any],
    *,
    sigma_n_S_m: float | None,
    n_snap: int,
    n_nodes: int,
) -> np.ndarray | None:
    if sigma_n_S_m is None or not np.isfinite(float(sigma_n_S_m)) or float(sigma_n_S_m) <= 0.0:
        return None
    if "edge_jn_snapshot_A_m2" not in history or "edge_jtot_snapshot_A_m2" not in history:
        return None
    edge_jn = np.asarray(history["edge_jn_snapshot_A_m2"], dtype=float)
    edge_jtot = np.asarray(history["edge_jtot_snapshot_A_m2"], dtype=float)
    if edge_jn.shape != edge_jtot.shape or edge_jn.ndim != 2 or edge_jn.shape[0] != n_snap:
        return None
    # Since j_n = sigma_n E on each edge projection, j_tot.E = j_tot*j_n/sigma_n.
    edge_joule = edge_jtot * edge_jn / float(sigma_n_S_m)
    return _edge_to_node_snapshots(edge_joule, history=history, n_nodes=n_nodes)


def _edge_to_node_snapshots(edge_values: np.ndarray, *, history: Mapping[str, Any], n_nodes: int) -> np.ndarray:
    edge_i = np.asarray(history.get("edge_i"), dtype=np.int64)
    edge_j = np.asarray(history.get("edge_j"), dtype=np.int64)
    if edge_i.ndim != 1 or edge_j.ndim != 1 or edge_i.size != edge_j.size:
        raise ValueError("history must include 1D edge_i and edge_j arrays for node projection.")
    if edge_values.shape[1] != edge_i.size:
        raise ValueError("edge snapshot array has inconsistent n_edges.")
    edge_length = np.asarray(history.get("edge_length_m", np.ones(edge_i.size)), dtype=float)
    dual = np.asarray(history.get("dual_face_length_m", np.ones(edge_i.size)), dtype=float)
    weights = np.maximum(dual / np.maximum(edge_length, 1.0e-300), 1.0e-300)
    out = np.zeros((edge_values.shape[0], n_nodes), dtype=float)
    wsum = np.zeros(n_nodes, dtype=float)
    np.add.at(wsum, edge_i, weights)
    np.add.at(wsum, edge_j, weights)
    wsum = np.maximum(wsum, 1.0e-300)
    for k in range(edge_values.shape[0]):
        acc = np.zeros(n_nodes, dtype=float)
        vals = weights * edge_values[k]
        np.add.at(acc, edge_i, vals)
        np.add.at(acc, edge_j, vals)
        out[k] = acc / wsum
    return out


def _nearest_indices(axis: np.ndarray, values: np.ndarray) -> np.ndarray:
    ax = np.asarray(axis, dtype=float)
    if ax.ndim != 1 or ax.size == 0:
        raise ValueError("lookup axis must be a non-empty 1D array.")
    vals = np.asarray(values, dtype=float)
    pos = np.searchsorted(ax, vals, side="left")
    pos = np.clip(pos, 0, ax.size - 1)
    left = np.clip(pos - 1, 0, ax.size - 1)
    choose_left = np.abs(vals - ax[left]) <= np.abs(vals - ax[pos])
    return np.where(choose_left, left, pos).astype(np.int64)


def _table_1d_or_nan(values: np.ndarray | None, indices: np.ndarray, *, shape: tuple[int, int]) -> np.ndarray:
    if values is None:
        return np.full(shape, np.nan, dtype=float)
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return np.full(shape, np.nan, dtype=float)
    idx = np.clip(indices, 0, arr.size - 1)
    return np.asarray(arr[idx], dtype=float)


def _is_snapshot_or_static_topology_key(key: str) -> bool:
    if "snapshot" in key:
        return True
    return key in {
        "delta0_meV",
        "javg_A_m2",
        "qref_m_inv",
        "edge_i",
        "edge_j",
        "edge_length_m",
        "edge_unit_x",
        "edge_unit_y",
        "dual_face_length_m",
    }


def _diagnostic_metadata_json(cat: Mapping[str, Any], power_table_npz: str | Path) -> str:
    import json

    metadata = cat.get("metadata", {})
    if not isinstance(metadata, dict):
        metadata = {}
    payload = {
        "source": "diagnostic_only_ss_snapshot_power_energy_maps",
        "power_table_npz": str(power_table_npz),
        "interpolation": "nearest_neighbor_on_PRE_power_table_axes",
        "thermal_policy": "does_not_couple_back_to_solver",
        "power_table_backend": metadata.get("backend"),
    }
    return json.dumps(payload, sort_keys=True)


__all__ = [
    "save_ss_snapshot_bundle_npz",
    "compute_ss_snapshot_power_diagnostics",
    "write_ss_snapshot_power_diagnostics",
]
