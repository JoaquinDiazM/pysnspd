"""Photon/phonon-bubble helpers for pySNSPD transient runs.

The first implementation follows the reduced two-temperature interpretation of
Vodolazov's phonon-bubble initial condition: a sudden local increase of phonon
energy density, followed by the existing Te/Tph runtime coupling.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

EV_J = 1.602176634e-19


@dataclass(frozen=True)
class PhotonBubbleParams:
    energy_eV: float = 0.0
    time_s: float = 0.0
    x_m: float | None = None
    y_m: float = 0.0
    sigma_m: float = 10.0e-9
    enabled: bool = True

    def as_dict(self) -> dict[str, float | bool | None]:
        return asdict(self)


def inject_phonon_bubble(
    *,
    mesh: Any,
    Tph_K: np.ndarray,
    power_table_npz: str | Path,
    thickness_m: float,
    params: PhotonBubbleParams,
) -> tuple[np.ndarray, dict[str, float | bool | None]]:
    """Inject photon energy into the phonon energy density.

    The update conserves the requested photon energy over the 2D mesh control
    volumes times film thickness:

    ``sum_i A_i d [u_ph(T_i+) - u_ph(T_i-)] = E_gamma``.

    For ``energy_eV == 0`` the state is returned unchanged, but metadata still
    records that the event was processed.
    """

    old_Tph = np.asarray(Tph_K, dtype=float).reshape(-1)
    if not bool(params.enabled):
        return old_Tph.copy(), {
            "enabled": False,
            "applied": False,
            "energy_eV": float(params.energy_eV),
            "reason": "photon disabled",
        }

    energy_J = float(params.energy_eV) * EV_J
    if not np.isfinite(energy_J) or energy_J < 0.0:
        raise ValueError("Photon energy must be non-negative and finite.")

    nodes = np.asarray(mesh.nodes, dtype=float)[:, :2]
    if old_Tph.size != nodes.shape[0]:
        raise ValueError(f"Tph_K has length {old_Tph.size}, expected {nodes.shape[0]}.")

    x0 = _resolve_x0(nodes, params.x_m)
    y0 = float(params.y_m)
    sigma = float(params.sigma_m)
    if not np.isfinite(sigma) or sigma <= 0.0:
        raise ValueError("Photon sigma must be positive and finite.")

    weights = np.exp(-0.5 * (((nodes[:, 0] - x0) / sigma) ** 2 + ((nodes[:, 1] - y0) / sigma) ** 2))
    areas = node_control_areas_m2(mesh)
    norm = float(np.sum(areas * weights))
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError("Photon bubble normalization failed; check mesh and sigma.")

    if energy_J == 0.0:
        return old_Tph.copy(), {
            "enabled": True,
            "applied": True,
            "energy_eV": 0.0,
            "energy_J": 0.0,
            "x_m": float(x0),
            "y_m": float(y0),
            "sigma_m": float(sigma),
            "normalization_area_m2": float(norm),
            "max_delta_u_ph_J_m3": 0.0,
            "energy_reconstructed_J": 0.0,
            "reason": "zero-energy photon bubble: state unchanged",
        }

    T_axis, u_axis = load_phonon_energy_table(power_table_npz)
    u_old = np.interp(old_Tph, T_axis, u_axis)
    delta_u = energy_J * weights / max(float(thickness_m) * norm, 1.0e-300)
    u_new = u_old + delta_u
    T_new = np.interp(u_new, u_axis, T_axis, left=float(T_axis[0]), right=float(T_axis[-1]))

    reconstructed = float(np.sum(areas * float(thickness_m) * (np.interp(T_new, T_axis, u_axis) - u_old)))

    return T_new, {
        "enabled": True,
        "applied": True,
        "energy_eV": float(params.energy_eV),
        "energy_J": float(energy_J),
        "x_m": float(x0),
        "y_m": float(y0),
        "sigma_m": float(sigma),
        "normalization_area_m2": float(norm),
        "max_delta_u_ph_J_m3": float(np.nanmax(delta_u)),
        "energy_reconstructed_J": reconstructed,
        "relative_energy_error": float((reconstructed - energy_J) / max(energy_J, 1.0e-300)),
        "reason": "phonon bubble applied",
    }


def load_phonon_energy_table(power_table_npz: str | Path) -> tuple[np.ndarray, np.ndarray]:
    path = Path(power_table_npz)
    if not path.exists():
        raise FileNotFoundError(f"Missing power table for phonon bubble: {path}")
    with np.load(path, allow_pickle=True) as data:
        if "Tph_values_K" not in data.files or "u_ph_J_m3" not in data.files:
            raise ValueError("power_table_catalog.npz must contain Tph_values_K and u_ph_J_m3.")
        T = np.asarray(data["Tph_values_K"], dtype=float).reshape(-1)
        u = np.asarray(data["u_ph_J_m3"], dtype=float).reshape(-1)

    if T.size != u.size or T.size < 2:
        raise ValueError("Invalid phonon energy table dimensions.")
    order = np.argsort(T)
    T = T[order]
    u = u[order]

    # Inverse interpolation needs monotone u.  Small numerical plateaus are OK.
    keep = np.concatenate(([True], np.diff(u) > 0.0))
    if np.count_nonzero(keep) >= 2:
        T = T[keep]
        u = u[keep]
    return T, u


def node_control_areas_m2(mesh: Any) -> np.ndarray:
    nodes = np.asarray(mesh.nodes, dtype=float)[:, :2]
    n = nodes.shape[0]
    tris = None
    if hasattr(mesh, "triangles"):
        tris = np.asarray(mesh.triangles, dtype=np.int64)
    elif hasattr(mesh, "elements"):
        tris = np.asarray(mesh.elements, dtype=np.int64)

    areas = np.zeros(n, dtype=float)
    if tris is not None and tris.size:
        p0 = nodes[tris[:, 0]]
        p1 = nodes[tris[:, 1]]
        p2 = nodes[tris[:, 2]]
        tri_area = 0.5 * np.abs(
            (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
            - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0])
        )
        share = tri_area / 3.0
        np.add.at(areas, tris[:, 0], share)
        np.add.at(areas, tris[:, 1], share)
        np.add.at(areas, tris[:, 2], share)

    if not np.any(areas > 0.0):
        total_area = float(getattr(mesh, "length_m", np.ptp(nodes[:, 0])) * getattr(mesh, "width_m", np.ptp(nodes[:, 1])))
        areas[:] = total_area / max(n, 1)
    else:
        positive = areas[areas > 0.0]
        fill = float(np.nanmedian(positive)) if positive.size else float(np.nanmean(areas))
        areas[areas <= 0.0] = fill

    return areas


def _resolve_x0(nodes_m: np.ndarray, x_m: float | None) -> float:
    if x_m is not None and np.isfinite(float(x_m)):
        return float(x_m)
    x = np.asarray(nodes_m[:, 0], dtype=float)
    return 0.5 * (float(np.nanmin(x)) + float(np.nanmax(x)))


__all__ = [
    "EV_J",
    "PhotonBubbleParams",
    "inject_phonon_bubble",
    "load_phonon_energy_table",
    "node_control_areas_m2",
]
