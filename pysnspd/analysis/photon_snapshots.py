"""Memory-bounded photon snapshot diagnostics for E3 plotting.

The transient already stores every physical state needed here.  This module
reconstructs runtime-consistent power, energy and heat-capacity fields without
advancing the solver.  All snapshots contribute to the reported color limits,
while node maps are retained only for the requested plotting times.
"""

from __future__ import annotations

from typing import Any, Callable, Mapping, Sequence

import numpy as np
from scipy import sparse

from pysnspd.thermal.evolution import PowerTableRuntimeInterpolator

MEV_J = 1.602176634e-22


def compute_photon_snapshot_plot_diagnostics(
    *,
    snapshots: Mapping[str, Any],
    selected_indices: Sequence[int],
    power_table_npz: str,
    ops: Any,
    sigma_n_S_m: float,
    thermal_active_mask: np.ndarray,
    thermal_bath_K: float,
    chunk_size: int = 64,
    global_limits_override: Mapping[str, Any] | None = None,
    progress_callback: Callable[[int, int], None] | None = None,
) -> dict[str, Any]:
    """Return selected maps and exact finite extrema over every snapshot."""

    real = _snapshot_matrix(snapshots, "psi_real_snapshot_J")
    imag = _snapshot_matrix(snapshots, "psi_imag_snapshot_J")
    delta_meV = _snapshot_matrix(snapshots, "delta_snapshot_meV")
    phi_V = _snapshot_matrix(snapshots, "phi_snapshot_V")
    Te_K = _snapshot_matrix(snapshots, "Te_snapshot_K")
    Tph_K = _snapshot_matrix(snapshots, "Tph_snapshot_K")
    shape = real.shape
    for name, values in (
        ("psi_imag_snapshot_J", imag),
        ("delta_snapshot_meV", delta_meV),
        ("phi_snapshot_V", phi_V),
        ("Te_snapshot_K", Te_K),
        ("Tph_snapshot_K", Tph_K),
    ):
        if values.shape != shape:
            raise ValueError(f"{name} must share snapshot shape {shape}; got {values.shape}.")

    n_snap, n_nodes = shape
    if int(getattr(ops, "n_nodes", n_nodes)) != n_nodes:
        raise ValueError("Finite-volume operators and photon snapshots disagree on n_nodes.")
    selected = np.asarray(selected_indices, dtype=np.int64).reshape(-1)
    if selected.size == 0 or np.any(selected < 0) or np.any(selected >= n_snap):
        raise ValueError("selected_indices must contain valid photon snapshot indices.")
    if np.unique(selected).size != selected.size:
        raise ValueError("selected_indices must be unique.")
    sigma = float(sigma_n_S_m)
    if not np.isfinite(sigma) or sigma <= 0.0:
        raise ValueError("sigma_n_S_m must be positive and finite.")
    active = np.asarray(thermal_active_mask, dtype=bool).reshape(-1)
    if active.size != n_nodes:
        raise ValueError("thermal_active_mask must contain one value per mesh node.")
    bath = float(thermal_bath_K)
    if not np.isfinite(bath) or bath <= 0.0:
        raise ValueError("thermal_bath_K must be positive and finite.")

    times_ps = _snapshot_times_ps(snapshots, n_snap=n_snap)
    interpolator = PowerTableRuntimeInterpolator(power_table_npz)
    projector = _SnapshotProjector(ops)
    selection_position = np.full(n_snap, -1, dtype=np.int64)
    selection_position[selected] = np.arange(selected.size, dtype=np.int64)

    keys = (
        "q_abs_snapshot_m_inv",
        "joule_snapshot_W_m3",
        "P_S_snapshot_W_m3",
        "P_R_snapshot_W_m3",
        "P_total_snapshot_W_m3",
        "P_diff_snapshot_W_m3",
        "P_esc_snapshot_W_m3",
        "u_e_snapshot_J_m3",
        "u_ph_snapshot_J_m3",
        "C_e_snapshot_J_m3_K",
        "C_ph_snapshot_J_m3_K",
        "kappa_s_snapshot_W_m_K",
    )
    selected_maps = {
        key: np.empty((selected.size, n_nodes), dtype=float) for key in keys
    }
    override = dict(global_limits_override or {})
    reuse_limits = all(
        key in override
        and np.asarray(override[key], dtype=float).size == 2
        and np.all(np.isfinite(np.asarray(override[key], dtype=float)))
        for key in keys
    )
    extrema = {
        key: (
            [float(value) for value in np.asarray(override[key], dtype=float).reshape(-1)]
            if reuse_limits
            else [float("inf"), float("-inf")]
        )
        for key in keys
    }

    stride = max(1, int(chunk_size))
    scan_indices = selected if reuse_limits else np.arange(n_snap, dtype=np.int64)
    for offset in range(0, scan_indices.size, stride):
        batch = scan_indices[offset : offset + stride]
        psi = real[batch] + 1j * imag[batch]
        edge_phase = np.angle(
            psi[:, projector.edge_j] * np.conjugate(psi[:, projector.edge_i])
        )
        q_plot = projector.phase_gradient_magnitude(edge_phase)
        q_thermal = projector.edge_absolute_to_node(
            np.abs(edge_phase) / projector.edge_length_m[None, :]
        )

        lookup = interpolator.evaluate(
            Te_K=Te_K[batch],
            Tph_K=Tph_K[batch],
            delta_J=np.maximum(delta_meV[batch], 0.0) * MEV_J,
            q_abs_m_inv=q_thermal,
        )
        joule = projector.joule_power_density(phi_V[batch], sigma_n_S_m=sigma)
        diffusion = projector.diffusion_power_density(
            Te_K[batch],
            lookup.kappa_s_W_m_K,
            active_mask=active,
            bath_K=bath,
        )
        chunk = {
            "q_abs_snapshot_m_inv": q_plot,
            "joule_snapshot_W_m3": joule,
            "P_S_snapshot_W_m3": lookup.P_S_W_m3,
            "P_R_snapshot_W_m3": lookup.P_R_W_m3,
            "P_total_snapshot_W_m3": lookup.P_total_W_m3,
            "P_diff_snapshot_W_m3": diffusion,
            "P_esc_snapshot_W_m3": lookup.P_esc_W_m3,
            "u_e_snapshot_J_m3": lookup.u_e_J_m3,
            "u_ph_snapshot_J_m3": lookup.u_ph_J_m3,
            "C_e_snapshot_J_m3_K": lookup.C_e_J_m3_K,
            "C_ph_snapshot_J_m3_K": lookup.C_ph_J_m3_K,
            "kappa_s_snapshot_W_m_K": lookup.kappa_s_W_m_K,
        }
        if not reuse_limits:
            for key, values in chunk.items():
                _update_extrema(extrema[key], values)

        keep = selection_position[batch] >= 0
        if np.any(keep):
            destinations = selection_position[batch[keep]]
            for key, values in chunk.items():
                selected_maps[key][destinations] = np.asarray(values, dtype=float)[keep]
        if progress_callback is not None:
            progress_callback(min(offset + batch.size, scan_indices.size), scan_indices.size)

    result: dict[str, Any] = {
        "snapshot_t_ps": times_ps[selected],
        "selected_indices": selected,
        **selected_maps,
        "snapshot_global_limits": {
            key: np.asarray(bounds, dtype=float) for key, bounds in extrema.items()
        },
        "global_scale_policy": "exact finite extrema over all persisted photon snapshots",
        "global_limits_reused": bool(reuse_limits),
    }
    return result


class _SnapshotProjector:
    """Sparse batched equivalents of the runtime edge/node operations."""

    def __init__(self, ops: Any) -> None:
        self.edge_i = np.asarray(ops.edge_i, dtype=np.int64)
        self.edge_j = np.asarray(ops.edge_j, dtype=np.int64)
        self.edge_length_m = np.asarray(ops.edge_length_m, dtype=float)
        self.dual_face_length_m = np.asarray(ops.dual_face_length_m, dtype=float)
        self.node_area_m2 = np.asarray(ops.node_area_m2, dtype=float)
        self.n_nodes = self.node_area_m2.size
        self.n_edges = self.edge_i.size
        if self.edge_j.size != self.n_edges or self.edge_length_m.size != self.n_edges:
            raise ValueError("Incomplete finite-volume edge geometry.")

        columns = np.concatenate(
            [np.arange(self.n_edges, dtype=np.int64)] * 2
        )
        rows = np.concatenate([self.edge_i, self.edge_j])
        dx = np.asarray(ops.edge_unit, dtype=float)[:, 0] * self.edge_length_m
        dy = np.asarray(ops.edge_unit, dtype=float)[:, 1] * self.edge_length_m
        self._phase_x = sparse.csr_matrix(
            (np.concatenate([dx, dx]), (rows, columns)),
            shape=(self.n_nodes, self.n_edges),
        )
        self._phase_y = sparse.csr_matrix(
            (np.concatenate([dy, dy]), (rows, columns)),
            shape=(self.n_nodes, self.n_edges),
        )
        self._Axx = np.bincount(rows, weights=np.concatenate([dx * dx, dx * dx]), minlength=self.n_nodes)
        self._Axy = np.bincount(rows, weights=np.concatenate([dx * dy, dx * dy]), minlength=self.n_nodes)
        self._Ayy = np.bincount(rows, weights=np.concatenate([dy * dy, dy * dy]), minlength=self.n_nodes)
        self._det = self._Axx * self._Ayy - self._Axy * self._Axy
        self._phase_good = np.isfinite(self._det) & (np.abs(self._det) > 1.0e-300)

        weights = self.dual_face_length_m / np.maximum(self.edge_length_m, 1.0e-300)
        weights = np.maximum(weights, 1.0e-300)
        weight_sum = np.bincount(rows, weights=np.concatenate([weights, weights]), minlength=self.n_nodes)
        normalized = np.concatenate([weights, weights]) / np.maximum(weight_sum[rows], 1.0e-300)
        self._edge_average = sparse.csr_matrix(
            (normalized, (rows, columns)),
            shape=(self.n_nodes, self.n_edges),
        )
        diffusion_values = np.concatenate(
            [np.ones(self.n_edges), -np.ones(self.n_edges)]
        )
        self._diffusion_incidence = sparse.csr_matrix(
            (diffusion_values, (rows, columns)),
            shape=(self.n_nodes, self.n_edges),
        )
        self._conductance = self.dual_face_length_m / np.maximum(
            self.edge_length_m, 1.0e-300
        )

    def phase_gradient_magnitude(self, edge_phase: np.ndarray) -> np.ndarray:
        phase = np.asarray(edge_phase, dtype=float)
        bx = np.asarray(self._phase_x @ phase.T).T
        by = np.asarray(self._phase_y @ phase.T).T
        qx = np.zeros((phase.shape[0], self.n_nodes), dtype=float)
        qy = np.zeros_like(qx)
        good = self._phase_good
        qx[:, good] = (
            self._Ayy[good][None, :] * bx[:, good]
            - self._Axy[good][None, :] * by[:, good]
        ) / self._det[good][None, :]
        qy[:, good] = (
            -self._Axy[good][None, :] * bx[:, good]
            + self._Axx[good][None, :] * by[:, good]
        ) / self._det[good][None, :]
        out = np.hypot(qx, qy)
        out[~np.isfinite(out)] = 0.0
        return out

    def edge_absolute_to_node(self, edge_values: np.ndarray) -> np.ndarray:
        return np.asarray(self._edge_average @ np.asarray(edge_values, dtype=float).T).T

    def joule_power_density(
        self,
        phi_V: np.ndarray,
        *,
        sigma_n_S_m: float,
    ) -> np.ndarray:
        phi = np.asarray(phi_V, dtype=float)
        gradient = (
            phi[:, self.edge_j] - phi[:, self.edge_i]
        ) / np.maximum(self.edge_length_m[None, :], 1.0e-300)
        edge_joule = float(sigma_n_S_m) * gradient * gradient
        return self.edge_absolute_to_node(edge_joule)

    def diffusion_power_density(
        self,
        Te_K: np.ndarray,
        kappa_W_m_K: np.ndarray,
        *,
        active_mask: np.ndarray,
        bath_K: float,
    ) -> np.ndarray:
        temperature = np.asarray(Te_K, dtype=float)
        kappa = np.asarray(kappa_W_m_K, dtype=float)
        active = np.asarray(active_mask, dtype=bool)
        Ti = np.where(active[self.edge_i][None, :], temperature[:, self.edge_i], float(bath_K))
        Tj = np.where(active[self.edge_j][None, :], temperature[:, self.edge_j], float(bath_K))
        ki = np.where(active[self.edge_i][None, :], kappa[:, self.edge_i], 0.0)
        kj = np.where(active[self.edge_j][None, :], kappa[:, self.edge_j], 0.0)
        kedge = 0.5 * (ki + kj)
        cross = active[self.edge_i] ^ active[self.edge_j]
        if np.any(cross):
            kedge[:, cross] = np.where(
                active[self.edge_i[cross]][None, :],
                kappa[:, self.edge_i[cross]],
                kappa[:, self.edge_j[cross]],
            )
        flux = kedge * self._conductance[None, :] * (Tj - Ti)
        masked_incidence = self._diffusion_incidence.multiply(
            active[:, None]
        )
        return np.asarray(masked_incidence @ flux.T).T / np.maximum(
            self.node_area_m2[None, :], 1.0e-300
        )


def _snapshot_matrix(snapshots: Mapping[str, Any], key: str) -> np.ndarray:
    values = np.asarray(snapshots.get(key, []), dtype=float)
    if values.ndim != 2:
        raise ValueError(f"Photon snapshots lack 2D field {key}.")
    return values


def _snapshot_times_ps(snapshots: Mapping[str, Any], *, n_snap: int) -> np.ndarray:
    if "snapshot_t_ps" in snapshots:
        values = np.asarray(snapshots["snapshot_t_ps"], dtype=float).reshape(-1)
    else:
        values = np.asarray(snapshots.get("snapshot_t_s", []), dtype=float).reshape(-1) / 1.0e-12
    if values.size != n_snap:
        raise ValueError("Photon snapshot times do not match the stored fields.")
    return values


def _update_extrema(bounds: list[float], values: np.ndarray) -> None:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size:
        bounds[0] = min(bounds[0], float(np.min(finite)))
        bounds[1] = max(bounds[1], float(np.max(finite)))


__all__ = ["compute_photon_snapshot_plot_diagnostics"]
