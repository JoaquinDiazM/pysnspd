"""Runtime two-temperature thermal coupling for the SS gTDGL solver.

The first coupled version intentionally keeps the thermal problem local and
conservative: only a configurable central strip is dynamic, external nodes are
ideal reservoirs at the bath temperature, and the Joule source is positive
``|j_n|^2/sigma_n``.  The PRE power table supplies the electron--phonon powers,
heat capacities, superconducting thermal conductivity and phonon escape term.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

MEV_J = 1.602176634e-22


@dataclass(frozen=True)
class ThermalRuntimeConfig:
    enabled: bool = False
    window_m: float = 100.0e-9
    start_time_s: float = 2.0e-12
    bath_K: float = 0.9
    min_K: float | None = None
    max_K: float | None = None
    max_step_K: float = 0.05
    max_substeps: int = 64
    stationarity_rate_K_per_ps: float = 1.0e-2


@dataclass(frozen=True)
class ThermalLookupResult:
    P_S_W_m3: np.ndarray
    P_R_W_m3: np.ndarray
    P_total_W_m3: np.ndarray
    u_e_J_m3: np.ndarray
    C_e_J_m3_K: np.ndarray
    kappa_s_W_m_K: np.ndarray
    u_ph_J_m3: np.ndarray
    C_ph_J_m3_K: np.ndarray
    P_esc_W_m3: np.ndarray


class PowerTableRuntimeInterpolator:
    """Small multilinear interpolator for PRE power-table catalogues."""

    def __init__(self, path: str | Path):
        self.path = Path(path)
        with np.load(self.path, allow_pickle=True) as data:
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
            self.Te_values_K = np.asarray(data["Te_values_K"], dtype=float)
            self.Tph_values_K = np.asarray(data["Tph_values_K"], dtype=float)
            self.delta_values_J = np.asarray(data["delta_values_J"], dtype=float)
            self.q_values_m_inv = np.asarray(data["q_values_m_inv"], dtype=float)
            self.P_S_W_m3 = np.asarray(data["P_S_W_m3"], dtype=float)
            self.P_R_W_m3 = np.asarray(data["P_R_W_m3"], dtype=float)
            self.P_total_W_m3 = np.asarray(data["P_total_W_m3"], dtype=float)
            self.u_e_J_m3 = np.asarray(data["u_e_J_m3"], dtype=float)
            self.C_e_J_m3_K = np.asarray(data["C_e_J_m3_K"], dtype=float)
            self.kappa_s_W_m_K = np.asarray(data.get("kappa_s_W_m_K", np.array([], dtype=float)), dtype=float)
            self.u_ph_J_m3 = np.asarray(data.get("u_ph_J_m3", np.array([], dtype=float)), dtype=float)
            self.C_ph_J_m3_K = np.asarray(data.get("C_ph_J_m3_K", np.array([], dtype=float)), dtype=float)
            self.P_esc_W_m3 = np.asarray(data.get("P_esc_W_m3", np.array([], dtype=float)), dtype=float)

        self.Te_min_K = float(np.nanmin(self.Te_values_K))
        self.Te_max_K = float(np.nanmax(self.Te_values_K))
        self.Tph_min_K = float(np.nanmin(self.Tph_values_K))
        self.Tph_max_K = float(np.nanmax(self.Tph_values_K))
        self.delta_min_J = float(np.nanmin(self.delta_values_J))
        self.delta_max_J = float(np.nanmax(self.delta_values_J))
        self.q_min_m_inv = float(np.nanmin(self.q_values_m_inv))
        self.q_max_m_inv = float(np.nanmax(self.q_values_m_inv))

    def evaluate(
        self,
        *,
        Te_K: np.ndarray,
        Tph_K: np.ndarray,
        delta_J: np.ndarray,
        q_abs_m_inv: np.ndarray,
    ) -> ThermalLookupResult:
        Te = np.asarray(Te_K, dtype=float)
        Tph = np.asarray(Tph_K, dtype=float)
        delta = np.asarray(delta_J, dtype=float)
        q_abs = np.asarray(q_abs_m_inv, dtype=float)
        Te = np.clip(Te, self.Te_min_K, self.Te_max_K)
        Tph = np.clip(Tph, self.Tph_min_K, self.Tph_max_K)
        delta = np.clip(delta, self.delta_min_J, self.delta_max_J)
        q_abs = np.clip(q_abs, self.q_min_m_inv, self.q_max_m_inv)

        axes4 = (self.Te_values_K, self.Tph_values_K, self.delta_values_J, self.q_values_m_inv)
        coords4 = (Te, Tph, delta, q_abs)
        axes3 = (self.Te_values_K, self.delta_values_J, self.q_values_m_inv)
        coords3 = (Te, delta, q_abs)
        axes2 = (self.Te_values_K, self.delta_values_J)
        coords2 = (Te, delta)

        P_S = _interp_nd(self.P_S_W_m3, axes4, coords4)
        P_R = _interp_nd(self.P_R_W_m3, axes4, coords4)
        P_total = _interp_nd(self.P_total_W_m3, axes4, coords4)
        u_e = _interp_nd(self.u_e_J_m3, axes3, coords3)
        C_e = _interp_nd(self.C_e_J_m3_K, axes3, coords3)
        if self.kappa_s_W_m_K.size:
            # Current PRE tables store kappa_s[Te,Delta].
            if self.kappa_s_W_m_K.ndim == 2:
                kappa = _interp_nd(self.kappa_s_W_m_K, axes2, coords2)
            else:
                kappa = _interp_nd(self.kappa_s_W_m_K, axes3, coords3)
        else:
            kappa = np.zeros_like(P_total)
        u_ph = _interp_1d_or_nan(self.u_ph_J_m3, self.Tph_values_K, Tph)
        C_ph = _interp_1d_or_nan(self.C_ph_J_m3_K, self.Tph_values_K, Tph)
        P_esc = _interp_1d_or_nan(self.P_esc_W_m3, self.Tph_values_K, Tph)
        return ThermalLookupResult(
            P_S_W_m3=np.asarray(P_S, dtype=float),
            P_R_W_m3=np.asarray(P_R, dtype=float),
            P_total_W_m3=np.asarray(P_total, dtype=float),
            u_e_J_m3=np.asarray(u_e, dtype=float),
            C_e_J_m3_K=np.asarray(C_e, dtype=float),
            kappa_s_W_m_K=np.asarray(kappa, dtype=float),
            u_ph_J_m3=np.asarray(u_ph, dtype=float),
            C_ph_J_m3_K=np.asarray(C_ph, dtype=float),
            P_esc_W_m3=np.asarray(P_esc, dtype=float),
        )


def build_central_thermal_mask(nodes_m: np.ndarray, *, window_m: float) -> np.ndarray:
    nodes = np.asarray(nodes_m, dtype=float)
    if nodes.ndim != 2 or nodes.shape[1] < 1:
        raise ValueError("nodes_m must have shape (n_nodes, >=1).")
    width = float(window_m)
    if width <= 0.0:
        raise ValueError("thermal window_m must be positive.")
    x = nodes[:, 0]
    xc = 0.5 * (float(np.nanmin(x)) + float(np.nanmax(x)))
    return np.abs(x - xc) <= 0.5 * width


class ThermalRuntimeController:
    """Explicit runtime thermal stepper coupled through mutable Te/Tph arrays."""

    def __init__(
        self,
        *,
        nodes_m: np.ndarray,
        ops: Any,
        material: Any,
        Te_K: np.ndarray,
        Tph_K: np.ndarray,
        power_table_npz: str | Path,
        config: ThermalRuntimeConfig,
    ):
        self.nodes_m = np.asarray(nodes_m, dtype=float)
        self.ops = ops
        self.material = material
        self.Te_K = np.asarray(Te_K, dtype=float)
        self.Tph_K = np.asarray(Tph_K, dtype=float)
        self.config = config
        self.interpolator = PowerTableRuntimeInterpolator(power_table_npz)
        self.mask = build_central_thermal_mask(self.nodes_m, window_m=float(config.window_m))
        self.bath_K = float(config.bath_K)
        self.min_K = float(config.min_K) if config.min_K is not None else self.bath_K
        cat_max = min(float(self.interpolator.Te_max_K), float(self.interpolator.Tph_max_K))
        self.max_K = float(config.max_K) if config.max_K is not None else cat_max
        self.max_step_K = max(float(config.max_step_K), 1.0e-12)
        self.max_substeps = max(1, int(config.max_substeps))
        self.Te_K[~self.mask] = self.bath_K
        self.Tph_K[~self.mask] = self.bath_K
        self.last: dict[str, float] = self._inactive_diag()

    def snapshot_payload(self) -> dict[str, np.ndarray]:
        return {
            "Te_snapshot_K": np.asarray(self.Te_K, dtype=float).copy(),
            "Tph_snapshot_K": np.asarray(self.Tph_K, dtype=float).copy(),
        }

    def step(
        self,
        *,
        time_s: float,
        dt_s: float,
        psi_dimensionless: np.ndarray,
        native_normal_current: np.ndarray,
        current_scale_A_m2: float,
    ) -> dict[str, float]:
        if (not self.config.enabled) or float(time_s) < float(self.config.start_time_s) or not np.any(self.mask):
            self.Te_K[~self.mask] = self.bath_K
            self.Tph_K[~self.mask] = self.bath_K
            self.last = self._inactive_diag(active_time_s=float(time_s))
            return self.last

        dt_total = max(float(dt_s), 0.0)
        if dt_total <= 0.0:
            self.last = self._inactive_diag(active_time_s=float(time_s))
            return self.last

        substeps = 1
        max_dTe = max_dTph = 0.0
        max_rate = 0.0
        max_PJ = max_Pep = max_Pesc = max_diff = 0.0
        for _ in range(self.max_substeps):
            deriv = self._derivatives(
                psi_dimensionless=psi_dimensionless,
                native_normal_current=native_normal_current,
                current_scale_A_m2=current_scale_A_m2,
            )
            active = self.mask
            dTe_trial = dt_total * deriv["dTe_dt_K_s"]
            dTph_trial = dt_total * deriv["dTph_dt_K_s"]
            largest = float(max(np.nanmax(np.abs(dTe_trial[active])), np.nanmax(np.abs(dTph_trial[active]))))
            if largest <= self.max_step_K or substeps >= self.max_substeps:
                break
            substeps = min(self.max_substeps, max(1, int(np.ceil(largest / self.max_step_K))))
            break

        dt_sub = dt_total / float(substeps)
        for _ in range(substeps):
            deriv = self._derivatives(
                psi_dimensionless=psi_dimensionless,
                native_normal_current=native_normal_current,
                current_scale_A_m2=current_scale_A_m2,
            )
            active = self.mask
            dTe = dt_sub * deriv["dTe_dt_K_s"]
            dTph = dt_sub * deriv["dTph_dt_K_s"]
            self.Te_K[active] = np.clip(self.Te_K[active] + dTe[active], self.min_K, self.max_K)
            self.Tph_K[active] = np.clip(self.Tph_K[active] + dTph[active], self.min_K, self.max_K)
            self.Te_K[~active] = self.bath_K
            self.Tph_K[~active] = self.bath_K
            max_dTe = max(max_dTe, float(np.nanmax(np.abs(dTe[active]))))
            max_dTph = max(max_dTph, float(np.nanmax(np.abs(dTph[active]))))
            max_rate = max(max_rate, float(np.nanmax(np.abs(deriv["dTe_dt_K_s"][active]))) / 1.0e12)
            max_rate = max(max_rate, float(np.nanmax(np.abs(deriv["dTph_dt_K_s"][active]))) / 1.0e12)
            max_PJ = max(max_PJ, float(np.nanmax(np.abs(deriv["P_J_W_m3"][active]))))
            max_Pep = max(max_Pep, float(np.nanmax(np.abs(deriv["P_ep_W_m3"][active]))))
            max_Pesc = max(max_Pesc, float(np.nanmax(np.abs(deriv["P_esc_W_m3"][active]))))
            max_diff = max(max_diff, float(np.nanmax(np.abs(deriv["P_diff_W_m3"][active]))))

        active = self.mask
        self.last = {
            "thermal_enabled": 1.0,
            "thermal_active": 1.0,
            "thermal_active_n_nodes": float(np.count_nonzero(active)),
            "thermal_substeps": float(substeps),
            "thermal_max_Te_K": float(np.nanmax(self.Te_K[active])),
            "thermal_mean_Te_K": float(np.nanmean(self.Te_K[active])),
            "thermal_max_Tph_K": float(np.nanmax(self.Tph_K[active])),
            "thermal_mean_Tph_K": float(np.nanmean(self.Tph_K[active])),
            "thermal_max_abs_dTe_K": max_dTe,
            "thermal_max_abs_dTph_K": max_dTph,
            "thermal_max_rate_K_per_ps": max_rate,
            "thermal_max_P_J_W_m3": max_PJ,
            "thermal_max_P_ep_W_m3": max_Pep,
            "thermal_max_P_esc_W_m3": max_Pesc,
            "thermal_max_P_diff_W_m3": max_diff,
        }
        return self.last

    def _derivatives(
        self,
        *,
        psi_dimensionless: np.ndarray,
        native_normal_current: np.ndarray,
        current_scale_A_m2: float,
    ) -> dict[str, np.ndarray]:
        active = self.mask
        delta_J = np.abs(np.asarray(psi_dimensionless, dtype=np.complex128)) * float(self.material.delta0_J)
        edge_q = np.abs(_edge_phase_gradient_from_psi(psi_dimensionless, self.ops))
        q_node = _edge_to_node_average(edge_q, self.ops, n_nodes=self.Te_K.size)
        edge_jn_A_m2 = np.asarray(native_normal_current, dtype=float) * float(current_scale_A_m2)
        edge_joule = edge_jn_A_m2 * edge_jn_A_m2 / max(float(self.material.sigma_n_S_m), 1.0e-300)
        P_J = _edge_to_node_average(edge_joule, self.ops, n_nodes=self.Te_K.size)

        lookup = self.interpolator.evaluate(
            Te_K=self.Te_K[active],
            Tph_K=self.Tph_K[active],
            delta_J=delta_J[active],
            q_abs_m_inv=q_node[active],
        )
        C_e = _positive_finite(lookup.C_e_J_m3_K, fallback=1.0, floor=1.0e-12)
        C_ph = _positive_finite(lookup.C_ph_J_m3_K, fallback=1.0, floor=1.0e-12)
        P_ep = np.zeros_like(self.Te_K)
        P_esc = np.zeros_like(self.Te_K)
        kappa = np.zeros_like(self.Te_K)
        P_ep[active] = np.nan_to_num(lookup.P_total_W_m3, nan=0.0, posinf=0.0, neginf=0.0)
        P_esc[active] = np.nan_to_num(lookup.P_esc_W_m3, nan=0.0, posinf=0.0, neginf=0.0)
        kappa[active] = _positive_finite(lookup.kappa_s_W_m_K, fallback=0.0, floor=0.0)
        P_diff = _thermal_diffusion_power_density(self.Te_K, kappa, active, bath_K=self.bath_K, ops=self.ops)

        dTe_dt = np.zeros_like(self.Te_K)
        dTph_dt = np.zeros_like(self.Tph_K)
        dTe_dt[active] = (P_diff[active] + P_J[active] - P_ep[active]) / C_e
        dTph_dt[active] = (P_ep[active] - P_esc[active]) / C_ph
        return {
            "dTe_dt_K_s": dTe_dt,
            "dTph_dt_K_s": dTph_dt,
            "P_J_W_m3": P_J,
            "P_ep_W_m3": P_ep,
            "P_esc_W_m3": P_esc,
            "P_diff_W_m3": P_diff,
        }

    def _inactive_diag(self, *, active_time_s: float = 0.0) -> dict[str, float]:
        active = self.mask if hasattr(self, "mask") else np.array([], dtype=bool)
        if active.size and np.any(active):
            max_te = float(np.nanmax(self.Te_K[active]))
            max_tph = float(np.nanmax(self.Tph_K[active]))
            mean_te = float(np.nanmean(self.Te_K[active]))
            mean_tph = float(np.nanmean(self.Tph_K[active]))
            n_active = float(np.count_nonzero(active))
        else:
            max_te = max_tph = mean_te = mean_tph = float(self.bath_K if hasattr(self, "bath_K") else 0.0)
            n_active = 0.0
        return {
            "thermal_enabled": float(bool(self.config.enabled)) if hasattr(self, "config") else 0.0,
            "thermal_active": 0.0,
            "thermal_active_n_nodes": n_active,
            "thermal_substeps": 0.0,
            "thermal_max_Te_K": max_te,
            "thermal_mean_Te_K": mean_te,
            "thermal_max_Tph_K": max_tph,
            "thermal_mean_Tph_K": mean_tph,
            "thermal_max_abs_dTe_K": 0.0,
            "thermal_max_abs_dTph_K": 0.0,
            "thermal_max_rate_K_per_ps": 0.0,
            "thermal_max_P_J_W_m3": 0.0,
            "thermal_max_P_ep_W_m3": 0.0,
            "thermal_max_P_esc_W_m3": 0.0,
            "thermal_max_P_diff_W_m3": 0.0,
        }


def thermal_stationarity_diagnostics(
    history: Mapping[str, Any],
    *,
    enabled: bool,
    start_time_s: float,
    requested_total_time_s: float,
    rate_tol_K_per_ps: float,
) -> dict[str, Any]:
    if not enabled:
        return {
            "enabled": False,
            "passes": True,
            "reason": "thermal coupling disabled",
            "rate_tol_K_per_ps": float(rate_tol_K_per_ps),
        }
    rate = np.asarray(history.get("thermal_max_rate_K_per_ps", []), dtype=float)
    t_s = np.asarray(history.get("t_s", []), dtype=float)
    if rate.size == 0 or t_s.size == 0 or float(requested_total_time_s) <= float(start_time_s):
        return {
            "enabled": True,
            "passes": False,
            "reason": "no active thermal relaxation interval was recorded",
            "rate_tol_K_per_ps": float(rate_tol_K_per_ps),
        }
    if rate.size != t_s.size:
        rate = np.resize(rate, t_s.size)
    active = t_s >= float(start_time_s)
    if not np.any(active):
        return {
            "enabled": True,
            "passes": False,
            "reason": "thermal start time was not reached",
            "rate_tol_K_per_ps": float(rate_tol_K_per_ps),
        }
    idx = np.flatnonzero(active)
    tail_n = max(4, int(np.ceil(0.10 * idx.size)))
    tail = rate[idx[-tail_n:]]
    tail_rate = float(np.nanmax(tail)) if tail.size else float("inf")
    return {
        "enabled": True,
        "passes": bool(np.isfinite(tail_rate) and tail_rate <= float(rate_tol_K_per_ps)),
        "reason": "tail thermal rate below tolerance" if np.isfinite(tail_rate) and tail_rate <= float(rate_tol_K_per_ps) else "tail thermal rate above tolerance",
        "rate_tol_K_per_ps": float(rate_tol_K_per_ps),
        "tail_max_rate_K_per_ps": tail_rate,
        "tail_window_steps": int(tail_n),
        "start_time_ps": float(start_time_s) / 1.0e-12,
    }


def _interp_axis(axis: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ax = np.asarray(axis, dtype=float)
    vals = np.asarray(values, dtype=float)
    if ax.ndim != 1 or ax.size == 0:
        raise ValueError("interpolation axis must be a non-empty 1D array.")
    if ax.size == 1:
        z = np.zeros(vals.shape, dtype=np.int64)
        return z, z, np.zeros(vals.shape, dtype=float)
    clipped = np.clip(vals, float(ax[0]), float(ax[-1]))
    hi = np.searchsorted(ax, clipped, side="right")
    hi = np.clip(hi, 1, ax.size - 1).astype(np.int64)
    lo = hi - 1
    denom = np.maximum(ax[hi] - ax[lo], 1.0e-300)
    frac = (clipped - ax[lo]) / denom
    return lo, hi, frac


def _interp_nd(table: np.ndarray, axes: tuple[np.ndarray, ...], coords: tuple[np.ndarray, ...]) -> np.ndarray:
    arr = np.asarray(table, dtype=float)
    if arr.ndim != len(axes):
        raise ValueError(f"table ndim {arr.ndim} does not match {len(axes)} axes")
    lo_hi_f = [_interp_axis(axis, coord) for axis, coord in zip(axes, coords)]
    out = np.zeros(np.asarray(coords[0]).shape, dtype=float)
    for bits in np.ndindex(*(2 for _ in axes)):
        weight = np.ones_like(out, dtype=float)
        index = []
        for bit, (lo, hi, frac) in zip(bits, lo_hi_f):
            if bit:
                index.append(hi)
                weight = weight * frac
            else:
                index.append(lo)
                weight = weight * (1.0 - frac)
        out = out + weight * arr[tuple(index)]
    return out


def _interp_1d_or_nan(table: np.ndarray, axis: np.ndarray, values: np.ndarray) -> np.ndarray:
    arr = np.asarray(table, dtype=float)
    vals = np.asarray(values, dtype=float)
    if arr.size == 0:
        return np.full(vals.shape, np.nan, dtype=float)
    return _interp_nd(arr, (axis,), (vals,))


def _positive_finite(values: np.ndarray, *, fallback: float, floor: float) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    finite_pos = arr[np.isfinite(arr) & (arr > floor)]
    fb = float(np.nanmedian(finite_pos)) if finite_pos.size else float(fallback)
    out = np.nan_to_num(arr, nan=fb, posinf=fb, neginf=fb)
    return np.maximum(out, float(floor))


def _edge_to_node_average(edge_values: np.ndarray, ops: Any, *, n_nodes: int) -> np.ndarray:
    vals = np.asarray(edge_values, dtype=float)
    edge_i = np.asarray(ops.edge_i, dtype=np.int64)
    edge_j = np.asarray(ops.edge_j, dtype=np.int64)
    weights = np.asarray(ops.dual_face_length_m, dtype=float) / np.maximum(np.asarray(ops.edge_length_m, dtype=float), 1.0e-300)
    weights = np.maximum(weights, 1.0e-300)
    out = np.zeros(int(n_nodes), dtype=float)
    wsum = np.zeros(int(n_nodes), dtype=float)
    np.add.at(out, edge_i, weights * vals)
    np.add.at(out, edge_j, weights * vals)
    np.add.at(wsum, edge_i, weights)
    np.add.at(wsum, edge_j, weights)
    return out / np.maximum(wsum, 1.0e-300)


def _edge_phase_gradient_from_psi(psi: np.ndarray, ops: Any) -> np.ndarray:
    z = np.asarray(psi, dtype=np.complex128)
    dtheta = np.angle(z[np.asarray(ops.edge_j, dtype=np.int64)] * np.conjugate(z[np.asarray(ops.edge_i, dtype=np.int64)]))
    return dtheta / np.maximum(np.asarray(ops.edge_length_m, dtype=float), 1.0e-300)


def _thermal_diffusion_power_density(
    Te_K: np.ndarray,
    kappa_s_W_m_K: np.ndarray,
    active_mask: np.ndarray,
    *,
    bath_K: float,
    ops: Any,
) -> np.ndarray:
    T = np.asarray(Te_K, dtype=float)
    kappa = np.asarray(kappa_s_W_m_K, dtype=float)
    active = np.asarray(active_mask, dtype=bool)
    edge_i = np.asarray(ops.edge_i, dtype=np.int64)
    edge_j = np.asarray(ops.edge_j, dtype=np.int64)
    conductance = np.asarray(ops.dual_face_length_m, dtype=float) / np.maximum(np.asarray(ops.edge_length_m, dtype=float), 1.0e-300)
    Ti = np.where(active[edge_i], T[edge_i], float(bath_K))
    Tj = np.where(active[edge_j], T[edge_j], float(bath_K))
    ki = np.where(active[edge_i], kappa[edge_i], 0.0)
    kj = np.where(active[edge_j], kappa[edge_j], 0.0)
    kedge = 0.5 * (ki + kj)
    # If exactly one side is active, use the active-side kappa and a bath ghost node.
    cross_ij = active[edge_i] ^ active[edge_j]
    kedge[cross_ij] = np.where(active[edge_i[cross_ij]], kappa[edge_i[cross_ij]], kappa[edge_j[cross_ij]])
    flux_i_to_j = kedge * conductance * (Tj - Ti)
    out = np.zeros(T.size, dtype=float)
    use_i = active[edge_i]
    use_j = active[edge_j]
    np.add.at(out, edge_i[use_i], flux_i_to_j[use_i])
    np.add.at(out, edge_j[use_j], -flux_i_to_j[use_j])
    return out / np.maximum(np.asarray(ops.node_area_m2, dtype=float), 1.0e-300)


__all__ = [
    "PowerTableRuntimeInterpolator",
    "ThermalLookupResult",
    "ThermalRuntimeConfig",
    "ThermalRuntimeController",
    "build_central_thermal_mask",
    "thermal_stationarity_diagnostics",
]
