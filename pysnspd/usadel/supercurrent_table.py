"""Matsubara Usadel supercurrent tables for the gTDGL backend.

The SS solver needs a local constitutive closure

    j_s = J(q, |Delta|, T_e) q_hat

instead of the older calibration-only one-dimensional relation J(q, T_bias).
This module builds that table from the same Matsubara dirty-limit equation used
by the critical-current calibration.  The table is deliberately independent of
BCS self-consistency: |Delta| is an axis supplied by the local gTDGL field.
"""
from __future__ import annotations

from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor
from typing import Any, Iterable, Mapping

import numpy as np

from pysnspd.usadel.calibration import (
    matsubara_energy_axis_J,
    solve_matsubara_s_values,
)
from pysnspd.usadel.parameters import E_CHARGE_C, HBAR_J_S, K_B_J_K


@dataclass(frozen=True)
class SupercurrentTable3D:
    """Container for a strict 3D Usadel/Matsubara supercurrent table.

    The canonical storage layout is

        js_T_delta_q_A_m2[iT, iDelta, iq].

    The q-axis is nonnegative.  Interpolators recover odd-in-q behavior by
    interpolating |q| and multiplying by sign(q).
    """

    Te_axis_K: np.ndarray
    delta_axis_J: np.ndarray
    q_axis_m_inv: np.ndarray
    js_T_delta_q_A_m2: np.ndarray
    metadata: dict[str, Any]

    @property
    def shape(self) -> tuple[int, int, int]:
        return tuple(int(v) for v in self.js_T_delta_q_A_m2.shape)


def build_matsubara_supercurrent_table_3d(
    *,
    Te_axis_K: np.ndarray,
    delta_axis_J: np.ndarray,
    q_axis_m_inv: np.ndarray,
    D_m2_s: float,
    sigma_n_S_m: float,
    n_matsubara: int,
    workers: int = 1,
) -> SupercurrentTable3D:
    """Build ``J(q, |Delta|, T)`` from the Matsubara Usadel equation.

    For every table point we solve, for all positive Matsubara energies,

        Delta sqrt(1 - s_n^2) = (eps_n + Gamma_q sqrt(1 - s_n^2)) s_n,

    then evaluate the same current density convention used by the calibration
    layer,

        j_s = (2*pi*k_B*T/|e|) sigma_n q sum_n s_n^2.

    This is intentionally a local closure table.  It does not solve the BCS
    self-consistency equation for Delta because in the gTDGL sector Delta is the
    evolved local field.
    """

    Te_axis = _clean_axis_1d(Te_axis_K, name="Te_axis_K", positive=True)
    delta_axis = _clean_axis_1d(delta_axis_J, name="delta_axis_J", nonnegative=True)
    q_axis = _clean_axis_1d(q_axis_m_inv, name="q_axis_m_inv", nonnegative=True)

    D = float(D_m2_s)
    sigma = float(sigma_n_S_m)
    n_m = int(n_matsubara)
    if D <= 0.0 or not np.isfinite(D):
        raise ValueError("D_m2_s must be positive and finite.")
    if sigma <= 0.0 or not np.isfinite(sigma):
        raise ValueError("sigma_n_S_m must be positive and finite.")
    if n_m <= 0:
        raise ValueError("n_matsubara must be positive.")

    tasks = [(float(T), delta_axis, q_axis, D, sigma, n_m) for T in Te_axis]
    n_workers = max(1, int(workers))
    if n_workers == 1 or len(tasks) <= 1:
        planes = [_compute_temperature_plane(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            planes = list(pool.map(_compute_temperature_plane, tasks))

    table = np.stack(planes, axis=0).astype(float, copy=False)
    table[~np.isfinite(table)] = 0.0

    metadata = {
        "backend": "matsubara_usadel_supercurrent_table_3d_v1",
        "layout": "js_T_delta_q_A_m2[Te, delta, q]",
        "current_relation": "j_s=(2*pi*k_B*T/|e|)*sigma_n*q*sum_n(s_n^2)",
        "gamma_definition": "Gamma_q=hbar*D*q^2/2",
        "self_consistency": "not imposed; |Delta| is an explicit gTDGL local-field axis",
        "n_Te": int(Te_axis.size),
        "n_delta": int(delta_axis.size),
        "n_q": int(q_axis.size),
        "n_matsubara": int(n_m),
        "workers": int(n_workers),
        "Te_min_K": float(np.min(Te_axis)),
        "Te_max_K": float(np.max(Te_axis)),
        "delta_min_J": float(np.min(delta_axis)),
        "delta_max_J": float(np.max(delta_axis)),
        "q_min_m_inv": float(np.min(q_axis)),
        "q_max_m_inv": float(np.max(q_axis)),
    }
    return SupercurrentTable3D(
        Te_axis_K=Te_axis,
        delta_axis_J=delta_axis,
        q_axis_m_inv=q_axis,
        js_T_delta_q_A_m2=table,
        metadata=metadata,
    )


def append_supercurrent_table_3d_to_npz(npz_path: str, table: SupercurrentTable3D) -> None:
    """Append the strict 3D supercurrent table to an existing PRE NPZ."""

    with np.load(npz_path, allow_pickle=True) as data:
        arrays = {key: data[key] for key in data.files}

    arrays["js_A_m2"] = np.asarray(table.js_T_delta_q_A_m2, dtype=float)
    arrays["j_s_A_m2"] = arrays["js_A_m2"]
    arrays["js_T_delta_q_A_m2"] = arrays["js_A_m2"]
    arrays["q_axis_m_inv"] = np.asarray(table.q_axis_m_inv, dtype=float)
    arrays["delta_axis_J"] = np.asarray(table.delta_axis_J, dtype=float)
    arrays["Te_axis_K"] = np.asarray(table.Te_axis_K, dtype=float)
    arrays["js_table_layout"] = np.array("Te,delta,q")
    arrays["js_table_backend"] = np.array(str(table.metadata["backend"]))
    arrays["js_table_n_matsubara"] = np.array(int(table.metadata["n_matsubara"]), dtype=np.int64)
    arrays["js_table_n_Te"] = np.array(int(table.metadata["n_Te"]), dtype=np.int64)
    arrays["js_table_n_delta"] = np.array(int(table.metadata["n_delta"]), dtype=np.int64)
    arrays["js_table_n_q"] = np.array(int(table.metadata["n_q"]), dtype=np.int64)

    np.savez_compressed(npz_path, **arrays)


def supercurrent_table_summary(table: SupercurrentTable3D) -> dict[str, Any]:
    """Return a manifest-friendly summary."""

    arr = np.asarray(table.js_T_delta_q_A_m2, dtype=float)
    finite = arr[np.isfinite(arr)]
    max_abs = float(np.max(np.abs(finite))) if finite.size else float("nan")
    return {
        **table.metadata,
        "table_key": "js_A_m2",
        "alias_keys": ["j_s_A_m2", "js_T_delta_q_A_m2"],
        "axis_keys": ["Te_axis_K", "delta_axis_J", "q_axis_m_inv"],
        "shape": list(table.shape),
        "js_max_abs_A_m2": max_abs,
        "strict_required_by_ss": True,
    }


def _compute_temperature_plane(task: tuple[float, np.ndarray, np.ndarray, float, float, int]) -> np.ndarray:
    T, delta_axis, q_axis, D, sigma, n_m = task
    eps = matsubara_energy_axis_J(T_K=float(T), n_matsubara=int(n_m))
    gamma_axis = 0.5 * HBAR_J_S * D * q_axis * q_axis
    plane = np.zeros((delta_axis.size, q_axis.size), dtype=float)
    prefactor = 2.0 * np.pi * K_B_J_K * float(T) * sigma / E_CHARGE_C
    for i, delta in enumerate(delta_axis):
        if delta <= 0.0:
            continue
        for j, (q, gamma) in enumerate(zip(q_axis, gamma_axis)):
            if q <= 0.0:
                continue
            s = solve_matsubara_s_values(
                delta_J=float(delta),
                gamma_J=float(gamma),
                eps_n_J=eps,
            )
            plane[i, j] = prefactor * float(q) * float(np.sum(s * s))
    return plane


def _clean_axis_1d(
    values: np.ndarray,
    *,
    name: str,
    positive: bool = False,
    nonnegative: bool = False,
) -> np.ndarray:
    arr = np.asarray(values, dtype=float).reshape(-1)
    if arr.size == 0:
        raise ValueError(f"{name} must not be empty.")
    if np.any(~np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values.")
    arr = np.unique(arr)
    arr.sort()
    if positive and np.any(arr <= 0.0):
        raise ValueError(f"{name} must contain positive values only.")
    if nonnegative and np.any(arr < 0.0):
        raise ValueError(f"{name} must contain nonnegative values only.")
    return arr


def temperature_axis_from_request(
    *,
    T_bias_K: float,
    Tc_K: float,
    n_Te: int,
    Te_min_K: float | None = None,
    Te_max_K: float | None = None,
) -> np.ndarray:
    """Build a compact default temperature axis for the PRE current table."""

    n = max(1, int(n_Te))
    T_bias = float(T_bias_K)
    Tc = float(Tc_K)
    if Te_min_K is None and Te_max_K is None and n == 1:
        return np.array([T_bias], dtype=float)
    lo = T_bias if Te_min_K is None else float(Te_min_K)
    hi = min(0.98 * Tc, max(T_bias, 0.98 * Tc)) if Te_max_K is None else float(Te_max_K)
    if n == 1:
        return np.array([lo], dtype=float)
    if hi < lo:
        raise ValueError(f"Te_max_K must be >= Te_min_K, got {hi} < {lo}.")
    return np.linspace(lo, hi, n)
