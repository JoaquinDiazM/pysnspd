"""Matsubara Usadel constitutive tables for the gTDGL backend.

The temporal solvers interpolate the finite stiffness

    kappa(T_e, |Delta|^2, |q|),        j_s = kappa |Delta|^2 q,

instead of interpolating ``j_s`` itself.  This preserves the dirty-limit
``j_s = O(|Delta|^2 q)`` asymptotic law below the first nonzero amplitude node.
The current-density table is stored alongside the stiffness because it is cheap
to form during PRE and remains useful for diagnostics and plotting.

All table rows ``(T_e, |Delta|)`` are independent and are evaluated in parallel
when ``workers > 1``.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np

from pysnspd.usadel.calibration import (
    matsubara_energy_axis_J,
    solve_matsubara_s_values,
)
from pysnspd.usadel.parameters import E_CHARGE_C, HBAR_J_S, K_B_J_K


@dataclass(frozen=True)
class SupercurrentTable3D:
    """Container for strict 3D Usadel/Matsubara constitutive tables.

    Canonical storage layout:

        js_stiffness_T_delta2_q_A_per_m_J2[iT, iDelta2, iq]

    The ``delta_axis_J`` and ``delta2_axis_J2`` entries have the same indices.
    The q-axis is nonnegative because the stiffness is even in q.
    """

    Te_axis_K: np.ndarray
    delta_axis_J: np.ndarray
    delta2_axis_J2: np.ndarray
    q_axis_m_inv: np.ndarray
    js_T_delta_q_A_m2: np.ndarray
    js_stiffness_T_delta2_q_A_per_m_J2: np.ndarray
    metadata: dict[str, Any]

    @property
    def shape(self) -> tuple[int, int, int]:
        return tuple(int(v) for v in self.js_stiffness_T_delta2_q_A_per_m_J2.shape)


def build_matsubara_supercurrent_table_3d(
    *,
    Te_axis_K: np.ndarray,
    delta_axis_J: np.ndarray,
    q_axis_m_inv: np.ndarray,
    D_m2_s: float,
    sigma_n_S_m: float,
    n_matsubara: int,
    workers: int = 1,
    backend: str = "process",
) -> SupercurrentTable3D:
    """Build ``kappa(T, |Delta|^2, |q|)`` and ``j_s`` from Matsubara Usadel.

    For every table point we solve, for all positive Matsubara energies,

        Delta sqrt(1 - s_n^2) = (eps_n + Gamma_q sqrt(1 - s_n^2)) s_n,

    then evaluate

        kappa = (2*pi*k_B*T/|e|) sigma_n sum_n s_n^2 / |Delta|^2,
        j_s = kappa |Delta|^2 q.

    At exact ``|Delta|=0`` the finite analytic limit

        kappa_0 = (2*pi*k_B*T/|e|) sigma_n
                  sum_n (epsilon_n + Gamma_q)^-2

    is used.  ``|Delta|`` is a local-field axis, not a self-consistent BCS
    solution.
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

    n_workers = max(1, int(workers))
    mode = str(backend or "process").lower()
    if mode not in {"process", "thread", "serial"}:
        mode = "process"
    if mode == "serial":
        n_workers = 1

    current = np.zeros((Te_axis.size, delta_axis.size, q_axis.size), dtype=float)
    stiffness = np.zeros_like(current)
    tasks = list(_current_row_tasks(Te_axis, delta_axis, q_axis, D, sigma, n_m))

    if n_workers == 1 or len(tasks) <= 1:
        rows = [_compute_current_row(task) for task in tasks]
    else:
        executor_cls = ThreadPoolExecutor if mode == "thread" else ProcessPoolExecutor
        with executor_cls(max_workers=n_workers) as pool:
            rows = list(pool.map(_compute_current_row, tasks))

    for iT, iD, current_row, stiffness_row in rows:
        current[iT, iD, :] = current_row
        stiffness[iT, iD, :] = stiffness_row
    if np.any(~np.isfinite(current)) or np.any(~np.isfinite(stiffness)):
        raise FloatingPointError("Matsubara constitutive table contains non-finite values.")
    if np.any(stiffness <= 0.0):
        raise FloatingPointError("Matsubara stiffness table must be strictly positive.")

    metadata = {
        "backend": "matsubara_usadel_stiffness_table_3d_v2",
        "layout": "js_stiffness_T_delta2_q_A_per_m_J2[Te, delta2, q]",
        "current_relation": "j_s=(2*pi*k_B*T/|e|)*sigma_n*q*sum_n(s_n^2)",
        "stiffness_relation": "kappa=j_s/(|Delta|^2*q), with the analytic |Delta|->0 limit",
        "gamma_definition": "Gamma_q=hbar*D*q^2/2",
        "self_consistency": "not imposed; |Delta| is an explicit gTDGL local-field axis",
        "n_Te": int(Te_axis.size),
        "n_delta": int(delta_axis.size),
        "n_q": int(q_axis.size),
        "n_matsubara": int(n_m),
        "workers": int(n_workers),
        "parallel_backend": mode,
        "parallel_tasks": int(len(tasks)),
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
        delta2_axis_J2=delta_axis * delta_axis,
        q_axis_m_inv=q_axis,
        js_T_delta_q_A_m2=current,
        js_stiffness_T_delta2_q_A_per_m_J2=stiffness,
        metadata=metadata,
    )


def append_supercurrent_table_3d_to_npz(npz_path: str, table: SupercurrentTable3D) -> None:
    """Append strict 3D stiffness/current resources to an existing PRE NPZ."""

    with np.load(npz_path, allow_pickle=True) as data:
        arrays = {key: data[key] for key in data.files}

    arrays["js_A_m2"] = np.asarray(table.js_T_delta_q_A_m2, dtype=float)
    arrays["js_stiffness_A_per_m_J2"] = np.asarray(
        table.js_stiffness_T_delta2_q_A_per_m_J2,
        dtype=float,
    )
    arrays["q_axis_m_inv"] = np.asarray(table.q_axis_m_inv, dtype=float)
    arrays["delta_axis_J"] = np.asarray(table.delta_axis_J, dtype=float)
    arrays["delta2_axis_J2"] = np.asarray(table.delta2_axis_J2, dtype=float)
    arrays["Te_axis_K"] = np.asarray(table.Te_axis_K, dtype=float)
    arrays["js_table_layout"] = np.array("Te,delta2,q")
    arrays["js_table_backend"] = np.array(str(table.metadata["backend"]))
    arrays["js_table_n_matsubara"] = np.array(int(table.metadata["n_matsubara"]), dtype=np.int64)
    arrays["js_table_n_Te"] = np.array(int(table.metadata["n_Te"]), dtype=np.int64)
    arrays["js_table_n_delta"] = np.array(int(table.metadata["n_delta"]), dtype=np.int64)
    arrays["js_table_n_q"] = np.array(int(table.metadata["n_q"]), dtype=np.int64)
    arrays["js_table_parallel_workers"] = np.array(int(table.metadata["workers"]), dtype=np.int64)
    arrays["js_table_parallel_tasks"] = np.array(int(table.metadata["parallel_tasks"]), dtype=np.int64)

    np.savez_compressed(npz_path, **arrays)


def supercurrent_table_summary(table: SupercurrentTable3D) -> dict[str, Any]:
    """Return a manifest-friendly summary."""

    current = np.asarray(table.js_T_delta_q_A_m2, dtype=float)
    stiffness = np.asarray(table.js_stiffness_T_delta2_q_A_per_m_J2, dtype=float)
    return {
        **table.metadata,
        "runtime_table_key": "js_stiffness_A_per_m_J2",
        "diagnostic_current_key": "js_A_m2",
        "axis_keys": ["Te_axis_K", "delta2_axis_J2", "q_axis_m_inv"],
        "shape": list(table.shape),
        "js_max_abs_A_m2": float(np.max(np.abs(current))),
        "stiffness_min_A_per_m_J2": float(np.min(stiffness)),
        "stiffness_max_A_per_m_J2": float(np.max(stiffness)),
        "strict_required_by_temporal_runs": True,
    }


def _current_row_tasks(
    Te_axis: np.ndarray,
    delta_axis: np.ndarray,
    q_axis: np.ndarray,
    D: float,
    sigma: float,
    n_m: int,
) -> Iterable[tuple[int, int, float, float, np.ndarray, float, float, int]]:
    for iT, T in enumerate(Te_axis):
        for iD, delta in enumerate(delta_axis):
            yield (int(iT), int(iD), float(T), float(delta), q_axis, float(D), float(sigma), int(n_m))


def _compute_current_row(
    task: tuple[int, int, float, float, np.ndarray, float, float, int],
) -> tuple[int, int, np.ndarray, np.ndarray]:
    iT, iD, T, delta, q_axis, D, sigma, n_m = task
    current = np.zeros(q_axis.size, dtype=float)
    stiffness = np.empty(q_axis.size, dtype=float)
    eps = matsubara_energy_axis_J(T_K=float(T), n_matsubara=int(n_m))
    gamma_axis = 0.5 * HBAR_J_S * D * q_axis * q_axis
    prefactor = 2.0 * np.pi * K_B_J_K * float(T) * sigma / E_CHARGE_C
    for j, (q, gamma) in enumerate(zip(q_axis, gamma_axis)):
        if delta == 0.0:
            kappa = prefactor * float(np.sum(1.0 / np.square(eps + float(gamma))))
        else:
            s = solve_matsubara_s_values(
                delta_J=float(delta),
                gamma_J=float(gamma),
                eps_n_J=eps,
            )
            kappa = prefactor * float(np.sum(s * s)) / (float(delta) * float(delta))
        stiffness[j] = kappa
        current[j] = kappa * float(delta) * float(delta) * float(q)
    return iT, iD, current, stiffness


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
