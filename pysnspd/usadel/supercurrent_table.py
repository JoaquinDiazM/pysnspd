"""Dirty-limit Usadel supercurrent-density table for PRE-run diagnostics.

This module belongs to the PRE stage.  It appends an explicit
``j_s^{Usadel}(Delta, q)`` table to ``usadel_dos_catalog.npz`` so SS diagnostics
interpolate a real PRE-computed Usadel quantity instead of trying to infer one
inside the time-stepper.

All public inputs/outputs are SI:

* ``Delta`` [J]
* ``q``     [m^-1], the phase-gradient wave number used by pySNSPD
* ``T``     [K]
* ``j_s``   [A m^-2]

The dirty-limit Matsubara current relation is evaluated as

    j_s = (2 pi k_B T / |e|) sigma_n q sum_n sin^2(theta_n),

where pySNSPD's ``q`` is the wave number.  This is equivalent to the common
``(2 pi k_B T / |e| hbar) sigma_n Q`` form with ``Q = hbar q``.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping
import sys

import numpy as np
import yaml

HBAR_J_S = 1.054571817e-34
KB_J_K = 1.380649e-23
QE_C = 1.602176634e-19


@dataclass(frozen=True)
class UsadelSupercurrentTableSummary:
    backend: str
    table_shape: tuple[int, int]
    T_K: float
    n_matsubara: int
    q_min_m_inv: float
    q_max_m_inv: float
    delta_min_J: float
    delta_max_J: float
    js_min_A_m2: float
    js_max_A_m2: float
    js_abs_max_A_m2: float
    workers: int = 1

    def as_dict(self) -> dict[str, Any]:
        return {
            "backend": self.backend,
            "table_shape": list(self.table_shape),
            "T_K": float(self.T_K),
            "n_matsubara": int(self.n_matsubara),
            "q_min_m_inv": float(self.q_min_m_inv),
            "q_max_m_inv": float(self.q_max_m_inv),
            "delta_min_J": float(self.delta_min_J),
            "delta_max_J": float(self.delta_max_J),
            "js_min_A_m2": float(self.js_min_A_m2),
            "js_max_A_m2": float(self.js_max_A_m2),
            "js_abs_max_A_m2": float(self.js_abs_max_A_m2),
            "workers": int(self.workers),
            "units_policy": "SI; q axis is wave number [m^-1], current density is [A m^-2].",
        }


def append_usadel_supercurrent_table_to_npz(
    usadel_npz: str | Path,
    *,
    config: Mapping[str, Any],
    summary_yaml: str | Path | None = None,
    output_summary_yaml: str | Path | None = None,
    n_matsubara: int | None = None,
    workers: int = 1,
    progress: bool = False,
) -> UsadelSupercurrentTableSummary:
    """Append ``js_A_m2`` to a PRE-run Usadel catalogue in-place.

    The existing Usadel DOS catalogue may contain object arrays used by older
    PRE metadata.  Repacking the catalogue therefore loads the existing archive
    with ``allow_pickle=True`` and writes it back with the new numeric current
    table.  The SS diagnostic loader remains strict and reads only the numeric
    arrays it needs.
    """
    usadel_npz = Path(usadel_npz)
    if not usadel_npz.exists():
        raise FileNotFoundError(usadel_npz)

    summary = _load_yaml(summary_yaml) if summary_yaml is not None else {}
    with np.load(usadel_npz, allow_pickle=True) as data:
        arrays = {key: data[key] for key in data.files}

    q_axis = _axis_from_arrays_or_summary(
        arrays,
        summary,
        names=("q_axis_m_inv", "q_values_m_inv", "q_grid_m_inv", "q_m_inv", "q_axis", "q_values"),
        min_key="q_min_m_inv",
        max_key="q_max_m_inv",
        n_key="n_q",
    )
    delta_axis = _axis_from_arrays_or_summary(
        arrays,
        summary,
        names=("delta_axis_J", "delta_values_J", "Delta_axis_J", "Delta_values_J", "delta_J", "Delta_J"),
        min_key="delta_min_J",
        max_key="delta_max_J",
        n_key="n_delta",
    )

    T_K = _bias_temperature_K(config, summary)
    D_m2_s = _material_value(config, summary, "D_m2_s")
    sigma_n_S_m = _material_value(config, summary, "sigma_n_S_m", aliases=("sigma_n", "sigma_S_m"))
    nmats = int(n_matsubara or _summary_int(summary, "n_matsubara", default=500))
    nworkers = max(1, int(workers or 1))

    js = compute_dirty_usadel_supercurrent_table(
        delta_axis_J=delta_axis,
        q_axis_m_inv=q_axis,
        T_K=T_K,
        D_m2_s=D_m2_s,
        sigma_n_S_m=sigma_n_S_m,
        n_matsubara=nmats,
        workers=nworkers,
        progress=progress,
    )

    arrays["q_axis_m_inv"] = np.asarray(q_axis, dtype=float)
    arrays["delta_axis_J"] = np.asarray(delta_axis, dtype=float)
    arrays["js_A_m2"] = np.asarray(js, dtype=float)
    arrays["js_table_T_K"] = np.array([float(T_K)], dtype=float)
    arrays["js_table_n_matsubara"] = np.array([int(nmats)], dtype=np.int64)
    arrays["js_table_workers"] = np.array([int(nworkers)], dtype=np.int64)
    arrays["js_table_axes"] = np.array(["delta_axis_J", "q_axis_m_inv"], dtype="U32")
    arrays["js_table_backend"] = np.array(["dirty_usadel_matsubara_current_relation_v1"], dtype="U64")

    np.savez_compressed(usadel_npz, **arrays)

    out = UsadelSupercurrentTableSummary(
        backend="dirty_usadel_matsubara_current_relation_v1",
        table_shape=tuple(int(v) for v in js.shape),
        T_K=T_K,
        n_matsubara=nmats,
        q_min_m_inv=float(np.nanmin(q_axis)),
        q_max_m_inv=float(np.nanmax(q_axis)),
        delta_min_J=float(np.nanmin(delta_axis)),
        delta_max_J=float(np.nanmax(delta_axis)),
        js_min_A_m2=float(np.nanmin(js)),
        js_max_A_m2=float(np.nanmax(js)),
        js_abs_max_A_m2=float(np.nanmax(np.abs(js))),
        workers=nworkers,
    )
    if output_summary_yaml is not None:
        output_summary_yaml = Path(output_summary_yaml)
        with output_summary_yaml.open("w", encoding="utf-8") as f:
            yaml.safe_dump(out.as_dict(), f, sort_keys=False)
    return out


def load_usadel_current_diagnostic_catalog(usadel_npz: str | Path) -> dict[str, np.ndarray]:
    """Load only the numeric arrays required by SS Usadel-current diagnostics."""
    usadel_npz = Path(usadel_npz)
    required = ("js_A_m2", "q_axis_m_inv", "delta_axis_J")
    with np.load(usadel_npz, allow_pickle=False) as data:
        missing = [key for key in required if key not in data.files]
        if missing:
            raise KeyError(
                "Usadel supercurrent diagnostics require PRE-computed arrays "
                f"{required}; missing {missing}. Rerun 01_prerun_pytdgl_like_template.py."
            )
        arrays = {key: np.asarray(data[key], dtype=float) for key in required}
        for optional in ("js_table_T_K", "js_table_n_matsubara", "js_table_workers"):
            if optional in data.files:
                arrays[optional] = data[optional]
    return arrays


def compute_dirty_usadel_supercurrent_table(
    *,
    delta_axis_J: np.ndarray,
    q_axis_m_inv: np.ndarray,
    T_K: float,
    D_m2_s: float,
    sigma_n_S_m: float,
    n_matsubara: int = 500,
    max_iter: int = 80,
    tol: float = 1.0e-13,
    workers: int = 1,
    progress: bool = False,
) -> np.ndarray:
    """Compute ``j_s^{Usadel}(Delta, q)`` on a tensor-product grid.

    Parallelization is over Delta rows.  This keeps the implementation
    deterministic and avoids changing the Usadel equations themselves.
    """
    delta_axis = np.asarray(delta_axis_J, dtype=float).reshape(-1)
    q_axis = np.asarray(q_axis_m_inv, dtype=float).reshape(-1)
    if delta_axis.size < 1:
        raise ValueError("delta_axis_J must contain at least one point.")
    if q_axis.size < 1:
        raise ValueError("q_axis_m_inv must contain at least one point.")
    if T_K <= 0:
        raise ValueError("T_K must be positive.")
    if D_m2_s <= 0 or sigma_n_S_m <= 0:
        raise ValueError("D_m2_s and sigma_n_S_m must be positive.")
    n_matsubara = int(n_matsubara)
    if n_matsubara <= 0:
        raise ValueError("n_matsubara must be positive.")

    table = np.empty((delta_axis.size, q_axis.size), dtype=float)
    nworkers = max(1, int(workers or 1))
    task_args = [
        (
            int(idelta),
            float(delta),
            q_axis,
            float(T_K),
            float(D_m2_s),
            float(sigma_n_S_m),
            int(n_matsubara),
            int(max_iter),
            float(tol),
        )
        for idelta, delta in enumerate(delta_axis)
    ]
    progress_bar = _SimpleProgress("Usadel supercurrent table", len(task_args), enabled=bool(progress))
    try:
        if nworkers == 1 or len(task_args) <= 1:
            for args in task_args:
                idelta, row = _compute_dirty_usadel_supercurrent_row(args)
                table[idelta, :] = row
                progress_bar.update()
        else:
            with ProcessPoolExecutor(max_workers=min(nworkers, len(task_args))) as executor:
                futures = [executor.submit(_compute_dirty_usadel_supercurrent_row, args) for args in task_args]
                for future in as_completed(futures):
                    idelta, row = future.result()
                    table[idelta, :] = row
                    progress_bar.update()
    finally:
        progress_bar.close()
    return table


def _compute_dirty_usadel_supercurrent_row(args: tuple[Any, ...]) -> tuple[int, np.ndarray]:
    (
        idelta,
        delta,
        q_axis,
        T_K,
        D_m2_s,
        sigma_n_S_m,
        n_matsubara,
        max_iter,
        tol,
    ) = args
    q_axis = np.asarray(q_axis, dtype=float).reshape(-1)
    row = np.empty(q_axis.size, dtype=float)
    delta = max(float(delta), 0.0)
    if delta == 0.0:
        row[:] = 0.0
        return int(idelta), row

    omega = (2 * np.arange(int(n_matsubara), dtype=float) + 1.0) * np.pi * KB_J_K * float(T_K)
    prefactor = 2.0 * np.pi * KB_J_K * float(T_K) * float(sigma_n_S_m) / QE_C
    for iq, q in enumerate(q_axis):
        q = float(q)
        gamma_J = 0.5 * HBAR_J_S * float(D_m2_s) * q * q
        sin_theta = _solve_usadel_sin_theta(delta, omega, gamma_J, max_iter=max_iter, tol=tol)
        row[iq] = prefactor * q * float(np.sum(sin_theta * sin_theta))
    return int(idelta), row


class _SimpleProgress:
    def __init__(self, label: str, total: int, *, enabled: bool) -> None:
        self.label = str(label)
        self.total = max(1, int(total))
        self.enabled = bool(enabled)
        self.count = 0
        self._use_tqdm = False
        self._bar = None
        if self.enabled:
            try:
                from tqdm import tqdm  # type: ignore

                self._bar = tqdm(total=self.total, desc=self.label, unit="Delta", dynamic_ncols=True)
                self._use_tqdm = True
            except Exception:
                self._print()

    def update(self, n: int = 1) -> None:
        if not self.enabled:
            return
        self.count = min(self.total, self.count + int(n))
        if self._use_tqdm and self._bar is not None:
            self._bar.update(int(n))
        else:
            self._print()

    def close(self) -> None:
        if not self.enabled:
            return
        if self._use_tqdm and self._bar is not None:
            self._bar.close()
        else:
            sys.stderr.write("\n")
            sys.stderr.flush()

    def _print(self) -> None:
        sys.stderr.write(f"\r{self.label}: {self.count}/{self.total}")
        sys.stderr.flush()


def _solve_usadel_sin_theta(
    delta_J: float,
    omega_J: np.ndarray,
    gamma_J: float,
    *,
    max_iter: int,
    tol: float,
) -> np.ndarray:
    """Solve the dirty Usadel theta equation for ``sin(theta_n)``.

    Equation used at fixed ``Delta`` and depairing energy ``Gamma``:

        Delta cos(theta_n) - (omega_n + Gamma cos(theta_n)) sin(theta_n) = 0.

    A bounded fixed-point iteration is stable for the PRE-grid regime used here.
    """
    omega = np.asarray(omega_J, dtype=float)
    s = delta_J / np.sqrt((omega + gamma_J) ** 2 + delta_J**2)
    s = np.clip(s, 0.0, 1.0)
    for _ in range(int(max_iter)):
        c = np.sqrt(np.maximum(1.0 - s * s, 0.0))
        new_s = delta_J / np.sqrt((omega + gamma_J * c) ** 2 + delta_J**2)
        new_s = np.clip(new_s, 0.0, 1.0)
        if float(np.max(np.abs(new_s - s))) < tol:
            s = new_s
            break
        s = new_s
    return s


def _load_yaml(path: str | Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    path = Path(path)
    if not path.exists():
        return {}
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    return data if isinstance(data, dict) else {}


def _axis_from_arrays_or_summary(
    arrays: Mapping[str, np.ndarray],
    summary: Mapping[str, Any],
    *,
    names: tuple[str, ...],
    min_key: str,
    max_key: str,
    n_key: str,
) -> np.ndarray:
    for name in names:
        if name in arrays:
            axis = np.asarray(arrays[name], dtype=float).reshape(-1)
            axis = axis[np.isfinite(axis)]
            if axis.size >= 1:
                return axis
    if min_key in summary and max_key in summary and n_key in summary:
        return np.linspace(float(summary[min_key]), float(summary[max_key]), int(summary[n_key]))
    raise KeyError(
        f"Could not resolve axis {names}; catalogue must contain one of these keys "
        f"or summary must contain {min_key}, {max_key}, {n_key}."
    )


def _bias_temperature_K(config: Mapping[str, Any], summary: Mapping[str, Any]) -> float:
    bias = config.get("bias", {}) if isinstance(config.get("bias", {}), Mapping) else {}
    if "T_bias_K" in bias:
        return float(bias["T_bias_K"])
    for key in ("T_bias_K", "Te_min_K"):
        if key in summary:
            return float(summary[key])
    raise KeyError("Could not resolve bias temperature T_bias_K for Usadel current table.")


def _material_value(
    config: Mapping[str, Any],
    summary: Mapping[str, Any],
    key: str,
    *,
    aliases: tuple[str, ...] = (),
) -> float:
    material = config.get("material", {}) if isinstance(config.get("material", {}), Mapping) else {}
    for name in (key, *aliases):
        if name in summary:
            return float(summary[name])
        if name in material:
            return float(material[name])
    raise KeyError(f"Could not resolve material parameter {key}.")


def _summary_int(summary: Mapping[str, Any], key: str, *, default: int) -> int:
    try:
        return int(summary.get(key, default))
    except Exception:
        return int(default)
