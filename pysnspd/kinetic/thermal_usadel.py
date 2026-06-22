"""Parallel thermal-Usadel grid helpers for PRE-run.

This module moves the expensive OE5 thermal Usadel audit into the PRE-run.
The heavy independent work is the recomputation of

    Delta_eq(Te, q)

over a temperature/q grid. This is naturally parallel over Te rows.

The resulting grid is saved under raw/<run_name>/pre/ so later SS/PHOTON
runs can load it instead of solving the Matsubara self-consistency equation
again.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np

from pysnspd.kinetic.powers import (
    HBAR_J_S,
    KB_J_K,
    E_CHARGE_C,
    MEV_J,
    ThermalUsadelGrid,
    build_thermal_usadel_grid,
)
from pysnspd.usadel.calibration import (
    matsubara_energy_axis_J,
    solve_gap_for_gamma_J,
    solve_matsubara_s_values,
)


def build_thermal_usadel_grid_parallel(
    usadel_catalog,
    Te_values_K: np.ndarray,
    *,
    n_q: int = 140,
    n_matsubara: int = 500,
    stable_lowT_branch_only: bool = True,
    workers: int = 1,
) -> ThermalUsadelGrid:
    """Build Delta_eq(Te,q), optionally in parallel.

    Parameters
    ----------
    workers:
        Number of process workers. If workers <= 1, this delegates to the
        serial implementation in ``pysnspd.kinetic.powers``.
    """
    workers = int(workers)
    if workers <= 1:
        return build_thermal_usadel_grid(
            usadel_catalog,
            Te_values_K,
            n_q=n_q,
            n_matsubara=n_matsubara,
            stable_lowT_branch_only=stable_lowT_branch_only,
        )

    prepared = _prepare_q_axis(
        usadel_catalog,
        n_q=n_q,
        stable_lowT_branch_only=stable_lowT_branch_only,
    )

    Te_values = np.asarray(Te_values_K, dtype=float)
    q_values = prepared["q_values_m_inv"]
    gamma_values = prepared["gamma_values_J"]

    metadata = dict(prepared["metadata"])
    width_m = float(metadata["width_m"])
    thickness_m = float(metadata["thickness_m"])
    area_m2 = width_m * thickness_m

    payloads = []
    for Te in Te_values:
        payloads.append(
            {
                "Te_K": float(Te),
                "q_values_m_inv": q_values,
                "gamma_values_J": gamma_values,
                "D_m2_s": float(metadata["D_m2_s"]),
                "sigma_n_S_m": float(metadata["sigma_n_S_m"]),
                "area_m2": area_m2,
                "Tc_K": float(metadata["Tc_K"]),
                "n_matsubara": int(n_matsubara),
            }
        )

    with ProcessPoolExecutor(max_workers=workers) as executor:
        rows = list(executor.map(_compute_thermal_usadel_row, payloads))

    delta_grid = np.vstack([row["delta_values_J"] for row in rows])
    current_grid = np.vstack([row["current_values_A"] for row in rows])
    current_density_grid = np.vstack([row["current_density_values_A_m2"] for row in rows])

    current_fraction_grid = np.zeros_like(current_grid)
    for i in range(current_grid.shape[0]):
        row_max = float(np.max(current_grid[i, :]))
        if row_max > 0.0:
            current_fraction_grid[i, :] = current_grid[i, :] / row_max

    metadata.update(
        {
            "backend": "thermal_usadel_grid_oe5_v5_parallel",
            "description": (
                "Delta_eq(Te,q) recomputed with the OE3 Matsubara Usadel "
                "self-consistency solver. This grid is generated during PRE-run "
                "and can be reused by SS/PHOTON runs."
            ),
            "n_Te": int(Te_values.size),
            "n_q": int(q_values.size),
            "n_matsubara": int(n_matsubara),
            "workers": int(workers),
            "parallel_axis": "Te rows",
            "stable_lowT_branch_only": bool(stable_lowT_branch_only),
        }
    )

    return ThermalUsadelGrid(
        Te_values_K=Te_values,
        q_values_m_inv=q_values,
        gamma_values_J=gamma_values,
        delta_eq_Tq_J=delta_grid,
        current_Tq_A=current_grid,
        current_density_Tq_A_m2=current_density_grid,
        current_fraction_Tq=current_fraction_grid,
        reference_current_fraction=prepared["reference_current_fraction"],
        metadata=metadata,
    )


def save_thermal_usadel_grid_npz(grid: ThermalUsadelGrid, path: str | Path) -> Path:
    """Save the thermal Usadel grid to a compressed NPZ file."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)

    np.savez_compressed(
        output,
        Te_values_K=grid.Te_values_K,
        q_values_m_inv=grid.q_values_m_inv,
        gamma_values_J=grid.gamma_values_J,
        delta_eq_Tq_J=grid.delta_eq_Tq_J,
        current_Tq_A=grid.current_Tq_A,
        current_density_Tq_A_m2=grid.current_density_Tq_A_m2,
        current_fraction_Tq=grid.current_fraction_Tq,
        reference_current_fraction=grid.reference_current_fraction,
        metadata=np.array(grid.metadata, dtype=object),
    )
    return output


def thermal_usadel_grid_summary(grid: ThermalUsadelGrid) -> dict[str, Any]:
    """Compact manifest-friendly summary."""
    return {
        "backend": str(grid.metadata.get("backend", "unknown")),
        "shape": [int(grid.Te_values_K.size), int(grid.q_values_m_inv.size)],
        "n_Te": int(grid.Te_values_K.size),
        "n_q": int(grid.q_values_m_inv.size),
        "workers": int(grid.metadata.get("workers", 1)),
        "Te_min_K": float(np.min(grid.Te_values_K)),
        "Te_max_K": float(np.max(grid.Te_values_K)),
        "q_min_m_inv": float(np.min(grid.q_values_m_inv)),
        "q_max_m_inv": float(np.max(grid.q_values_m_inv)),
        "gamma_min_meV": float(np.min(grid.gamma_values_J) / MEV_J),
        "gamma_max_meV": float(np.max(grid.gamma_values_J) / MEV_J),
        "delta_min_meV": float(np.min(grid.delta_eq_Tq_J) / MEV_J),
        "delta_max_meV": float(np.max(grid.delta_eq_Tq_J) / MEV_J),
        "delta_is_finite": bool(np.all(np.isfinite(grid.delta_eq_Tq_J))),
        "current_is_finite": bool(np.all(np.isfinite(grid.current_Tq_A))),
        "Tc_K": float(grid.metadata.get("Tc_K", np.nan)),
        "T_bias_K": float(grid.metadata.get("T_bias_K", np.nan)),
        "Ic_bias_A": float(grid.metadata.get("Ic_bias_A", np.nan)),
        "q_c_bias_m_inv": float(grid.metadata.get("q_c_bias_m_inv", np.nan)),
    }


def _prepare_q_axis(
    usadel_catalog,
    *,
    n_q: int,
    stable_lowT_branch_only: bool,
) -> dict[str, Any]:
    D_m2_s = float(usadel_catalog.metadata["D_m2_s"])
    sigma_n = float(usadel_catalog.metadata["sigma_n_S_m"])
    width_m = float(usadel_catalog.metadata["width_m"])
    thickness_m = float(usadel_catalog.metadata["thickness_m"])
    Tc_K = float(usadel_catalog.metadata["Tc_K"])
    T_bias_K = float(usadel_catalog.metadata["T_bias_K"])

    q_cal = np.asarray(usadel_catalog.calibration_q_values_m_inv, dtype=float)
    current_cal = np.asarray(usadel_catalog.calibration_current_values_A, dtype=float)

    valid = np.isfinite(q_cal) & np.isfinite(current_cal)
    valid &= q_cal >= 0.0
    valid &= current_cal >= 0.0
    if np.sum(valid) < 3:
        raise ValueError("Invalid OE3 calibration branch.")

    q_cal = q_cal[valid]
    current_cal = current_cal[valid]
    order = np.argsort(q_cal)
    q_cal = q_cal[order]
    current_cal = current_cal[order]

    idx_ic = int(np.argmax(current_cal))
    Ic_bias_A = float(current_cal[idx_ic])
    q_c_bias = float(q_cal[idx_ic])

    if stable_lowT_branch_only:
        q_min = 0.0
        q_max = q_c_bias
    else:
        q_min = float(np.min(q_cal))
        q_max = float(np.max(q_cal))

    q_values = np.linspace(q_min, q_max, int(n_q))
    gamma_values = 0.5 * HBAR_J_S * D_m2_s * q_values * q_values

    current_reference = np.interp(q_values, q_cal, current_cal)
    reference_fraction = current_reference / Ic_bias_A if Ic_bias_A > 0.0 else np.zeros_like(q_values)

    metadata = {
        "Tc_K": Tc_K,
        "T_bias_K": T_bias_K,
        "D_m2_s": D_m2_s,
        "sigma_n_S_m": sigma_n,
        "width_m": width_m,
        "thickness_m": thickness_m,
        "Ic_bias_A": Ic_bias_A,
        "q_c_bias_m_inv": q_c_bias,
    }

    return {
        "q_values_m_inv": q_values,
        "gamma_values_J": gamma_values,
        "reference_current_fraction": reference_fraction,
        "metadata": metadata,
    }


def _compute_thermal_usadel_row(payload: dict[str, Any]) -> dict[str, np.ndarray]:
    Te_K = float(payload["Te_K"])
    q_values = np.asarray(payload["q_values_m_inv"], dtype=float)
    gamma_values = np.asarray(payload["gamma_values_J"], dtype=float)

    sigma_n = float(payload["sigma_n_S_m"])
    area_m2 = float(payload["area_m2"])
    Tc_K = float(payload["Tc_K"])
    n_matsubara = int(payload["n_matsubara"])

    eps_n = matsubara_energy_axis_J(T_K=Te_K, n_matsubara=n_matsubara)

    delta_values = np.zeros_like(q_values)
    current_values = np.zeros_like(q_values)
    current_density_values = np.zeros_like(q_values)

    for i, (q, gamma) in enumerate(zip(q_values, gamma_values, strict=True)):
        delta = solve_gap_for_gamma_J(
            gamma_J=float(gamma),
            T_K=Te_K,
            Tc_K=Tc_K,
            eps_n_J=eps_n,
        )
        delta_values[i] = delta

        if delta > 0.0 and q > 0.0:
            s = solve_matsubara_s_values(
                delta_J=float(delta),
                gamma_J=float(gamma),
                eps_n_J=eps_n,
            )
            sum_s2 = float(np.sum(s * s))
            current = (
                area_m2
                * (2.0 * np.pi * KB_J_K * Te_K / E_CHARGE_C)
                * sigma_n
                * float(q)
                * sum_s2
            )
        else:
            current = 0.0

        current_values[i] = current
        current_density_values[i] = current / area_m2

    return {
        "delta_values_J": delta_values,
        "current_values_A": current_values,
        "current_density_values_A_m2": current_density_values,
    }