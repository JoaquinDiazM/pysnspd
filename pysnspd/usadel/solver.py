"""Uniform dirty-limit Usadel real-axis spectral solver.

OE3 v2:
- Solve the real-axis dirty-limit quartic Usadel equation.
- Keep the catalogue interface used by downstream modules.
- Allow the independent ``(|Delta|, Gamma_q)`` catalogue cells to be
  evaluated in parallel for official PRE-runs.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from typing import Iterable

import numpy as np


def usadel_quartic_residual(
    c: complex,
    *,
    z_J: complex,
    delta_J: float,
    gamma_J: float,
) -> complex:
    """Residual of the real-axis uniform dirty-limit Usadel quartic."""

    a = gamma_J * c - 1j * z_J
    return a * a * (1.0 - c * c) - delta_J * delta_J * c * c


def usadel_quartic_derivative(
    c: complex,
    *,
    z_J: complex,
    delta_J: float,
    gamma_J: float,
) -> complex:
    """Derivative of :func:`usadel_quartic_residual` with respect to ``c``."""

    a = gamma_J * c - 1j * z_J
    return (
        2.0 * gamma_J * a * (1.0 - c * c)
        - 2.0 * c * a * a
        - 2.0 * delta_J * delta_J * c
    )


def bcs_complex_cos_theta(
    energy_J: np.ndarray,
    *,
    delta_J: float,
    eta_J: float,
) -> np.ndarray:
    """Analytic BCS/Dynes branch for the ``Gamma=0`` limit."""

    energy = np.asarray(energy_J, dtype=float)
    if delta_J < 0.0:
        raise ValueError("delta_J must be nonnegative.")
    if eta_J < 0.0:
        raise ValueError("eta_J must be nonnegative.")
    if delta_J == 0.0:
        return np.ones_like(energy, dtype=complex)

    z = energy + 1j * eta_J
    c = z / np.sqrt(z * z - delta_J * delta_J)
    if np.real(c[-1]) < 0.0:
        c = -c
    return c


def solve_usadel_cos_theta_branch(
    energy_J: np.ndarray,
    *,
    delta_J: float,
    gamma_J: float,
    eta_J: float,
    newton_tol: float = 1.0e-11,
    max_iter: int = 40,
) -> np.ndarray:
    """Solve the physical branch ``c(E)=cos(theta(E))`` by continuation."""

    energy = np.asarray(energy_J, dtype=float)
    if energy.ndim != 1:
        raise ValueError("energy_J must be one-dimensional.")
    if energy.size < 2:
        raise ValueError("energy_J must contain at least two points.")
    if np.any(np.diff(energy) < 0.0):
        raise ValueError("energy_J must be sorted in ascending order.")
    if delta_J < 0.0:
        raise ValueError("delta_J must be nonnegative.")
    if gamma_J < 0.0:
        raise ValueError("gamma_J must be nonnegative.")
    if eta_J < 0.0:
        raise ValueError("eta_J must be nonnegative.")
    if delta_J == 0.0:
        return np.ones_like(energy, dtype=complex)
    if gamma_J == 0.0:
        return bcs_complex_cos_theta(energy, delta_J=delta_J, eta_J=eta_J)

    out = np.empty_like(energy, dtype=complex)
    bcs_branch = bcs_complex_cos_theta(energy, delta_J=delta_J, eta_J=eta_J)
    c_prev = bcs_branch[-1]
    if abs(c_prev - 1.0) > 0.5:
        c_prev = 1.0 + 0.0j

    for index in range(energy.size - 1, -1, -1):
        z = complex(energy[index], eta_J)
        c = c_prev
        converged = False
        for _ in range(max_iter):
            f = usadel_quartic_residual(c, z_J=z, delta_J=delta_J, gamma_J=gamma_J)
            scale = max(delta_J * delta_J, abs(z) * abs(z), 1.0e-300)
            if abs(f) < newton_tol * scale:
                converged = True
                break
            df = usadel_quartic_derivative(c, z_J=z, delta_J=delta_J, gamma_J=gamma_J)
            if abs(df) < 1.0e-300:
                break
            step = f / df
            damping = 1.0
            c_candidate = c - damping * step
            while not np.isfinite(c_candidate.real) or not np.isfinite(c_candidate.imag):
                damping *= 0.5
                c_candidate = c - damping * step
                if damping < 1.0e-6:
                    break
            c = c_candidate
            if abs(step) < newton_tol * max(1.0, abs(c)):
                converged = True
                break
        if not converged:
            c = _fallback_quartic_root(
                z_J=z,
                delta_J=delta_J,
                gamma_J=gamma_J,
                previous_root=c_prev,
            )
        if np.real(c) < 0.0:
            c = -c
        out[index] = c
        c_prev = c
    return out


def _fallback_quartic_root(
    *,
    z_J: complex,
    delta_J: float,
    gamma_J: float,
    previous_root: complex,
) -> complex:
    """Robust fallback using explicit quartic roots."""

    g = gamma_J
    z = z_J
    d = delta_J
    coeffs = np.array(
        [
            g * g,
            -2.0j * g * z,
            d * d - g * g - z * z,
            2.0j * g * z,
            z * z,
        ],
        dtype=complex,
    )
    roots = np.roots(coeffs)
    finite = roots[np.isfinite(roots)]
    if finite.size == 0:
        return previous_root
    distances = np.abs(finite - previous_root)
    return complex(finite[int(np.argmin(distances))])


def usadel_dos(
    energy_J: np.ndarray,
    *,
    delta_J: float,
    gamma_J: float = 0.0,
    eta_J: float = 0.0,
) -> np.ndarray:
    """Compute normalized DOS ``rho(E)=Re[cos(theta(E))]``."""

    c = solve_usadel_cos_theta_branch(
        energy_J,
        delta_J=delta_J,
        gamma_J=gamma_J,
        eta_J=eta_J,
    )
    rho = np.real(c)
    rho = np.nan_to_num(rho, nan=0.0, posinf=0.0, neginf=0.0)
    return np.maximum(rho, 0.0).astype(float)


def usadel_anomalous_abs(
    energy_J: np.ndarray,
    *,
    delta_J: float,
    gamma_J: float = 0.0,
    eta_J: float = 0.0,
) -> np.ndarray:
    """Return ``|s(E)|`` for the physical Usadel branch."""

    energy = np.asarray(energy_J, dtype=float)
    if delta_J == 0.0:
        return np.zeros_like(energy, dtype=float)
    c = solve_usadel_cos_theta_branch(
        energy,
        delta_J=delta_J,
        gamma_J=gamma_J,
        eta_J=eta_J,
    )
    z = energy + 1j * eta_J
    denominator = 1j * z - gamma_J * c
    with np.errstate(divide="ignore", invalid="ignore"):
        s = -delta_J * c / denominator
    value = np.abs(s)
    value = np.nan_to_num(value, nan=0.0, posinf=0.0, neginf=0.0)
    return value.astype(float)


def compute_dos_grid(
    energy_values_J: np.ndarray,
    delta_values_J: np.ndarray,
    gamma_values_J: np.ndarray,
    *,
    eta_J: float,
    workers: int = 1,
    backend: str = "process",
) -> tuple[np.ndarray, np.ndarray]:
    """Compute DOS and anomalous grids.

    Returns arrays with shape ``(n_delta, n_gamma, n_energy)``.  Catalogue
    cells for different ``|Delta|`` rows are independent, so official PRE-runs
    can evaluate the rows in parallel using the same worker policy as the
    current table builder.
    """

    energy = np.asarray(energy_values_J, dtype=float)
    deltas = np.asarray(delta_values_J, dtype=float)
    gammas = np.asarray(gamma_values_J, dtype=float)
    if energy.ndim != 1 or deltas.ndim != 1 or gammas.ndim != 1:
        raise ValueError("energy_values_J, delta_values_J and gamma_values_J must be 1D arrays.")

    rho = np.empty((deltas.size, gammas.size, energy.size), dtype=float)
    anomalous = np.empty_like(rho)

    n_workers = max(1, int(workers))
    mode = str(backend or "process").lower()
    if mode == "serial":
        n_workers = 1
    if n_workers == 1 or deltas.size <= 1:
        for task in _dos_row_tasks(energy, deltas, gammas, float(eta_J)):
            i, row_rho, row_anom = _compute_dos_row(task)
            rho[i, :, :] = row_rho
            anomalous[i, :, :] = row_anom
        return rho, anomalous

    executor_cls = ThreadPoolExecutor if mode == "thread" else ProcessPoolExecutor
    with executor_cls(max_workers=n_workers) as pool:
        for i, row_rho, row_anom in pool.map(_compute_dos_row, _dos_row_tasks(energy, deltas, gammas, float(eta_J))):
            rho[i, :, :] = row_rho
            anomalous[i, :, :] = row_anom
    return rho, anomalous


def _dos_row_tasks(
    energy: np.ndarray,
    deltas: np.ndarray,
    gammas: np.ndarray,
    eta_J: float,
) -> Iterable[tuple[int, float, np.ndarray, np.ndarray, float]]:
    for i, delta in enumerate(deltas):
        yield (int(i), float(delta), energy, gammas, float(eta_J))


def _compute_dos_row(task: tuple[int, float, np.ndarray, np.ndarray, float]) -> tuple[int, np.ndarray, np.ndarray]:
    i, delta, energy, gammas, eta_J = task
    row_rho = np.empty((gammas.size, energy.size), dtype=float)
    row_anom = np.empty_like(row_rho)
    for j, gamma in enumerate(gammas):
        row_rho[j, :] = usadel_dos(
            energy,
            delta_J=float(delta),
            gamma_J=float(gamma),
            eta_J=float(eta_J),
        )
        row_anom[j, :] = usadel_anomalous_abs(
            energy,
            delta_J=float(delta),
            gamma_J=float(gamma),
            eta_J=float(eta_J),
        )
    return i, row_rho, row_anom


def dos_diagnostics(
    rho_delta_gamma_E: np.ndarray,
    energy_values_J: np.ndarray,
) -> dict:
    """Return simple diagnostics for a DOS catalogue."""

    rho = np.asarray(rho_delta_gamma_E, dtype=float)
    energy = np.asarray(energy_values_J, dtype=float)
    if rho.ndim != 3:
        raise ValueError("rho_delta_gamma_E must have shape (n_delta, n_gamma, n_energy).")
    if energy.ndim != 1:
        raise ValueError("energy_values_J must be one-dimensional.")
    if rho.shape[-1] != energy.size:
        raise ValueError("DOS energy dimension does not match energy axis.")
    return {
        "rho_min": float(np.min(rho)),
        "rho_max": float(np.max(rho)),
        "rho_mean": float(np.mean(rho)),
        "rho_is_finite": bool(np.all(np.isfinite(rho))),
        "energy_min_J": float(np.min(energy)),
        "energy_max_J": float(np.max(energy)),
    }


# Backward-compatible names from OE3 v1.
