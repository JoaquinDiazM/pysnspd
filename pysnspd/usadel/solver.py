"""
Uniform dirty-limit Usadel real-axis spectral solver.

OE3 v2:
- Replace the initial Dynes proxy by the real-axis quartic Usadel equation.
- Keep the same catalogue interface used by downstream modules.

We solve the uniform current/depairing spectral equation

    (i z - Gamma c) s + Delta c = 0,
    c^2 + s^2 = 1,

where

    c(E) = cos(theta(E)),
    s(E) = sin(theta(E)),
    z = E + i eta.

Eliminating s gives the quartic equation

    (Gamma c - i z)^2 (1 - c^2) - Delta^2 c^2 = 0.

The physical branch is selected by continuation from high energy, where
c -> 1 and rho(E) = Re[c(E)] -> 1.
"""

from __future__ import annotations

import numpy as np


def usadel_quartic_residual(
    c: complex,
    *,
    z_J: complex,
    delta_J: float,
    gamma_J: float,
) -> complex:
    """
    Residual of the quartic Usadel equation.

    Parameters
    ----------
    c:
        Complex value of cos(theta).
    z_J:
        Complex energy ``E + i eta`` in Joules.
    delta_J:
        Gap magnitude in Joules.
    gamma_J:
        Depairing energy in Joules.

    Returns
    -------
    complex
        Quartic residual.
    """
    a = gamma_J * c - 1j * z_J
    return a * a * (1.0 - c * c) - delta_J * delta_J * c * c


def usadel_quartic_derivative(
    c: complex,
    *,
    z_J: complex,
    delta_J: float,
    gamma_J: float,
) -> complex:
    """
    Derivative of the quartic residual with respect to ``c``.
    """
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
    """
    Analytic BCS/Dynes branch for Gamma = 0.

    This is the exact Gamma=0 limit of the Usadel equation, not the current
    broadened backend.
    """
    energy = np.asarray(energy_J, dtype=float)

    if delta_J < 0.0:
        raise ValueError("delta_J must be nonnegative.")

    if eta_J < 0.0:
        raise ValueError("eta_J must be nonnegative.")

    if delta_J == 0.0:
        return np.ones_like(energy, dtype=complex)

    z = energy + 1j * eta_J
    c = z / np.sqrt(z * z - delta_J * delta_J)

    # Choose the branch approaching +1 at high energy.
    high = c[-1]
    if np.real(high) < 0.0:
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
    """
    Solve the physical Usadel branch ``c(E)=cos(theta(E))``.

    The branch is followed from high to low energy. At high energy the
    physical solution satisfies ``c -> 1``. The previous energy solution is
    used as the initial guess for the next energy value.

    Parameters
    ----------
    energy_J:
        Ascending real energy axis.
    delta_J:
        Gap magnitude in Joules.
    gamma_J:
        Depairing energy in Joules.
    eta_J:
        Positive numerical broadening.
    newton_tol:
        Absolute Newton residual tolerance.
    max_iter:
        Maximum Newton iterations per energy point.

    Returns
    -------
    numpy.ndarray
        Complex array ``c(E)`` with the same shape as ``energy_J``.
    """
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
        return bcs_complex_cos_theta(
            energy,
            delta_J=delta_J,
            eta_J=eta_J,
        )

    out = np.empty_like(energy, dtype=complex)

    # Start at high energy from the Gamma=0 branch. This is close to +1.
    bcs_branch = bcs_complex_cos_theta(
        energy,
        delta_J=delta_J,
        eta_J=eta_J,
    )

    c_prev = bcs_branch[-1]
    if abs(c_prev - 1.0) > 0.5:
        c_prev = 1.0 + 0.0j

    for index in range(energy.size - 1, -1, -1):
        z = complex(energy[index], eta_J)
        c = c_prev

        converged = False
        for _ in range(max_iter):
            f = usadel_quartic_residual(
                c,
                z_J=z,
                delta_J=delta_J,
                gamma_J=gamma_J,
            )
            df = usadel_quartic_derivative(
                c,
                z_J=z,
                delta_J=delta_J,
                gamma_J=gamma_J,
            )

            if abs(f) < newton_tol * max(delta_J * delta_J, abs(z) * abs(z), 1.0e-300):
                converged = True
                break

            if abs(df) < 1.0e-300:
                break

            step = f / df

            # Mild damping avoids occasional branch jumps near the smeared gap edge.
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

        # Enforce the DOS branch with nonnegative real part when possible.
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
    """
    Robust fallback using explicit quartic roots.

    The selected root is the one closest to the previous physical branch value.
    """
    g = gamma_J
    z = z_J
    d = delta_J

    # Expanded polynomial:
    # g^2 c^4 - 2 i g z c^3 + (d^2 - g^2 - z^2)c^2
    # + 2 i g z c + z^2 = 0
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
    root = finite[int(np.argmin(distances))]

    return complex(root)


def usadel_dos(
    energy_J: np.ndarray,
    *,
    delta_J: float,
    gamma_J: float = 0.0,
    eta_J: float = 0.0,
) -> np.ndarray:
    """
    Compute the normalized DOS from the quartic Usadel branch.

    Returns

        rho(E) = Re[c(E)].
    """
    c = solve_usadel_cos_theta_branch(
        energy_J,
        delta_J=delta_J,
        gamma_J=gamma_J,
        eta_J=eta_J,
    )

    rho = np.real(c)
    rho = np.nan_to_num(rho, nan=0.0, posinf=0.0, neginf=0.0)
    rho = np.maximum(rho, 0.0)

    return rho.astype(float)


def usadel_anomalous_abs(
    energy_J: np.ndarray,
    *,
    delta_J: float,
    gamma_J: float = 0.0,
    eta_J: float = 0.0,
) -> np.ndarray:
    """
    Return ``abs(s(E))`` for the physical Usadel branch.

    This is stored as an anomalous spectral proxy for later phase-space
    kernels. The final kernel implementation may replace this by the exact
    spectral factors used in the appendix.
    """
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
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute DOS and anomalous grids.

    Returns arrays with shape

        ``(n_delta, n_gamma, n_energy)``.
    """
    energy = np.asarray(energy_values_J, dtype=float)
    deltas = np.asarray(delta_values_J, dtype=float)
    gammas = np.asarray(gamma_values_J, dtype=float)

    rho = np.empty((deltas.size, gammas.size, energy.size), dtype=float)
    anomalous = np.empty_like(rho)

    for i, delta in enumerate(deltas):
        for j, gamma in enumerate(gammas):
            rho[i, j, :] = usadel_dos(
                energy,
                delta_J=float(delta),
                gamma_J=float(gamma),
                eta_J=float(eta_J),
            )
            anomalous[i, j, :] = usadel_anomalous_abs(
                energy,
                delta_J=float(delta),
                gamma_J=float(gamma),
                eta_J=float(eta_J),
            )

    return rho, anomalous


def dos_diagnostics(
    rho_delta_gamma_E: np.ndarray,
    energy_values_J: np.ndarray,
) -> dict:
    """
    Return simple diagnostics for a DOS catalogue.
    """
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
def dynes_bcs_dos(
    energy_J: np.ndarray,
    *,
    delta_J: float,
    gamma_J: float = 0.0,
    eta_J: float = 0.0,
) -> np.ndarray:
    """
    Backward-compatible alias.

    In OE3 v2 this function now calls the quartic Usadel solver instead of the
    initial Dynes proxy.
    """
    return usadel_dos(
        energy_J,
        delta_J=delta_J,
        gamma_J=gamma_J,
        eta_J=eta_J,
    )


def anomalous_proxy(
    energy_J: np.ndarray,
    *,
    delta_J: float,
    gamma_J: float = 0.0,
    eta_J: float = 0.0,
) -> np.ndarray:
    """
    Backward-compatible alias for the anomalous spectral proxy.
    """
    return usadel_anomalous_abs(
        energy_J,
        delta_J=delta_J,
        gamma_J=gamma_J,
        eta_J=eta_J,
    )