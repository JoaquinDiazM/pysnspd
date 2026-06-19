"""
First spectral solver for the pySNSPD Usadel block.

OE3 first implementation:
- Provide a robust BCS/Dynes-like DOS with current-broadening proxy.
- Keep the API compatible with a future full dirty-limit Usadel solver.

The full real-axis quartic/self-consistent Usadel machinery will replace this
backend later without changing downstream catalogue consumers.
"""

from __future__ import annotations

import numpy as np


def dynes_bcs_dos(
    energy_J: np.ndarray,
    *,
    delta_J: float,
    gamma_J: float = 0.0,
    eta_J: float = 0.0,
) -> np.ndarray:
    """
    Compute a regularized normalized quasiparticle DOS.

    Parameters
    ----------
    energy_J:
        Energy axis in Joules.
    delta_J:
        Local superconducting gap magnitude in Joules.
    gamma_J:
        Depairing/current-broadening proxy in Joules.
    eta_J:
        Small numerical broadening in Joules.

    Returns
    -------
    numpy.ndarray
        Normalized DOS ``rho(E) = N(E)/N(0)``.

    Notes
    -----
    For ``gamma_J = eta_J = 0`` this approaches the BCS DOS. A finite
    ``gamma_J + eta_J`` smooths the gap edge and produces subgap tails. This is
    a first OE3 proxy for the current-smeared Usadel DOS.
    """
    E = np.asarray(energy_J, dtype=float)

    if delta_J < 0.0:
        raise ValueError("delta_J must be nonnegative.")

    if gamma_J < 0.0:
        raise ValueError("gamma_J must be nonnegative.")

    if eta_J < 0.0:
        raise ValueError("eta_J must be nonnegative.")

    if delta_J == 0.0:
        return np.ones_like(E, dtype=float)

    z = E + 1j * (eta_J + gamma_J)
    denominator = np.sqrt(z * z - delta_J * delta_J)

    rho = np.real(z / denominator)
    rho = np.nan_to_num(rho, nan=0.0, posinf=0.0, neginf=0.0)
    rho = np.maximum(rho, 0.0)

    return rho.astype(float)


def anomalous_proxy(
    energy_J: np.ndarray,
    *,
    delta_J: float,
    gamma_J: float = 0.0,
    eta_J: float = 0.0,
) -> np.ndarray:
    """
    Compute a regularized anomalous spectral proxy.

    This is stored for future phase-space catalogues. It is not yet the final
    Vodolazov/Simon anomalous spectral function.
    """
    E = np.asarray(energy_J, dtype=float)

    if delta_J <= 0.0:
        return np.zeros_like(E, dtype=float)

    z = E + 1j * (eta_J + gamma_J)
    denominator = np.sqrt(z * z - delta_J * delta_J)

    value = np.abs(delta_J / denominator)
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
    Compute DOS and anomalous-proxy grids.

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
            rho[i, j, :] = dynes_bcs_dos(
                energy,
                delta_J=float(delta),
                gamma_J=float(gamma),
                eta_J=float(eta_J),
            )
            anomalous[i, j, :] = anomalous_proxy(
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