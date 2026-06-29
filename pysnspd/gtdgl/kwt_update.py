"""Local KWT update rules for OE7 stationary relaxation."""
from __future__ import annotations

import numpy as np

from pysnspd.gtdgl.material import E_CHARGE_C, HBAR_J_S, GTDGLMaterial
from pysnspd.gtdgl.state import FormulaFields

VALID_PHI_PHASE_POLICIES = {"plus", "none", "minus"}


def kwt_delta_update_attempt(
    *,
    psi_J: np.ndarray,
    phi_V: np.ndarray,
    defs: FormulaFields,
    dt_s: float,
    material: GTDGLMaterial,
    phi_phase_policy: str = "plus",
) -> tuple[np.ndarray | None, float]:
    """Notebook local semi-implicit KWT update attempt.

    Solves locally
        Delta^{n+1} + z |Delta^{n+1}|^2 = w
    with a selectable temporal gauge link. The default ``plus`` preserves the
    previous OE7 convention,

        U = exp(+i 2e varphi dt / hbar),

    and the local algebra uses ``Uinv = conj(U)``. The diagnostic policies are
    ``none`` for no electrostatic phase link and ``minus`` for the opposite
    sign convention.
    """
    tau0 = material.tau0_GL_s
    Delta_n = np.asarray(psi_J, dtype=np.complex128)
    phi = np.asarray(phi_V, dtype=float)
    amp2_n = np.abs(Delta_n) ** 2

    phi_phase_policy = str(phi_phase_policy)
    if phi_phase_policy not in VALID_PHI_PHASE_POLICIES:
        raise ValueError(
            "phi_phase_policy must be one of "
            f"{sorted(VALID_PHI_PHASE_POLICIES)}, got {phi_phase_policy!r}."
        )

    phase = (2.0 * E_CHARGE_C / HBAR_J_S) * phi * float(dt_s)
    if phi_phase_policy == "plus":
        Uinv = np.exp(-1j * phase)
    elif phi_phase_policy == "minus":
        Uinv = np.exp(+1j * phase)
    else:  # phi_phase_policy == "none"
        Uinv = np.ones_like(phi, dtype=np.complex128)

    alpha = defs.alpha_kwt_J_inv2
    z = alpha * Uinv * Delta_n
    w = Uinv * (
        Delta_n
        + alpha * Delta_n * amp2_n
        + (float(dt_s) / tau0) * defs.rho * defs.forcing_J
    )

    ccoef = np.real(w * np.conjugate(z))
    absz2 = np.abs(z) ** 2
    absw2 = np.abs(w) ** 2
    B = 1.0 + 2.0 * ccoef
    discr = B**2 - 4.0 * absz2 * absw2
    discr_min = float(np.nanmin(discr))
    if discr_min < -1.0e-14:
        return None, discr_min

    discr = np.maximum(discr, 0.0)
    denom = B + np.sqrt(discr)
    denom = np.where(np.abs(denom) < 1.0e-300, np.inf, denom)
    amp2_new = 2.0 * absw2 / denom
    amp2_new = np.maximum(amp2_new, 0.0)
    Delta_new = w - z * amp2_new

    if not (
        np.all(np.isfinite(np.real(Delta_new)))
        and np.all(np.isfinite(np.imag(Delta_new)))
    ):
        return None, discr_min
    return Delta_new, discr_min


def kwt_local_update(
    *,
    psi_J: np.ndarray,
    phi_V: np.ndarray,
    Te_K: np.ndarray,
    forcing_J: np.ndarray,
    dt_s: float,
    material: GTDGLMaterial,
    max_phase_step_rad: float = 0.25,
    use_phi_phase: bool = False,
) -> tuple[np.ndarray, bool, float]:
    """Backward-compatible local KWT update used by legacy tests.

    The production OE7 notebook-port path uses ``kwt_delta_update_attempt``
    because it needs the full notebook formula-field bundle. This wrapper keeps
    the old public API available for tests and scripts that only pass an
    externally computed forcing. With zero forcing and ``use_phi_phase=False``
    it preserves the input state exactly, matching the old smoke-test contract.
    """
    psi = np.asarray(psi_J, dtype=np.complex128)
    phi = np.asarray(phi_V, dtype=float)
    Te = np.asarray(Te_K, dtype=float)
    forcing = np.asarray(forcing_J, dtype=np.complex128)

    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive.")
    if psi.shape != forcing.shape:
        raise ValueError("psi_J and forcing_J must have the same shape.")

    rho = np.maximum(material.rho_kwt(Te, np.abs(psi)), 1.0e-30)
    psi_new = psi + (float(dt_s) / material.tau0_GL_s) * forcing / rho

    max_abs_phase = 0.0
    if use_phi_phase:
        phase_step = (2.0 * E_CHARGE_C / HBAR_J_S) * phi * float(dt_s)
        max_abs_phase = float(np.max(np.abs(phase_step))) if phase_step.size else 0.0
        if max_abs_phase > max_phase_step_rad:
            return psi.copy(), False, max_abs_phase
        psi_new = np.exp(1j * phase_step) * psi_new

    if not (
        np.all(np.isfinite(np.real(psi_new)))
        and np.all(np.isfinite(np.imag(psi_new)))
    ):
        return psi.copy(), False, max_abs_phase

    return psi_new, True, max_abs_phase

