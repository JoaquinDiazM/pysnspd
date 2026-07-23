"""Runtime power and energy catalogues for the coupled thermal layer.

The phase-space catalogue stores the expensive superconducting kernels
``J_S(Omega; Te, |Delta|, q)`` and ``J_R(Omega; Te, |Delta|, q)``.  Those
kernels are still inconvenient inside a gTDGL time loop because the thermal
solver would have to integrate over the phonon-energy axis at every node and
at every time step.

This module performs the next PRE-run reduction: it contracts the phase-space
kernels with the Eliashberg spectrum and the Bose imbalance on a configurable
``Tph`` axis.  The resulting arrays are local lookup tables for
``P_ep^S(Te,Tph,|Delta|,q)`` and ``P_ep^R(Te,Tph,|Delta|,q)``.  It also stores
an electronic-energy table consistent with the Appendix-A energy functional,
so later OE6 code can advance the temperature through an energy variable rather
than redoing microscopic integrals.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from pysnspd.config import validate_config
from pysnspd.kinetic.eliashberg import EliashbergSpectrum, j_to_mev
from pysnspd.kinetic.phase_space import PhaseSpaceCatalog, fermi_positive_energy
from pysnspd.usadel.catalog import UsadelCatalog, J_to_meV
from pysnspd.usadel.parameters import HBAR_J_S, K_B_J_K

MEV_J = 1.602176634e-22
E_CHARGE_C = 1.602176634e-19
H_J_S = 2.0 * np.pi * HBAR_J_S


def electronic_density_of_states_from_sigma_D(
    sigma_n_S_m: float,
    D_m2_s: float,
) -> float:
    """Return the single-spin density of states from the Einstein relation."""
    if sigma_n_S_m <= 0.0:
        raise ValueError("sigma_n_S_m must be positive.")
    if D_m2_s <= 0.0:
        raise ValueError("D_m2_s must be positive.")
    return float(sigma_n_S_m / (2.0 * E_CHARGE_C**2 * D_m2_s))


def bose_positive_energy(omega_J: np.ndarray, T_K: float) -> np.ndarray:
    """Return the Bose-Einstein occupation for positive phonon energy."""
    if T_K <= 0.0:
        raise ValueError("Temperature must be positive.")

    omega = np.asarray(omega_J, dtype=float)
    out = np.zeros_like(omega, dtype=float)
    positive = omega > 0.0
    x = omega[positive] / (K_B_J_K * float(T_K))
    vals = np.zeros_like(x)
    small = x < 1.0e-6
    vals[small] = 1.0 / x[small] - 0.5 + x[small] / 12.0
    vals[~small] = 1.0 / np.expm1(np.clip(x[~small], 0.0, 700.0))
    out[positive] = vals
    return out


@dataclass(frozen=True)
class PowerTableCatalog:
    """Projected power and local-energy catalogue.

    Array layout
    ------------
    ``P_*_WTdq`` use shape ``(n_Te, n_Tph, n_delta, n_q)``.
    ``u_e_Tdq_J_m3`` and ``C_e_Tdq_J_m3_K`` use shape
    ``(n_Te, n_delta, n_q)``.
    """

    Te_values_K: np.ndarray
    Tph_values_K: np.ndarray
    delta_values_J: np.ndarray
    gamma_values_J: np.ndarray
    q_values_m_inv: np.ndarray
    omega_values_J: np.ndarray
    alpha2F: np.ndarray
    phdos_states_per_THz: np.ndarray
    P_S_W_m3: np.ndarray
    P_R_W_m3: np.ndarray
    P_total_W_m3: np.ndarray
    u_e_J_m3: np.ndarray
    C_e_J_m3_K: np.ndarray
    kappa_s_W_m_K: np.ndarray
    u_ph_J_m3: np.ndarray
    C_ph_J_m3_K: np.ndarray
    P_esc_W_m3: np.ndarray
    u_ph_weighted_J: np.ndarray
    C_ph_weighted_J_K: np.ndarray
    metadata: dict[str, Any]


def build_power_table_catalog(
    *,
    phase_space_catalog: PhaseSpaceCatalog,
    usadel_catalog: UsadelCatalog,
    spectrum: EliashbergSpectrum,
    config: Mapping[str, Any],
    n_Tph: int | None = None,
    Tph_min_K: float | None = None,
    Tph_max_K: float | None = None,
    omega_max_meV: float | None = None,
    workers: int | None = None,
    parallel_backend: str | None = None,
    progress: bool = True,
) -> PowerTableCatalog:
    """Build PRE-run projected powers and energy tables.

    The expensive energy integral over the quasiparticle axis has already been
    reduced to ``J_S`` and ``J_R`` by the phase-space catalogue.  This builder
    performs the remaining Omega contraction for all ``Tph`` values and stores
    the result as lookup tables for the future OE6 thermal time step.
    """
    cfg = validate_config(config, require_big_data_root_exists=False)
    n_workers, backend = _resolve_parallel(cfg, workers, parallel_backend)

    phase_cfg = cfg.get("catalogs", {}).get("phase_space", {})
    if n_Tph is None:
        n_Tph = int(phase_cfg.get("n_Tph", phase_space_catalog.Te_values_K.size))
    if Tph_min_K is None:
        Tph_min_K = float(phase_cfg.get("Tph_min_K", cfg["bias"]["T_bias_K"]))
    if Tph_max_K is None:
        Tph_max_K = float(phase_cfg.get("Tph_max_K", float(np.max(phase_space_catalog.Te_values_K))))
    if omega_max_meV is None:
        omega_max_meV = float(j_to_mev(np.max(phase_space_catalog.omega_values_J)))

    if int(n_Tph) <= 0:
        raise ValueError("n_Tph must be positive.")
    if float(Tph_min_K) <= 0.0 or float(Tph_max_K) < float(Tph_min_K):
        raise ValueError("Require 0 < Tph_min_K <= Tph_max_K.")
    if float(omega_max_meV) <= 0.0:
        raise ValueError("omega_max_meV must be positive.")

    Te_values = np.asarray(phase_space_catalog.Te_values_K, dtype=float)
    Tph_values = np.linspace(float(Tph_min_K), float(Tph_max_K), int(n_Tph))
    delta_values = np.asarray(phase_space_catalog.delta_values_J, dtype=float)
    gamma_values = np.asarray(phase_space_catalog.gamma_values_J, dtype=float)
    q_values = np.asarray(phase_space_catalog.q_values_m_inv, dtype=float)
    omega_all = np.asarray(phase_space_catalog.omega_values_J, dtype=float)

    omega_max_J = float(omega_max_meV) * MEV_J
    omega_mask = omega_all <= omega_max_J
    if np.count_nonzero(omega_mask) < 2:
        raise ValueError("Power table requires at least two Omega samples.")
    omega = omega_all[omega_mask]
    alpha = np.asarray(spectrum.alpha2F_on_omega_J(omega), dtype=float)
    phdos = np.asarray(spectrum.phdos_on_omega_J(omega), dtype=float)

    sigma_n = _metadata_float(usadel_catalog.metadata, "sigma_n_S_m")
    D_m2_s = _metadata_float(usadel_catalog.metadata, "D_m2_s")
    N0 = electronic_density_of_states_from_sigma_D(sigma_n, D_m2_s)
    delta0_J = _delta0_from_catalog(usadel_catalog)
    ion_density_m3, ion_density_source = _ion_density_from_config(cfg)
    tau_esc_s, tau_esc_source = _tau_esc_from_config(cfg)

    shape_power = (Te_values.size, Tph_values.size, delta_values.size, q_values.size)
    P_S = np.empty(shape_power, dtype=float)
    P_R = np.empty(shape_power, dtype=float)
    u_e = np.empty((Te_values.size, delta_values.size, q_values.size), dtype=float)

    delta_indices = np.asarray(phase_space_catalog.delta_indices, dtype=np.int64)
    q_indices = np.asarray(phase_space_catalog.q_indices, dtype=np.int64)
    energy = np.asarray(usadel_catalog.energy_values_J, dtype=float)

    tasks: list[tuple[Any, ...]] = []
    for iT, Te in enumerate(Te_values):
        for id_local, id_parent in enumerate(delta_indices):
            tasks.append(
                (
                    int(iT),
                    float(Te),
                    int(id_local),
                    float(delta_values[id_local]),
                    omega,
                    alpha,
                    Tph_values,
                    np.asarray(phase_space_catalog.J_S_TdqO_J[iT, id_local, :, :][:, omega_mask], dtype=float),
                    np.asarray(phase_space_catalog.J_R_TdqO_J[iT, id_local, :, :][:, omega_mask], dtype=float),
                    energy,
                    np.asarray(usadel_catalog.rho_delta_gamma_E[id_parent, q_indices, :], dtype=float),
                    float(N0),
                    float(delta0_J),
                )
            )

    progress_bar = _PowerTableProgress(
        total_chunks=len(tasks),
        states_per_chunk=Tph_values.size * q_values.size,
        workers=n_workers,
        backend=backend,
        enabled=bool(progress),
    )
    progress_bar.begin()

    if n_workers <= 1 or backend == "serial":
        for task in tasks:
            iT, id_local, PS_block, PR_block, ue_block = _power_T_delta_task(task)
            P_S[iT, :, id_local, :] = PS_block
            P_R[iT, :, id_local, :] = PR_block
            u_e[iT, id_local, :] = ue_block
            progress_bar.update()
    else:
        executor_cls = ThreadPoolExecutor if backend == "thread" else ProcessPoolExecutor
        with executor_cls(max_workers=n_workers) as executor:
            futures = [executor.submit(_power_T_delta_task, task) for task in tasks]
            for future in as_completed(futures):
                iT, id_local, PS_block, PR_block, ue_block = future.result()
                P_S[iT, :, id_local, :] = PS_block
                P_R[iT, :, id_local, :] = PR_block
                u_e[iT, id_local, :] = ue_block
                progress_bar.update()

    progress_bar.done()

    P_total = P_S + P_R
    C_e = _temperature_gradient(u_e, Te_values, axis=0)
    kappa_s = _bardeen_kappa_s_table(
        Te_values_K=Te_values,
        delta_values_J=delta_values,
        D_m2_s=float(D_m2_s),
        N0_J_m3=float(N0),
    )
    u_ph = _phonon_energy_density_simon(
        omega,
        phdos,
        Tph_values,
        ion_density_m3=float(ion_density_m3),
    )
    C_ph = _temperature_gradient(u_ph, Tph_values, axis=0)
    T_bath_K = float(cfg.get("bias", {}).get("T_bias_K", Tph_values[0]))
    u_ph_bath = float(
        _phonon_energy_density_simon(
            omega,
            phdos,
            np.asarray([T_bath_K], dtype=float),
            ion_density_m3=float(ion_density_m3),
        )[0]
    )
    P_esc = (u_ph - u_ph_bath) / float(tau_esc_s)

    metadata = {
        "backend": "projected_power_table_from_phase_space_oe6_pre_v2",
        "description": (
            "Runtime-oriented PRE catalogue. It contracts J_S/J_R with alpha2F, "
            "Omega and the Bose imbalance on a Tph axis, avoiding Omega integrals "
            "inside the coupled gTDGL thermal loop."
        ),
        "sign_convention": "Positive P_S/P_R means energy leaves electrons and enters phonons.",
        "units": {
            "Te_values_K": "K",
            "Tph_values_K": "K",
            "delta_values_J": "J",
            "gamma_values_J": "J",
            "q_values_m_inv": "m^-1",
            "omega_values_J": "J",
            "P_S_W_m3": "W m^-3",
            "P_R_W_m3": "W m^-3",
            "P_total_W_m3": "W m^-3",
            "u_e_J_m3": "J m^-3",
            "C_e_J_m3_K": "J m^-3 K^-1",
            "kappa_s_W_m_K": "W m^-1 K^-1",
            "u_ph_J_m3": "J m^-3",
            "C_ph_J_m3_K": "J m^-3 K^-1",
            "P_esc_W_m3": "W m^-3, positive means phonons lose energy to bath",
            "u_ph_weighted_J": "legacy diagnostic weighted phonon integral",
            "C_ph_weighted_J_K": "legacy diagnostic weighted phonon integral derivative",
        },
        "N0_J_m3": float(N0),
        "D_m2_s": float(D_m2_s),
        "sigma_n_S_m": float(sigma_n),
        "delta0_J": float(delta0_J),
        "delta0_meV": float(J_to_meV(delta0_J)),
        "ion_density_m3": float(ion_density_m3),
        "ion_density_nm3": float(ion_density_m3) * 1.0e-27,
        "ion_density_source": str(ion_density_source),
        "tau_esc_s": float(tau_esc_s),
        "tau_esc_ps": float(tau_esc_s) / 1.0e-12,
        "tau_esc_source": str(tau_esc_source),
        "T_bath_K": float(T_bath_K),
        "u_ph_bath_J_m3": float(u_ph_bath),
        "kappa_s_formula": "Bardeen/Allmaras Eq. 3.9: kappa_s=(2*pi^2*D*k_B^2*N0/3)*Te*[1-(6/pi^2)*int_0^(Delta/kBT) x^2 exp(x)/(exp(x)+1)^2 dx]",
        "phonon_energy_formula": "u_ph=N*int dnu_THz F(nu) Omega(nu) n_B(Omega,Tph), with N from material.ion_density_* and PhDOS in states/THz",
        "phonon_escape_formula": "P_esc=(u_ph(Tph)-u_ph(T_bath))/tau_esc",
        "omega_max_meV_requested": float(omega_max_meV),
        "omega_max_meV_used": float(j_to_mev(float(omega[-1]))),
        "alpha2F_source": spectrum.metadata.get("source", ""),
        "alpha2F_path": spectrum.metadata.get("path", ""),
        "alpha2F_policy": spectrum.metadata.get("alpha2F_policy", ""),
        "n_Te": int(Te_values.size),
        "n_Tph": int(Tph_values.size),
        "n_delta": int(delta_values.size),
        "n_q": int(q_values.size),
        "n_omega_used": int(omega.size),
        "parallel_workers": int(n_workers),
        "parallel_backend": str(backend),
        "parallel_tasks": int(len(tasks)),
        "parallel_task_layout": "one task per (Te, Delta), each task computes all Tph and q states",
        "parallel_states": int(Te_values.size * Tph_values.size * delta_values.size * q_values.size),
        "phase_space_backend": phase_space_catalog.metadata.get("backend", ""),
        "phase_space_shape": list(phase_space_catalog.shape),
        "energy_functional_policy": (
            "u_e follows Appendix A / Simon A2 under the thermal closure, with |Delta| "
            "and q treated as local state variables."
        ),
        "phonon_energy_policy": (
            "u_ph and C_ph are volume-normalized Simon phonon energy tables using "
            "u_ph=N*int F(Omega) Omega n_B dOmega, with the tabulated PhDOS interpreted "
            "as states per THz and N supplied by the material configuration."
        ),
        "thermal_conductivity_policy": (
            "kappa_s_W_m_K is tabulated on (Te, Delta) using the Bardeen-type "
            "superconducting electronic thermal conductivity quoted by Allmaras Eq. 3.9."
        ),
    }

    return PowerTableCatalog(
        Te_values_K=Te_values,
        Tph_values_K=Tph_values,
        delta_values_J=delta_values,
        gamma_values_J=gamma_values,
        q_values_m_inv=q_values,
        omega_values_J=omega,
        alpha2F=alpha,
        phdos_states_per_THz=phdos,
        P_S_W_m3=P_S,
        P_R_W_m3=P_R,
        P_total_W_m3=P_total,
        u_e_J_m3=u_e,
        C_e_J_m3_K=C_e,
        kappa_s_W_m_K=kappa_s,
        u_ph_J_m3=u_ph,
        C_ph_J_m3_K=C_ph,
        P_esc_W_m3=P_esc,
        u_ph_weighted_J=u_ph,
        C_ph_weighted_J_K=C_ph,
        metadata=metadata,
    )


def _power_T_delta_task(task: tuple[Any, ...]) -> tuple[int, int, np.ndarray, np.ndarray, np.ndarray]:
    (
        iT,
        Te_K,
        id_local,
        delta_J,
        omega,
        alpha,
        Tph_values,
        JS_block,
        JR_block,
        energy,
        rho_block,
        N0,
        delta0_J,
    ) = task

    omega = np.asarray(omega, dtype=float)
    alpha = np.asarray(alpha, dtype=float)
    Tph_values = np.asarray(Tph_values, dtype=float)
    JS_block = _as_q_omega_block(JS_block, omega_size=omega.size, name="JS_block")
    JR_block = _as_q_omega_block(JR_block, omega_size=omega.size, name="JR_block")

    n_e = bose_positive_energy(omega, float(Te_K))
    n_ph = np.vstack([bose_positive_energy(omega, float(Tph)) for Tph in Tph_values])
    bose_diff = n_e[None, :] - n_ph
    weight = alpha[None, :] * omega[None, :] * bose_diff

    integrand_S = weight[:, None, :] * JS_block[None, :, :]
    integrand_R = weight[:, None, :] * JR_block[None, :, :]
    P_S = (8.0 * np.pi * float(N0) / HBAR_J_S) * np.trapezoid(integrand_S, omega, axis=2)
    P_R = (4.0 * np.pi * float(N0) / HBAR_J_S) * np.trapezoid(integrand_R, omega, axis=2)

    u_e = _electron_energy_density_block(
        energy,
        rho_block,
        Te_K=float(Te_K),
        delta_J=float(delta_J),
        delta0_J=float(delta0_J),
        N0_J_m3=float(N0),
    )
    return int(iT), int(id_local), np.asarray(P_S, dtype=float), np.asarray(P_R, dtype=float), u_e


def _as_q_omega_block(block: np.ndarray, *, omega_size: int, name: str) -> np.ndarray:
    """Return a phase-space block with shape ``(n_q, n_omega)``.

    Numpy advanced indexing can accidentally transpose ``(:, omega_mask)`` into
    ``(n_omega, n_q)`` when the boolean mask is applied directly.  The PRE-run
    builder now slices in two steps, but this runtime guard keeps the table
    contraction robust for smoke objects and older catalogues.
    """
    arr = np.asarray(block, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"{name} must be two-dimensional, got shape {arr.shape}.")
    if arr.shape[1] == int(omega_size):
        return arr
    if arr.shape[0] == int(omega_size):
        return arr.T
    raise ValueError(
        f"{name} must have one axis with n_omega={omega_size}; got shape {arr.shape}."
    )


def _electron_energy_density_block(
    energy_J: np.ndarray,
    rho_block: np.ndarray,
    *,
    Te_K: float,
    delta_J: float,
    delta0_J: float,
    N0_J_m3: float,
) -> np.ndarray:
    E = np.asarray(energy_J, dtype=float)
    rho = np.asarray(rho_block, dtype=float)
    if rho.ndim != 2 or rho.shape[1] != E.size:
        raise ValueError("rho_block must have shape (n_q, n_energy).")
    f = fermi_positive_energy(E, Te_K)
    qp = np.trapezoid(rho * (E[None, :] * f[None, :]), E, axis=1)
    cond = 0.0
    if delta_J > 0.0 and delta0_J > 0.0:
        cond = 0.25 * delta_J * delta_J * (0.5 + np.log(delta0_J / max(delta_J, 1.0e-300)))
    out = 4.0 * N0_J_m3 * (qp - cond)
    return np.asarray(out, dtype=float)


def _bardeen_kappa_s_table(
    *,
    Te_values_K: np.ndarray,
    delta_values_J: np.ndarray,
    D_m2_s: float,
    N0_J_m3: float,
) -> np.ndarray:
    """Bardeen/Allmaras superconducting electronic thermal conductivity.

    Returns a table with shape ``(n_Te, n_delta)`` in W m^-1 K^-1.
    The normal-state prefactor is reduced by the standard Bardeen integral,
    so ``Delta=0`` recovers the Wiedemann-like dirty-limit normal value.
    """
    Te_values = np.asarray(Te_values_K, dtype=float)
    delta_values = np.asarray(delta_values_J, dtype=float)
    out = np.empty((Te_values.size, delta_values.size), dtype=float)
    prefactor = 2.0 * np.pi * np.pi * float(D_m2_s) * K_B_J_K * K_B_J_K * float(N0_J_m3) / 3.0
    for iT, Te in enumerate(Te_values):
        if Te <= 0.0 or not np.isfinite(Te):
            out[iT, :] = np.nan
            continue
        kappa_n = prefactor * float(Te)
        for idelta, delta in enumerate(delta_values):
            upper = float(delta) / (K_B_J_K * float(Te)) if delta > 0.0 else 0.0
            reduction = _bardeen_reduction_integral(upper)
            factor = 1.0 - (6.0 / (np.pi * np.pi)) * reduction
            out[iT, idelta] = kappa_n * float(np.clip(factor, 0.0, 1.0))
    return np.asarray(out, dtype=float)


def _bardeen_reduction_integral(upper: float) -> float:
    if not np.isfinite(upper) or upper <= 0.0:
        return 0.0
    # The integral is essentially saturated by x~50.  The sech form is stable
    # and equal to exp(x)/(exp(x)+1)^2.
    xmax = min(float(upper), 80.0)
    n = max(256, min(4096, int(64 * xmax) + 1))
    x = np.linspace(0.0, xmax, n, dtype=float)
    integrand = 0.25 * x * x / np.cosh(0.5 * np.minimum(x, 80.0)) ** 2
    return float(np.trapezoid(integrand, x))


def _phonon_energy_density_simon(
    omega_J: np.ndarray,
    phdos_states_per_THz: np.ndarray,
    Tph_values_K: np.ndarray,
    *,
    ion_density_m3: float,
) -> np.ndarray:
    """Volume-normalized Simon phonon energy density.

    The Simon/MIT PhDOS column is tabulated per THz, while the PRE catalogue
    uses a phonon-energy axis in joules.  We therefore integrate over the
    corresponding frequency in THz and multiply by the ion density ``N``.
    """
    omega = np.asarray(omega_J, dtype=float)
    density = np.asarray(phdos_states_per_THz, dtype=float)
    freq_THz = omega / (H_J_S * 1.0e12)
    out = np.zeros_like(np.asarray(Tph_values_K, dtype=float), dtype=float)
    for i, Tph in enumerate(Tph_values_K):
        n = bose_positive_energy(omega, float(Tph))
        integrand = omega * density * n
        integrand = np.nan_to_num(integrand, nan=0.0, posinf=0.0, neginf=0.0)
        out[i] = float(ion_density_m3) * float(np.trapezoid(integrand, freq_THz))
    return out


def _temperature_gradient(values: np.ndarray, axis_values_K: np.ndarray, *, axis: int) -> np.ndarray:
    vals = np.asarray(values, dtype=float)
    x = np.asarray(axis_values_K, dtype=float)
    if x.size < 2:
        return np.zeros_like(vals)
    return np.asarray(np.gradient(vals, x, axis=axis, edge_order=1), dtype=float)


def save_power_table_catalog_npz(catalog: PowerTableCatalog, path: str | Path) -> Path:
    """Save a power-table catalogue to compressed ``.npz``."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        Te_values_K=catalog.Te_values_K,
        Tph_values_K=catalog.Tph_values_K,
        delta_values_J=catalog.delta_values_J,
        gamma_values_J=catalog.gamma_values_J,
        q_values_m_inv=catalog.q_values_m_inv,
        omega_values_J=catalog.omega_values_J,
        alpha2F=catalog.alpha2F,
        phdos_states_per_THz=catalog.phdos_states_per_THz,
        P_S_W_m3=catalog.P_S_W_m3,
        P_R_W_m3=catalog.P_R_W_m3,
        P_total_W_m3=catalog.P_total_W_m3,
        u_e_J_m3=catalog.u_e_J_m3,
        C_e_J_m3_K=catalog.C_e_J_m3_K,
        kappa_s_W_m_K=catalog.kappa_s_W_m_K,
        u_ph_J_m3=catalog.u_ph_J_m3,
        C_ph_J_m3_K=catalog.C_ph_J_m3_K,
        P_esc_W_m3=catalog.P_esc_W_m3,
        u_ph_weighted_J=catalog.u_ph_weighted_J,
        C_ph_weighted_J_K=catalog.C_ph_weighted_J_K,
        metadata=np.array(catalog.metadata, dtype=object),
    )
    return output


def power_table_summary(catalog: PowerTableCatalog) -> dict[str, Any]:
    """Return a compact diagnostic summary for manifests and terminal output."""
    P_S = np.asarray(catalog.P_S_W_m3, dtype=float)
    P_R = np.asarray(catalog.P_R_W_m3, dtype=float)
    P_total = np.asarray(catalog.P_total_W_m3, dtype=float)
    ue = np.asarray(catalog.u_e_J_m3, dtype=float)
    Ce = np.asarray(catalog.C_e_J_m3_K, dtype=float)
    kappa = np.asarray(catalog.kappa_s_W_m_K, dtype=float)
    uph = np.asarray(catalog.u_ph_J_m3, dtype=float)
    Cph = np.asarray(catalog.C_ph_J_m3_K, dtype=float)
    Pesc = np.asarray(catalog.P_esc_W_m3, dtype=float)

    equal_power = []
    for iT, Te in enumerate(catalog.Te_values_K):
        iTph = int(np.argmin(np.abs(catalog.Tph_values_K - Te)))
        if abs(float(catalog.Tph_values_K[iTph]) - float(Te)) <= 1.0e-10:
            equal_power.append(np.asarray(P_total[iT, iTph, :, :], dtype=float))
    if equal_power:
        max_abs_equal = float(np.max(np.abs(np.stack(equal_power, axis=0))))
    else:
        max_abs_equal = float("nan")

    return {
        "backend": str(catalog.metadata.get("backend", "")),
        "n_Te": int(catalog.Te_values_K.size),
        "n_Tph": int(catalog.Tph_values_K.size),
        "n_delta": int(catalog.delta_values_J.size),
        "n_q": int(catalog.q_values_m_inv.size),
        "n_omega_used": int(catalog.omega_values_J.size),
        "Te_min_K": float(np.min(catalog.Te_values_K)),
        "Te_max_K": float(np.max(catalog.Te_values_K)),
        "Tph_min_K": float(np.min(catalog.Tph_values_K)),
        "Tph_max_K": float(np.max(catalog.Tph_values_K)),
        "delta_min_meV": float(J_to_meV(np.min(catalog.delta_values_J))),
        "delta_max_meV": float(J_to_meV(np.max(catalog.delta_values_J))),
        "q_min_m_inv": float(np.min(catalog.q_values_m_inv)),
        "q_max_m_inv": float(np.max(catalog.q_values_m_inv)),
        "omega_max_meV_used": float(j_to_mev(np.max(catalog.omega_values_J))),
        "P_S_min_W_m3": float(np.min(P_S)),
        "P_S_max_W_m3": float(np.max(P_S)),
        "P_R_min_W_m3": float(np.min(P_R)),
        "P_R_max_W_m3": float(np.max(P_R)),
        "P_total_min_W_m3": float(np.min(P_total)),
        "P_total_max_W_m3": float(np.max(P_total)),
        "P_total_is_finite": bool(np.all(np.isfinite(P_total))),
        "max_abs_equal_T_power_W_m3": max_abs_equal,
        "u_e_min_J_m3": float(np.min(ue)),
        "u_e_max_J_m3": float(np.max(ue)),
        "C_e_min_J_m3_K": float(np.min(Ce)),
        "C_e_max_J_m3_K": float(np.max(Ce)),
        "u_e_is_finite": bool(np.all(np.isfinite(ue))),
        "C_e_is_finite": bool(np.all(np.isfinite(Ce))),
        "kappa_s_min_W_m_K": float(np.nanmin(kappa)),
        "kappa_s_max_W_m_K": float(np.nanmax(kappa)),
        "kappa_s_is_finite": bool(np.all(np.isfinite(kappa))),
        "u_ph_min_J_m3": float(np.nanmin(uph)),
        "u_ph_max_J_m3": float(np.nanmax(uph)),
        "C_ph_min_J_m3_K": float(np.nanmin(Cph)),
        "C_ph_max_J_m3_K": float(np.nanmax(Cph)),
        "P_esc_min_W_m3": float(np.nanmin(Pesc)),
        "P_esc_max_W_m3": float(np.nanmax(Pesc)),
        "P_esc_is_finite": bool(np.all(np.isfinite(Pesc))),
        "ion_density_nm3": float(catalog.metadata.get("ion_density_nm3", np.nan)),
        "tau_esc_ps": float(catalog.metadata.get("tau_esc_ps", np.nan)),
        "parallel_workers": int(catalog.metadata.get("parallel_workers", 1)),
        "parallel_backend": str(catalog.metadata.get("parallel_backend", "serial")),
        "parallel_tasks": int(catalog.metadata.get("parallel_tasks", 1)),
        "parallel_states": int(catalog.metadata.get("parallel_states", 1)),
        "sign_convention": str(catalog.metadata.get("sign_convention", "")),
        "energy_functional_policy": str(catalog.metadata.get("energy_functional_policy", "")),
        "phonon_energy_policy": str(catalog.metadata.get("phonon_energy_policy", "")),
    }


def _ion_density_from_config(cfg: Mapping[str, Any]) -> tuple[float, str]:
    material = cfg.get("material", {}) if isinstance(cfg, Mapping) else {}
    for key in ("ion_density_m3", "ion_density_per_m3", "N_m3", "N_per_m3"):
        if key in material:
            value = float(material[key])
            if np.isfinite(value) and value > 0.0:
                return value, f"material.{key}"
    for key in ("ion_density_nm3", "ion_density_per_nm3", "N_nm3", "N_per_nm3"):
        if key in material:
            value = float(material[key]) * 1.0e27
            if np.isfinite(value) and value > 0.0:
                return value, f"material.{key}"
    # Simon NbN table value; kept as fallback so older smoke configs still run.
    return 48.0e27, "fallback: Simon NbN N=48 nm^-3"


def _tau_esc_from_config(cfg: Mapping[str, Any]) -> tuple[float, str]:
    material = cfg.get("material", {}) if isinstance(cfg, Mapping) else {}
    for key, scale in (("tau_esc_s", 1.0), ("tau_escape_s", 1.0), ("tau_esc_ps", 1.0e-12), ("tau_escape_ps", 1.0e-12)):
        if key in material:
            value = float(material[key]) * scale
            if np.isfinite(value) and value > 0.0:
                return value, f"material.{key}"
    return 20.0e-12, "fallback: 20 ps"


def _resolve_parallel(
    cfg: Mapping[str, Any],
    workers: int | None,
    backend: str | None,
) -> tuple[int, str]:
    parallel = cfg.get("parallel", {}) if isinstance(cfg, Mapping) else {}
    resolved_backend = str(backend or parallel.get("backend", "process")).lower()
    if resolved_backend not in {"process", "thread", "serial"}:
        resolved_backend = "process"
    if workers is None:
        if bool(parallel.get("enabled", False)):
            resolved_workers = int(parallel.get("workers", 1))
        else:
            resolved_workers = 1
    else:
        resolved_workers = int(workers)
    if resolved_backend == "serial":
        resolved_workers = 1
    return max(1, resolved_workers), resolved_backend


def _metadata_float(metadata: Mapping[str, Any], key: str) -> float:
    try:
        return float(metadata[key])
    except Exception as exc:
        raise KeyError(f"Required metadata key is missing or non-numeric: {key}") from exc


def _delta0_from_catalog(usadel_catalog: UsadelCatalog) -> float:
    md = getattr(usadel_catalog, "metadata", {})
    for key in ("delta0_J", "Delta0_J", "gap0_J"):
        if key in md:
            val = float(md[key])
            if np.isfinite(val) and val > 0.0:
                return val
    for key in ("delta0_meV", "Delta0_meV", "gap0_meV"):
        if key in md:
            val = float(md[key]) * MEV_J
            if np.isfinite(val) and val > 0.0:
                return val
    if hasattr(usadel_catalog, "calibration_delta_eq_values_J"):
        vals = np.asarray(usadel_catalog.calibration_delta_eq_values_J, dtype=float)
        vals = vals[np.isfinite(vals) & (vals > 0.0)]
        if vals.size:
            return float(np.max(vals))
    vals = np.asarray(usadel_catalog.delta_values_J, dtype=float)
    vals = vals[np.isfinite(vals) & (vals > 0.0)]
    if vals.size:
        return float(np.max(vals))
    raise ValueError("Could not infer a positive Delta0 from the Usadel catalogue.")


class _PowerTableProgress:
    """Single-line dependency-free progress bar for projected-power chunks."""

    def __init__(
        self,
        *,
        total_chunks: int,
        states_per_chunk: int,
        workers: int,
        backend: str,
        enabled: bool,
        width: int = 34,
    ) -> None:
        self.total_chunks = max(1, int(total_chunks))
        self.states_per_chunk = max(1, int(states_per_chunk))
        self.total_states = self.total_chunks * self.states_per_chunk
        self.workers = int(workers)
        self.backend = str(backend)
        self.enabled = bool(enabled)
        self.width = int(width)
        self.done_chunks = 0
        self._last_percent = -1
        self._prefix = (
            "Power-table "
            f"({self.total_chunks} chunks, {self.total_states} states, "
            f"workers={self.workers}, backend={self.backend})"
        )

    def begin(self) -> None:
        if not self.enabled:
            return
        self._print(force=True, final=False)

    def update(self) -> None:
        if not self.enabled:
            return
        self.done_chunks = min(self.done_chunks + 1, self.total_chunks)
        self._print(force=False, final=False)

    def done(self) -> None:
        if not self.enabled:
            return
        self.done_chunks = self.total_chunks
        self._print(force=True, final=True)

    def _print(self, *, force: bool, final: bool) -> None:
        percent = int(round(100.0 * self.done_chunks / self.total_chunks))
        if not force and percent == self._last_percent and self.done_chunks != self.total_chunks:
            return
        self._last_percent = percent
        frac = self.done_chunks / self.total_chunks
        filled = int(round(self.width * frac))
        bar = "=" * filled + "-" * (self.width - filled)
        states_done = self.done_chunks * self.states_per_chunk
        line = (
            f"\r{self._prefix}: [{bar}] {percent:3d}% "
            f"chunks={self.done_chunks}/{self.total_chunks} "
            f"states={states_done}/{self.total_states}"
        )
        print(line, end="\n" if final else "", flush=True)


__all__ = [
    "PowerTableCatalog",
    "bose_positive_energy",
    "build_power_table_catalog",
    "electronic_density_of_states_from_sigma_D",
    "power_table_summary",
    "save_power_table_catalog_npz",
]
