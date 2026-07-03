"""
Phase-space catalogues for the pySNSPD kinetic block.

This module builds tabulated phase-space integrals J_S and J_R from the Usadel
DOS catalogue. The expensive catalogue cells are independent in (Te, |Delta|,
q), so the builder supports the same worker/backend policy used by the PRE-run
for the Usadel/DOS and strict Matsubara supercurrent tables.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping

import numpy as np

from pysnspd.config import validate_config
from pysnspd.usadel.catalog import UsadelCatalog, J_to_meV
from pysnspd.usadel.parameters import K_B_J_K


@dataclass(frozen=True)
class PhaseSpaceCatalog:
    """Container for phase-space integrals.

    Arrays use the shape ``(n_Te, n_delta, n_q, n_omega)``.
    """

    Te_values_K: np.ndarray
    omega_values_J: np.ndarray
    delta_values_J: np.ndarray
    gamma_values_J: np.ndarray
    q_values_m_inv: np.ndarray
    J_S_TdqO_J: np.ndarray
    J_R_TdqO_J: np.ndarray
    delta_indices: np.ndarray
    q_indices: np.ndarray
    metadata: dict[str, Any]

    @property
    def shape(self) -> tuple[int, int, int, int]:
        """Return catalogue shape ``(n_Te, n_delta, n_q, n_omega)``."""
        return tuple(int(v) for v in self.J_S_TdqO_J.shape)


def build_phase_space_catalog_from_usadel_catalog(
    usadel_catalog: UsadelCatalog,
    config: Mapping[str, Any],
    *,
    n_Te: int | None = None,
    n_delta: int | None = None,
    n_q: int | None = None,
    n_omega: int | None = None,
    Te_min_K: float | None = None,
    Te_max_K: float | None = None,
    omega_max_meV: float | None = None,
    workers: int | None = None,
    parallel_backend: str | None = None,
    progress: bool = True,
) -> PhaseSpaceCatalog:
    """Build a phase-space catalogue from a Usadel DOS catalogue.

    OE4-v1.4 keeps the finite-energy-window diagnostics from OE4-v1.3 and adds
    parallel construction over independent ``(Te, Delta)`` blocks. Each block
    computes all selected q slices for one temperature and one gap value.

    Parameters
    ----------
    workers, parallel_backend:
        Worker count and backend policy. ``parallel_backend`` accepts
        ``"process"``, ``"thread"`` or ``"serial"``. If omitted, values are
        resolved from the project config.
    progress:
        Print a dedicated phase-space progress bar, separate from the coarse
        PRE-run stage progress bar.
    """
    cfg = validate_config(config, require_big_data_root_exists=False)
    phase_cfg = cfg["catalogs"]["phase_space"]
    n_workers, backend = _resolve_parallel(cfg, workers, parallel_backend)

    if n_Te is None:
        n_Te = int(phase_cfg["n_Te"])
    if n_delta is None:
        n_delta = int(phase_cfg["n_delta"])
    if n_q is None:
        n_q = int(phase_cfg["n_q"])
    if n_omega is None:
        n_omega = int(phase_cfg["n_omega"])
    if Te_min_K is None:
        Te_min_K = float(phase_cfg.get("Te_min_K", cfg["bias"]["T_bias_K"]))
    if Te_max_K is None:
        Te_max_K = float(phase_cfg.get("Te_max_K", max(4.0 * cfg["material"]["Tc_K"], 30.0)))
    if omega_max_meV is None and "omega_max_meV" in phase_cfg:
        omega_max_meV = float(phase_cfg["omega_max_meV"])

    if n_Te <= 0 or n_delta <= 0 or n_q <= 0 or n_omega <= 1:
        raise ValueError("Phase-space grid sizes must be positive and n_omega > 1.")
    if Te_min_K <= 0.0 or Te_max_K <= Te_min_K:
        raise ValueError("Require 0 < Te_min_K < Te_max_K.")

    meV_J = 1.602176634e-22
    energy = np.asarray(usadel_catalog.energy_values_J, dtype=float)
    energy_max_J = float(np.max(energy))
    energy_max_meV = J_to_meV(energy_max_J)

    if omega_max_meV is None:
        omega_max_J = energy_max_J
        omega_axis_source = "parent_dos_energy_max"
    else:
        if omega_max_meV <= 0.0:
            raise ValueError("omega_max_meV must be positive when provided.")
        omega_max_J = float(omega_max_meV) * meV_J
        omega_axis_source = "user_requested_meV"
        if omega_max_J > energy_max_J:
            raise ValueError(
                "The requested phase-space Omega range exceeds the parent DOS range: "
                f"omega_max_meV={J_to_meV(omega_max_J):.6g}, "
                f"energy_max_meV={energy_max_meV:.6g}. "
                "Increase --energy-max-factor or reduce --phase-omega-max-meV."
            )

    delta_indices = _select_axis_indices(usadel_catalog.delta_values_J.size, int(n_delta))
    q_indices = _select_axis_indices(usadel_catalog.q_values_m_inv.size, int(n_q))
    Te_values_K = np.linspace(float(Te_min_K), float(Te_max_K), int(n_Te))
    omega_values_J = np.linspace(0.0, omega_max_J, int(n_omega))
    delta_values_J = usadel_catalog.delta_values_J[delta_indices]
    gamma_values_J = usadel_catalog.gamma_values_J[q_indices]
    q_values_m_inv = usadel_catalog.q_values_m_inv[q_indices]

    shape = (
        Te_values_K.size,
        delta_values_J.size,
        q_values_m_inv.size,
        omega_values_J.size,
    )
    J_S = np.empty(shape, dtype=float)
    J_R = np.empty(shape, dtype=float)

    tasks = [
        (
            int(iT),
            float(Te_K),
            int(id_local),
            int(id_parent),
            float(usadel_catalog.delta_values_J[id_parent]),
            np.asarray(usadel_catalog.rho_delta_gamma_E[id_parent, q_indices, :], dtype=float),
            energy,
            omega_values_J,
        )
        for iT, Te_K in enumerate(Te_values_K)
        for id_local, id_parent in enumerate(delta_indices)
    ]

    progress_bar = _PhaseSpaceProgress(
        total_chunks=len(tasks),
        cells_per_chunk=q_indices.size,
        workers=n_workers,
        backend=backend,
        enabled=bool(progress),
    )
    progress_bar.begin()

    if n_workers <= 1 or backend == "serial":
        for task in tasks:
            iT, id_local, JS_block, JR_block = _phase_space_T_delta_task(task)
            J_S[iT, id_local, :, :] = JS_block
            J_R[iT, id_local, :, :] = JR_block
            progress_bar.update()
    else:
        executor_cls = ThreadPoolExecutor if backend == "thread" else ProcessPoolExecutor
        with executor_cls(max_workers=n_workers) as executor:
            futures = [executor.submit(_phase_space_T_delta_task, task) for task in tasks]
            for future in as_completed(futures):
                iT, id_local, JS_block, JR_block = future.result()
                J_S[iT, id_local, :, :] = JS_block
                J_R[iT, id_local, :, :] = JR_block
                progress_bar.update()

    progress_bar.done()

    energy_floor_J = max(float(energy[1]) * 0.5, 1.0e-300)
    js_lower_by_delta_J = np.maximum(delta_values_J, energy_floor_J)
    js_hard_cutoff_by_delta_J = np.maximum(0.0, energy_max_J - js_lower_by_delta_J)
    jr_threshold_by_delta_J = np.maximum(0.0, 2.0 * delta_values_J)
    scattering_window_margin_J = float(np.min(js_hard_cutoff_by_delta_J)) - omega_max_J

    metadata = {
        "backend": "phase_space_from_usadel_dos_oe4_v1_4_parallel",
        "description": (
            "OE4 phase-space catalogue. It tabulates J_S and J_R from the "
            "Usadel DOS catalogue. It does not yet include alpha^2F(Omega), "
            "phonon DOS, T_ph, escape, or final power integrals."
        ),
        "units": {
            "Te_values_K": "K",
            "omega_values_J": "J",
            "delta_values_J": "J",
            "gamma_values_J": "J",
            "q_values_m_inv": "m^-1",
            "J_S_TdqO_J": "J",
            "J_R_TdqO_J": "J",
        },
        "parent_usadel_backend": str(usadel_catalog.metadata.get("backend", "unknown")),
        "parent_usadel_shape": list(usadel_catalog.shape),
        "delta_indices": delta_indices.tolist(),
        "q_indices": q_indices.tolist(),
        "grid_is_downsampled": bool(
            delta_indices.size < usadel_catalog.delta_values_J.size
            or q_indices.size < usadel_catalog.q_values_m_inv.size
        ),
        "Te_min_K": float(Te_values_K[0]),
        "Te_max_K": float(Te_values_K[-1]),
        "omega_axis_source": omega_axis_source,
        "omega_max_J": float(omega_values_J[-1]),
        "omega_max_meV": J_to_meV(float(omega_values_J[-1])),
        "energy_max_J": energy_max_J,
        "energy_max_meV": J_to_meV(energy_max_J),
        "energy_floor_J": energy_floor_J,
        "phase_space_parallel_workers": int(n_workers),
        "phase_space_parallel_backend": str(backend),
        "phase_space_parallel_tasks": int(len(tasks)),
        "phase_space_parallel_task_layout": "one task per (Te, Delta), each task computes all selected q slices",
        "phase_space_parallel_cells": int(Te_values_K.size * delta_values_J.size * q_values_m_inv.size),
        "js_hard_cutoff_by_delta_J": js_hard_cutoff_by_delta_J.tolist(),
        "js_hard_cutoff_by_delta_meV": [J_to_meV(float(v)) for v in js_hard_cutoff_by_delta_J],
        "jr_threshold_by_delta_J": jr_threshold_by_delta_J.tolist(),
        "jr_threshold_by_delta_meV": [J_to_meV(float(v)) for v in jr_threshold_by_delta_J],
        "scattering_window_margin_J": scattering_window_margin_J,
        "scattering_window_margin_meV": J_to_meV(scattering_window_margin_J),
        "scattering_window_is_truncated": bool(scattering_window_margin_J < 0.0),
        "energy_window_policy": (
            "PDF Appendix A writes the scattering energy moment with an upper "
            "limit extending to infinity. OE4 uses a finite Usadel DOS grid "
            "and a separately configurable Omega axis. J_S is reliable only "
            "where both E and E+Omega lie inside the parent DOS catalogue."
        ),
        "coherence_factor_policy": (
            "PDF Appendix A notes that the most general dirty-limit structure "
            "would use N1(E)N1(E') +/- R2(E)R2(E'). OE4 follows the Simon/BCS "
            "reduced coherence factors rho(E)rho(E')[1 +/- Delta^2/(EE')], "
            "with rho supplied by the Usadel DOS catalogue."
        ),
        "threshold_policy": (
            "J_R(Omega)=0 for Omega <= 2 Delta, following the Simon/BCS "
            "integration limits in Appendix A. In a strongly depaired Usadel "
            "spectrum, the physical spectral onset may be closer to 2E_g."
        ),
        "normal_limit_policy": (
            "J_R is set to zero for Delta <= 0 because it is treated as a "
            "distinct superconducting recombination/pair-breaking channel."
        ),
        "thermal_closure_policy": (
            "OE4 assumes f(E)->f_FD(E,Te). The full kinetic equations evolve "
            "nonthermal f(E,t) and n(Omega,t)."
        ),
    }

    return PhaseSpaceCatalog(
        Te_values_K=Te_values_K,
        omega_values_J=omega_values_J,
        delta_values_J=delta_values_J,
        gamma_values_J=gamma_values_J,
        q_values_m_inv=q_values_m_inv,
        J_S_TdqO_J=J_S,
        J_R_TdqO_J=J_R,
        delta_indices=delta_indices,
        q_indices=q_indices,
        metadata=metadata,
    )


def _phase_space_T_delta_task(task: tuple[Any, ...]) -> tuple[int, int, np.ndarray, np.ndarray]:
    """Compute all q slices for one (Te, Delta) phase-space task."""
    (
        iT,
        Te_K,
        id_local,
        _id_parent,
        delta_J,
        rho_q_E,
        energy,
        omega_values_J,
    ) = task
    rho_q_E = np.asarray(rho_q_E, dtype=float)
    n_q = int(rho_q_E.shape[0])
    n_omega = int(np.asarray(omega_values_J).size)
    JS_block = np.empty((n_q, n_omega), dtype=float)
    JR_block = np.empty((n_q, n_omega), dtype=float)
    for iq in range(n_q):
        rho_E = rho_q_E[iq, :]
        JS_block[iq, :] = scattering_phase_space_spectrum(
            energy,
            rho_E,
            omega_values_J,
            Te_K=float(Te_K),
            delta_J=float(delta_J),
        )
        JR_block[iq, :] = recombination_phase_space_spectrum(
            energy,
            rho_E,
            omega_values_J,
            Te_K=float(Te_K),
            delta_J=float(delta_J),
        )
    return int(iT), int(id_local), JS_block, JR_block


def scattering_phase_space_spectrum(
    energy_values_J: np.ndarray,
    rho_E: np.ndarray,
    omega_values_J: np.ndarray,
    *,
    Te_K: float,
    delta_J: float,
) -> np.ndarray:
    """Compute J_S(Omega) for one (Te, Delta, q) slice."""
    E = np.asarray(energy_values_J, dtype=float)
    rho = np.asarray(rho_E, dtype=float)
    omega = np.asarray(omega_values_J, dtype=float)
    if E.ndim != 1 or rho.shape != E.shape:
        raise ValueError("energy_values_J and rho_E must be one-dimensional arrays with the same shape.")

    f_E = fermi_positive_energy(E, Te_K)
    out = np.zeros_like(omega, dtype=float)

    E_max = float(E[-1])
    E_floor = max(float(E[1]) * 0.5, 1.0e-300)
    lower = max(delta_J, E_floor) if delta_J > 0.0 else E_floor

    for i, Om in enumerate(omega):
        Ep = E + Om
        mask = (E >= lower) & (Ep <= E_max)
        if delta_J > 0.0:
            mask &= Ep >= delta_J
        if np.count_nonzero(mask) < 2:
            out[i] = 0.0
            continue
        E_m = E[mask]
        Ep_m = Ep[mask]
        rho_E_m = rho[mask]
        rho_Ep = np.interp(Ep_m, E, rho, left=0.0, right=0.0)
        denom = np.maximum(E_m * Ep_m, E_floor * E_floor)
        coherence = 1.0 - delta_J * delta_J / denom
        f_Ep = fermi_positive_energy(Ep_m, Te_K)
        integrand = rho_E_m * rho_Ep * coherence * (f_E[mask] - f_Ep)
        integrand = np.nan_to_num(integrand, nan=0.0, posinf=0.0, neginf=0.0)
        out[i] = float(np.trapezoid(integrand, E_m))
    return np.maximum(out, 0.0)


def recombination_phase_space_spectrum(
    energy_values_J: np.ndarray,
    rho_E: np.ndarray,
    omega_values_J: np.ndarray,
    *,
    Te_K: float,
    delta_J: float,
) -> np.ndarray:
    """Compute J_R(Omega) for one (Te, Delta, q) slice."""
    E = np.asarray(energy_values_J, dtype=float)
    rho = np.asarray(rho_E, dtype=float)
    omega = np.asarray(omega_values_J, dtype=float)
    if E.ndim != 1 or rho.shape != E.shape:
        raise ValueError("energy_values_J and rho_E must be one-dimensional arrays with the same shape.")

    out = np.zeros_like(omega, dtype=float)
    if delta_J <= 0.0:
        return out

    E_max = float(E[-1])
    E_floor = max(float(E[1]) * 0.5, 1.0e-300)
    for i, Om in enumerate(omega):
        if Om <= 2.0 * delta_J:
            out[i] = 0.0
            continue
        Ep = Om - E
        mask = (
            (E >= delta_J)
            & (E <= Om - delta_J)
            & (Ep >= delta_J)
            & (Ep <= E_max)
        )
        if np.count_nonzero(mask) < 2:
            out[i] = 0.0
            continue
        E_m = E[mask]
        Ep_m = Ep[mask]
        rho_E_m = rho[mask]
        rho_Ep = np.interp(Ep_m, E, rho, left=0.0, right=0.0)
        denom = np.maximum(E_m * Ep_m, E_floor * E_floor)
        coherence = 1.0 + delta_J * delta_J / denom
        thermal_factor = pair_recombination_thermal_factor(E_m, Ep_m, Te_K)
        integrand = rho_E_m * rho_Ep * coherence * thermal_factor
        integrand = np.nan_to_num(integrand, nan=0.0, posinf=0.0, neginf=0.0)
        out[i] = float(np.trapezoid(integrand, E_m))
    return np.maximum(out, 0.0)


def fermi_positive_energy(energy_J: np.ndarray, T_K: float) -> np.ndarray:
    """Fermi function for positive quasiparticle energies."""
    if T_K <= 0.0:
        raise ValueError("T_K must be positive.")
    x = np.asarray(energy_J, dtype=float) / (K_B_J_K * T_K)
    x = np.clip(x, 0.0, 700.0)
    return 1.0 / (np.exp(x) + 1.0)


def pair_recombination_thermal_factor(
    E_J: np.ndarray,
    Ep_J: np.ndarray,
    T_K: float,
) -> np.ndarray:
    """Stable version of f(E) f(E') [exp((E+E')/kBT) - 1]."""
    if T_K <= 0.0:
        raise ValueError("T_K must be positive.")
    a = np.asarray(E_J, dtype=float) / (K_B_J_K * T_K)
    b = np.asarray(Ep_J, dtype=float) / (K_B_J_K * T_K)
    exp_minus_a = np.exp(-np.minimum(a, 700.0))
    exp_minus_b = np.exp(-np.minimum(b, 700.0))
    exp_minus_sum = np.exp(-np.minimum(a + b, 700.0))
    numerator = 1.0 - exp_minus_sum
    denominator = (1.0 + exp_minus_a) * (1.0 + exp_minus_b)
    return numerator / denominator


def phase_space_summary(catalog: PhaseSpaceCatalog) -> dict[str, Any]:
    """Return a compact summary dictionary."""
    JS = catalog.J_S_TdqO_J
    JR = catalog.J_R_TdqO_J
    window = phase_space_energy_window_summary(catalog)
    summary = {
        "backend": str(catalog.metadata.get("backend", "unknown")),
        "shape": list(catalog.shape),
        "n_Te": int(catalog.Te_values_K.size),
        "n_delta": int(catalog.delta_values_J.size),
        "n_q": int(catalog.q_values_m_inv.size),
        "n_omega": int(catalog.omega_values_J.size),
        "Te_min_K": float(np.min(catalog.Te_values_K)),
        "Te_max_K": float(np.max(catalog.Te_values_K)),
        "omega_max_J": float(np.max(catalog.omega_values_J)),
        "omega_max_meV": J_to_meV(float(np.max(catalog.omega_values_J))),
        "omega_axis_source": str(catalog.metadata.get("omega_axis_source", "")),
        "scattering_window_margin_J": float(catalog.metadata.get("scattering_window_margin_J", float("nan"))),
        "scattering_window_margin_meV": float(catalog.metadata.get("scattering_window_margin_meV", float("nan"))),
        "delta_min_meV": J_to_meV(float(np.min(catalog.delta_values_J))),
        "delta_max_meV": J_to_meV(float(np.max(catalog.delta_values_J))),
        "gamma_min_meV": J_to_meV(float(np.min(catalog.gamma_values_J))),
        "gamma_max_meV": J_to_meV(float(np.max(catalog.gamma_values_J))),
        "q_min_m_inv": float(np.min(catalog.q_values_m_inv)),
        "q_max_m_inv": float(np.max(catalog.q_values_m_inv)),
        "J_S_min_J": float(np.min(JS)),
        "J_S_max_J": float(np.max(JS)),
        "J_S_is_finite": bool(np.all(np.isfinite(JS))),
        "J_R_min_J": float(np.min(JR)),
        "J_R_max_J": float(np.max(JR)),
        "J_R_is_finite": bool(np.all(np.isfinite(JR))),
        "grid_is_downsampled": bool(catalog.metadata.get("grid_is_downsampled", False)),
        "parallel_workers": int(catalog.metadata.get("phase_space_parallel_workers", 1)),
        "parallel_backend": str(catalog.metadata.get("phase_space_parallel_backend", "serial")),
        "parallel_tasks": int(catalog.metadata.get("phase_space_parallel_tasks", 1)),
        "parallel_cells": int(catalog.metadata.get("phase_space_parallel_cells", 1)),
        "parallel_task_layout": str(catalog.metadata.get("phase_space_parallel_task_layout", "")),
        "normal_limit_policy": str(catalog.metadata.get("normal_limit_policy", "")),
        "threshold_policy": str(catalog.metadata.get("threshold_policy", "")),
        "energy_window_policy": str(catalog.metadata.get("energy_window_policy", "")),
        "coherence_factor_policy": str(catalog.metadata.get("coherence_factor_policy", "")),
        "thermal_closure_policy": str(catalog.metadata.get("thermal_closure_policy", "")),
    }
    summary.update(window)
    return summary


def phase_space_energy_window_summary(catalog: PhaseSpaceCatalog) -> dict[str, Any]:
    """Summarize finite-energy-window limitations of the phase-space catalogue."""
    omega_max_J = float(np.max(catalog.omega_values_J))
    energy_max_J = float(catalog.metadata.get("energy_max_J", omega_max_J))
    js_cutoff = np.asarray(catalog.metadata.get("js_hard_cutoff_by_delta_J", []), dtype=float)
    if js_cutoff.size != catalog.delta_values_J.size:
        energy_floor_J = float(catalog.metadata.get("energy_floor_J", 1.0e-300))
        lower = np.maximum(catalog.delta_values_J, energy_floor_J)
        js_cutoff = np.maximum(0.0, energy_max_J - lower)
    jr_threshold = np.asarray(catalog.metadata.get("jr_threshold_by_delta_J", []), dtype=float)
    if jr_threshold.size != catalog.delta_values_J.size:
        jr_threshold = np.maximum(0.0, 2.0 * catalog.delta_values_J)
    finite_js = js_cutoff[np.isfinite(js_cutoff)]
    finite_jr = jr_threshold[np.isfinite(jr_threshold)]
    js_min = float(np.min(finite_js)) if finite_js.size else float("nan")
    js_max = float(np.max(finite_js)) if finite_js.size else float("nan")
    jr_min = float(np.min(finite_jr)) if finite_jr.size else float("nan")
    jr_max = float(np.max(finite_jr)) if finite_jr.size else float("nan")
    return {
        "energy_max_J": energy_max_J,
        "energy_max_meV": J_to_meV(energy_max_J),
        "J_S_hard_cutoff_min_J": js_min,
        "J_S_hard_cutoff_min_meV": J_to_meV(js_min) if np.isfinite(js_min) else float("nan"),
        "J_S_hard_cutoff_max_J": js_max,
        "J_S_hard_cutoff_max_meV": J_to_meV(js_max) if np.isfinite(js_max) else float("nan"),
        "J_R_threshold_min_J": jr_min,
        "J_R_threshold_min_meV": J_to_meV(jr_min) if np.isfinite(jr_min) else float("nan"),
        "J_R_threshold_max_J": jr_max,
        "J_R_threshold_max_meV": J_to_meV(jr_max) if np.isfinite(jr_max) else float("nan"),
        "scattering_window_is_truncated": bool(omega_max_J > js_min) if np.isfinite(js_min) else False,
    }


def save_phase_space_catalog_npz(catalog: PhaseSpaceCatalog, path: str | Path) -> Path:
    """Save phase-space catalogue to compressed ``.npz``."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        Te_values_K=catalog.Te_values_K,
        omega_values_J=catalog.omega_values_J,
        delta_values_J=catalog.delta_values_J,
        gamma_values_J=catalog.gamma_values_J,
        q_values_m_inv=catalog.q_values_m_inv,
        J_S_TdqO_J=catalog.J_S_TdqO_J,
        J_R_TdqO_J=catalog.J_R_TdqO_J,
        delta_indices=catalog.delta_indices,
        q_indices=catalog.q_indices,
        metadata=np.array(catalog.metadata, dtype=object),
    )
    return output


def load_phase_space_catalog_npz(path: str | Path) -> PhaseSpaceCatalog:
    """Load a phase-space catalogue saved by :func:`save_phase_space_catalog_npz`."""
    source = Path(path)
    with np.load(source, allow_pickle=True) as data:
        metadata = data["metadata"].item()
        return PhaseSpaceCatalog(
            Te_values_K=np.asarray(data["Te_values_K"], dtype=float),
            omega_values_J=np.asarray(data["omega_values_J"], dtype=float),
            delta_values_J=np.asarray(data["delta_values_J"], dtype=float),
            gamma_values_J=np.asarray(data["gamma_values_J"], dtype=float),
            q_values_m_inv=np.asarray(data["q_values_m_inv"], dtype=float),
            J_S_TdqO_J=np.asarray(data["J_S_TdqO_J"], dtype=float),
            J_R_TdqO_J=np.asarray(data["J_R_TdqO_J"], dtype=float),
            delta_indices=np.asarray(data["delta_indices"], dtype=np.int64),
            q_indices=np.asarray(data["q_indices"], dtype=np.int64),
            metadata=dict(metadata),
        )


def _select_axis_indices(n_total: int, n_requested: int) -> np.ndarray:
    """Select approximately evenly spaced indices from an axis."""
    if n_total <= 0:
        raise ValueError("n_total must be positive.")
    if n_requested <= 0:
        raise ValueError("n_requested must be positive.")
    if n_requested >= n_total:
        return np.arange(n_total, dtype=np.int64)
    return np.unique(np.round(np.linspace(0, n_total - 1, n_requested)).astype(np.int64))


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


class _PhaseSpaceProgress:
    """Single-line dependency-free progress bar for phase-space chunks."""

    def __init__(
        self,
        *,
        total_chunks: int,
        cells_per_chunk: int,
        workers: int,
        backend: str,
        enabled: bool,
        width: int = 34,
    ) -> None:
        self.total_chunks = max(1, int(total_chunks))
        self.cells_per_chunk = max(1, int(cells_per_chunk))
        self.total_cells = self.total_chunks * self.cells_per_chunk
        self.workers = int(workers)
        self.backend = str(backend)
        self.enabled = bool(enabled)
        self.width = int(width)
        self.done_chunks = 0
        self._last_percent = -1
        self._prefix = (
            "Phase-space "
            f"({self.total_chunks} chunks, {self.total_cells} cells, "
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
        cells_done = self.done_chunks * self.cells_per_chunk
        line = (
            f"\r{self._prefix}: [{bar}] {percent:3d}% "
            f"chunks={self.done_chunks}/{self.total_chunks} "
            f"cells={cells_done}/{self.total_cells}"
        )
        print(line, end="\n" if final else "", flush=True)


__all__ = [
    "PhaseSpaceCatalog",
    "build_phase_space_catalog_from_usadel_catalog",
    "scattering_phase_space_spectrum",
    "recombination_phase_space_spectrum",
    "fermi_positive_energy",
    "pair_recombination_thermal_factor",
    "phase_space_summary",
    "phase_space_energy_window_summary",
    "save_phase_space_catalog_npz",
    "load_phase_space_catalog_npz",
]
