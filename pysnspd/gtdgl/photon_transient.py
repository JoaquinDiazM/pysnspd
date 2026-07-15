"""Conservative circuit + photon transient driver for pySNSPD.

The implementation is intentionally an orchestration layer around the already
validated ``solve_stationary_pytdgl_like`` adapter.  It runs short mesoscopic
chunks, updates the lumped circuit, optionally applies one phonon-bubble event,
and saves a compact transient history.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Mapping

import numpy as np
import yaml

from pysnspd.circuit.readout import (
    CircuitParams,
    CircuitState,
    central_tdgl_voltage_V,
    circuit_observables,
    initialize_circuit_from_ss,
    step_circuit_rk2,
)
from pysnspd.gtdgl import solve_stationary_pytdgl_like
from pysnspd.gtdgl.photon import PhotonBubbleParams, inject_phonon_bubble
from pysnspd.gtdgl.snapshot_diagnostics import save_ss_snapshot_bundle_npz


@dataclass(frozen=True)
class CoupledTransientConfig:
    total_time_s: float
    mesoscopic_dt_s: float
    chunk_time_s: float
    n_snapshots: int = 8
    center_voltage_width_m: float = 100.0e-9
    center_voltage_probe_band_m: float | None = None
    thermal_enabled: bool = True
    thermal_window_m: float = 100.0e-9
    thermal_max_step_K: float = 0.05
    thermal_max_substeps: int = 64
    terminal_psi: float = 0.0
    terminal_healing_xi: float | None = None
    terminal_healing_fraction: float = 0.95
    supercurrent_law: str = "usadel_poisson"
    allmaras_phase_direct_amplitude_fraction: float = 1.0e-2
    allmaras_phase_convergence_tol: float = 1.0e-3
    allmaras_phase_convergence_max_iterations: int = 64
    progress: bool = False

    def validated(self) -> "CoupledTransientConfig":
        for key in ("total_time_s", "mesoscopic_dt_s", "chunk_time_s"):
            value = float(getattr(self, key))
            if not np.isfinite(value) or value <= 0.0:
                raise ValueError(f"{key} must be positive and finite.")
        if int(self.n_snapshots) <= 0:
            raise ValueError("n_snapshots must be positive.")
        if not (0.0 < float(self.allmaras_phase_direct_amplitude_fraction) < 1.0):
            raise ValueError("allmaras_phase_direct_amplitude_fraction must lie in (0, 1).")
        if float(self.allmaras_phase_convergence_tol) <= 0.0:
            raise ValueError("allmaras_phase_convergence_tol must be positive.")
        if int(self.allmaras_phase_convergence_max_iterations) < 1:
            raise ValueError("allmaras_phase_convergence_max_iterations must be at least one.")
        return self


def run_coupled_transient(
    *,
    mesh: Any,
    edge_data: Any,
    ops: Any,
    material: Any,
    initial_state_npz: str | Path,
    initial_current_A: float,
    usadel_catalog: Any,
    power_table_npz: str | Path | None,
    output_dir: str | Path,
    config: CoupledTransientConfig,
    circuit_params: CircuitParams,
    photon_params: PhotonBubbleParams,
) -> dict[str, Any]:
    cfg = config.validated()
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    current = _load_initial_state_as_seed(initial_state_npz)
    nodes_m = np.asarray(mesh.nodes, dtype=float)[:, :2]

    V_tdgl0 = central_tdgl_voltage_V(
        nodes_m=nodes_m,
        phi_V=current.node_phi_electric_V,
        center_width_m=float(cfg.center_voltage_width_m),
        probe_band_m=cfg.center_voltage_probe_band_m,
    )

    circuit, circuit_params = initialize_circuit_from_ss(
        I_ss_A=float(initial_current_A),
        V_tdgl_ss_V=float(V_tdgl0),
        params=circuit_params,
    )

    snapshot_times = np.linspace(0.0, float(cfg.total_time_s), max(2, int(cfg.n_snapshots)))
    next_snapshot = 0
    snapshots: dict[str, list[np.ndarray | float]] = {
        "snapshot_t_s": [],
        "psi_real_snapshot_J": [],
        "psi_imag_snapshot_J": [],
        "delta_snapshot_meV": [],
        "phi_snapshot_V": [],
        "Te_snapshot_K": [],
        "Tph_snapshot_K": [],
    }

    history: dict[str, list[float | int | bool]] = {
        "t_s": [],
        "I_b_A": [],
        "I_s_A": [],
        "I_rf_A": [],
        "V_out_V": [],
        "v_c_V": [],
        "V_tdgl_center_V": [],
        "chunk_dt_s": [],
        "photon_applied": [],
        "min_delta_over_delta0": [],
        "mean_delta_over_delta0": [],
        "max_pairbreaking_ratio": [],
        "max_Te_K": [],
        "max_Tph_K": [],
    }

    photon_applied = False
    photon_metadata: dict[str, Any] = {
        "enabled": bool(photon_params.enabled),
        "applied": False,
        "energy_eV": float(photon_params.energy_eV),
    }

    t = 0.0
    _append_snapshot_if_due(
        t,
        current,
        material=material,
        snapshot_times=snapshot_times,
        next_snapshot_ref=[next_snapshot],
        snapshots=snapshots,
    )
    next_snapshot = len(snapshots["snapshot_t_s"])

    progress = _ProgressReporter(
        enabled=bool(cfg.progress),
        total_ps=float(cfg.total_time_s) / 1.0e-12,
        description="Pipeline 03 transient",
    )

    max_loop_count = int(np.ceil(float(cfg.total_time_s) / max(float(cfg.chunk_time_s), 1.0e-30))) + 8
    loop_count = 0

    while t < float(cfg.total_time_s) - 1.0e-30:
        loop_count += 1
        if loop_count > max_loop_count:
            progress.close()
            raise RuntimeError(
                "Pipeline 03 transient loop exceeded the expected number of chunks; "
                "check chunk_time_s, photon_time_s and total_time_s."
            )

        remaining = float(cfg.total_time_s) - t
        if remaining <= 1.0e-30:
            break

        dt_chunk = min(float(cfg.chunk_time_s), remaining)

        # Split exactly at photon time, so the state change happens between two
        # valid mesoscopic chunks.  With energy_eV=0 this tests the machinery
        # without changing Tph.
        if bool(photon_params.enabled) and not photon_applied:
            t_ph = float(photon_params.time_s)
            if t < t_ph < t + dt_chunk:
                dt_chunk = t_ph - t

        if dt_chunk > 0.0:
            result = solve_stationary_pytdgl_like(
                mesh=mesh,
                edge_data=edge_data,
                seed=current,
                material=material,
                ops=ops,
                total_time_s=float(dt_chunk),
                dt_s=float(cfg.mesoscopic_dt_s),
                target_current_A=float(circuit.I_s_A),
                usadel_catalog=usadel_catalog,
                terminal_psi=float(cfg.terminal_psi),
                adaptive=True,
                n_snapshots=2,
                progress=False,
                supercurrent_law=str(cfg.supercurrent_law),
                terminal_healing_xi=cfg.terminal_healing_xi,
                terminal_healing_fraction=float(cfg.terminal_healing_fraction),
                stationarity_eta=1.0e-5,
                convergence_min_steps=500,
                stop_on_convergence=False,
                allmaras_phase_direct_amplitude_fraction=float(cfg.allmaras_phase_direct_amplitude_fraction),
                allmaras_phase_convergence_tol=float(cfg.allmaras_phase_convergence_tol),
                allmaras_phase_convergence_max_iterations=int(cfg.allmaras_phase_convergence_max_iterations),
                thermal_enabled=bool(cfg.thermal_enabled and power_table_npz is not None),
                thermal_power_table_npz=(None if power_table_npz is None else str(power_table_npz)),
                thermal_window_m=float(cfg.thermal_window_m),
                thermal_start_time_s=0.0,
                thermal_bath_K=float(np.nanmedian(current.node_Tph_K)),
                thermal_min_K=float(np.nanmin(current.node_Tph_K)),
                thermal_max_K=None,
                thermal_max_step_K=float(cfg.thermal_max_step_K),
                thermal_max_substeps=int(cfg.thermal_max_substeps),
            )

            state = result.state
            current = _state_to_seed(state)
            t += float(dt_chunk)

            V_tdgl = central_tdgl_voltage_V(
                nodes_m=nodes_m,
                phi_V=current.node_phi_electric_V,
                center_width_m=float(cfg.center_voltage_width_m),
                probe_band_m=cfg.center_voltage_probe_band_m,
            )
            circuit = step_circuit_rk2(
                circuit,
                V_tdgl_V=float(V_tdgl),
                dt_s=float(dt_chunk),
                params=circuit_params,
            )
            obs = circuit_observables(circuit, params=circuit_params, V_tdgl_V=V_tdgl)
            _append_history(history, t=t, obs=obs, state=state, material=material, dt_chunk=dt_chunk, photon_applied=False)
            progress.update(
                dt_ps=float(dt_chunk) / 1.0e-12,
                I_s_uA=float(obs["I_s_A"]) * 1.0e6,
                V_out_uV=float(obs["V_out_V"]) * 1.0e6,
                Vtdgl_uV=float(obs["V_tdgl_center_V"]) * 1.0e6,
                photon=photon_applied,
            )
            _append_snapshot_if_due(
                t,
                current,
                material=material,
                snapshot_times=snapshot_times,
                next_snapshot_ref=[next_snapshot],
                snapshots=snapshots,
            )
            next_snapshot = len(snapshots["snapshot_t_s"])

        if bool(photon_params.enabled) and not photon_applied and t >= float(photon_params.time_s) - 1.0e-30:
            if power_table_npz is None:
                raise FileNotFoundError("Photon bubble requires raw_pre/power_table_catalog.npz.")
            Tph_new, photon_metadata = inject_phonon_bubble(
                mesh=mesh,
                Tph_K=current.node_Tph_K,
                power_table_npz=str(power_table_npz),
                thickness_m=float(material.thickness_m),
                params=photon_params,
            )
            current = _copy_seed_with_Tph(current, Tph_new)
            photon_applied = True
            # Record a zero-duration marker after the state mutation.
            V_tdgl = central_tdgl_voltage_V(
                nodes_m=nodes_m,
                phi_V=current.node_phi_electric_V,
                center_width_m=float(cfg.center_voltage_width_m),
                probe_band_m=cfg.center_voltage_probe_band_m,
            )
            obs = circuit_observables(circuit, params=circuit_params, V_tdgl_V=V_tdgl)
            history["t_s"].append(float(t))
            history["I_b_A"].append(float(obs["I_b_A"]))
            history["I_s_A"].append(float(obs["I_s_A"]))
            history["I_rf_A"].append(float(obs["I_rf_A"]))
            history["V_out_V"].append(float(obs["V_out_V"]))
            history["v_c_V"].append(float(obs["v_c_V"]))
            history["V_tdgl_center_V"].append(float(obs["V_tdgl_center_V"]))
            history["chunk_dt_s"].append(0.0)
            history["photon_applied"].append(True)
            abspsi = np.abs(current.node_psi_real_J + 1j * current.node_psi_imag_J)
            history["min_delta_over_delta0"].append(float(np.nanmin(abspsi) / material.delta0_J))
            history["mean_delta_over_delta0"].append(float(np.nanmean(abspsi) / material.delta0_J))
            history["max_pairbreaking_ratio"].append(float("nan"))
            history["max_Te_K"].append(float(np.nanmax(current.node_Te_K)))
            history["max_Tph_K"].append(float(np.nanmax(current.node_Tph_K)))
            progress.mark_photon(float(t) / 1.0e-12, photon_metadata)

    progress.close()

    # Always append final state.
    _append_snapshot_force(float(cfg.total_time_s), current, material=material, snapshots=snapshots)

    final_state_npz = _save_final_state_like_stationary(current, out / "final_state.npz")
    history_npz = _save_transient_history(history, out / "transient_history.npz")
    snapshots_npz = _save_transient_snapshots(snapshots, ops=ops, output_path=out / "transient_snapshots.npz")

    summary = {
        "backend": "pipeline03_split_coupled_circuit_photon_v1",
        "initial_state_npz": str(initial_state_npz),
        "initial_current_A": float(initial_current_A),
        "initial_V_tdgl_center_V": float(V_tdgl0),
        "final_time_ps": float(cfg.total_time_s / 1.0e-12),
        "n_history_rows": int(len(history["t_s"])),
        "n_snapshots": int(len(snapshots["snapshot_t_s"])),
        "config": asdict(cfg),
        "circuit": {
            "params": circuit_params.as_dict(),
            "final_state": circuit.as_dict(),
        },
        "photon": photon_metadata,
        "outputs": {
            "final_state_npz": str(final_state_npz),
            "transient_history_npz": str(history_npz),
            "transient_snapshots_npz": str(snapshots_npz),
        },
    }
    summary_path = out / "photon_summary.yaml"
    with summary_path.open("w", encoding="utf-8") as f:
        yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)

    return summary


def _load_initial_state_as_seed(path: str | Path) -> SimpleNamespace:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Missing SS stationary state: {p}")
    with np.load(p, allow_pickle=True) as data:
        psi = np.asarray(data["psi_real_J"], dtype=float) + 1j * np.asarray(data["psi_imag_J"], dtype=float)
        return SimpleNamespace(
            node_psi_real_J=np.real(psi).copy(),
            node_psi_imag_J=np.imag(psi).copy(),
            node_phi_electric_V=np.asarray(data["phi_V"], dtype=float).copy(),
            node_Te_K=np.asarray(data.get("Te_K", np.zeros_like(np.real(psi))), dtype=float).copy(),
            node_Tph_K=np.asarray(data.get("Tph_K", np.zeros_like(np.real(psi))), dtype=float).copy(),
        )


def _state_to_seed(state: Any) -> SimpleNamespace:
    return SimpleNamespace(
        node_psi_real_J=np.real(state.psi_J).copy(),
        node_psi_imag_J=np.imag(state.psi_J).copy(),
        node_phi_electric_V=np.asarray(state.phi_V, dtype=float).copy(),
        node_Te_K=np.asarray(state.Te_K, dtype=float).copy(),
        node_Tph_K=np.asarray(state.Tph_K, dtype=float).copy(),
        currents=state.currents,
        metadata=getattr(state, "metadata", {}),
    )


def _copy_seed_with_Tph(seed: SimpleNamespace, Tph_K: np.ndarray) -> SimpleNamespace:
    return SimpleNamespace(
        node_psi_real_J=np.asarray(seed.node_psi_real_J, dtype=float).copy(),
        node_psi_imag_J=np.asarray(seed.node_psi_imag_J, dtype=float).copy(),
        node_phi_electric_V=np.asarray(seed.node_phi_electric_V, dtype=float).copy(),
        node_Te_K=np.asarray(seed.node_Te_K, dtype=float).copy(),
        node_Tph_K=np.asarray(Tph_K, dtype=float).copy(),
        currents=getattr(seed, "currents", None),
        metadata=getattr(seed, "metadata", {}),
    )


def _append_history(
    history: dict[str, list],
    *,
    t: float,
    obs: Mapping[str, float],
    state: Any,
    material: Any,
    dt_chunk: float,
    photon_applied: bool,
) -> None:
    history["t_s"].append(float(t))
    history["I_b_A"].append(float(obs["I_b_A"]))
    history["I_s_A"].append(float(obs["I_s_A"]))
    history["I_rf_A"].append(float(obs["I_rf_A"]))
    history["V_out_V"].append(float(obs["V_out_V"]))
    history["v_c_V"].append(float(obs["v_c_V"]))
    history["V_tdgl_center_V"].append(float(obs["V_tdgl_center_V"]))
    history["chunk_dt_s"].append(float(dt_chunk))
    history["photon_applied"].append(bool(photon_applied))
    abspsi = np.abs(np.asarray(state.psi_J, dtype=np.complex128))
    history["min_delta_over_delta0"].append(float(np.nanmin(abspsi) / material.delta0_J))
    history["mean_delta_over_delta0"].append(float(np.nanmean(abspsi) / material.delta0_J))
    try:
        history["max_pairbreaking_ratio"].append(float(np.nanmax(state.currents.node_pairbreaking_ratio)))
    except Exception:
        history["max_pairbreaking_ratio"].append(float("nan"))
    history["max_Te_K"].append(float(np.nanmax(state.Te_K)))
    history["max_Tph_K"].append(float(np.nanmax(state.Tph_K)))


def _append_snapshot_if_due(
    t: float,
    seed: SimpleNamespace,
    *,
    material: Any,
    snapshot_times: np.ndarray,
    next_snapshot_ref: list[int],
    snapshots: dict[str, list],
) -> None:
    idx = int(next_snapshot_ref[0])
    while idx < snapshot_times.size and t >= float(snapshot_times[idx]) - 1.0e-30:
        _append_snapshot_force(float(t), seed, material=material, snapshots=snapshots)
        idx += 1
    next_snapshot_ref[0] = idx


def _append_snapshot_force(
    t: float,
    seed: SimpleNamespace,
    *,
    material: Any,
    snapshots: dict[str, list],
) -> None:
    psi = np.asarray(seed.node_psi_real_J, dtype=float) + 1j * np.asarray(seed.node_psi_imag_J, dtype=float)
    snapshots["snapshot_t_s"].append(float(t))
    snapshots["psi_real_snapshot_J"].append(np.real(psi).copy())
    snapshots["psi_imag_snapshot_J"].append(np.imag(psi).copy())
    snapshots["delta_snapshot_meV"].append(np.abs(psi) / 1.602176634e-22)
    snapshots["phi_snapshot_V"].append(np.asarray(seed.node_phi_electric_V, dtype=float).copy())
    snapshots["Te_snapshot_K"].append(np.asarray(seed.node_Te_K, dtype=float).copy())
    snapshots["Tph_snapshot_K"].append(np.asarray(seed.node_Tph_K, dtype=float).copy())


def _save_final_state_like_stationary(seed: SimpleNamespace, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_path,
        psi_real_J=np.asarray(seed.node_psi_real_J, dtype=float),
        psi_imag_J=np.asarray(seed.node_psi_imag_J, dtype=float),
        phi_V=np.asarray(seed.node_phi_electric_V, dtype=float),
        Te_K=np.asarray(seed.node_Te_K, dtype=float),
        Tph_K=np.asarray(seed.node_Tph_K, dtype=float),
    )
    return output_path


def _save_transient_history(history: dict[str, list], output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {key: np.asarray(value) for key, value in history.items()}
    if "t_s" in arrays:
        arrays["t_ps"] = np.asarray(arrays["t_s"], dtype=float) / 1.0e-12
    np.savez_compressed(output_path, **arrays)
    return output_path


def _save_transient_snapshots(snapshots: dict[str, list], *, ops: Any, output_path: Path) -> Path:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {}
    for key, values in snapshots.items():
        arrays[key] = np.asarray(values)
    if "snapshot_t_s" in arrays:
        arrays["snapshot_t_ps"] = np.asarray(arrays["snapshot_t_s"], dtype=float) / 1.0e-12

    # Keep the same static FV topology keys used by SS snapshot diagnostics.
    for key in ("edge_i", "edge_j", "edge_length_m", "edge_unit_x", "edge_unit_y", "dual_face_length_m"):
        if hasattr(ops, key):
            arrays[key] = np.asarray(getattr(ops, key))

    np.savez_compressed(output_path, **arrays)
    return output_path



class _ProgressReporter:
    """Optional progress wrapper for pipeline 03.

    Uses tqdm when available.  If tqdm is unavailable, it falls back to sparse
    text updates, so the transient still runs on minimal environments.
    """

    def __init__(self, *, enabled: bool, total_ps: float, description: str) -> None:
        self.enabled = bool(enabled)
        self.total_ps = float(total_ps)
        self.current_ps = 0.0
        self._bar = None
        self._last_print_fraction = -1.0

        if not self.enabled:
            return

        try:
            from tqdm.auto import tqdm  # type: ignore

            self._bar = tqdm(
                total=self.total_ps,
                desc=description,
                unit="ps",
                dynamic_ncols=True,
                leave=True,
                mininterval=0.5,
            )
        except Exception:
            self._bar = None
            print(f"{description}: 0.000 / {self.total_ps:.3f} ps")

    def update(
        self,
        *,
        dt_ps: float,
        I_s_uA: float,
        V_out_uV: float,
        Vtdgl_uV: float,
        photon: bool,
    ) -> None:
        if not self.enabled:
            return

        dt = max(float(dt_ps), 0.0)
        self.current_ps = min(self.current_ps + dt, self.total_ps)
        postfix = {
            "I_s[uA]": f"{float(I_s_uA):.6g}",
            "Vout[uV]": f"{float(V_out_uV):.3g}",
            "Vtdgl[uV]": f"{float(Vtdgl_uV):.3g}",
            "photon": "yes" if bool(photon) else "no",
        }

        if self._bar is not None:
            self._bar.update(dt)
            self._bar.set_postfix(postfix)
            return

        fraction = self.current_ps / max(self.total_ps, 1.0e-300)
        if fraction >= self._last_print_fraction + 0.10 or fraction >= 1.0:
            self._last_print_fraction = fraction
            print(
                f"Pipeline 03 transient: {self.current_ps:.3f}/{self.total_ps:.3f} ps "
                f"({100.0 * fraction:5.1f}%) "
                f"I_s={float(I_s_uA):.6g} uA "
                f"Vout={float(V_out_uV):.3g} uV "
                f"Vtdgl={float(Vtdgl_uV):.3g} uV "
                f"photon={'yes' if bool(photon) else 'no'}"
            )

    def mark_photon(self, t_ps: float, metadata: Mapping[str, Any]) -> None:
        if not self.enabled:
            return
        energy = metadata.get("energy_eV", float("nan"))
        reason = metadata.get("reason", "photon event")
        msg = f"photon event at t={float(t_ps):.6g} ps, E={energy} eV: {reason}"
        if self._bar is not None:
            self._bar.write(msg)
        else:
            print(msg)

    def close(self) -> None:
        if self._bar is not None:
            self._bar.close()


__all__ = [
    "CoupledTransientConfig",
    "run_coupled_transient",
]
