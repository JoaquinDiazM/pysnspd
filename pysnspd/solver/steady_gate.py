"""Online total-stationarity gate for long SS runs."""

from __future__ import annotations

from typing import Any

import numpy as np

from pysnspd.circuit.readout import circuit_stationarity_diagnostics
from pysnspd.gtdgl.currents import native_edge_currents_to_current_fields
from pysnspd.solver.callbacks import _terminal_edge_mask_from_device
from pysnspd.solver.targets import (
    contact_recovery_diagnostics,
    continuity_diagnostics,
    dynamic_stationarity_diagnostics,
    stationarity_diagnostics,
)
from pysnspd.thermal.evolution import thermal_stationarity_diagnostics


def classify_photon_readiness(
    *,
    strict_stationarity_passes: bool,
    dynamic_stationarity_passes: bool,
    contact_recovery_passes: bool,
    continuity_passes: bool,
    thermal_stationarity_passes: bool,
    circuit_stationarity_passes: bool,
    phase_drive_converged: bool,
) -> dict[str, bool | str]:
    """Combine hard validity gates with alternative mesoscopic SS modes."""

    mesoscopic_passes = bool(
        strict_stationarity_passes or dynamic_stationarity_passes
    )
    mesoscopic_mode = (
        "strict_fixed_point"
        if strict_stationarity_passes
        else (
            "weak_dynamic_attractor"
            if dynamic_stationarity_passes
            else "none"
        )
    )
    passes = bool(
        mesoscopic_passes
        and contact_recovery_passes
        and continuity_passes
        and thermal_stationarity_passes
        and circuit_stationarity_passes
        and phase_drive_converged
    )
    return {
        "passes": passes,
        "mesoscopic_passes": mesoscopic_passes,
        "mesoscopic_mode": mesoscopic_mode,
    }


def build_total_stationarity_callback(
    *,
    mesh,
    edge_data,
    ops,
    material,
    device,
    target_current_A: float,
    circuit_controller,
    circuit_runtime,
    thermal_config,
    Te_K: np.ndarray,
    Tph_K: np.ndarray,
    terminal_healing_fraction: float,
    recovery_min_xi: float,
    recovery_max_xi: float,
    stationarity_phase_gradient_rel: float,
    stationarity_phi_gradient_rel: float,
    stationarity_q_abs_m_inv: float,
    stationarity_phi_gradient_abs_V_m: float,
    stationarity_edge_active_threshold: float,
    stationarity_bulk_exclusion_xi: float,
    dynamic_stationarity_minimum_tail_ps: float,
    dynamic_stationarity_profile_rel: float,
    dynamic_stationarity_voltage_rel: float,
    dynamic_stationarity_psl_threshold: float,
    continuity_rms_tol: float,
    continuity_max_tol: float,
    continuity_poisson_tol: float,
    thermal_stationarity_relative_tolerance: float,
    thermal_stationarity_p99_tolerance_K: float,
    thermal_stationarity_projection_horizon_ps: float,
    thermal_stationarity_projection_relative_tolerance: float,
    photon_ready_consecutive_evaluations: int,
    tau0: float,
    state: dict[str, Any],
):
    """Build the snapshot-rate gate used for SS early termination."""

    nodes = np.asarray(mesh.nodes, dtype=float)[:, :2]
    edge_i = np.asarray(ops.edge_i, dtype=np.int64)
    edge_j = np.asarray(ops.edge_j, dtype=np.int64)
    edge_length = np.maximum(
        np.linalg.norm(nodes[edge_j] - nodes[edge_i], axis=1),
        1.0e-300,
    )
    edge_center_x = 0.5 * (nodes[edge_i, 0] + nodes[edge_j, 0])
    xmin = float(np.nanmin(nodes[:, 0]))
    xmax = float(np.nanmax(nodes[:, 0]))
    edge_distance = np.minimum(
        np.maximum(edge_center_x - xmin, 0.0),
        np.maximum(xmax - edge_center_x, 0.0),
    )
    blocked_edges = _terminal_edge_mask_from_device(device, ops)
    frames: dict[str, list[Any]] = {
        "t_s": [],
        "psi_real": [],
        "psi_imag": [],
        "phi": [],
        "edge_q": [],
        "terminal_voltage": [],
        "Te": [],
        "Tph": [],
    }
    state["last_evaluation_time_ps"] = None
    state["passes"] = False
    state["instantaneous_passes"] = False
    state["consecutive_passes"] = 0
    state["required_consecutive_passes"] = max(
        1,
        int(photon_ready_consecutive_evaluations),
    )

    def callback(*, time: float, frame, running_state) -> str | None:
        time_s = float(time) * float(tau0)
        phi_V = np.asarray(frame.mu, dtype=float) * device.voltage_scale_V
        phi_V = phi_V - float(np.mean(phi_V))
        current_A = (
            float(circuit_controller.state.I_s_A)
            if circuit_controller is not None
            else float(target_current_A)
        )
        currents, _ = native_edge_currents_to_current_fields(
            psi_dimensionless=np.asarray(frame.psi, dtype=np.complex128),
            native_supercurrent=np.asarray(frame.supercurrent, dtype=float),
            native_normal_current=np.asarray(frame.normal_current, dtype=float),
            device=device,
            mesh=mesh,
            edge_data=edge_data,
            ops=ops,
            material=material,
            Te_K=Te_K,
            target_current_A=current_A,
        )
        frames["t_s"].append(time_s)
        frames["psi_real"].append(np.real(frame.psi) * material.delta0_J)
        frames["psi_imag"].append(np.imag(frame.psi) * material.delta0_J)
        frames["phi"].append(phi_V)
        frames["edge_q"].append(np.asarray(currents.edge_Q_m_inv, dtype=float))
        frames["terminal_voltage"].append(float(np.ptp(phi_V)))
        frames["Te"].append(np.asarray(Te_K, dtype=float).copy())
        frames["Tph"].append(np.asarray(Tph_K, dtype=float).copy())
        retain_window_s = max(
            float(dynamic_stationarity_minimum_tail_ps) * 1.0e-12,
            float(circuit_runtime.hold_time_s),
        ) + 1.0e-12
        while (
            len(frames["t_s"]) > 3
            and float(frames["t_s"][1]) < time_s - retain_window_s
        ):
            for values in frames.values():
                values.pop(0)

        minimum_time_s = max(
            float(circuit_runtime.start_time_s),
            float(thermal_config.start_time_s),
        ) + max(
            float(circuit_runtime.hold_time_s),
            float(dynamic_stationarity_minimum_tail_ps) * 1.0e-12,
            0.5e-12,
        )
        last_eval = state.get("last_evaluation_time_s")
        if time_s < minimum_time_s:
            return None
        if last_eval is not None and time_s - float(last_eval) < 0.5e-12:
            return None
        state["last_evaluation_time_s"] = time_s
        state["last_evaluation_time_ps"] = time_s / 1.0e-12

        xi2 = np.asarray(material.xi_mod_squared_m2(Te_K), dtype=float)
        xi_m = float(np.sqrt(np.nanmedian(np.maximum(xi2, 1.0e-300))))
        raw = running_state.data
        online_history = {
            "snapshot_t_s": np.asarray(frames["t_s"], dtype=float),
            "psi_snapshot_real_J": np.asarray(frames["psi_real"], dtype=float),
            "psi_snapshot_imag_J": np.asarray(frames["psi_imag"], dtype=float),
            "phi_snapshot_V": np.asarray(frames["phi"], dtype=float),
            "edge_phase_gradient_snapshot_m_inv": np.asarray(
                frames["edge_q"],
                dtype=float,
            ),
            "terminal_voltage_V": np.asarray(frames["terminal_voltage"], dtype=float),
            "Te_snapshot_K": np.asarray(frames["Te"], dtype=float),
            "Tph_snapshot_K": np.asarray(frames["Tph"], dtype=float),
            "edge_i": edge_i,
            "edge_j": edge_j,
            "edge_length_m": edge_length,
            "edge_distance_from_contact_m": edge_distance,
            "stationarity_xi_m": np.asarray([xi_m]),
            "normal_terminal_edge_mask": blocked_edges,
            "eta_R": np.asarray(raw.get("max_d_abs_sq_psi", []), dtype=float),
            "t_s": np.asarray(frames["t_s"], dtype=float),
            "pytdgl_like_poisson_residual_rel": np.asarray(
                raw.get("poisson_residual_rel", []),
                dtype=float,
            ),
        }
        stationarity = stationarity_diagnostics(
            history=online_history,
            material=material,
            phase_gradient_rel_tol=float(stationarity_phase_gradient_rel),
            phi_gradient_rel_tol=float(stationarity_phi_gradient_rel),
            phase_gradient_abs_tol_m_inv=float(stationarity_q_abs_m_inv),
            phi_gradient_abs_tol_V_m=float(stationarity_phi_gradient_abs_V_m),
            edge_active_threshold=float(stationarity_edge_active_threshold),
            bulk_exclusion_xi=float(stationarity_bulk_exclusion_xi),
        )
        dynamic_stationarity = dynamic_stationarity_diagnostics(
            history=online_history,
            nodes_m=nodes,
            delta0_J=float(material.delta0_J),
            minimum_tail_duration_ps=float(dynamic_stationarity_minimum_tail_ps),
            profile_relative_tolerance=float(dynamic_stationarity_profile_rel),
            voltage_relative_tolerance=float(dynamic_stationarity_voltage_rel),
            psl_threshold_over_delta0=float(dynamic_stationarity_psl_threshold),
            bulk_exclusion_xi=float(stationarity_bulk_exclusion_xi),
        )
        contact = contact_recovery_diagnostics(
            psi_dimensionless=np.asarray(frame.psi, dtype=np.complex128),
            nodes_m=nodes,
            material=material,
            Te_K=Te_K,
            threshold_fraction=float(terminal_healing_fraction),
            min_allowed_xi=float(recovery_min_xi),
            max_allowed_xi=float(recovery_max_xi),
            bin_width_m=float(getattr(mesh, "target_spacing_m", 0.0) or 0.0),
        )
        continuity = continuity_diagnostics(
            currents=currents,
            mesh=mesh,
            material=material,
            target_current_A=current_A,
            history=online_history,
            rms_tol=float(continuity_rms_tol),
            max_tol=float(continuity_max_tol),
            poisson_tol=float(continuity_poisson_tol),
        )
        thermal = thermal_stationarity_diagnostics(
            online_history,
            enabled=bool(thermal_config.enabled),
            start_time_s=float(thermal_config.start_time_s),
            bath_K=float(thermal_config.bath_K),
            minimum_tail_duration_ps=float(dynamic_stationarity_minimum_tail_ps),
            relative_tolerance=float(thermal_stationarity_relative_tolerance),
            p99_tolerance_K=float(thermal_stationarity_p99_tolerance_K),
            projection_horizon_ps=float(thermal_stationarity_projection_horizon_ps),
            projection_relative_tolerance=float(
                thermal_stationarity_projection_relative_tolerance
            ),
        )
        circuit = (
            circuit_stationarity_diagnostics(raw, config=circuit_runtime)
            if circuit_controller is not None
            else {"enabled": False, "passes": True}
        )
        phase_convergence = np.asarray(
            raw.get("allmaras_phase_convergence_converged", []),
            dtype=bool,
        ).reshape(-1)
        phase_converged = bool(
            phase_convergence.size and phase_convergence[-1]
        )
        readiness = classify_photon_readiness(
            strict_stationarity_passes=stationarity.passes,
            dynamic_stationarity_passes=dynamic_stationarity.passes,
            contact_recovery_passes=contact.passes,
            continuity_passes=continuity.passes,
            thermal_stationarity_passes=bool(thermal.get("passes", False)),
            circuit_stationarity_passes=bool(circuit.get("passes", False)),
            phase_drive_converged=phase_converged,
        )
        mesoscopic_passes = bool(readiness["mesoscopic_passes"])
        instantaneous_passes = bool(readiness["passes"])
        consecutive = (
            int(state.get("consecutive_passes", 0)) + 1
            if instantaneous_passes
            else 0
        )
        required = int(state["required_consecutive_passes"])
        passes = bool(instantaneous_passes and consecutive >= required)
        state.update(
            {
                "passes": passes,
                "instantaneous_passes": instantaneous_passes,
                "consecutive_passes": consecutive,
                "phase_drive_converged": phase_converged,
                "mesoscopic_passes": mesoscopic_passes,
                "mesoscopic_mode": str(readiness["mesoscopic_mode"]),
                "stationarity": stationarity.as_dict(),
                "dynamic_stationarity": dynamic_stationarity.as_dict(),
                "contact_recovery": contact.as_dict(),
                "continuity": continuity.as_dict(),
                "thermal_stationarity": thermal,
                "circuit_stationarity": circuit,
            }
        )
        return "total_stationarity_stop" if passes else None

    return callback


__all__ = ["build_total_stationarity_callback", "classify_photon_readiness"]
