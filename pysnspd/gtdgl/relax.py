"""Stationary gTDGL/Poisson relaxation driver for OE7.

The low-level current, Poisson, boundary-condition, update, diagnostic and I/O
helpers live in smaller sibling modules. This module keeps the public API used
by ``pipelines/02_ss_run_template.py`` and older tests, while the main driver
remains here.
"""
from __future__ import annotations

from typing import Any

import numpy as np

try:
    from tqdm.auto import trange
except Exception:  # pragma: no cover
    trange = None

from pysnspd.gtdgl.material import GTDGLMaterial
from pysnspd.gtdgl.operators import FVOperators, boundary_currents_from_node_vectors, terminal_voltage
from pysnspd.gtdgl.state import (
    MEV_J,
    CurrentFields,
    FormulaFields,
    PoissonResult,
    StepInfo,
    GTDGLStationaryState,
    RelaxationResult,
    PhiBoundaryConditions,
    _PoissonOperator,
)
from pysnspd.gtdgl.fields import (
    compute_current_fields,
    compute_formula_fields,
    edge_supercurrent_usadel,
    edge_supercurrent_gl,
    safe_abs_delta,
    pairbreaking_ratio_edges,
    edge_to_node_weighted_average,
)
from pysnspd.gtdgl.poisson_projection import (
    VALID_POISSON_TERMINAL_POLICIES,
    build_poisson_operator,
    _build_constrained_poisson_operator,
    solve_varphi_poisson,
    solve_poisson_potential,
    target_terminal_boundary_accum_A_m,
    build_phi_boundary_conditions,
    _edge_indices_and_signs_for_pairs,
)
from pysnspd.gtdgl.kwt_update import (
    VALID_PHI_PHASE_POLICIES,
    kwt_delta_update_attempt,
    kwt_local_update,
)
from pysnspd.gtdgl.stationary_boundary import (
    VALID_DELTA_BOUNDARY_POLICIES,
    apply_delta_boundary_policy,
    apply_stationary_boundary_conditions,
    clip_gap_amplitude,
    boundary_node_masks,
    terminal_inner_node_pairs,
    nearest_inward_boundary_pairs,
    _boundary_temperature_nodes,
    _edge_aware_boundary_pairs,
    _node_adjacency_from_ops_or_mesh,
)
from pysnspd.gtdgl.diagnostics import (
    current_residual,
    max_current_residual,
    normal_current_fraction_rms,
    current_density_maxima_A_m2,
    normal_current_fraction_max,
    seed_target_current_A,
    seed_q_bias_m_inv,
    seed_delta_bias_J,
    target_current_density_A_m2,
    suggest_next_dt,
    stationary_trial_rejection_reason,
)
from pysnspd.gtdgl.state_io import save_stationary_state_npz, save_relaxation_history_npz


def relax_stationary_gtdgl(
    *,
    mesh,
    edge_data,
    seed,
    material: GTDGLMaterial,
    ops: FVOperators,
    steps: int = 2000,
    dt_s: float = 1.0e-17,
    min_steps: int = 10,
    tolerance_eta: float = 1.0e-9,
    tolerance_current_residual: float = 1.0e-6,
    eta_reject: float = 5.0e-4,
    adapt_dt: bool = True,
    dt_min_s: float = 1.0e-22,
    dt_max_s: float = 1.0e-13,
    max_pairbreaking_accept: float = 0.95,
    max_js_over_javg_accept: float = 2.0,
    max_jtot_over_javg_accept: float = 2.5,
    spike_shrink_factor: float = 0.35,
    max_trial_retries: int = 40,
    lock_terminals: bool = True,
    delta_boundary_policy: str = "current_inversion",
    poisson_terminal_policy: str = "target_flux",
    target_current_A: float | None = None,
    progress: bool = False,
    n_phi_snapshots: int = 6,
    phi_phase_policy: str = "plus",
    use_phi_phase: bool = True,
) -> RelaxationResult:
    """Relax the OE6 seed with frozen temperatures and notebook solver ordering."""
    if not use_phi_phase and str(phi_phase_policy) == "plus":
        # Backward-compatible legacy knob. New diagnostics should prefer the
        # explicit ``phi_phase_policy`` argument.
        phi_phase_policy = "none"

    if phi_phase_policy not in VALID_PHI_PHASE_POLICIES:
        raise ValueError(
            "phi_phase_policy must be one of "
            f"{sorted(VALID_PHI_PHASE_POLICIES)}, got {phi_phase_policy!r}."
        )

    if delta_boundary_policy not in VALID_DELTA_BOUNDARY_POLICIES:
        raise ValueError(
            "delta_boundary_policy must be one of "
            f"{sorted(VALID_DELTA_BOUNDARY_POLICIES)}, got {delta_boundary_policy!r}."
        )
    if poisson_terminal_policy not in VALID_POISSON_TERMINAL_POLICIES:
        raise ValueError(
            "poisson_terminal_policy must be one of "
            f"{sorted(VALID_POISSON_TERMINAL_POLICIES)}, got {poisson_terminal_policy!r}."
        )
    if not lock_terminals:
        delta_boundary_policy = "none"


    if steps <= 0:
        raise ValueError("steps must be positive.")
    if dt_s <= 0.0:
        raise ValueError("dt_s must be positive.")
    if min_steps < 0:
        raise ValueError("min_steps must be non-negative.")

    if target_current_A is None:
        target_current_A = seed_target_current_A(seed)
    target_current_A = float(target_current_A)
    q_bias = seed_q_bias_m_inv(seed, target_current_A=target_current_A)
    javg = target_current_density_A_m2(material, target_current_A)
    q_ref = abs(q_bias) if abs(q_bias) > 0.0 else 1.0

    psi0 = (
        np.asarray(seed.node_psi_real_J, dtype=float)
        + 1j * np.asarray(seed.node_psi_imag_J, dtype=float)
    )
    phi0 = np.asarray(seed.node_phi_electric_V, dtype=float).copy()
    Te = np.asarray(seed.node_Te_K, dtype=float).copy()
    Tph = np.asarray(seed.node_Tph_K, dtype=float).copy()

    psi = apply_delta_boundary_policy(
        psi_trial_J=psi0,
        mesh=mesh,
        seed=seed,
        q_bias_m_inv=q_bias,
        material=material,
        ops=ops,
        Te_K=Te,
        target_current_A=target_current_A,
        policy=delta_boundary_policy,
    )
    phi = phi0 - float(np.mean(phi0))

    if poisson_terminal_policy == "target_flux":
        boundary_accum = target_terminal_boundary_accum_A_m(
            edge_data=edge_data,
            ops=ops,
            material=material,
            target_current_A=target_current_A,
        )
    else:
        boundary_accum = np.zeros(ops.n_nodes, dtype=float)
    # pyTDGL Eq. (17) Poisson projection: no first-edge electric constraints.
    # The terminal/outward flux enters only through ``boundary_accum`` in the
    # RHS, and the mean-zero Neumann matrix is factorized once by sparse LU.
    phi_bc = None
    poisson_op = build_poisson_operator(material=material, ops=ops, phi_bc=None)

    # Notebook initial projection: compute defs, then Poisson, then recompute fields.
    defs0 = compute_formula_fields(psi_J=psi, Te_K=Te, material=material, ops=ops)
    poisson0 = solve_varphi_poisson(
        edge_js_us_A_m2=defs0.edge_js_us_A_m2,
        material=material,
        ops=ops,
        poisson_op=poisson_op,
        boundary_accum_A_m=boundary_accum,
        phi_bc=phi_bc,
    )
    phi = poisson0.phi_V
    currents = compute_current_fields(
        psi_J=psi,
        phi_V=phi,
        Te_K=Te,
        material=material,
        ops=ops,
        boundary_accum_A_m=boundary_accum,
    )

    t_s = 0.0
    accepted = 0
    rejected = 0
    spike_rejected = 0
    last_rejection_reason = ""
    converged = False

    hist_keys = [
        "t_s",
        "dt_s",
        "retries",
        "discr_min",
        "eta_R",
        "max_amp2_change_rel",
        "current_residual",
        "current_residual_max",
        "terminal_voltage_V",
        "pairbreaking_max",
        "delta_min_over_delta0",
        "delta_max_over_delta0",
        "normal_current_fraction_rms",
        "normal_current_fraction_max",
        "normal_current_max_A_m2",
        "total_current_max_A_m2",
        "median_Q_m_inv",
        "p95_Q_m_inv",
        "max_Q_m_inv",
        "max_js_A_m2",
        "max_j_A_m2",
        "edge_pairbreaking_max",
        "edge_js_over_javg",
        "edge_jtot_over_javg",
    ]
    hist: dict[str, list[float]] = {key: [] for key in hist_keys}

    n_phi_snapshots = max(2, int(n_phi_snapshots))
    snapshot_steps = set(np.linspace(0, int(steps) - 1, n_phi_snapshots, dtype=int).tolist())
    snapshots: dict[str, list[np.ndarray] | list[float]] = {
        "snapshot_t_s": [],
        "psi_snapshot_real_J": [],
        "psi_snapshot_imag_J": [],
        "delta_snapshot_meV": [],
        "phi_snapshot_V": [],
        "current_density_snapshot_A_m2": [],
        "current_density_snapshot_x_A_m2": [],
        "current_density_snapshot_y_A_m2": [],
        "supercurrent_density_snapshot_A_m2": [],
        "supercurrent_density_snapshot_x_A_m2": [],
        "supercurrent_density_snapshot_y_A_m2": [],
        "normal_current_density_snapshot_A_m2": [],
        "normal_current_density_snapshot_x_A_m2": [],
        "normal_current_density_snapshot_y_A_m2": [],
        "divergence_snapshot_A_m3": [],
        "pairbreaking_ratio_snapshot": [],
        "edge_Q_snapshot_m_inv": [],
        "edge_js_us_snapshot_A_m2": [],
        "edge_jn_snapshot_A_m2": [],
        "edge_jtot_snapshot_A_m2": [],
    }

    def append_snapshot() -> None:
        jtot_mag = np.sqrt(currents.node_jtot_x_A_m2**2 + currents.node_jtot_y_A_m2**2)
        js_mag = np.sqrt(currents.node_js_us_x_A_m2**2 + currents.node_js_us_y_A_m2**2)
        jn_mag = np.sqrt(currents.node_jn_x_A_m2**2 + currents.node_jn_y_A_m2**2)
        snapshots["snapshot_t_s"].append(float(t_s))
        snapshots["psi_snapshot_real_J"].append(np.real(psi).copy())
        snapshots["psi_snapshot_imag_J"].append(np.imag(psi).copy())
        snapshots["delta_snapshot_meV"].append(np.abs(psi).copy() / MEV_J)
        snapshots["phi_snapshot_V"].append(phi.copy())
        snapshots["current_density_snapshot_A_m2"].append(jtot_mag.copy())
        snapshots["current_density_snapshot_x_A_m2"].append(currents.node_jtot_x_A_m2.copy())
        snapshots["current_density_snapshot_y_A_m2"].append(currents.node_jtot_y_A_m2.copy())
        snapshots["supercurrent_density_snapshot_A_m2"].append(js_mag.copy())
        snapshots["supercurrent_density_snapshot_x_A_m2"].append(currents.node_js_us_x_A_m2.copy())
        snapshots["supercurrent_density_snapshot_y_A_m2"].append(currents.node_js_us_y_A_m2.copy())
        snapshots["normal_current_density_snapshot_A_m2"].append(jn_mag.copy())
        snapshots["normal_current_density_snapshot_x_A_m2"].append(currents.node_jn_x_A_m2.copy())
        snapshots["normal_current_density_snapshot_y_A_m2"].append(currents.node_jn_y_A_m2.copy())
        snapshots["divergence_snapshot_A_m3"].append(currents.node_div_jtot_A_m3.copy())
        snapshots["pairbreaking_ratio_snapshot"].append(currents.node_pairbreaking_ratio.copy())
        snapshots["edge_Q_snapshot_m_inv"].append(currents.edge_Q_m_inv.copy())
        snapshots["edge_js_us_snapshot_A_m2"].append(currents.edge_js_us_A_m2.copy())
        snapshots["edge_jn_snapshot_A_m2"].append(currents.edge_jn_A_m2.copy())
        snapshots["edge_jtot_snapshot_A_m2"].append(currents.edge_jtot_A_m2.copy())

    append_snapshot()

    iterator = range(int(steps))
    if progress and trange is not None:
        iterator = trange(int(steps), desc="OE7 notebook KWT", leave=True)

    for n in iterator:
        retries = 0
        dt_eff = float(dt_s)
        while True:
            defs_n = compute_formula_fields(psi_J=psi, Te_K=Te, material=material, ops=ops)
            psi_new, discr_min = kwt_delta_update_attempt(
                psi_J=psi,
                phi_V=phi,
                defs=defs_n,
                dt_s=dt_eff,
                material=material,
                phi_phase_policy=phi_phase_policy,
            )
            if psi_new is not None:
                break
            retries += 1
            rejected += 1
            if retries > 30:
                raise RuntimeError(
                    "Failed KWT update: negative discriminant after "
                    f"{retries} retries. Last min={discr_min:.3e}"
                )
            dt_eff = max(dt_min_s, 0.5 * dt_eff)

        psi_trial = apply_delta_boundary_policy(
            psi_trial_J=psi_new,
            mesh=mesh,
            seed=seed,
            q_bias_m_inv=q_bias,
            material=material,
            ops=ops,
            Te_K=Te,
            target_current_A=target_current_A,
            policy=delta_boundary_policy,
        )
        defs_pre = compute_formula_fields(psi_J=psi_trial, Te_K=Te, material=material, ops=ops)
        poisson = solve_varphi_poisson(
            edge_js_us_A_m2=defs_pre.edge_js_us_A_m2,
            material=material,
            ops=ops,
            poisson_op=poisson_op,
            boundary_accum_A_m=boundary_accum,
            phi_bc=phi_bc,
        )
        phi_trial = poisson.phi_V
        trial_currents = compute_current_fields(
            psi_J=psi_trial,
            phi_V=phi_trial,
            Te_K=Te,
            material=material,
            ops=ops,
            boundary_accum_A_m=boundary_accum,
        )

        amp2_change_rel = float(
            np.nanmax(np.abs(np.abs(psi_trial) ** 2 - np.abs(psi) ** 2))
            / material.delta0_J**2
        )

        trial_edge_pb_max = float(np.nanmax(trial_currents.edge_pairbreaking_ratio))
        trial_js_max_A_m2 = float(np.nanmax(np.abs(trial_currents.edge_js_us_A_m2)))
        trial_jtot_max_A_m2 = float(np.nanmax(np.abs(trial_currents.edge_jtot_A_m2)))
        javg_abs = max(abs(float(javg)), 1.0e-300)
        trial_js_over_javg = trial_js_max_A_m2 / javg_abs
        trial_jtot_over_javg = trial_jtot_max_A_m2 / javg_abs

        rejection_reason = stationary_trial_rejection_reason(
            amp2_change_rel=amp2_change_rel,
            edge_pairbreaking_max=trial_edge_pb_max,
            edge_js_over_javg=trial_js_over_javg,
            edge_jtot_over_javg=trial_jtot_over_javg,
            eta_reject=float(eta_reject),
            max_pairbreaking_accept=float(max_pairbreaking_accept),
            max_js_over_javg_accept=float(max_js_over_javg_accept),
            max_jtot_over_javg_accept=float(max_jtot_over_javg_accept),
        )
        if adapt_dt and rejection_reason is not None and dt_eff > 1.01 * float(dt_min_s):
            retries += 1
            rejected += 1
            spike_rejected += 1
            last_rejection_reason = str(rejection_reason)
            if retries > int(max_trial_retries):
                raise RuntimeError(
                    "Failed OE7 stationary step after adaptive trial rejections. "
                    f"last_reason={last_rejection_reason}, "
                    f"dt_eff_s={dt_eff:.3e}, "
                    f"amp2_change_rel={amp2_change_rel:.3e}, "
                    f"edge_pairbreaking_max={trial_edge_pb_max:.3e}, "
                    f"edge_js_over_javg={trial_js_over_javg:.3e}, "
                    f"edge_jtot_over_javg={trial_jtot_over_javg:.3e}"
                )
            dt_eff = max(float(dt_min_s), float(spike_shrink_factor) * dt_eff)
            continue

        psi = psi_trial
        phi = phi_trial
        currents = trial_currents
        t_s += dt_eff
        accepted += 1

        residual = current_residual(currents, mesh, material, target_current_A)
        residual_max = max_current_residual(currents, mesh, material, target_current_A)
        voltage = terminal_voltage(np.asarray(mesh.nodes, dtype=float), phi, length_m=float(mesh.length_m))
        pb_max = float(np.nanmax(currents.node_pairbreaking_ratio))
        edge_pb_max = float(np.nanmax(currents.edge_pairbreaking_ratio))
        delta_min_ratio = float(np.nanmin(np.abs(psi)) / material.delta0_J)
        delta_max_ratio = float(np.nanmax(np.abs(psi)) / material.delta0_J)
        normal_frac = normal_current_fraction_rms(currents)
        normal_max_frac = normal_current_fraction_max(currents)
        jn_max_A_m2, jt_max_A_m2 = current_density_maxima_A_m2(currents)
        Qabs = np.abs(currents.edge_Q_m_inv)
        js_max_A_m2 = float(np.nanmax(np.abs(currents.edge_js_us_A_m2)))
        jtot_max_A_m2 = float(np.nanmax(np.abs(currents.edge_jtot_A_m2)))
        js_over_javg = js_max_A_m2 / max(abs(float(javg)), 1.0e-300)
        jtot_over_javg = jtot_max_A_m2 / max(abs(float(javg)), 1.0e-300)

        values = {
            "t_s": t_s,
            "dt_s": dt_eff,
            "retries": float(retries),
            "discr_min": float(discr_min),
            "eta_R": amp2_change_rel,
            "max_amp2_change_rel": amp2_change_rel,
            "current_residual": residual,
            "current_residual_max": residual_max,
            "terminal_voltage_V": voltage,
            "pairbreaking_max": pb_max,
            "delta_min_over_delta0": delta_min_ratio,
            "delta_max_over_delta0": delta_max_ratio,
            "normal_current_fraction_rms": normal_frac,
            "normal_current_fraction_max": normal_max_frac,
            "normal_current_max_A_m2": jn_max_A_m2,
            "total_current_max_A_m2": jt_max_A_m2,
            "median_Q_m_inv": float(np.nanmedian(Qabs)),
            "p95_Q_m_inv": float(np.nanpercentile(Qabs, 95.0)),
            "max_Q_m_inv": float(np.nanmax(Qabs)),
            "max_js_A_m2": js_max_A_m2,
            "max_j_A_m2": jtot_max_A_m2,
            "edge_pairbreaking_max": edge_pb_max,
            "edge_js_over_javg": js_over_javg,
            "edge_jtot_over_javg": jtot_over_javg,
        }
        for key in hist_keys:
            hist[key].append(values[key])

        if n in snapshot_steps:
            append_snapshot()

        if progress and hasattr(iterator, "set_postfix") and accepted % 10 == 0:
            iterator.set_postfix(
                dA2=f"{amp2_change_rel:.2e}",
                eps=f"{residual:.2e}",
                V=f"{voltage:.2e}",
                chi=f"{pb_max:.3g}",
                dt_fs=f"{dt_eff / 1.0e-15:.3g}",
            )

        if accepted >= min_steps and amp2_change_rel < tolerance_eta and residual < tolerance_current_residual:
            converged = True
            break

        dt_s = suggest_next_dt(
            dt_s=dt_eff,
            max_amp2_change_rel=amp2_change_rel,
            retries=retries,
            adaptive=adapt_dt,
            target=float(eta_reject),
            shrink_factor=0.55,
            grow_factor=1.03,
            dt_min_s=dt_min_s,
            dt_max_s=dt_max_s,
        )

        if not np.all(np.isfinite(psi)) or not np.all(np.isfinite(phi)):
            raise FloatingPointError(f"Stopped: non-finite state at accepted step {accepted}.")

    if len(snapshots["snapshot_t_s"]) == 0 or snapshots["snapshot_t_s"][-1] != t_s:
        append_snapshot()

    # Keep exactly n_phi_snapshots, preserving first and final snapshots.
    n_snap = len(snapshots["snapshot_t_s"])
    if n_snap > n_phi_snapshots:
        keep = np.unique(np.rint(np.linspace(0, n_snap - 1, n_phi_snapshots)).astype(int))
        if keep[-1] != n_snap - 1:
            keep[-1] = n_snap - 1
        for key, seq in list(snapshots.items()):
            snapshots[key] = [seq[int(i)] for i in keep]

    boundary = boundary_currents_from_node_vectors(
        mesh=mesh,
        edge_data=edge_data,
        jx_A_m2=currents.node_jtot_x_A_m2,
        jy_A_m2=currents.node_jtot_y_A_m2,
        thickness_m=material.thickness_m,
    )
    voltage = terminal_voltage(np.asarray(mesh.nodes, dtype=float), phi, length_m=float(mesh.length_m))
    normal_ohmic_voltage = (
        float(target_current_A)
        * float(mesh.length_m)
        / (material.sigma_n_S_m * material.width_m * material.thickness_m)
    )
    normal_max_A_m2, total_max_A_m2 = current_density_maxima_A_m2(currents)

    summary = {
        "backend": "oe7_notebook_order_kwt_poisson_v1",
        "gauge_policy": "selectable_temporal_gauge_link_in_kwt",
        "phi_phase_policy": str(phi_phase_policy),
        "converged": bool(converged),
        "accepted_steps": int(accepted),
        "rejected_steps": int(rejected),
        "spike_rejected_steps": int(spike_rejected),
        "last_rejection_reason": str(last_rejection_reason),
        "stationary_acceptance_policy": {
            "eta_reject": float(eta_reject),
            "max_pairbreaking_accept": float(max_pairbreaking_accept),
            "max_js_over_javg_accept": float(max_js_over_javg_accept),
            "max_jtot_over_javg_accept": float(max_jtot_over_javg_accept),
            "spike_shrink_factor": float(spike_shrink_factor),
            "max_trial_retries": int(max_trial_retries),
        },
        "final_time_ps": float(t_s / 1.0e-12),
        "tau_scale": float(material.tau_scale),
        "tau_ee_Tc_effective_ps": float(material.tau_scale * material.tau_ee_Tc_s / 1.0e-12),
        "tau_ep_Tc_effective_ps": float(material.tau_scale * material.tau_ep_Tc_s / 1.0e-12),
        "target_current_A": float(target_current_A),
        "delta_boundary_policy": str(delta_boundary_policy),
        "poisson_terminal_policy": str(poisson_terminal_policy),
        "target_q_bias_m_inv": float(q_bias),
        "target_j_bias_A_m2": float(javg),
        "terminal_voltage_V": float(voltage),
        "normal_ohmic_voltage_V": float(normal_ohmic_voltage),
        "terminal_voltage_over_normal": float(
            voltage / normal_ohmic_voltage if normal_ohmic_voltage != 0.0 else float("nan")
        ),
        "normal_current_fraction_rms": float(normal_current_fraction_rms(currents)),
        "normal_current_fraction_max": float(normal_max_A_m2 / max(total_max_A_m2, 1.0e-300)),
        "normal_current_max_A_m2": float(normal_max_A_m2),
        "total_current_max_A_m2": float(total_max_A_m2),
        "current_residual": float(current_residual(currents, mesh, material, target_current_A)),
        "eta_R_final": float(hist["eta_R"][-1]) if hist["eta_R"] else float("nan"),
        "divergence_rms_A_m3": float(np.sqrt(np.nanmean(currents.node_div_jtot_A_m3**2))),
        "min_delta_over_delta0": float(np.nanmin(np.abs(psi)) / material.delta0_J),
        "mean_delta_over_delta0": float(np.nanmean(np.abs(psi)) / material.delta0_J),
        "max_pairbreaking_ratio": float(np.nanmax(currents.node_pairbreaking_ratio)),
        "edge_pairbreaking_ratio_max": float(np.nanmax(currents.edge_pairbreaking_ratio)),
        "p99_pairbreaking_ratio": float(np.nanpercentile(currents.node_pairbreaking_ratio, 99.0)),
        "edge_Q_max_m_inv": float(np.nanmax(np.abs(currents.edge_Q_m_inv))),
        "edge_js_over_javg_max": float(np.nanmax(np.abs(currents.edge_js_us_A_m2)) / max(abs(float(javg)), 1.0e-300)),
        "edge_jtot_over_javg_max": float(np.nanmax(np.abs(currents.edge_jtot_A_m2)) / max(abs(float(javg)), 1.0e-300)),
        "boundary_currents_A": boundary,
    }

    metadata = {
        "backend": summary["backend"],
        "description": "Notebook-order frozen-temperature gTDGL/Poisson relaxation.",
        "thermal_policy": "frozen_Te_Tph",
        "circuit_policy": "inactive",
        "boundary_policy": str(delta_boundary_policy),
        "poisson_terminal_policy": str(poisson_terminal_policy),
        "phi_phase_policy": str(phi_phase_policy),
        "poisson_policy": "notebook_conservative_FV_mean_zero_gauge",
        "pairbreaking_ratio": "xi^2 Q^2 / (1 - T/Tc)",
        "adaptive_rejection_policy": "reject accepted-looking SS trials if eta, pairbreaking, js/javg, or jtot/javg exceed stationary ceilings",
    }

    state = GTDGLStationaryState(
        psi_J=psi,
        phi_V=phi,
        Te_K=Te,
        Tph_K=Tph,
        currents=currents,
        metadata=metadata,
    )

    history: dict[str, np.ndarray] = {key: np.asarray(val, dtype=float) for key, val in hist.items()}
    history["spike_rejected_steps"] = np.asarray([spike_rejected], dtype=float)
    snapshot_t_s = np.asarray(snapshots["snapshot_t_s"], dtype=float)
    history.update(
        {
            "delta0_meV": np.asarray([material.delta0_J / MEV_J], dtype=float),
            "javg_A_m2": np.asarray([javg], dtype=float),
            "qref_m_inv": np.asarray([q_ref], dtype=float),
            "snapshot_t_s": snapshot_t_s,
            "phi_snapshot_t_s": snapshot_t_s,
            "phi_snapshot_V": np.asarray(snapshots["phi_snapshot_V"], dtype=float),
            "psi_snapshot_t_s": snapshot_t_s,
            "psi_snapshot_real_J": np.asarray(snapshots["psi_snapshot_real_J"], dtype=float),
            "psi_snapshot_imag_J": np.asarray(snapshots["psi_snapshot_imag_J"], dtype=float),
            "delta_snapshot_t_s": snapshot_t_s,
            "delta_snapshot_meV": np.asarray(snapshots["delta_snapshot_meV"], dtype=float),
            "current_snapshot_t_s": snapshot_t_s,
            "current_density_snapshot_A_m2": np.asarray(snapshots["current_density_snapshot_A_m2"], dtype=float),
            "current_density_snapshot_x_A_m2": np.asarray(snapshots["current_density_snapshot_x_A_m2"], dtype=float),
            "current_density_snapshot_y_A_m2": np.asarray(snapshots["current_density_snapshot_y_A_m2"], dtype=float),
            "jtot_snapshot_t_s": snapshot_t_s,
            "jtot_snapshot_mag_A_m2": np.asarray(snapshots["current_density_snapshot_A_m2"], dtype=float),
            "jtot_snapshot_x_A_m2": np.asarray(snapshots["current_density_snapshot_x_A_m2"], dtype=float),
            "jtot_snapshot_y_A_m2": np.asarray(snapshots["current_density_snapshot_y_A_m2"], dtype=float),
            "supercurrent_snapshot_t_s": snapshot_t_s,
            "supercurrent_density_snapshot_A_m2": np.asarray(snapshots["supercurrent_density_snapshot_A_m2"], dtype=float),
            "supercurrent_density_snapshot_x_A_m2": np.asarray(snapshots["supercurrent_density_snapshot_x_A_m2"], dtype=float),
            "supercurrent_density_snapshot_y_A_m2": np.asarray(snapshots["supercurrent_density_snapshot_y_A_m2"], dtype=float),
            "js_us_snapshot_t_s": snapshot_t_s,
            "js_us_snapshot_mag_A_m2": np.asarray(snapshots["supercurrent_density_snapshot_A_m2"], dtype=float),
            "js_us_snapshot_x_A_m2": np.asarray(snapshots["supercurrent_density_snapshot_x_A_m2"], dtype=float),
            "js_us_snapshot_y_A_m2": np.asarray(snapshots["supercurrent_density_snapshot_y_A_m2"], dtype=float),
            "normal_current_snapshot_t_s": snapshot_t_s,
            "normal_current_density_snapshot_A_m2": np.asarray(snapshots["normal_current_density_snapshot_A_m2"], dtype=float),
            "normal_current_density_snapshot_x_A_m2": np.asarray(snapshots["normal_current_density_snapshot_x_A_m2"], dtype=float),
            "normal_current_density_snapshot_y_A_m2": np.asarray(snapshots["normal_current_density_snapshot_y_A_m2"], dtype=float),
            "jn_snapshot_t_s": snapshot_t_s,
            "jn_snapshot_mag_A_m2": np.asarray(snapshots["normal_current_density_snapshot_A_m2"], dtype=float),
            "jn_snapshot_x_A_m2": np.asarray(snapshots["normal_current_density_snapshot_x_A_m2"], dtype=float),
            "jn_snapshot_y_A_m2": np.asarray(snapshots["normal_current_density_snapshot_y_A_m2"], dtype=float),
            "divergence_snapshot_t_s": snapshot_t_s,
            "divergence_snapshot_A_m3": np.asarray(snapshots["divergence_snapshot_A_m3"], dtype=float),
            "pairbreaking_snapshot_t_s": snapshot_t_s,
            "pairbreaking_ratio_snapshot": np.asarray(snapshots["pairbreaking_ratio_snapshot"], dtype=float),

            # Exact edge topology used by the FV solver. These static arrays let
            # plotting diagnostics inspect node-to-node edge currents without
            # rebuilding/reordering edges from the triangulation.
            "edge_i": np.asarray(ops.edge_i, dtype=np.int64),
            "edge_j": np.asarray(ops.edge_j, dtype=np.int64),
            "edge_length_m": np.asarray(ops.edge_length_m, dtype=float),
            "edge_unit_x": np.asarray(ops.edge_unit[:, 0], dtype=float),
            "edge_unit_y": np.asarray(ops.edge_unit[:, 1], dtype=float),
            "dual_face_length_m": np.asarray(ops.dual_face_length_m, dtype=float),

            # Edge-current snapshots. These are the literal edge projections
            # that feed the node-vector reconstructions and the Poisson balance.
            "edge_snapshot_t_s": snapshot_t_s,
            "edge_Q_snapshot_m_inv": np.asarray(snapshots["edge_Q_snapshot_m_inv"], dtype=float),
            "edge_js_us_snapshot_A_m2": np.asarray(snapshots["edge_js_us_snapshot_A_m2"], dtype=float),
            "edge_jn_snapshot_A_m2": np.asarray(snapshots["edge_jn_snapshot_A_m2"], dtype=float),
            "edge_jtot_snapshot_A_m2": np.asarray(snapshots["edge_jtot_snapshot_A_m2"], dtype=float),
        }
    )

    return RelaxationResult(state=state, history=history, summary=summary)

