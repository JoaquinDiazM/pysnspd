"""
Basic plotting utilities for pySNSPD.

OE2 scope:
- Mesh geometry plot.
- Boundary tag diagnostic plot.
- Clear visual distinction between interior nodes, contacts and boundaries.

These functions save figures to disk and do not display interactive windows.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg", force=True)

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from pysnspd.mesh.delaunay import MeshData
from pysnspd.mesh.edges import EdgeData


_COLORS = {
    "interior_node": "#4b5563",
    "interior_edge": "#cbd5e1",
    "left": "#d62728",
    "right": "#9467bd",
    "top": "#1f77b4",
    "bottom": "#2ca02c",
    "boundary_unknown": "#111827",
}


def plot_mesh_geometry(
    mesh: MeshData,
    edge_data: EdgeData,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """
    Plot nodes, interior edges and boundary edges with distinct styling.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    nodes_nm = mesh.nodes * 1.0e9
    node_tags = _classify_nodes(mesh)

    fig, ax = plt.subplots(figsize=(9.0, 3.2))

    _draw_edge_group(
        ax,
        nodes_nm,
        edge_data,
        edge_data.tags == "interior",
        color=_COLORS["interior_edge"],
        linewidth=0.35,
        alpha=0.55,
        label=None,
    )

    for tag in ["top", "bottom", "left", "right", "boundary_unknown"]:
        mask = edge_data.tags == tag
        if np.any(mask):
            _draw_edge_group(
                ax,
                nodes_nm,
                edge_data,
                mask,
                color=_COLORS[tag],
                linewidth=0.45 if tag in {"top", "bottom"} else 0.5,
                alpha=0.95,
                #label=f"{tag} edges",
                label=None,
            )

    _draw_node_group(
        ax,
        nodes_nm,
        node_tags == "interior",
        color=_COLORS["interior_node"],
        size=1,
        alpha=0.65,
        label="interior nodes",
    )

    for tag in ["top", "bottom", "left", "right", "boundary_unknown"]:
        mask = node_tags == tag
        if np.any(mask):
            _draw_node_group(
                ax,
                nodes_nm,
                mask,
                color=_COLORS[tag],
                size=2 if tag in {"left", "right"} else 1.5,
                alpha=0.95,
                label=f"{tag} nodes",
            )

    ax.set_title("Nanowire mesh")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    _set_mesh_limits(ax, mesh)
    ax.grid(False)
    _legend_below(ax, ncol=3, fontsize=7, y_offset=-0.25)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output


def plot_boundary_tags(
    mesh: MeshData,
    edge_data: EdgeData,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """
    Plot only the boundary-tag diagnostic with a light interior background.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    nodes_nm = mesh.nodes * 1.0e9

    fig, ax = plt.subplots(figsize=(9.0, 3.2))

    _draw_edge_group(
        ax,
        nodes_nm,
        edge_data,
        edge_data.tags == "interior",
        color=_COLORS["interior_edge"],
        linewidth=0.25,
        alpha=0.30,
        label=None,
    )

    for tag in ["left", "right", "top", "bottom", "boundary_unknown"]:
        mask = edge_data.tags == tag
        if not np.any(mask):
            continue

        _draw_edge_group(
            ax,
            nodes_nm,
            edge_data,
            mask,
            color=_COLORS[tag],
            linewidth=0.45 if tag in {"left", "right"} else 0.5,
            alpha=0.98,
            label=f"{tag} ({int(np.count_nonzero(mask))})",
        )

    ax.set_title("Boundary tags")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    _set_mesh_limits(ax, mesh)
    ax.grid(False)
    _legend_below(ax, ncol=2, fontsize=8, y_offset=-0.25)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output


def _draw_edge_group(
    ax,
    nodes_nm: np.ndarray,
    edge_data: EdgeData,
    mask: np.ndarray,
    *,
    color: str,
    linewidth: float,
    alpha: float,
    label: str | None,
) -> None:
    """
    Draw a masked edge group.
    """
    if not np.any(mask):
        return

    segments = nodes_nm[edge_data.edges[mask]]
    collection = LineCollection(
        segments,
        colors=color,
        linewidths=linewidth,
        alpha=alpha,
        label=label,
    )
    ax.add_collection(collection)


def _draw_node_group(
    ax,
    nodes_nm: np.ndarray,
    mask: np.ndarray,
    *,
    color: str,
    size: float,
    alpha: float,
    label: str,
) -> None:
    """
    Draw a masked node group.
    """
    if not np.any(mask):
        return

    ax.scatter(
        nodes_nm[mask, 0],
        nodes_nm[mask, 1],
        s=size,
        c=color,
        alpha=alpha,
        linewidths=0.0,
        label=label,
        zorder=3,
    )


def _classify_nodes(mesh: MeshData) -> np.ndarray:
    """
    Classify nodes as interior, left, right, top, bottom or boundary_unknown.

    Corners are assigned to left/right contacts. This is intentional because
    contact treatment will be more important than top/bottom insulation at
    the longitudinal terminals.
    """
    nodes = mesh.nodes
    tags = np.full(mesh.n_nodes, "interior", dtype="<U32")

    atol = max(1.0e-15, mesh.target_spacing_m * 1.0e-6)

    left = np.isclose(nodes[:, 0], 0.0, atol=atol, rtol=0.0)
    right = np.isclose(nodes[:, 0], mesh.length_m, atol=atol, rtol=0.0)
    bottom = np.isclose(nodes[:, 1], -0.5 * mesh.width_m, atol=atol, rtol=0.0)
    top = np.isclose(nodes[:, 1], 0.5 * mesh.width_m, atol=atol, rtol=0.0)

    tags[top] = "top"
    tags[bottom] = "bottom"
    tags[left] = "left"
    tags[right] = "right"

    boundary = left | right | top | bottom
    unknown = boundary & (tags == "interior")
    tags[unknown] = "boundary_unknown"

    return tags


def _set_mesh_limits(ax, mesh: MeshData) -> None:
    """
    Set plot limits with a small padding.
    """
    length_nm = mesh.length_m * 1.0e9
    width_nm = mesh.width_m * 1.0e9

    ax.set_xlim(-0.04 * length_nm, 1.04 * length_nm)
    ax.set_ylim(-0.68 * width_nm, 0.68 * width_nm)

def _legend_below(ax, *, ncol: int, fontsize: int = 7, y_offset: float = -0.2) -> None:
    """
    Place the legend below the axes so it does not cover the mesh.
    """
    handles, labels = ax.get_legend_handles_labels()
    if not handles:
        return

    ax.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, y_offset),
        ncol=ncol,
        fontsize=fontsize,
        frameon=True,
    )


def plot_usadel_dos_slices(
    catalog,
    output_path: str | Path,
    *,
    dpi: int = 220,
) -> Path:
    """
    Plot representative DOS slices from a Usadel/Dynes catalogue.

    The plot shows the largest-Delta slice for several depairing values.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    energy_meV = catalog.energy_values_J / 1.602176634e-22
    rho = catalog.rho_delta_gamma_E

    delta_index = rho.shape[0] - 1
    n_gamma = rho.shape[1]

    if n_gamma <= 4:
        gamma_indices = list(range(n_gamma))
    else:
        gamma_indices = sorted(set([0, n_gamma // 3, 2 * n_gamma // 3, n_gamma - 1]))

    fig, ax = plt.subplots(figsize=(7.0, 4.2))

    for idx in gamma_indices:
        gamma_meV = catalog.gamma_values_J[idx] / 1.602176634e-22
        ax.plot(
            energy_meV,
            rho[delta_index, idx, :],
            linewidth=1.2,
            label=rf"$\Gamma_q={gamma_meV:.3f}$ meV",
        )

    delta_meV = catalog.delta_values_J[delta_index] / 1.602176634e-22
    ax.axvline(delta_meV, linestyle="--", linewidth=0.9, alpha=0.7)

    ax.set_title("DOS catalogue diagnostic")
    ax.set_xlabel("E [meV]")
    ax.set_ylabel(r"$\rho(E;|\Delta|,\Gamma_q)$")
    ax.grid(True, linewidth=0.25, alpha=0.35)
    ax.legend(loc="best", fontsize=8, frameon=True)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output

def plot_usadel_calibration_sweep(
    catalog,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """
    Plot the Usadel critical-current calibration sweep.

    The figure shows I(q) and Delta_eq(q). The calibrated critical point is
    marked explicitly.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    q = catalog.calibration_q_values_m_inv
    current_uA = catalog.calibration_current_values_A * 1.0e6
    delta_meV = catalog.calibration_delta_eq_values_J / 1.602176634e-22

    metadata = catalog.metadata
    cal = metadata["calibration"]

    q_c = float(cal["q_critical_m_inv"])
    Ic_uA = float(cal["Ic_model_A"]) * 1.0e6
    delta_c_meV = float(cal["delta_critical_meV"])
    D_m2_s = float(metadata["D_m2_s"])
    jc_A_m2 = float(cal["jc_target_A_m2"])

    # Color convention
    current_color = "tab:blue"
    current_light = "cornflowerblue"
    gap_color = "tab:red"
    gap_light = "lightcoral"
    qcrit_color = "0.35"

    fig, ax1 = plt.subplots(figsize=(7.2, 4.4))

    # Current branch
    ax1.plot(
        q * 1.0e-6,
        current_uA,
        color=current_color,
        linewidth=1.5,
        label=r"$I_s(q)$",
    )
    ax1.axvline(
        q_c * 1.0e-6,
        color=qcrit_color,
        linestyle="--",
        linewidth=0.9,
        alpha=0.8,
        label=r"$q_c$",
    )
    ax1.axhline(
        Ic_uA,
        color=current_light,
        linestyle=":",
        linewidth=0.9,
        alpha=0.9,
    )
    ax1.scatter(
        [q_c * 1.0e-6],
        [Ic_uA],
        color=current_color,
        edgecolor="white",
        linewidth=0.6,
        s=28,
        zorder=4,
        label=rf"$I_c={Ic_uA:.3f}\,\mu\mathrm{{A}}$",
    )

    ax1.set_xlabel(r"$q$ [$10^6\,\mathrm{m}^{-1}$]")
    ax1.set_ylabel(r"$I_s$ [$\mu$A]", color=current_color)
    ax1.tick_params(axis="y", colors=current_color)
    ax1.spines["left"].set_color(current_color)
    ax1.grid(True, linewidth=0.25, alpha=0.35)

    # Gap branch
    ax2 = ax1.twinx()
    ax2.plot(
        q * 1.0e-6,
        delta_meV,
        color=gap_color,
        linewidth=1.2,
        linestyle="--",
        label=r"$\Delta_{\rm eq}(q)$",
    )
    ax2.scatter(
        [q_c * 1.0e-6],
        [delta_c_meV],
        color=gap_color,
        edgecolor="white",
        linewidth=0.6,
        s=28,
        zorder=4,
        label=rf"$\Delta(q_c)={delta_c_meV:.3f}$ meV",
    )
    ax2.set_ylabel(r"$\Delta_{\rm eq}$ [meV]", color=gap_color)
    ax2.tick_params(axis="y", colors=gap_color)
    ax2.spines["right"].set_color(gap_color)

    title = (
        "Usadel calibration sweep\n"
        rf"$D={D_m2_s:.3e}\,\mathrm{{m^2/s}}$, "
        rf"$j_c={jc_A_m2:.3e}\,\mathrm{{A/m^2}}$"
    )
    ax1.set_title(title)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(
        h1 + h2,
        l1 + l2,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.20),
        ncol=2,
        fontsize=8,
        frameon=True,
    )

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output

def plot_phase_space_slices(
    catalog,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot representative phase-space slices J_S and J_R.

    The diagnostic uses the largest Delta slice, an intermediate q slice and a
    few Te values.

    NOTE [PDF Appendix A mismatch]:
    Appendix A writes the scattering energy moment with an infinite upper
    energy limit. The plotted OE4 catalogue uses a finite Usadel DOS grid.
    Therefore the vertical line at E_max - Delta and the shaded region mark a
    numerical window limitation of J_S, not a physical suppression mechanism.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    meV_J = 1.602176634e-22
    omega_meV = catalog.omega_values_J / meV_J

    delta_index = catalog.delta_values_J.size - 1
    q_index = catalog.q_values_m_inv.size // 2

    n_Te = catalog.Te_values_K.size
    if n_Te <= 3:
        Te_indices = list(range(n_Te))
    else:
        Te_indices = sorted(set([0, n_Te // 2, n_Te - 1]))

    fig, ax = plt.subplots(figsize=(7.8, 4.8))

    for iT in Te_indices:
        Te = catalog.Te_values_K[iT]
        JS_meV = catalog.J_S_TdqO_J[iT, delta_index, q_index, :] / meV_J
        JR_meV = catalog.J_R_TdqO_J[iT, delta_index, q_index, :] / meV_J

        ax.plot(
            omega_meV,
            JS_meV,
            linewidth=1.2,
            label=rf"$\mathcal{{J}}_S$, $T_e={Te:.1f}$ K",
        )
        ax.plot(
            omega_meV,
            JR_meV,
            linewidth=1.2,
            linestyle="--",
            label=rf"$\mathcal{{J}}_R$, $T_e={Te:.1f}$ K",
        )

    delta_meV = catalog.delta_values_J[delta_index] / meV_J
    gamma_meV = catalog.gamma_values_J[q_index] / meV_J

    # Recombination/pair-breaking threshold from the Appendix-A Simon/BCS
    # integration interval E in [Delta, Omega-Delta].
    jr_threshold_meV = 2.0 * delta_meV
    ax.axvline(
        jr_threshold_meV,
        linestyle=":",
        linewidth=1.0,
        alpha=0.85,
        label=rf"$2\Delta={jr_threshold_meV:.2f}$ meV",
    )

    # Finite-DOS support limit for scattering: E' = E + Omega must remain
    # inside the parent Usadel catalogue. This line is numerical, not physical.
    js_cutoffs = catalog.metadata.get("js_hard_cutoff_by_delta_J", None)
    if js_cutoffs is not None and len(js_cutoffs) == catalog.delta_values_J.size:
        js_cutoff_meV = float(js_cutoffs[delta_index]) / meV_J
    else:
        energy_max_J = float(
            catalog.metadata.get("energy_max_J", np.max(catalog.omega_values_J))
        )
        js_cutoff_meV = max(0.0, energy_max_J - catalog.delta_values_J[delta_index]) / meV_J

    if np.isfinite(js_cutoff_meV) and js_cutoff_meV > 0.0:
        ax.axvline(
            js_cutoff_meV,
            linestyle="-.",
            linewidth=1.0,
            alpha=0.85,
            label=rf"$E_{{\max}}-\Delta={js_cutoff_meV:.2f}$ meV",
        )
        if js_cutoff_meV < float(np.max(omega_meV)):
            ax.axvspan(
                js_cutoff_meV,
                float(np.max(omega_meV)),
                alpha=0.08,
                label=r"$\mathcal{J}_S$ finite-window zone",
            )

    ax.set_title(
        "Phase-space catalogue diagnostic\n"
        rf"$\Delta={delta_meV:.3f}$ meV, "
        rf"$\Gamma_q={gamma_meV:.3f}$ meV"
    )
    ax.set_xlabel(r"$\Omega$ [meV]")
    ax.set_ylabel(r"$\mathcal{J}_{S,R}$ [meV]")
    ax.grid(True, linewidth=0.25, alpha=0.35)
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        ncol=2,
        fontsize=8,
        frameon=True,
    )

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def plot_stationary_state(config, run_name):
    """Plot stationary gTDGL fields and current-conservation diagnostics."""
    return 0


def plot_photon_transient(config, run_name):
    """Plot photon-run snapshots, histories, and circuit observables."""
    return 0


def plot_catalog_diagnostics(config, run_name):
    """Plot Usadel and phase-space catalog sanity checks."""
    return 0
