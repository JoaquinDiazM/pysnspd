"""Diagnostic plots for OE7 stationary gTDGL/Poisson runs."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from pysnspd.gtdgl.operators import unwrap_phase_graph

MEV_J = 1.602176634e-22


def plot_ss_state_delta(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot the relaxed order-parameter amplitude in meV."""
    return _plot_node_scalar(
        mesh,
        np.abs(state.psi_J) / MEV_J,
        output_path,
        title="OE7 SS: relaxed Δ",
        label="Δ [meV]",
        dpi=dpi,
    )


def plot_ss_state_phase(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot an x-sorted unwrapped phase diagnostic."""
    theta = _unwrap_phase_by_x(mesh, np.angle(state.psi_J))
    return _plot_node_scalar(
        mesh,
        theta,
        output_path,
        title="OE7 SS: unwrapped phase θ",
        label="θ [rad]",
        dpi=dpi,
    )


def plot_ss_state_phi(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot electrostatic potential."""
    return _plot_node_scalar(
        mesh,
        state.phi_V,
        output_path,
        title="OE7 SS: electrostatic potential φ",
        label="φ [V]",
        dpi=dpi,
    )


def plot_ss_state_divergence(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot finite-volume current divergence."""
    div = np.asarray(state.currents.node_div_jtot_A_m3, dtype=float)
    vmax = float(np.max(np.abs(div))) if div.size else 1.0
    vmax = max(vmax, 1.0e-30)
    return _plot_node_scalar(
        mesh,
        div,
        output_path,
        title="OE7 SS: finite-volume div(j)",
        label="div(j) [A m$^{-3}$]",
        vmin=-vmax,
        vmax=vmax,
        dpi=dpi,
    )


def plot_ss_state_current_density(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot total current-density magnitude and node-averaged vectors."""
    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    triangles = np.asarray(mesh.triangles, dtype=np.int64)

    jx = np.asarray(state.currents.node_jtot_x_A_m2, dtype=float)
    jy = np.asarray(state.currents.node_jtot_y_A_m2, dtype=float)
    mag = np.sqrt(jx**2 + jy**2)

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.0, 3.2), constrained_layout=True)
    tri = mtri.Triangulation(x_nm, y_nm, triangles)
    im = ax.tripcolor(tri, mag, shading="gouraud")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$|\vec{j}|$ [A m$^{-2}$]")

    n = max(1, mag.size // 120)
    scale = np.nanmax(mag)
    if np.isfinite(scale) and scale > 0.0:
        ax.quiver(
            x_nm[::n],
            y_nm[::n],
            jx[::n] / scale,
            jy[::n] / scale,
            angles="xy",
            scale_units="xy",
            scale=0.040,
            width=0.0025,
        )

    ax.set_title("OE7 SS: total current density")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_ss_boundary_currents(summary: dict, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot integrated terminal and transverse boundary currents."""
    boundary = dict(summary["boundary_currents_A"])
    labels = ["left", "right", "bottom", "top"]
    values = [boundary.get(f"{name}_A", 0.0) for name in labels]

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6.2, 3.2), constrained_layout=True)
    ax.bar(labels, values)
    ax.axhline(0.0, linewidth=0.8)
    ax.set_title("OE7 SS: integrated boundary currents")
    ax.set_ylabel("current [A]")
    ax.grid(False)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def plot_ss_relaxation_history(history: dict, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot stationary relaxation residual history."""
    t_ps = np.asarray(history["t_s"], dtype=float) / 1.0e-12
    eta = np.asarray(history["eta_R"], dtype=float)
    residual = np.asarray(history["current_residual"], dtype=float)
    voltage = np.abs(np.asarray(history["terminal_voltage_V"], dtype=float))

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6.4, 3.4), constrained_layout=True)
    if t_ps.size:
        ax.semilogy(t_ps, np.maximum(eta, 1.0e-300), label=r"$\eta_R$")
        ax.semilogy(t_ps, np.maximum(residual, 1.0e-300), label=r"$\epsilon_{\nabla\cdot j}$")
        ax.semilogy(t_ps, np.maximum(voltage, 1.0e-300), label=r"$|V_{\rm TDGL}|$ [V]")
    ax.set_title("OE7 SS: relaxation diagnostics")
    ax.set_xlabel("t [ps]")
    ax.set_ylabel("diagnostic value")
    ax.grid(False)
    ax.legend(frameon=False)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def _plot_node_scalar(
    mesh,
    values,
    output_path: str | Path,
    *,
    title: str,
    label: str,
    vmin=None,
    vmax=None,
    dpi: int = 480,
) -> Path:
    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    triangles = np.asarray(mesh.triangles, dtype=np.int64)
    z = np.asarray(values, dtype=float).reshape(-1)

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.0, 3.2), constrained_layout=True)
    tri = mtri.Triangulation(x_nm, y_nm, triangles)
    im = ax.tripcolor(tri, z, shading="gouraud", vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(label)
    ax.set_title(title)
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output


def _unwrap_phase_by_x(mesh, theta_wrapped: np.ndarray) -> np.ndarray:
    psi = np.exp(1j * np.asarray(theta_wrapped, dtype=float))
    edges = _edges_from_triangles(mesh.triangles)
    return unwrap_phase_graph(psi, edges)


def _edges_from_triangles(triangles: np.ndarray) -> np.ndarray:
    triangles = np.asarray(triangles, dtype=np.int64)
    pairs = np.vstack(
        [
            triangles[:, [0, 1]],
            triangles[:, [1, 2]],
            triangles[:, [2, 0]],
        ]
    )
    pairs.sort(axis=1)
    return np.unique(pairs, axis=0)