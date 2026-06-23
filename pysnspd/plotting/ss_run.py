"""Diagnostic plots for OE7 stationary gTDGL/Poisson runs."""
from __future__ import annotations

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from pysnspd.gtdgl.operators import unwrap_phase_graph

from pysnspd.gtdgl.operators import (
    boundary_currents_from_edge_scalar_least_squares,
    boundary_currents_from_node_vectors,
    edge_scalar_to_node_vector_least_squares,
    strip_transport_current_profile_from_node_vectors,
)

MEV_J = 1.602176634e-22


def plot_ss_state_delta(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot the relaxed order-parameter amplitude in meV.

    The color scale starts at zero to avoid visually amplifying tiny numerical
    variations around the nearly uniform superconducting gap.
    """
    delta_meV = np.abs(state.psi_J) / MEV_J
    vmax = max(float(np.nanmax(delta_meV)), 1.0e-30)

    return _plot_node_scalar(
        mesh,
        delta_meV,
        output_path,
        title="OE7 SS: relaxed Δ",
        label="Δ [meV]",
        vmin=0.0,
        vmax=vmax,
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

def plot_ss_phi_snapshots(
    mesh,
    history: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
    ncols: int = 3,
) -> Path:
    """Plot electrostatic-potential snapshots during OE7 relaxation.

    Expects history keys:
        phi_snapshot_t_s : shape [n_snapshots]
        phi_snapshot_V   : shape [n_snapshots, n_nodes]
    """
    if "phi_snapshot_V" not in history or "phi_snapshot_t_s" not in history:
        raise KeyError("history must contain phi_snapshot_V and phi_snapshot_t_s.")

    phi = np.asarray(history["phi_snapshot_V"], dtype=float)
    t_s = np.asarray(history["phi_snapshot_t_s"], dtype=float)

    if phi.ndim != 2:
        raise ValueError(f"phi_snapshot_V must be 2D, got shape {phi.shape}.")
    if t_s.ndim != 1 or t_s.size != phi.shape[0]:
        raise ValueError(
            "phi_snapshot_t_s must be 1D and match the number of phi snapshots."
        )

    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    triangles = np.asarray(mesh.triangles, dtype=np.int64)

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    n_snap = int(phi.shape[0])
    ncols = max(1, int(ncols))
    nrows = int(np.ceil(n_snap / ncols))

    vmax = float(np.nanmax(np.abs(phi))) if phi.size else 1.0
    vmax = max(vmax, 1.0e-30)
    vmin = -vmax

    tri = mtri.Triangulation(x_nm, y_nm, triangles)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3.4 * ncols, 2.7 * nrows),
        constrained_layout=True,
        squeeze=False,
    )

    last_im = None
    for k in range(nrows * ncols):
        ax = axes.flat[k]
        if k >= n_snap:
            ax.axis("off")
            continue

        last_im = ax.tripcolor(
            tri,
            phi[k],
            shading="gouraud",
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_title(f"t = {t_s[k] / 1.0e-12:.3g} ps")
        ax.set_xlabel("x [nm]")
        ax.set_ylabel("y [nm]")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(False)

    if last_im is not None:
        cbar = fig.colorbar(last_im, ax=axes.ravel().tolist())
        cbar.set_label("φ [V]")

    fig.suptitle("OE7 SS: electrostatic potential φ snapshots")
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output

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


def plot_ss_state_current_density(
    mesh,
    state,
    output_path: str | Path,
    *,
    ops=None,
    dpi: int = 480,
) -> Path:
    """Plot total current-density magnitude and vectors.

    If FV operators are provided, reconstruct the node vector field from edge
    current projections using local least squares. This matches the diagnostic
    philosophy of the older notebook better than the simple node average.
    """
    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    triangles = np.asarray(mesh.triangles, dtype=np.int64)

    if ops is not None:
        jx, jy = edge_scalar_to_node_vector_least_squares(
            state.currents.edge_jtot_A_m2,
            ops,
        )
        title = "OE7 SS: total current density"
    else:
        jx = np.asarray(state.currents.node_jtot_x_A_m2, dtype=float)
        jy = np.asarray(state.currents.node_jtot_y_A_m2, dtype=float)
        title = "OE7 SS: total current density"

    mag = np.sqrt(jx**2 + jy**2)
    vmax = max(float(np.nanmax(mag)), 1.0e-30)

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.0, 3.2), constrained_layout=True)
    tri = mtri.Triangulation(x_nm, y_nm, triangles)

    im = ax.tripcolor(tri, mag, shading="gouraud", vmin=0.0, vmax=vmax)
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
            scale=0.030,
            width=0.002,
        )

    ax.set_title(title)
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

def plot_ss_boundary_current_reconstruction_comparison(
    *,
    mesh,
    edge_data,
    ops,
    state,
    output_path: str | Path,
    target_current_A: float | None = None,
    thickness_m: float,
    dpi: int = 480,
) -> Path:
    """Compare terminal currents from different diagnostic reconstructions."""
    node_avg = boundary_currents_from_node_vectors(
        mesh=mesh,
        edge_data=edge_data,
        jx_A_m2=state.currents.node_jtot_x_A_m2,
        jy_A_m2=state.currents.node_jtot_y_A_m2,
        thickness_m=thickness_m,
    )

    ls = boundary_currents_from_edge_scalar_least_squares(
        mesh=mesh,
        edge_data=edge_data,
        ops=ops,
        edge_current_i_to_j=state.currents.edge_jtot_A_m2,
        thickness_m=thickness_m,
    )

    labels = [
        "left\nnode avg",
        "left\nLS",
        "right\nnode avg",
        "right\nLS",
    ]
    values = [
        node_avg.get("left_A", 0.0),
        ls.get("left_A", 0.0),
        node_avg.get("right_A", 0.0),
        ls.get("right_A", 0.0),
    ]

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6.8, 3.4), constrained_layout=True)
    ax.bar(labels, values)
    ax.axhline(0.0, linewidth=0.8)

    if target_current_A is not None:
        I = float(target_current_A)
        ax.axhline(+I, linestyle="--", linewidth=0.9, label=r"$+I_{\rm target}$")
        ax.axhline(-I, linestyle="--", linewidth=0.9, label=r"$-I_{\rm target}$")
        ax.legend(frameon=False)

    ax.set_title("OE7 SS: boundary-current reconstruction")
    ax.set_ylabel("current [A]")
    ax.grid(False)
    fig.savefig(output, dpi=dpi)
    plt.close(fig)
    return output

def _strip_transport_current_profile_from_edges(
    *,
    mesh,
    ops,
    edge_current_i_to_j_A_m2: np.ndarray,
    thickness_m: float,
    target_current_A: float | None = None,
    n_cuts: int = 41,
) -> tuple[np.ndarray, np.ndarray]:
    """Notebook-style longitudinal current profile from edge currents.

    For each vertical cut x=x_c, sum the x-projected edge current carried by
    edges crossing that cut:

        I(x_c) ~= d * sum_edges_crossing_cut j_e e_x s_e.

    Endpoints are set to target_current_A when available, because the terminal
    current is imposed as a boundary condition.
    """
    nodes = np.asarray(mesh.nodes, dtype=float)
    current = np.asarray(edge_current_i_to_j_A_m2, dtype=float)

    if current.shape != (ops.n_edges,):
        raise ValueError(
            f"edge_current_i_to_j_A_m2 must have shape ({ops.n_edges},), "
            f"got {current.shape}."
        )

    if n_cuts < 2:
        raise ValueError("n_cuts must be at least 2.")

    x = nodes[:, 0]
    xmin = float(np.min(x))
    xmax = float(np.max(x))

    xs = np.linspace(xmin, xmax, int(n_cuts))
    I_A = np.zeros_like(xs)

    if target_current_A is not None and np.isfinite(float(target_current_A)):
        I_A[0] = float(target_current_A)
        I_A[-1] = float(target_current_A)
    else:
        I_A[0] = np.nan
        I_A[-1] = np.nan

    xi = x[ops.edge_i]
    xj = x[ops.edge_j]

    jx_edge = current * ops.edge_unit[:, 0]
    face = ops.dual_face_length_m

    eps = max(1.0e-30, 1.0e-12 * max(xmax - xmin, 1.0e-300))

    for k, xc in enumerate(xs[1:-1], start=1):
        crosses = ((xi - xc) * (xj - xc) <= 0.0) & (np.abs(xj - xi) > eps)

        if np.any(crosses):
            I_A[k] = float(
                thickness_m * np.sum(jx_edge[crosses] * face[crosses])
            )
        else:
            I_A[k] = np.nan

    return xs, I_A


def plot_ss_transport_current_profile(
    *,
    mesh,
    ops,
    state,
    output_path: str | Path,
    target_current_A: float | None = None,
    thickness_m: float,
    n_bins: int = 41,
    dpi: int = 480,
) -> Path:
    """Plot notebook-style longitudinal transport-current profile."""
    x_m, I_A = _strip_transport_current_profile_from_edges(
        mesh=mesh,
        ops=ops,
        edge_current_i_to_j_A_m2=state.currents.edge_jtot_A_m2,
        thickness_m=thickness_m,
        target_current_A=target_current_A,
        n_cuts=n_bins,
    )

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(6.6, 3.4), constrained_layout=True)

    ax.plot(
        x_m * 1.0e9,
        I_A,
        marker="o",
        markersize=2.5,
        linewidth=1.0,
        label="edge-cut profile",
    )

    if target_current_A is not None:
        ax.axhline(
            float(target_current_A),
            linestyle="--",
            linewidth=0.9,
            label=r"$I_{\rm target}$",
        )

    ax.legend(frameon=False)
    ax.set_title("OE7 SS: transport-current profile")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("I(x) [A]")
    ax.grid(False)

    fig.savefig(output, dpi=dpi)
    plt.close(fig)

    return output