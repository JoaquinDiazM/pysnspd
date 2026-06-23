"""Diagnostic plots for OE7 stationary gTDGL/Poisson runs."""

from __future__ import annotations

from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

from pysnspd.gtdgl.operators import (
    boundary_currents_from_node_vectors,
    strip_transport_current_profile_from_node_vectors,
    unwrap_phase_graph,
)

MEV_J = 1.602176634e-22


def plot_ss_state_delta(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot the relaxed order-parameter amplitude in meV."""

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
    """Plot graph-unwrapped phase."""

    theta = unwrap_phase_graph(
        np.asarray(state.psi_J, dtype=np.complex128),
        np.asarray(_mesh_edges_from_triangles(mesh), dtype=np.int64),
        seed_index=_center_node_index(mesh),
        subtract_mean=False,
    )

    return _plot_node_scalar(
        mesh,
        theta,
        output_path,
        title="OE7 SS: unwrapped phase θ",
        label="θ [rad]",
        dpi=dpi,
    )


def plot_ss_state_phi(mesh, state, output_path: str | Path, *, dpi: int = 480) -> Path:
    """Plot final electrostatic potential."""

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
    """Plot electrostatic-potential snapshots during OE7 relaxation."""

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
        figsize=(3.7 * ncols, 2.7 * nrows),
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
        ax.set_title(f"t = {t_s[k] / 1.0e-12:.4g} ps")
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
    vmax = float(np.nanmax(np.abs(div))) if div.size else 1.0
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
    """Plot total current-density magnitude and sparse vectors.

    This intentionally uses the node-vector fields stored in the state, which
    are direct FV edge-to-node averages. It does not use the old LS
    reconstruction, because that diagnostic was visually doubling the current
    scale on this mesh.
    """

    del ops

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    triangles = np.asarray(mesh.triangles, dtype=np.int64)

    tri = mtri.Triangulation(x_nm, y_nm, triangles)

    jx = np.asarray(state.currents.node_jtot_x_A_m2, dtype=float)
    jy = np.asarray(state.currents.node_jtot_y_A_m2, dtype=float)
    jmag = np.sqrt(jx * jx + jy * jy)

    vmax = max(float(np.nanmax(jmag)), 1.0e-30)

    fig, ax = plt.subplots(figsize=(8.0, 3.2))

    im = ax.tripcolor(
        tri,
        jmag,
        shading="gouraud",
        vmin=0.0,
        vmax=vmax,
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$|\vec{j}|$ [A m$^{-2}$]")

    step = max(1, nodes.shape[0] // 150)
    ax.quiver(
        x_nm[::step],
        y_nm[::step],
        jx[::step],
        jy[::step],
        angles="xy",
        scale_units="xy",
        scale=None,
        width=0.002,
    )

    ax.set_title("OE7 SS: total current density")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output


def plot_ss_pairbreaking_ratio(
    mesh,
    state,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot chi_pb = xi^2 Q^2 / (1 - T/Tc).

    chi_pb = 1 is the local GL pairbreaking threshold where the stationary
    amplitude predicted by the local GL term goes to zero.
    """

    chi = np.asarray(state.currents.node_pairbreaking_ratio, dtype=float)

    vmax = float(np.nanpercentile(chi[np.isfinite(chi)], 99.5)) if np.any(np.isfinite(chi)) else 1.0
    vmax = max(vmax, 1.0)

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9
    tri = mtri.Triangulation(
        x_nm,
        y_nm,
        np.asarray(mesh.triangles, dtype=np.int64),
    )

    fig, ax = plt.subplots(figsize=(8.0, 3.2))

    im = ax.tripcolor(
        tri,
        chi,
        shading="gouraud",
        vmin=0.0,
        vmax=vmax,
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(r"$\chi_{\rm pb}=\xi^2Q^2/(1-T/T_c)$")

    if np.nanmin(chi) <= 1.0 <= np.nanmax(chi):
        ax.tricontour(tri, chi, levels=[1.0], linewidths=1.0)

    ax.set_title("OE7 SS: pairbreaking diagnostic")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output


def plot_ss_boundary_currents(
    summary: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot integrated boundary currents from the final state summary."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    currents = summary["boundary_currents_A"]
    labels = ["left", "right", "bottom", "top"]
    values = [currents[f"{label}_A"] for label in labels]

    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    ax.bar(labels, values)
    ax.axhline(0.0, linewidth=0.8)

    ax.set_title("OE7 SS: integrated boundary currents")
    ax.set_ylabel("current [A]")
    ax.grid(False)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output


def plot_ss_transport_current_profile(
    *,
    mesh,
    ops,
    state,
    output_path: str | Path,
    target_current_A: float,
    thickness_m: float,
    dpi: int = 480,
    n_bins: int = 41,
) -> Path:
    """Plot longitudinal transport-current profile.

    Uses the stored FV node-averaged jx field. This is a visualization
    diagnostic, not a solver ingredient.
    """

    del ops

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    x_m, current_A = strip_transport_current_profile_from_node_vectors(
        mesh=mesh,
        jx_A_m2=np.asarray(state.currents.node_jtot_x_A_m2, dtype=float),
        thickness_m=float(thickness_m),
        n_bins=int(n_bins),
    )

    fig, ax = plt.subplots(figsize=(8.0, 3.2))

    ax.plot(x_m * 1.0e9, current_A, marker="o", label="node-avg profile")
    ax.axhline(float(target_current_A), linestyle="--", label=r"$I_{\rm target}$")

    ax.set_title("OE7 SS: transport-current profile")
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("I(x) [A]")
    ax.grid(False)
    ax.legend()

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output


def plot_ss_relaxation_history(
    history: dict,
    output_path: str | Path,
    *,
    dpi: int = 480,
) -> Path:
    """Plot compact relaxation diagnostics."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    t_ps = np.asarray(history.get("t_s", []), dtype=float) / 1.0e-12

    fig, ax = plt.subplots(figsize=(8.0, 3.6))

    if t_ps.size:
        if "eta_R" in history:
            ax.semilogy(t_ps, np.asarray(history["eta_R"], dtype=float), label=r"$\eta_R$")
        if "current_residual" in history:
            ax.semilogy(
                t_ps,
                np.asarray(history["current_residual"], dtype=float),
                label=r"$\epsilon_{\nabla\cdot j}$",
            )
        if "terminal_voltage_V" in history:
            ax.semilogy(
                t_ps,
                np.maximum(np.abs(np.asarray(history["terminal_voltage_V"], dtype=float)), 1.0e-300),
                label=r"$|V_{\rm TDGL}|$ [V]",
            )
        if "pairbreaking_max" in history:
            ax.semilogy(
                t_ps,
                np.maximum(np.asarray(history["pairbreaking_max"], dtype=float), 1.0e-300),
                label=r"$\max \chi_{\rm pb}$",
            )
        if "delta_min_over_delta0" in history:
            ax.semilogy(
                t_ps,
                np.maximum(np.asarray(history["delta_min_over_delta0"], dtype=float), 1.0e-300),
                label=r"$\min |\Delta|/\Delta_0$",
            )

    ax.set_title("OE7 SS: relaxation diagnostics")
    ax.set_xlabel("t [ps]")
    ax.set_ylabel("diagnostic value")
    ax.grid(False)
    ax.legend()

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
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
    """Common triangular node-scalar plot."""

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    nodes = np.asarray(mesh.nodes, dtype=float)
    x_nm = nodes[:, 0] * 1.0e9
    y_nm = nodes[:, 1] * 1.0e9

    tri = mtri.Triangulation(
        x_nm,
        y_nm,
        np.asarray(mesh.triangles, dtype=np.int64),
    )

    z = np.asarray(values, dtype=float).reshape(-1)

    fig, ax = plt.subplots(figsize=(8.0, 3.2))

    im = ax.tripcolor(
        tri,
        z,
        shading="gouraud",
        vmin=vmin,
        vmax=vmax,
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(label)

    ax.set_title(title)
    ax.set_xlabel("x [nm]")
    ax.set_ylabel("y [nm]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)

    fig.tight_layout()
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)

    return output


def _center_node_index(mesh) -> int:
    """Return node closest to the geometric center."""

    nodes = np.asarray(mesh.nodes, dtype=float)
    center = np.array(
        [
            0.5 * (float(np.min(nodes[:, 0])) + float(np.max(nodes[:, 0]))),
            0.5 * (float(np.min(nodes[:, 1])) + float(np.max(nodes[:, 1]))),
        ]
    )
    dist2 = np.sum((nodes[:, :2] - center[None, :]) ** 2, axis=1)
    return int(np.argmin(dist2))


def _mesh_edges_from_triangles(mesh) -> np.ndarray:
    """Build unique undirected edges from mesh triangles for phase unwrapping."""

    tri = np.asarray(mesh.triangles, dtype=np.int64)
    edges = np.vstack(
        [
            tri[:, [0, 1]],
            tri[:, [1, 2]],
            tri[:, [2, 0]],
        ]
    )
    edges = np.sort(edges, axis=1)
    edges = np.unique(edges, axis=0)
    return edges