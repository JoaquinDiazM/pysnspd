"""Small pyTDGL-compatible device/mesh adapters for pySNSPD OE7 data."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace
from typing import Sequence

import numpy as np

from pysnspd.gtdgl.material import GTDGLMaterial
from pysnspd.mesh.operators import FVOperators


@dataclass(frozen=True)
class TerminalInfo:
    """Minimal terminal metadata with the attributes used by pyTDGL's solver."""

    name: str
    site_indices: np.ndarray
    boundary_edge_indices: np.ndarray
    length: float


@dataclass
class EdgeMeshAdapter:
    edges: np.ndarray
    directions: np.ndarray
    normalized_directions: np.ndarray
    edge_lengths: np.ndarray
    dual_edge_lengths: np.ndarray
    centers: np.ndarray
    boundary_edge_indices: np.ndarray


@dataclass
class MeshAdapter:
    sites: np.ndarray
    areas: np.ndarray
    edge_mesh: EdgeMeshAdapter
    original_mesh: object
    original_ops: FVOperators


@dataclass
class PySNSPDTDGLDevice:
    """Minimal device accepted by ``TDGLSolver``.

    The class intentionally exposes the pyTDGL names used by the solver:
    ``mesh``, ``probe_point_indices``, ``layer.u``, ``layer.gamma`` and
    ``terminal_info``.
    """

    mesh: MeshAdapter
    material: GTDGLMaterial
    original_mesh: object
    edge_data: object
    length_scale_m: float
    current_scale_A: float
    voltage_scale_V: float
    terminal_infos: Sequence[TerminalInfo]
    probe_point_indices: np.ndarray | None = None
    layer: object | None = None

    def __post_init__(self) -> None:
        if self.layer is None:
            object.__setattr__(
                self,
                "layer",
                SimpleNamespace(u=1.0, gamma=0.0, z0=0.0),
            )

    def terminal_info(self) -> Sequence[TerminalInfo]:
        return list(self.terminal_infos)

    @property
    def terminal_neumann_current_unit_A(self) -> float:
        """Physical current corresponding to one dimensionless terminal flux.

        The pyTDGL-shaped algebra uses dimensionless coordinates ``x' = x/L0``
        and dimensionless potential ``mu = phi/V0`` internally.  pySNSPD still
        passes terminal currents in amperes.  For a film of thickness ``d``,
        the physical Neumann condition

            n · grad(phi) = -I / (sigma_n d L_terminal)

        becomes, in the internal coordinate,

            n · grad'(mu) = -I / (sigma_n d V0 L'_terminal).

        Therefore ``sigma_n * d * V0`` is the SI conversion factor from an
        integrated terminal current in amperes to the dimensionless Neumann
        value used by the pyTDGL-style Poisson operator.  This is only an
        internal operator conversion; user-facing inputs and outputs remain SI.
        """

        return float(self.material.sigma_n_S_m * self.material.thickness_m * self.voltage_scale_V)

    def terminal_mu_boundary_value(
        self,
        *,
        terminal: TerminalInfo,
        terminal_currents_A: dict[str, float],
    ) -> float:
        """Convert SI terminal currents to the internal Neumann value for mu.

        ``terminal_currents_A`` is always interpreted in amperes.  The returned
        value is the dimensionless boundary derivative required by the
        pyTDGL-like sparse Poisson operator.
        """

        current_A_per_dimless_length = (-1.0 / terminal.length) * sum(
            float(terminal_currents_A.get(name, 0.0))
            for name in terminal_currents_A
            if name != terminal.name
        )
        unit_A = max(self.terminal_neumann_current_unit_A, 1.0e-300)
        return float(current_A_per_dimless_length / unit_A)


def build_pytdgl_like_device(
    *,
    mesh,
    edge_data,
    material: GTDGLMaterial,
    ops: FVOperators,
    Te_K: np.ndarray,
    target_current_A: float,
    length_scale_m: float | None = None,
    voltage_scale_V: float | None = None,
    current_scale_A: float | None = None,
) -> PySNSPDTDGLDevice:
    """Build a pyTDGL-shaped device from pySNSPD mesh/operator objects."""

    if length_scale_m is None:
        # Use a physically meaningful scale when possible; fall back to median h.
        xi2 = material.xi_mod_squared_m2(np.asarray(Te_K, dtype=float))
        xi = float(np.sqrt(np.nanmedian(np.maximum(xi2, 1.0e-300))))
        if not np.isfinite(xi) or xi <= 0:
            xi = float(getattr(ops, "xi_mesh_m", 0.0) or np.sqrt(np.nanmedian(ops.node_area_m2)))
        length_scale_m = xi
    length_scale_m = float(length_scale_m)
    if length_scale_m <= 0:
        raise ValueError("length_scale_m must be positive.")

    if voltage_scale_V is None:
        # Josephson temporal link scale: dimensionless mu = (2e/hbar) phi * tau.
        # The first comparison backend uses tau0_GL_s as the TDGL time scale.
        from pysnspd.gtdgl.material import E_CHARGE_C, HBAR_J_S

        tau0 = float(material.tau0_GL_s)
        voltage_scale_V = HBAR_J_S / (2.0 * E_CHARGE_C * tau0)
    voltage_scale_V = float(voltage_scale_V)
    if voltage_scale_V <= 0:
        raise ValueError("voltage_scale_V must be positive.")

    if current_scale_A is None:
        # This is not a user-facing normalization of the bias current.  It is
        # the physical SI conversion factor that maps one dimensionless
        # pyTDGL-style terminal Neumann value to an integrated terminal current
        # in amperes: I_unit = sigma_n * d * V0.  The adapter still passes
        # terminal currents in amperes.
        current_scale_A = material.sigma_n_S_m * material.thickness_m * voltage_scale_V
    current_scale_A = float(current_scale_A)
    if current_scale_A <= 0:
        raise ValueError("current_scale_A must be positive.")

    nodes = np.asarray(mesh.nodes, dtype=float)[:, :2]
    sites = nodes / length_scale_m
    edge_vec = ops.edge_vec_m / length_scale_m
    edge_lengths = ops.edge_length_m / length_scale_m
    dual_lengths = ops.dual_face_length_m / length_scale_m
    areas = ops.node_area_m2 / (length_scale_m**2)

    tags = np.asarray(edge_data.tags).astype(str) if hasattr(edge_data, "tags") else np.full(ops.n_edges, "interior")
    boundary_mask = tags != "interior"
    if not np.any(boundary_mask):
        boundary_mask = np.isin(tags, ["left", "right", "top", "bottom"])
    boundary_edge_indices = np.flatnonzero(boundary_mask).astype(np.int64)

    centers = 0.5 * (sites[ops.edge_i] + sites[ops.edge_j])
    edge_mesh = EdgeMeshAdapter(
        edges=np.asarray(ops.edges, dtype=np.int64),
        directions=edge_vec,
        normalized_directions=np.asarray(ops.edge_unit, dtype=float),
        edge_lengths=edge_lengths,
        dual_edge_lengths=dual_lengths,
        centers=centers,
        boundary_edge_indices=boundary_edge_indices,
    )
    mesh_adapter = MeshAdapter(
        sites=sites,
        areas=areas,
        edge_mesh=edge_mesh,
        original_mesh=mesh,
        original_ops=ops,
    )

    terminal_infos = _terminal_infos_from_tags(tags, ops, edge_lengths)
    return PySNSPDTDGLDevice(
        mesh=mesh_adapter,
        material=material,
        original_mesh=mesh,
        edge_data=edge_data,
        length_scale_m=length_scale_m,
        current_scale_A=current_scale_A,
        voltage_scale_V=voltage_scale_V,
        terminal_infos=terminal_infos,
    )


def _terminal_infos_from_tags(tags: np.ndarray, ops: FVOperators, edge_lengths_dimless: np.ndarray) -> list[TerminalInfo]:
    infos: list[TerminalInfo] = []
    for name in ("left", "right"):
        edge_idx = np.flatnonzero(tags == name).astype(np.int64)
        if edge_idx.size:
            sites = np.unique(ops.edges[edge_idx].ravel()).astype(np.int64)
            length = float(np.sum(edge_lengths_dimless[edge_idx]))
        else:
            if name == "left" and ops.left_nodes is not None:
                sites = np.asarray(ops.left_nodes, dtype=np.int64)
                length = float(np.sum(ops.left_boundary_measure_m) / np.nanmean(ops.edge_length_m)) if ops.left_boundary_measure_m is not None else float(len(sites))
            elif name == "right" and ops.right_nodes is not None:
                sites = np.asarray(ops.right_nodes, dtype=np.int64)
                length = float(np.sum(ops.right_boundary_measure_m) / np.nanmean(ops.edge_length_m)) if ops.right_boundary_measure_m is not None else float(len(sites))
            else:
                sites = np.array([], dtype=np.int64)
                length = 0.0
        infos.append(
            TerminalInfo(
                name=name,
                site_indices=sites,
                boundary_edge_indices=edge_idx,
                length=max(length, 1.0e-300),
            )
        )
    return infos
