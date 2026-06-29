"""pyTDGL-like finite-volume operators adapted to pySNSPD meshes.

The public function and method names mirror ``tdgl.finite_volume.operators``:
``build_divergence``, ``build_gradient``, ``build_laplacian``,
``build_neumann_boundary_laplacian`` and ``MeshOperators``.  Internally the
module consumes the small mesh adapter defined in ``device.py``.
"""
from __future__ import annotations

from typing import Callable, Tuple, Union

import numpy as np
import scipy.sparse as sp

from .options import SparseSolver


def build_divergence(mesh) -> sp.csr_array:
    """Build the divergence matrix that maps edge fields onto sites."""

    edge_mesh = mesh.edge_mesh
    edge_indices = np.arange(len(edge_mesh.edges))
    weights = edge_mesh.dual_edge_lengths
    edges0 = edge_mesh.edges[:, 0]
    edges1 = edge_mesh.edges[:, 1]
    rows = np.concatenate([edges0, edges1])
    cols = np.concatenate([edge_indices, edge_indices])
    values = np.concatenate(
        [weights / mesh.areas[edges0], -weights / mesh.areas[edges1]]
    )
    return sp.csr_array((values, (rows, cols)), shape=(len(mesh.sites), len(edge_mesh.edges)))


def build_gradient(
    mesh,
    link_exponents: Union[np.ndarray, None] = None,
    weights: Union[np.ndarray, None] = None,
) -> sp.csr_array:
    """Build the covariant edge gradient for node fields."""

    edge_mesh = mesh.edge_mesh
    edge_indices = np.arange(len(edge_mesh.edges))
    if weights is None:
        weights = 1 / edge_mesh.edge_lengths
    if link_exponents is None:
        link_variable_weights = np.ones(len(weights), dtype=np.complex128)
    else:
        link_variable_weights = np.exp(
            -1j * np.einsum("ij, ij -> i", link_exponents, edge_mesh.directions)
        )
    rows = np.concatenate([edge_indices, edge_indices])
    cols = np.concatenate([edge_mesh.edges[:, 1], edge_mesh.edges[:, 0]])
    values = np.concatenate([link_variable_weights * weights, -weights])
    return sp.csr_array((values, (rows, cols)), shape=(len(edge_mesh.edges), len(mesh.sites)))


def build_laplacian(
    mesh,
    link_exponents: Union[np.ndarray, None] = None,
    fixed_sites: Union[np.ndarray, None] = None,
    free_rows: Union[np.ndarray, None] = None,
    fixed_sites_eigenvalues: float = 1,
    weights: Union[np.ndarray, None] = None,
) -> Tuple[sp.csc_array, np.ndarray]:
    """Build the pyTDGL-style finite-volume Laplacian matrix."""

    if fixed_sites is None:
        fixed_sites = np.array([], dtype=int)
    fixed_sites = np.asarray(fixed_sites, dtype=int)

    edge_mesh = mesh.edge_mesh
    if weights is None:
        weights = edge_mesh.dual_edge_lengths / edge_mesh.edge_lengths
    if link_exponents is None:
        link_variable_weights = np.ones(len(weights), dtype=np.complex128)
    else:
        link_variable_weights = np.exp(
            -1j * np.einsum("ij, ij -> i", link_exponents, edge_mesh.directions)
        )

    edges0 = edge_mesh.edges[:, 0]
    edges1 = edge_mesh.edges[:, 1]
    rows = np.concatenate([edges0, edges1, edges0, edges1])
    cols = np.concatenate([edges1, edges0, edges0, edges1])
    areas0 = mesh.areas[edges0]
    areas1 = mesh.areas[edges1]
    values = np.concatenate(
        [
            weights * link_variable_weights / areas0,
            weights * link_variable_weights.conjugate() / areas1,
            -weights / areas0,
            -weights / areas1,
        ]
    )

    if free_rows is None:
        free_rows = np.isin(rows, fixed_sites, invert=True)
    rows = rows[free_rows]
    cols = cols[free_rows]
    values = values[free_rows]

    rows = np.concatenate([rows, fixed_sites])
    cols = np.concatenate([cols, fixed_sites])
    values = np.concatenate(
        [values, fixed_sites_eigenvalues * np.ones(len(fixed_sites), dtype=float)]
    )
    laplacian = sp.csc_array((values, (rows, cols)), shape=(len(mesh.sites), len(mesh.sites)))
    return laplacian, free_rows


def build_neumann_boundary_laplacian(
    mesh,
    fixed_sites: Union[np.ndarray, None] = None,
) -> sp.csr_array:
    """Build the matrix multiplying non-homogeneous Neumann boundary values."""

    edge_mesh = mesh.edge_mesh
    boundary_index = np.arange(len(edge_mesh.boundary_edge_indices))
    boundary_edges = edge_mesh.edges[edge_mesh.boundary_edge_indices]
    boundary_edges_length = edge_mesh.edge_lengths[edge_mesh.boundary_edge_indices]
    rows = np.concatenate([boundary_edges[:, 0], boundary_edges[:, 1]])
    cols = np.concatenate([boundary_index, boundary_index])
    values = np.concatenate(
        [
            boundary_edges_length / (2 * mesh.areas[boundary_edges[:, 0]]),
            boundary_edges_length / (2 * mesh.areas[boundary_edges[:, 1]]),
        ]
    )
    neumann_laplacian = sp.csr_array(
        (values, (rows, cols)), shape=(len(mesh.sites), len(boundary_index))
    )
    if fixed_sites is not None:
        neumann_laplacian = neumann_laplacian.tolil()
        neumann_laplacian[np.asarray(fixed_sites, dtype=int), :] = 0
    return neumann_laplacian.tocsr(copy=False)


class MeshOperators:
    """A container for the finite volume operators for a given mesh."""

    def __init__(
        self,
        mesh,
        sparse_solver: SparseSolver,
        use_cupy: bool = False,
        fixed_sites: Union[np.ndarray, None] = None,
        fix_psi: bool = True,
    ):
        if use_cupy:
            raise NotImplementedError("pySNSPD pytdgl_like implements CPU/SuperLU only.")
        self.mesh = mesh
        self.areas = mesh.areas
        edge_mesh = mesh.edge_mesh
        self.edges = edge_mesh.edges
        self.edge_directions = edge_mesh.directions
        self.use_cupy = False
        self.sparse_solver = sparse_solver
        self.fixed_sites = np.array([], dtype=int) if fixed_sites is None else np.asarray(fixed_sites, dtype=int)
        self.fix_psi = bool(fix_psi)
        self.laplacian_free_rows: Union[np.ndarray, None] = None
        self.divergence: Union[sp.spmatrix, None] = None
        self.mu_laplacian: Union[sp.spmatrix, None] = None
        self.mu_boundary_laplacian: Union[sp.spmatrix, None] = None
        self.mu_laplacian_lu: Union[Callable, None] = None
        self.mu_gradient: Union[sp.spmatrix, None] = None
        self.psi_gradient: Union[sp.spmatrix, None] = None
        self.psi_laplacian: Union[sp.spmatrix, None] = None
        self.link_exponents: Union[np.ndarray, None] = None
        self.gradient_weights = 1 / edge_mesh.edge_lengths
        self.laplacian_weights = edge_mesh.dual_edge_lengths / edge_mesh.edge_lengths

    def build_operators(self) -> None:
        """Construct vector-potential-independent operators."""

        mesh = self.mesh
        self.mu_laplacian, _ = build_laplacian(mesh, weights=self.laplacian_weights)
        self.mu_boundary_laplacian = build_neumann_boundary_laplacian(mesh)
        self.mu_gradient = build_gradient(mesh, weights=self.gradient_weights)
        self.divergence = build_divergence(mesh)
        self.mu_laplacian_lu = sp.linalg.factorized(self.mu_laplacian)

    def set_link_exponents(self, link_exponents: np.ndarray) -> None:
        """Set link variables and construct covariant gradient/laplacian for psi."""

        mesh = self.mesh
        self.link_exponents = np.asarray(link_exponents, dtype=float)
        self.psi_gradient = build_gradient(
            mesh,
            link_exponents=self.link_exponents,
            weights=self.gradient_weights,
        )
        fixed_sites = self.fixed_sites if self.fix_psi else None
        self.psi_laplacian, self.laplacian_free_rows = build_laplacian(
            mesh,
            link_exponents=self.link_exponents,
            fixed_sites=fixed_sites,
            weights=self.laplacian_weights,
        )

    def get_supercurrent(self, psi: np.ndarray):
        """Compute the pyTDGL dimensionless supercurrent on the edges."""

        return (psi.conjugate()[self.edges[:, 0]] * (self.psi_gradient @ psi)).imag
