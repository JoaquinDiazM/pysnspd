"""Static full-2D tensor-GLL common energy with independent boundary nodes.

This opt-in geometry reuses the mixed functional's energy, Cartesian gradient,
gauge-link current and local 2D principal symbol. Its prolongation is identity:
there are no 1D leads, transverse trace identifications or implicit reservoirs.
Boundary values are free input data; this object neither fixes a boundary
condition nor evolves a transient. Dense matrices target small diagnostics.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix, eye

from .mixed_spatial import MixedSpatialFunctional
from .spatial_functional import PeriodicSpatialFunctional
from .spatial_open import lobatto_operators


class RectangularSpatialFunctional(MixedSpatialFunctional):
    """Rectangle x in [0,L], y in [-W/2,W/2], with one field per GLL node.

    The quadrature mass integrates area divided by W/ell0. Thus its sum is
    L/ell0, and the common energy unit remains N0*Delta0**2*(W*d)*ell0.
    All four sides retain their transverse variations. ``boundary_dofs``
    identifies sides without choosing a terminal injection or constraint.
    """

    scheme = "full2D_GLL_rectangle_independent_boundaries_v1"

    def __init__(self, catalog, length_m, width_m, thickness_m, *,
                 elements_x, elements_y, degree=4, delta_regularizer_bar=.1):
        for value in (length_m, width_m, thickness_m):
            if not np.isfinite(value) or value <= 0:
                raise ValueError("all geometric lengths must be finite and positive")
        PeriodicSpatialFunctional.__init__(
            self, catalog, length_m, width_m*thickness_m, 3,
            delta_regularizer_bar=delta_regularizer_bar)
        self.width_m = float(width_m)
        self.thickness_m = float(thickness_m)
        self.rectangle_length_m = float(length_m)
        self.elements_x, self.elements_y, self.degree = elements_x, elements_y, degree
        ell = self.ell0_m
        width_bar = width_m/ell
        x, mx, dx, kx = lobatto_operators(elements_x, length_m/ell, degree)
        y, my, dy, ky = lobatto_operators(elements_y, width_bar, degree)
        y -= width_bar/2
        nx, ny = len(x), len(y)
        self.rectangle_shape = (nx, ny)
        self.cells = self.quadrature_size = nx*ny
        self.h_bar = None
        self.quadrature_to_dof = np.arange(self.cells)
        self.rectangle_quadrature = self.quadrature_to_dof.reshape(nx, ny)
        self.sector_slices = {"rectangle": slice(0, self.cells)}
        self.interface_dofs = {}
        self.boundary_dofs = dict(left=self.rectangle_quadrature[0].copy(),
                                  right=self.rectangle_quadrature[-1].copy(),
                                  bottom=self.rectangle_quadrature[:, 0].copy(),
                                  top=self.rectangle_quadrature[:, -1].copy())
        self.prolongation = eye(self.cells, format="csr")
        self.quadrature_coordinates_bar = np.column_stack((np.repeat(x, ny), np.tile(y, nx)))
        self.dof_coordinates_bar = self.quadrature_coordinates_bar.copy()
        self.quadrature_mass_bar = np.kron(mx, my/width_bar)
        self.mass_bar = self.quadrature_mass_bar.copy()
        self.quadrature_volumes_m3 = self.cross_section_m2*ell*self.quadrature_mass_bar
        self.node_volumes_m3 = self.quadrature_volumes_m3.copy()
        self.spatial_dimensions = np.full(self.cells, 2, dtype=int)
        self._derivatives = np.array([np.kron(dx, np.eye(ny)), np.kron(np.eye(nx), dy)])
        self._stiffnesses = np.array([np.kron(kx, np.diag(my/width_bar)),
                                     np.kron(np.diag(mx/width_bar), ky)])
        self._lines = []
        edges, lengths, areas = [], [], []

        def line(indices, axis, area):
            indices = np.asarray(indices, dtype=int)
            edge_ids = np.arange(len(edges), len(edges)+len(indices)-1)
            edges.extend(zip(indices[:-1], indices[1:]))
            lengths.extend(np.diff(self.quadrature_coordinates_bar[indices, axis])*ell)
            areas.extend([area]*(len(indices)-1))
            self._lines.append((indices, axis, edge_ids))

        for j in range(ny):
            line(self.rectangle_quadrature[:, j], 0, thickness_m*ell*my[j])
        for i in range(nx):
            line(self.rectangle_quadrature[i], 1, thickness_m*ell*mx[i])
        self.graph_edges = np.asarray(edges, dtype=int)
        self.edge_lengths_m = np.asarray(lengths)
        self.face_areas_m2 = np.asarray(areas)
        nedge = len(edges)
        self.incidence = coo_matrix((np.r_[np.ones(nedge), -np.ones(nedge)],
            (np.r_[self.graph_edges[:, 0], self.graph_edges[:, 1]],
             np.tile(np.arange(nedge), 2))), shape=(self.cells, nedge)).tocsr()

    def _interface_reactions(self, raw_gradient):
        """No identified interface exists in the full rectangle."""
        return {}
