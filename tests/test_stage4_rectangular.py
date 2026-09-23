"""Manufactured full-2D identities; no transient or material admission."""
from types import SimpleNamespace

import numpy as np
import pytest

from pysnspd.experimental.cell_closures import CellScales, KWTMobility
from pysnspd.experimental.energy_catalog import BCS_GAP_RATIO, HBAR_J_S, K_B_J_K
from pysnspd.experimental.rectangular_spatial import RectangularSpatialFunctional
from pysnspd.experimental.spatial_dynamics import kwt_spatial_response


class SyntheticCatalog:
    def __init__(self):
        self.count_nodes = np.array([.1, .5, 1.2])
        self.count_weights = np.array([.2, .3, .5])
        self.vacuum = SimpleNamespace(metadata={"Tc_K": 8.65}, D_m2_s=5e-5,
            N0_per_J_m3=2.1e47, delta0_J=BCS_GAP_RATIO*K_B_J_K*8.65,
            gamma_axis=np.array([0., 2.]), delta_axis=np.array([0., 2.]),
            evaluate=self.vacuum_values)

    @staticmethod
    def vacuum_values(a, gamma):
        if not 0 <= a <= 2 or not 0 <= gamma <= 2:
            raise ValueError("synthetic support")
        return (.5*(a*a-1)**2+.3*gamma*a*a+.07*gamma**2,
                2*a*(a*a-1)+.6*gamma*a, .3*a*a+.14*gamma)

    def energy_kernel(self, a, gamma):
        self.vacuum_values(a, gamma)
        root = np.sqrt(self.count_nodes**2+a*a)
        alpha = .13+.03*self.count_nodes
        return (root+alpha*gamma*(1+.2*a*a), a/root+.4*alpha*gamma*a,
                alpha*(1+.2*a*a))


def gradient(result):
    return result.gradient_cartesian_bar[:, 0]+1j*result.gradient_cartesian_bar[:, 1]


def fd5(function, h=2e-5):
    return (function(-2*h)-8*function(-h)+8*function(h)-function(2*h))/(12*h)


@pytest.fixture
def case():
    catalog = SyntheticCatalog()
    ell = np.sqrt(HBAR_J_S*catalog.vacuum.D_m2_s/(2*K_B_J_K*8.65))
    model = RectangularSpatialFunctional(catalog, 8*ell, 4*ell, 7e-9,
        elements_x=2, elements_y=1, degree=3, delta_regularizer_bar=.2)
    x, y = model.dof_coordinates_bar.T
    z = (.88+.004*x+.007*y)*np.exp(1j*(.02*x+.014*y))
    links = .004*np.sin(np.arange(len(model.graph_edges))*.7)
    p = np.array([(.045+.0002*i)*np.exp(-catalog.count_nodes)
                  for i in range(model.quadrature_size)])
    return model, z, links, p


def test_independent_boundaries_and_physical_volume(case):
    m, z, _, _ = case
    np.testing.assert_array_equal(m.quadrature_to_dof, np.arange(m.cells))
    np.testing.assert_array_equal(m.prolongation.toarray(), np.eye(m.cells))
    assert m.interface_dofs == {}
    assert m.cells == m.quadrature_size == np.prod(m.rectangle_shape)
    assert m.mass_bar.sum() == pytest.approx(m.length_bar)
    assert m.node_volumes_m3.sum() == pytest.approx(m.length_m*m.width_m*m.thickness_m)
    assert np.all(m.mass_bar > 0) and np.all(m.face_areas_m2 > 0)
    assert np.all(m.edge_lengths_m > 0)
    assert len(set(m.boundary_dofs["left"])) == m.rectangle_shape[1]
    assert np.ptp(z[m.boundary_dofs["left"]].imag) > 0  # Retains transverse variation.
    assert len(m.graph_edges) == (m.rectangle_shape[0]-1)*m.rectangle_shape[1]+(m.rectangle_shape[1]-1)*m.rectangle_shape[0]


def test_tensor_polynomial_derivatives_and_positive_remainder(case):
    m, z, links, p = case
    x, y = m.dof_coordinates_bar.T
    field = .9+.001*x*x*y+.003*y*y
    fields = m.sample_fields(field)
    np.testing.assert_allclose(fields.derivative_quadrature_bar[:, 0], .002*x*y, atol=3e-15)
    np.testing.assert_allclose(fields.derivative_quadrature_bar[:, 1], .001*x*x+.006*y, atol=3e-15)
    _, _, derivative, stiffness = m._operators(z, links)
    for d, k in zip(derivative, stiffness):
        compatible = k-d.conj().T@(m.mass_bar[:, None]*d)
        assert np.linalg.eigvalsh(compatible)[0] > -1e-12
    result = m.evaluate(z, p, link_phases=links, require_stability=False)
    assert result.gradient_remainder_bar >= 0
    assert result.interface_reactions_bar == {}


def test_cartesian_and_current_variations_include_all_boundaries(case):
    m, z, links, p = case
    result = m.evaluate(z, p, link_phases=links, require_stability=False)
    dz = .07*np.cos(np.arange(m.cells))+.09j*np.sin(np.arange(m.cells))
    dl = .05*np.cos(np.arange(len(links)))
    def varied(t):
        return m.evaluate(z+t*dz, p, link_phases=links+t*dl, require_stability=False)
    predicted = np.real(np.vdot(gradient(result), dz))+np.dot(result.link_current_bar, dl)
    assert fd5(lambda t: varied(t).energy_bar) == pytest.approx(predicted, abs=3e-9)
    # Independent endpoint variations check that transverse boundary DOFs are not collapsed.
    for node in (m.boundary_dofs["left"][1], m.boundary_dofs["top"][-2]):
        direction = np.zeros(m.cells, complex)
        direction[node] = 1j
        measured = fd5(lambda t: m.evaluate(z+t*direction, p, link_phases=links,
                                           require_stability=False).energy_bar)
        assert measured == pytest.approx(gradient(result)[node].imag, abs=3e-9)


def test_gauge_noether_and_positive_full2D_symbol(case):
    m, z, links, p = case
    result = m.evaluate(z, p, link_phases=links)
    chi = .17*np.cos(np.arange(m.cells)*.9)
    newlinks = links+chi[m.graph_edges[:, 0]]-chi[m.graph_edges[:, 1]]
    changed = m.evaluate(z*np.exp(1j*chi), p, link_phases=newlinks)
    assert changed.energy_bar == pytest.approx(result.energy_bar, abs=5e-13)
    np.testing.assert_allclose(gradient(changed), gradient(result)*np.exp(1j*chi), atol=5e-13)
    np.testing.assert_allclose(changed.link_current_bar, result.link_current_bar, atol=5e-13)
    np.testing.assert_allclose(result.noether_residual_bar, 0., atol=5e-14)
    np.testing.assert_allclose(result.line_phase_residuals_bar, 0., atol=5e-14)
    assert len(result.principal_symbols) == m.cells
    assert all(s.matrix.shape == (4, 4) and s.stable for s in result.principal_symbols)


def test_thermal_envelope_and_one_hot_KWT_dissipation(case):
    m, z, links, _ = case
    thermal = m.evaluate_thermal(z, .25, link_phases=links, require_stability=False)
    fixed = m.evaluate(z, thermal.p_quadrature, link_phases=links, require_stability=False)
    np.testing.assert_allclose(thermal.gradient_cartesian_bar, fixed.gradient_cartesian_bar, atol=1e-14)
    dz = .1*np.cos(np.arange(m.cells))+.1j*np.sin(np.arange(m.cells))
    numerical = fd5(lambda t: m.evaluate_thermal(z+t*dz, .25, link_phases=links,
                                                require_stability=False).free_energy_bar)
    assert numerical == pytest.approx(np.real(np.vdot(gradient(thermal), dz)), abs=3e-9)
    scales = CellScales(m.catalog.vacuum.delta0_J, m.catalog.vacuum.N0_per_J_m3,
                        m.Tc_K, .9, 1.)
    response = kwt_spatial_response(z, fixed.gradient_cartesian_bar, m.quadrature_mass_bar,
        KWTMobility(scales), np.full(m.cells, .25), np.zeros(m.cells),
        quadrature_to_node=m.quadrature_to_dof)
    assert np.all(response.heat_density_bar >= 0)
    assert response.condensate_heat_rate_bar > 0
    assert abs(response.identity_residual_bar) < 2e-13
    assert response.field_energy_rate_bar == pytest.approx(-response.condensate_heat_rate_bar)


@pytest.mark.parametrize("options", [dict(elements_x=0), dict(elements_y=0), dict(degree=1),
                                     dict(width_m=0), dict(delta_regularizer_bar=0)])
def test_invalid_geometry_and_core_are_rejected(options):
    kwargs = dict(length_m=160e-9, width_m=80e-9, thickness_m=7e-9,
                  elements_x=1, elements_y=1, degree=3)
    kwargs.update(options)
    with pytest.raises(ValueError):
        RectangularSpatialFunctional(SyntheticCatalog(), **kwargs)
