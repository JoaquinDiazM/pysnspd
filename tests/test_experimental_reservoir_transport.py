"""A fixed reservoir reuses native fixed-energy events and their exact ledgers."""
from copy import copy
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.refined_cells import refined_count_catalog
from pysnspd.experimental.reservoir_transport import (
    ReservoirFaceGeometry, ElectronicReservoirTransport)


@pytest.fixture(scope="module")
def cells():
    base = OccupationEnergyCatalog.load(Path(__file__).resolve().parents[1]/
        "docs/implementation/stage1_r2/catalogs/occupation_catalog.npz")
    catalog = refined_count_catalog(base)
    return ElectronicCell(catalog, .43, .02), ElectronicCell(catalog, .91, .06)


def geometry(**changes):
    data = dict(cell_volume_m3=120e-9*7e-9*15e-9, face_area_m2=120e-9*7e-9,
                center_to_reservoir_m=7.5e-9, time_scale_s=1e-12)
    data.update(changes)
    return ReservoirFaceGeometry(**data)


def test_common_energy_FD_equilibrium_with_different_spectra(cells):
    active, bath = cells
    theta = .3
    boundary = ElectronicReservoirTransport(active, bath, theta, geometry())
    p = active.fermi_dirac(theta)
    # Equal count labels are not equal physical energies. Directly subtracting
    # these vectors would produce a false boundary flux at thermal equilibrium.
    assert np.max(abs(p-boundary.reservoir_population)) > .05
    result = boundary.evaluate(p)
    assert 4*np.dot(active.weights, abs(result.cell_rhs)) < 1e-12
    assert abs(result.number_rate_bar_into_cell) < 1e-12
    assert abs(result.energy_rate_bar_into_cell) < 1e-12
    for side, cell in enumerate((bath, active)):
        represented = np.sum(cell.energies[boundary.events.indices[side]]*
                             boundary.events.barycentric[side], axis=1)
        np.testing.assert_allclose(represented, boundary.events.energies, atol=3e-15, rtol=2e-15)
    assert boundary.events.energies.min() > max(active.energies[0], bath.energies[0])
    assert boundary.events.energies.max() < min(active.energies[-1], bath.energies[-1])


def test_rhs_and_opposite_number_energy_ledgers_share_actual_native_moments(cells):
    active, bath = cells
    boundary = ElectronicReservoirTransport(active, bath, .12, geometry())
    initial_bath = boundary.reservoir_population.copy()
    hot = active.fermi_dirac(.4)
    result = boundary.evaluate(hot)
    expected_rhs = boundary.events.rhs((boundary.reservoir_population, hot))[1]
    np.testing.assert_array_equal(result.cell_rhs, expected_rhs)
    number = 4*np.dot(active.weights, result.cell_rhs)
    energy = 4*np.dot(active.weights*active.energies, result.cell_rhs)
    assert result.number_rate_bar_into_cell == number
    assert result.energy_rate_bar_into_cell == energy
    assert energy < 0  # A hot active cell delivers energy to the colder bath.
    assert result.reservoir_energy_rate_bar+result.energy_rate_bar_into_cell == 0
    assert result.reservoir_number_rate_bar+result.number_rate_bar_into_cell == 0
    assert result.cell_power_W+result.reservoir_power_W == 0
    assert result.cell_quasiparticle_rate_per_s+result.reservoir_quasiparticle_rate_per_s == 0
    assert abs(result.event_number_moment_residual) < 2e-14
    assert abs(result.event_energy_moment_residual) < 2e-14
    vacuum = active.catalog.vacuum
    rate_scale = vacuum.N0_per_J_m3*vacuum.delta0_J*boundary.geometry.cell_volume_m3/boundary.geometry.time_scale_s
    assert result.cell_quasiparticle_rate_per_s == pytest.approx(number*rate_scale)
    assert result.cell_power_W == pytest.approx(energy*rate_scale*vacuum.delta0_J)
    np.testing.assert_array_equal(boundary.reservoir_population, initial_bath)
    assert not boundary.reservoir_population.flags.writeable


def test_geometry_changes_density_rhs_but_not_invented_bath_volume(cells):
    active, bath = cells
    first_geometry = geometry()
    first = ElectronicReservoirTransport(active, bath, .12, first_geometry)
    second = ElectronicReservoirTransport(active, bath, .12,
        geometry(cell_volume_m3=2*first_geometry.cell_volume_m3))
    p = active.fermi_dirac(.4)
    a, b = first.evaluate(p), second.evaluate(p)
    np.testing.assert_allclose(b.cell_rhs, a.cell_rhs/2, rtol=3e-15, atol=0)
    assert a.cell_power_W == pytest.approx(b.cell_power_W, rel=3e-15)
    assert a.cell_quasiparticle_rate_per_s == pytest.approx(b.cell_quasiparticle_rate_per_s, rel=3e-15)
    expected = (active.catalog.vacuum.D_m2_s*first_geometry.time_scale_s*first_geometry.face_area_m2/
                (first_geometry.cell_volume_m3*first_geometry.center_to_reservoir_m))
    assert first.coefficient_per_reference_time == expected
    third = ElectronicReservoirTransport(active, bath, .12,
        geometry(time_scale_s=2*first_geometry.time_scale_s))
    c = third.evaluate(p)
    np.testing.assert_allclose(c.cell_rhs, 2*a.cell_rhs, rtol=3e-15, atol=0)
    assert c.cell_power_W == pytest.approx(a.cell_power_W, rel=3e-15)


def test_pauli_faces_remain_inward_and_fixed_vacuum_bath_is_allowed(cells):
    active, bath = cells
    boundary = ElectronicReservoirTransport(active, bath, .2, geometry())
    p = np.full(len(active.energies), .3)
    p[::5] = 0.; p[1::5] = 1.
    saved = p.copy()
    result = boundary.evaluate(p)
    assert np.all(result.cell_rhs[p == 0] >= 0)
    assert np.all(result.cell_rhs[p == 1] <= 0)
    np.testing.assert_array_equal(p, saved)
    zero_bath = ElectronicReservoirTransport(active, bath, 0., geometry())
    equilibrium = zero_bath.evaluate(np.zeros_like(p))
    np.testing.assert_array_equal(equilibrium.cell_rhs, np.zeros_like(p))
    assert equilibrium.cell_power_W == 0


@pytest.mark.parametrize("name", ["cell_volume_m3", "face_area_m2", "center_to_reservoir_m", "time_scale_s"])
def test_geometry_has_no_zero_or_implicit_scale(name):
    for value in (0., -1., np.nan, np.inf):
        with pytest.raises(ValueError, match=name):
            geometry(**{name: value})


@pytest.mark.parametrize("bad", [-1e-16, 1+1e-15, np.nan, np.inf])
def test_population_is_not_clipped(cells, bad):
    active, bath = cells
    boundary = ElectronicReservoirTransport(active, bath, .2, geometry())
    p = active.fermi_dirac(.2); p[0] = bad
    with pytest.raises(ValueError, match="occupation"):
        boundary.evaluate(p)


def test_mismatched_material_units_or_native_quadratures_are_rejected(cells):
    active, bath = cells
    for field, value in (("delta0_J", bath.catalog.vacuum.delta0_J*2),
                         ("N0_per_J_m3", bath.catalog.vacuum.N0_per_J_m3*2),
                         ("D_m2_s", bath.catalog.vacuum.D_m2_s*2)):
        bad = copy(bath)
        metadata = {name: getattr(bath.catalog.vacuum, name)
                    for name in ("delta0_J", "N0_per_J_m3", "D_m2_s")}
        metadata[field] = value
        proxy = SimpleNamespace(vacuum=SimpleNamespace(**metadata), eta=bath.catalog.eta,
            count_nodes=bath.catalog.count_nodes, count_weights=bath.catalog.count_weights)
        object.__setattr__(bad, "catalog", proxy)
        with pytest.raises(ValueError, match=field):
            ElectronicReservoirTransport(active, bad, .2, geometry())
    bad = copy(bath)
    proxy = SimpleNamespace(vacuum=bath.catalog.vacuum, eta=bath.catalog.eta,
        count_nodes=bath.catalog.count_nodes, count_weights=2*bath.catalog.count_weights)
    object.__setattr__(bad, "catalog", proxy)
    with pytest.raises(ValueError, match="quadrature"):
        ElectronicReservoirTransport(active, bad, .2, geometry())
