from pathlib import Path
import numpy as np
import pytest
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.cell_transport import NativeTransportEvents


@pytest.fixture(scope="module")
def cells():
    path = Path(__file__).parents[1]/"docs/implementation/stage1_r2/catalogs/occupation_catalog.npz"
    catalog = OccupationEnergyCatalog.load(path)
    return ElectronicCell(catalog, .43, .02), ElectronicCell(catalog, 1.2, .2)


def test_events_conserve_native_energy_and_number(cells):
    event = NativeTransportEvents(*cells)
    p = np.asarray([c.fermi_dirac(t) for c, t in zip(cells, (.4, .12))])
    dp = event.rhs(p)
    assert abs(sum(4*np.dot(c.weights, v) for c, v in zip(cells, dp))) < 1e-13
    assert abs(sum(4*np.dot(c.weights*c.energies, v) for c, v in zip(cells, dp))) < 1e-13
    assert np.max(abs(event.event_energy_residual())) < 1e-13
    assert event.transferred_power(p) > 0


def test_thermal_equilibrium_with_different_gaps(cells):
    event = NativeTransportEvents(*cells)
    p = np.asarray([c.fermi_dirac(.3) for c in cells])
    assert np.max(abs(event.rhs(p))) < 1e-11


def test_pauli_faces_and_entropy(cells):
    event = NativeTransportEvents(*cells)
    random = np.random.default_rng(744)
    p = random.uniform(.05, .9, (2, len(cells[0].energies)))
    p[:, ::7] = 0; p[:, 1::7] = 1
    dp = event.rhs(p)
    assert np.all(dp[p == 0] >= 0)
    assert np.all(dp[p == 1] <= 0)
    q = random.uniform(.01, .99, p.shape)
    dq = event.rhs(q)
    entropy_rate = sum(-4*np.dot(c.weights*np.log(f/(1-f)), v)
                       for c, f, v in zip(cells, q, dq))
    assert entropy_rate >= 0


def test_private_states_no_extrapolation(cells):
    event = NativeTransportEvents(*cells)
    for side, cell in enumerate(cells):
        represented = np.sum(cell.energies[event.indices[side]]*event.barycentric[side], axis=1)
        np.testing.assert_allclose(represented, event.energies, atol=1e-14)
    with pytest.raises(ValueError):
        NativeTransportEvents(*cells, diffusion_over_length_squared=-1)
