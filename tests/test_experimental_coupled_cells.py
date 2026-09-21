from pathlib import Path
import importlib.util
import numpy as np
import pytest

ROOT = Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location('stage2_runner', ROOT/'sandbox/stage2_cells/run_coupled.py')
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


@pytest.mark.parametrize('case', ['one', 'two'])
def test_common_native_energy_balance(case):
    system, state, _ = runner.setup(case, phonon_nodes=33)
    assert abs(system.instantaneous_energy_residual(state)) < 1e-11
    assert abs(system.ledgers(system.rhs(0, state))[0][:, 1].sum()-.01) < 1e-14


def test_complete_equilibrium_at_unequal_gaps():
    system, state, _ = runner.setup('two', phonon_nodes=33, heating=0, scenario='equilibrium')
    amplitude, _, _ = system.unpack(state)
    assert abs(amplitude[0]-amplitude[1]) > .01
    derivative = system.rhs(0, state)
    assert np.max(abs(derivative)) < 1e-10


@pytest.mark.parametrize('scenario', ['phonon_vacuum', 'sparse_electrons'])
def test_initial_population_boundary_points_inward(scenario):
    system, state, _ = runner.setup('one', phonon_nodes=33, heating=0, scenario=scenario)
    _, p, n = system.unpack(state)
    rate = system.rhs(0, state)
    dp = rate[1:1+system.electron_size][None, :]
    dn = rate[1+system.electron_size:system.block_size][None, :]
    assert np.all(dp[p == 0] >= 0)
    assert np.all(dp[p == 1] <= 0)
    assert np.all(dn[n == 0] >= 0)


def test_invalid_population_rejected_without_modification():
    system, state, _ = runner.setup('one', phonon_nodes=33)
    state[1] = -1e-18
    before = state.copy()
    with pytest.raises(ValueError):
        system.rhs(0, state)
    np.testing.assert_array_equal(state, before)
