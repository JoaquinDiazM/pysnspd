"""SI graph-port and prescribed-input circuit checks; no detector transient."""
import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from pysnspd.experimental.electrical_ports import (
    ThesisCircuitParameters, ThesisCircuitState, solve_potential,
)


def incidence(nodes, edges):
    matrix = np.zeros((nodes, len(edges)))
    for k, (i, j) in enumerate(edges):
        matrix[i, k] = 1
        matrix[j, k] = -1
    return matrix


def test_uniform_normal_wire_SI_and_passive_sign():
    sigma = 2.5e5
    area = 120e-9*7e-9
    length = 360e-9
    segments = 3
    current = 30e-6
    conductance = np.full(segments, sigma*area/(length/segments))
    solution = solve_potential([[0, 1], [1, 2], [2, 3]], conductance,
        np.zeros(segments), [current, 0., 0., -current], left_node=0, right_node=3)
    resistance = length/(sigma*area)
    assert solution.Vdev_V == pytest.approx(current*resistance, rel=2e-15)
    np.testing.assert_allclose(solution.phi_V, current*resistance*np.linspace(1, 0, 4), rtol=2e-15)
    np.testing.assert_allclose(solution.I_normal_A, current, rtol=2e-15)
    np.testing.assert_allclose(solution.I_total_A, current, rtol=2e-15)
    assert solution.port_power_W == pytest.approx(current**2*resistance, rel=2e-15)
    assert solution.normal_joule_W == pytest.approx(solution.port_power_W, rel=2e-15)
    assert solution.super_power_W == 0
    assert max(abs(solution.current_residual_A)) < 1e-19
    assert abs(solution.power_residual_W) < 1e-20
    assert solution.reference_node == 3 and solution.phi_V[3] == 0


def test_graph_cycles_parallel_edges_and_reference_gauge():
    edges = np.array([[0, 1], [1, 2], [2, 3], [0, 3], [0, 1]])
    conductance = np.array([.2, .5, .3, .7, .6])
    supercurrent = np.array([.003, -.002, .001, .004, -.005])
    exact_phi = np.array([.04, .01, -.02, 0.])
    B = incidence(4, edges)
    expected_total = supercurrent+conductance*(B.T@exact_phi)
    injection = B@expected_total
    first = solve_potential(edges, conductance, supercurrent, injection, left_node=0, right_node=3)
    shifted = solve_potential(edges, conductance, supercurrent, injection,
        left_node=0, right_node=3, reference_node=1, reference_potential_V=2.)
    np.testing.assert_allclose(first.phi_V, exact_phi, atol=3e-17)
    np.testing.assert_allclose(shifted.phi_V, exact_phi-exact_phi[1]+2., atol=1e-15)
    np.testing.assert_allclose(first.I_total_A, expected_total, atol=2e-17)
    np.testing.assert_allclose(first.I_total_A, shifted.I_total_A, atol=2e-17)
    assert first.Vdev_V == pytest.approx(shifted.Vdev_V, abs=1e-16)
    assert first.port_power_W == pytest.approx(shifted.port_power_W, abs=5e-17)
    assert first.port_power_W == pytest.approx(first.face_power_W, abs=1e-17)
    assert first.face_power_W == pytest.approx(first.normal_joule_W+first.super_power_W, abs=1e-17)


def test_reversing_edge_orientation_is_not_reversing_the_physical_current():
    edges = np.array([[0, 1], [1, 2]])
    original = solve_potential(edges, [2., 3.], [.4, -.2], [.3, 0., -.3], left_node=0, right_node=2)
    reversed_ = solve_potential(edges[:, ::-1], [2., 3.], [-.4, .2], [.3, 0., -.3], left_node=0, right_node=2)
    np.testing.assert_allclose(reversed_.phi_V, original.phi_V, atol=1e-16)
    np.testing.assert_allclose(reversed_.I_total_A, -original.I_total_A, atol=1e-16)
    assert reversed_.Vdev_V == pytest.approx(original.Vdev_V)
    assert reversed_.port_power_W == pytest.approx(original.port_power_W)


def test_supercurrent_work_can_be_negative_without_absolute_value_or_repair():
    result = solve_potential([[0, 1]], [.01], [60e-6], [30e-6, -30e-6], left_node=0, right_node=1)
    assert result.Vdev_V == pytest.approx(-.003)
    assert result.I_normal_A[0] == pytest.approx(-30e-6)
    assert result.I_total_A[0] == pytest.approx(30e-6)
    assert result.normal_joule_W > 0
    assert result.super_power_W < 0 and result.port_power_W < 0
    assert abs(result.power_residual_W) < 1e-21
    assert abs(result.decomposition_residual_W) < 1e-21


@pytest.mark.parametrize("changes,match", [
    ({"edges": [[0, 0]]}, "self edges"),
    ({"edges": [[0, 2]]}, "outside"),
    ({"edges": [[0., 1.]]}, "integer"),
    ({"conductance_S": [0.]}, "strictly positive"),
    ({"conductance_S": [-1.]}, "strictly positive"),
    ({"conductance_S": [np.inf]}, "finite"),
    ({"supercurrent_A": [np.nan]}, "finite"),
    ({"supercurrent_A": [1j]}, "real"),
    ({"supercurrent_A": [0., 0.]}, "per edge"),
    ({"injection_A": [1e-9, -1.001e-9]}, "do not balance"),
    ({"right_node": 0}, "distinct"),
    ({"reference_node": True}, "integer"),
    ({"reference_potential_V": np.nan}, "finite"),
])
def test_graph_rejects_invalid_contract(changes, match):
    args = dict(edges=[[0, 1]], conductance_S=[1.], supercurrent_A=[0.],
                injection_A=[1e-9, -1e-9], left_node=0, right_node=1)
    args.update(changes)
    with pytest.raises(ValueError, match=match):
        solve_potential(**args)


def test_disconnected_graph_is_rejected_before_solving():
    with pytest.raises(ValueError, match="connected"):
        solve_potential([[0, 1], [2, 3]], [1., 1.], [0., 0.], [1., -1., 0., 0.], left_node=0, right_node=3)


def test_external_inductance_is_explicit_and_passive_parameters_positive():
    with pytest.raises(TypeError):
        ThesisCircuitParameters()
    for name in ("Lk_ext_H", "R_bias_ohm", "L_bias_H", "R_load_ohm", "C_couple_F"):
        for invalid in (0., -1., np.inf):
            args = dict(Lk_ext_H=7e-9)
            args[name] = invalid
            with pytest.raises(ValueError):
                ThesisCircuitParameters(**args)
    with pytest.raises(ValueError):
        ThesisCircuitParameters(7e-9, V_bias_V=np.nan)
    assert ThesisCircuitParameters(7e-9).Lk_ext_H == 7e-9


def test_CM4_outputs_energy_and_CM8_for_both_voltage_signs():
    circuit = ThesisCircuitParameters(Lk_ext_H=7e-9)
    state = ThesisCircuitState(30e-6, 24e-6, 80e-6)
    for voltage in (-.002, 0., .003):
        output = circuit.outputs(state, voltage)
        assert output.I_RF_A == pytest.approx(6e-6)
        assert output.Vout_V == pytest.approx(300e-6)
        assert output.Vd_V == pytest.approx(380e-6)
        assert output.Vd_V != output.Vout_V
        np.testing.assert_allclose(circuit.rhs(state, voltage),
            [-380., (380e-6-voltage)/7e-9, 6e-6/100e-12], rtol=2e-13)
        expected = .5*(1e-6*(30e-6)**2+7e-9*(24e-6)**2+100e-12*(80e-6)**2)
        assert circuit.energy(state) == pytest.approx(expected, rel=2e-15)
        balance = circuit.power_balance(state, voltage)
        assert balance.device_power_W == pytest.approx(24e-6*voltage)
        assert balance.source_power_W == pytest.approx(.3*30e-6)
        assert abs(balance.residual_W) < 5e-21
        # A deliberately incorrect derivative cannot earn a zero power residual.
        wrong = circuit.power_balance(state, voltage, derivative=np.zeros(3))
        assert abs(wrong.residual_W) > 1e-10


@pytest.mark.parametrize("voltage", [-.001, 0., .001])
def test_prescribed_voltage_stationary_oracle(voltage):
    circuit = ThesisCircuitParameters(Lk_ext_H=7e-9)
    state = circuit.stationary_for_prescribed_voltage(voltage)
    assert state[0] == state[1]
    assert state[2] == voltage
    np.testing.assert_allclose(circuit.rhs(state, voltage), 0., atol=1e-10)
    assert circuit.outputs(state, voltage).Vout_V == 0


@pytest.mark.parametrize("resistance", [0., 2., 1000.])
def test_prescribed_resistance_stationary_oracle(resistance):
    circuit = ThesisCircuitParameters(Lk_ext_H=7e-9)
    state = circuit.stationary_for_prescribed_resistance(resistance)
    assert state[0] == pytest.approx(.3/(1e4+resistance))
    np.testing.assert_allclose(circuit.rhs(state, resistance*state[1]), 0., atol=1e-9)
    assert state[2] == pytest.approx(resistance*state[1])
    with pytest.raises(ValueError):
        circuit.stationary_for_prescribed_resistance(-1.)


@pytest.mark.parametrize("resistance", [None, 200.])
def test_prescribed_circuit_transient_matches_independent_matrix_exponential(resistance):
    circuit = ThesisCircuitParameters(Lk_ext_H=7e-9)
    # Matrix assembled independently from the three circuit connections. This
    # test evolves only three circuit variables with a prescribed port law.
    Rb, Lb, Rl, C, Le = 1e4, 1e-6, 50., 100e-12, 7e-9
    A = np.array([[-(Rb+Rl)/Lb, Rl/Lb, -1/Lb],
                  [Rl/Le, -(Rl+(resistance or 0.))/Le, 1/Le],
                  [1/C, -1/C, 0.]])
    voltage = .002 if resistance is None else 0.
    forcing = np.array([.3/Lb, -voltage/Le, 0.])
    augmented = np.zeros((4, 4))
    augmented[:3, :3] = A
    augmented[:3, 3] = forcing
    initial = np.array([30e-6, 27e-6, 50e-6])
    times = np.linspace(0., 2e-9, 9)
    reference = np.array([(expm(augmented*t)@np.r_[initial, 1.])[:3] for t in times])
    rhs = lambda t, y: circuit.rhs(y, voltage if resistance is None else resistance*y[1])
    measured = solve_ivp(rhs, (times[0], times[-1]), initial, t_eval=times,
                         method="DOP853", rtol=2e-11, atol=1e-15)
    assert measured.success
    np.testing.assert_allclose(measured.y.T, reference, rtol=2e-9, atol=3e-14)
    for state in reference:
        vdev = voltage if resistance is None else resistance*state[1]
        assert abs(circuit.power_balance(state, vdev).residual_W) < 5e-20


def test_instantaneous_graph_and_circuit_port_work_cancel():
    circuit = ThesisCircuitParameters(Lk_ext_H=7e-9)
    state = np.array([31e-6, 23e-6, 1e-4])
    current = state[1]
    port = solve_potential([[0, 1], [1, 2]], [.03, .02], [20e-6, 19e-6],
        [current, 0., -current], left_node=0, right_node=2)
    balance = circuit.power_balance(state, port.Vdev_V)
    assert port.port_power_W == pytest.approx(balance.device_power_W, abs=1e-22)
    assert abs(balance.stored_energy_rate_W+port.port_power_W+balance.bias_joule_W+
               balance.load_joule_W-balance.source_power_W) < 5e-21
    # This algebraic interface test makes no time-coupled detector claim.


@pytest.mark.parametrize("state", [[1., 2.], [1., 2., np.nan], [1j, 2., 3.]])
def test_circuit_rejects_bad_state(state):
    circuit = ThesisCircuitParameters(Lk_ext_H=7e-9)
    with pytest.raises(ValueError):
        circuit.rhs(state, 0.)
    with pytest.raises(ValueError):
        circuit.energy(state)
