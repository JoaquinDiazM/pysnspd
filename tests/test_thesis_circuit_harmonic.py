import unittest
import numpy as np

from pysnspd.experimental.electrical_ports import ThesisCircuitParameters
from pysnspd.experimental.energy_catalog import HBAR_J_S, E_CHARGE_C
from pysnspd.experimental.thesis_circuit_harmonic import partition_inductance, solve_harmonic


class ThesisCircuitHarmonicTest(unittest.TestCase):
    def test_full_circuit_frequency_response_matches_existing_time_rhs_and_power(self):
        partition = partition_inductance(10e-9, .2e-9*(2*E_CHARGE_C)/HBAR_J_S)
        self.assertAlmostEqual(partition.external_H, 9.8e-9, delta=1e-22)
        parameters = ThesisCircuitParameters(partition.external_H, V_bias_V=0.)
        for omega in (1e6, 1e9, 1e12):
            impedance = 20.+1j*omega*.2e-9
            result = solve_harmonic(parameters, omega, impedance, bias_voltage_peak_V=.001+.0003j)
            y = result.state_peak
            actual_rhs = (parameters.rhs(y.real, result.device_voltage_peak_V.real)
                          +1j*parameters.rhs(y.imag, result.device_voltage_peak_V.imag))
            actual_rhs[0] += result.bias_voltage_peak_V/parameters.L_bias_H
            np.testing.assert_allclose(actual_rhs, 1j*omega*y, rtol=1e-11, atol=1e-12)
            self.assertLess(abs(result.complex_power_residual_VA), abs(result.source_complex_power_VA)*1e-12)
            self.assertLess(abs(result.average_power_residual_W), abs(result.source_complex_power_VA)*1e-12)
            # Independent phase sampling: stored-energy rate averages to zero,
            # and source/device/resistor powers recover the phasor identity.
            angles = np.arange(128)*2*np.pi/128
            wave = np.exp(1j*angles)
            state = np.real(y[:, None]*wave)
            derivative = np.real((1j*omega*y)[:, None]*wave)
            drive = np.real(result.bias_voltage_peak_V*wave)
            voltage = np.real(result.device_voltage_peak_V*wave)
            ib, detector, vc = state
            storage = (parameters.L_bias_H*ib*derivative[0]
                       +parameters.Lk_ext_H*detector*derivative[1]
                       +parameters.C_couple_F*vc*derivative[2])
            residual = (storage+voltage*detector+parameters.R_bias_ohm*ib**2
                        +parameters.R_load_ohm*(ib-detector)**2-drive*ib)
            self.assertLess(np.max(abs(residual)), abs(result.source_complex_power_VA)*1e-10)
            self.assertLess(abs(np.mean(storage)), abs(result.source_complex_power_VA)*1e-12)

    def test_invalid_partition_is_not_clipped_and_active_port_is_not_relabelled_passive(self):
        with self.assertRaises(ValueError):
            partition_inductance(10e-9, 12e-9*(2*E_CHARGE_C)/HBAR_J_S)
        p = ThesisCircuitParameters(9e-9)
        result = solve_harmonic(p, 1e9, -10.+20j)
        self.assertLess(result.device_complex_power_VA.real, 0.)
        self.assertLess(abs(result.complex_power_residual_VA), abs(result.source_complex_power_VA)*1e-12)
        with self.assertRaises(ValueError):
            solve_harmonic(p, 0., 10.)

    def test_opposite_fourier_conventions_require_opposite_impedance(self):
        p = ThesisCircuitParameters(9e-9)
        omega = 1e10
        positive = solve_harmonic(p, omega, 13.+20j, bias_voltage_peak_V=.001+.002j)
        negative = solve_harmonic(p, omega, 13.-20j, bias_voltage_peak_V=.001-.002j,
                                  time_convention='exp(-iwt)')
        np.testing.assert_allclose(negative.state_peak, positive.state_peak.conj(), rtol=1e-14)
        self.assertAlmostEqual(negative.device_complex_power_VA.real,
                               positive.device_complex_power_VA.real)
        self.assertAlmostEqual(negative.external_reactive_power_var,
                               -positive.external_reactive_power_var)
        with self.assertRaises(ValueError):
            solve_harmonic(p, omega, 13., time_convention='automatic')


if __name__ == '__main__':
    unittest.main()
