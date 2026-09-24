import unittest
import numpy as np

from sandbox.stage4_core.coupled_response import energy_rule, kernel, film_port_coefficients
from sandbox.stage4_core.harmonic_uniform_reference import uniform_reference, korzh_scales


class ReciprocalEnergyTailTest(unittest.TestCase):
    def test_tail_measure_preserves_inverse_square_integral_and_finite_rule(self):
        finite_E, finite_w = energy_rule(1.76, .2, .01, 8, 32.)
        energy, weight = energy_rule(1.76, .2, .01, 8, 32., 12)
        np.testing.assert_array_equal(energy[12:-12], finite_E)
        np.testing.assert_array_equal(weight[12:-12], finite_w)
        # The unchanged panel rule can contain equal floating-point nodes when
        # two edge offsets differ below roundoff. Their positive weights remain.
        self.assertTrue(np.all(np.diff(energy)>=0))
        self.assertTrue(np.all(weight>0))
        self.assertAlmostEqual(np.sum(weight[:12]/energy[:12]**2), 1/32., delta=1e-15)
        self.assertAlmostEqual(np.sum(weight[-12:]/energy[-12:]**2), 1/32., delta=1e-15)

    def test_same_integrand_tail_removes_phase_charge_cutoff_defect(self):
        reference=uniform_reference();scales=korzh_scales()
        base=dict(d=reference.gap_bar,t=reference.temperature_ratio,omega=6.,eta=.002*reference.gap_bar,
            order=24,matsubara_count=256,mode=0,**{'lambda':.008508906458670626},
            tau_ee_Tc_ps=.5,tau_ep_Tc_ps=2.47,tD_ps=scales['tD_s']*1e12,
            graph_conductance=.5,sheet_resistance_ohm=608.)
        results=[]
        for cutoff in (32.,64.):
            row=kernel(dict(base,cutoff=cutoff,tail_order=12))
            results.append(np.asarray(row['candidate_response']['real'])+1j*np.asarray(row['candidate_response']['imag']))
            port=film_port_coefficients(dict(base,cutoff=cutoff,tail_order=12))
            self.assertLess(abs(port['neutral']), 2e-9)
        self.assertLess(np.linalg.norm(results[0]-results[1])/np.linalg.norm(results[1]), 2e-7)


if __name__=='__main__':
    unittest.main()
