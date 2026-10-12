"""Geometric/phase diagnostics only: no physical evolution or model admission."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from sandbox.stage5_effective_validation import vortex_diagnostics as vd

class VortexGeometryTests(unittest.TestCase):
    def setUp(self):
        self.xy=np.array([[0.,0.],[1.,0.],[0.,1.]])
        self.tri=np.array([[0,1,2]])
        self.center=np.array([.25,.25])
        self.z=(self.xy[:,0]-.25)+1j*(self.xy[:,1]-.25)

    def test_ccw_sign_and_orientation(self):
        q,m,rows,_=vd.triangle_topology(self.z,self.xy,self.tri)
        self.assertEqual(q.tolist(),[1]);self.assertEqual(m.tolist(),[0.])
        np.testing.assert_allclose(rows[0]["linear_complex_zero_nm"],self.center*1e9)
        reverse,count=vd.ccw_triangles(self.xy,self.tri[:,[0,2,1]])
        self.assertEqual(count,1);np.testing.assert_array_equal(reverse,self.tri)
        q,*_=vd.triangle_topology(np.conj(self.z),self.xy,self.tri)
        self.assertEqual(q.tolist(),[-1])

    def test_edge_minimum_is_not_nodal_minimum(self):
        z=np.array([1+1j,1-1j,3+0j])
        minimum,inside,_=vd.complex_triangle_minimum(z,self.tri)
        np.testing.assert_allclose(minimum,[1.])
        self.assertFalse(inside[0]);self.assertGreater(abs(z).min(),minimum[0])

    def test_interpolated_edge_zero(self):
        z=np.array([0+1j,0-1j,2+1j])
        minimum,_,_=vd.complex_triangle_minimum(z,self.tri)
        np.testing.assert_allclose(minimum,[0.])

    def test_nodal_zero_phase_is_undefined(self):
        with self.assertRaisesRegex(ValueError,"undefined"):
            vd.triangle_topology(np.array([0.,1j,-1.]),self.xy,self.tri)

    def test_boundary_charge_agrees(self):
        q,*_=vd.triangle_topology(self.z,self.xy,self.tri)
        self.assertEqual(int(q.sum()),vd.boundary_charge(self.z,self.tri)[0])

    def test_clipped_weights_preserve_physical_area(self):
        xy=np.array([[-80,-40],[80,-40],[80,40],[-80,40]])*1e-9
        tri=np.array([[0,1,2],[0,2,3]])
        weights=vd.central_linear_weights(xy,tri,[0,0])
        self.assertAlmostEqual(weights.sum()*1e18,9600.,places=10)
        np.testing.assert_allclose(weights@np.ones(4),9600e-18,rtol=1e-14,atol=0)
        shifted=vd.central_linear_weights(xy+[123e-9,0],tri,[123e-9,0])
        np.testing.assert_allclose(shifted,weights,rtol=1e-13,atol=0)

    def test_L2_uses_physical_position_and_complex_interpolant(self):
        xy=np.array([[-80,-40],[80,-40],[80,40],[-80,40]])*1e-9
        tri=np.array([[0,1,2],[0,2,3]])
        z=np.array([[1+1j,2+1j,2+2j,1+2j]])
        a=dict(coordinates_m=xy,triangles=tri,delta_J=z)
        b=dict(coordinates_m=xy+[100e-9,0],triangles=tri,delta_J=z*np.exp(1j*.45))
        same=vd.common_field_l2(a,b,0,0,[0,0],[100e-9,0],1.)
        self.assertLess(same["relative_L2"],1e-14)
        # Identical nodal magnitudes can have different interpolated magnitudes.
        a["delta_J"]=np.array([[1+0j,0+1j,-1+0j,0-1j]])
        b["coordinates_m"]=xy
        b["delta_J"]=np.ones((1,4),complex)
        different=vd.common_field_l2(a,b,0,0,[0,0],[0,0],1.)
        self.assertGreater(different["relative_L2"],.1)

class PhaseSamplingTests(unittest.TestCase):
    def test_certified_small_increments(self):
        exact=np.array([0.,.2,.5,.9])
        phase,info=vd.phase_unwrap_diagnostic(np.angle(np.exp(1j*exact)),exact*vd.PHI0/(2*np.pi))
        self.assertTrue(info["certified"]);np.testing.assert_allclose(phase,exact)

    def test_branch_cut_is_unwrapped_when_sampling_is_certified(self):
        exact=np.array([0.,2.,4.,6.])
        phase,info=vd.phase_unwrap_diagnostic(np.angle(np.exp(1j*exact)),exact*vd.PHI0/(2*np.pi))
        self.assertTrue(info["certified"]);np.testing.assert_allclose(phase,exact,atol=1e-15)

    def test_hidden_full_turn_is_explicitly_ambiguous(self):
        phase,info=vd.phase_unwrap_diagnostic(np.array([0.,0.]),np.array([0.,vd.PHI0]))
        self.assertFalse(info["certified"]);self.assertEqual(info["ambiguous_interval_indices"],[0])
        self.assertAlmostEqual(float(phase[-1]),0.) # No turn repaired/invented from CM.

    def test_pi_increment_is_not_certified(self):
        _,info=vd.phase_unwrap_diagnostic(np.array([0.,np.pi]),np.array([0.,vd.PHI0/2]))
        self.assertFalse(info["certified"])

    def test_nonfinite_phase_rejected(self):
        with self.assertRaises(ValueError):
            vd.phase_unwrap_diagnostic(np.array([0.,np.nan]),np.array([0.,vd.PHI0]))

class ReportGateTests(unittest.TestCase):
    def test_optional_analysis_cannot_bypass_global_gate(self):
        root=Path(__file__).resolve().parents[1]
        with tempfile.TemporaryDirectory(dir="/home/jdiaz/scratch/tmp/pysnspd") as temp:
            temp=Path(temp);admission=temp/"no_admission.json";admission.write_text('{"stage5":false}')
            output=temp/"report"
            proc=subprocess.run([sys.executable,str(root/"sandbox/stage5_effective_validation/build_final_report.py"),
                "--admission",str(admission),"--build-root",str(output),"--prepare-only"],
                capture_output=True,text=True,timeout=15)
            self.assertNotEqual(proc.returncode,0)
            self.assertIn("stage5: true",proc.stderr)
            self.assertFalse(output.exists())

if __name__=="__main__":unittest.main()
