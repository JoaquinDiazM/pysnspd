"""Temporal comparison must not hide a common error behind baseline subtraction."""
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'sandbox/stage4_core'))
from thermal_nonlinear_time import compare_exact
from pysnspd.experimental.thermal_spatial_usadel import ThermalGraph


class TestNonlinearTimeComparison(unittest.TestCase):
    def test_common_baseline_error_is_detected_when_differences_agree(self):
        graph=ThermalGraph(np.ones(3),np.array([[0,1],[1,2]]),np.ones(2),
                           np.array([[0.,0.],[1.,0.],[2.,0.]]),np.array([0,2]))
        shape=np.array([0.,1.,0.]);d0=np.ones(3,dtype=complex)
        def fields(common):
            values={'d0':d0,'G0':np.zeros(3,dtype=complex)}
            for name,displacement in (('baseline',common*shape),
                ('amplitude',(common+.1)*shape),('angular_phase',(common+.1j)*shape)):
                values[name+'_displacement']=displacement
                values[name+'_gap_gradient_increment']=displacement.copy()
                values[name+'_current_increment']=np.array([displacement[1].real,-displacement[1].real])
            for name in ('amplitude','angular_phase'):
                values[name+'_difference']=values[name+'_displacement']-values['baseline_displacement']
            return values
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            for name,common in (('initial',0.),('coarse',.02),('fine',0.)):
                np.savez(root/(name+'.npz'),**fields(common))
            first={'observations':[{'time_ps':0.,'fields_path':'initial.npz'},
                                   {'time_ps':1.,'fields_path':'coarse.npz'}]}
            second={'observations':[{'time_ps':0.,'fields_path':'initial.npz'},
                                    {'time_ps':1.,'fields_path':'fine.npz'}]}
            result=compare_exact({'refinement_absolute_norm':1e-10,
                'refinement_relative_limit':.01},graph,root,first,second)
        baseline=next(row for row in result['records'] if row['time_ps']==1.
            and row['probe']=='baseline' and row['observable']=='displacement')
        self.assertFalse(baseline['admitted'])
        self.assertAlmostEqual(baseline['relative_to_initial'],.2)
        self.assertFalse(result['all_comparisons_met'])
        for row in result['records']:
            if row['probe']!='baseline' and row['observable']=='displacement':
                self.assertTrue(row['admitted'])


if __name__=='__main__':unittest.main()
