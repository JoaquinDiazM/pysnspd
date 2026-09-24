import unittest
import json
import numpy as np
from scipy.linalg import expm
from pysnspd.experimental.affine_krylov import affine_action,advance,exponential_step


class AffineKrylovTests(unittest.TestCase):
    def test_affine_nonnormal_system_matches_dense_exponential(self):
        A=np.array([[-4.,3.,0.],[0.,-1.,2.],[0.,0.,-.3]])
        b=np.array([.1,-.2,.3]);x=np.array([.001,0.,-.002,1.])
        matrix=np.zeros((4,4));matrix[:3,:3]=A;matrix[:3,3]=b
        result,history=advance(affine_action(lambda v:A@v,b),x,.7,rtol=1e-4,atol=1e-8,max_dimension=8)
        np.testing.assert_allclose(result,expm(.7*matrix)@x,atol=1e-13)
        self.assertTrue(all(r['accepted'] for r in history));self.assertAlmostEqual(result[-1],1.)

    def test_full_space_defect_matches_arnoldi_expression(self):
        rng=np.random.default_rng(4);A=rng.normal(size=(30,30))-.5*np.eye(30);x=rng.normal(size=30)
        result=exponential_step(lambda v:A@v,x,.2,max_dimension=8,rtol=1e-8,atol=1e-10)
        self.assertAlmostEqual(np.linalg.norm(A@result.state-result.derivative),result.endpoint_defect,delta=2e-13)
        self.assertLess(np.linalg.norm(result.state-expm(.2*A)@x),2*result.integrated_defect)

    def test_restart_controls_stiff_full_space_evolution(self):
        rates=np.linspace(1,80,30);x=np.ones(30)
        result,history=advance(lambda v:-rates*v,x,.4,max_dimension=8,rtol=1e-3,atol=1e-7)
        self.assertTrue(any(not r['accepted'] for r in history))
        json.dumps(history)  # Progress records must be JSON-native, including Boolean decisions.
        np.testing.assert_allclose(result,np.exp(-rates*.4),atol=3e-4)

    def test_baseline_subtraction_preserves_affine_forcing(self):
        A=np.diag([-1.,-2.]);b=np.array([2.,-1.]);dx=np.array([.001,-.001])
        base,_=advance(affine_action(lambda v:A@v,b),np.array([0.,0.,1.]),.5,max_dimension=4)
        delta,_=advance(lambda v:A@v,dx,.5,max_dimension=4)
        pert,_=advance(affine_action(lambda v:A@v,b),np.r_[dx,1.],.5,max_dimension=4)
        np.testing.assert_allclose(base[:2]+delta,pert[:2],atol=1e-14)

    def test_affine_constant_is_excluded_from_accuracy_scale(self):
        result=exponential_step(affine_action(lambda v:-v,np.array([1e-5])),np.array([0.,1.]),.1,
            rtol=1e-3,atol=1e-9,max_dimension=2,active_size=1)
        self.assertLess(result.tolerance,2e-9)


if __name__=='__main__':unittest.main()
