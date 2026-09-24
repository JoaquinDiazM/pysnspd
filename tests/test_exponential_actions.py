"""Independent ODE references for physical-coordinate exponential actions."""
import dataclasses
import json
import unittest
import numpy as np
from scipy.integrate import solve_ivp
from pysnspd.experimental.exponential_actions import phi_action,frozen_affine_step,etd2_trial


class TestExponentialActions(unittest.TestCase):
    def test_phi_actions_match_independent_polynomially_forced_ode(self):
        A=np.array([[-3.,20.,0.],[0.,-1.,.5],[0.,0.,0.]])
        b=np.array([.2,-.3,.1]);h=.7
        for order in (1,2):
            reference=solve_ivp(lambda t,y:A@y+(t/h)**(order-1)*b,
                (0,h),np.zeros(3),rtol=2e-12,atol=2e-14,method='DOP853').y[:,-1]
            result=phi_action(lambda x:A@x,b,h,order=order,rtol=1e-10,atol=1e-12)
            np.testing.assert_allclose(result.value,reference,rtol=2e-11,atol=2e-12)
            self.assertTrue(result.accepted)

    def test_zero_eigenvalue_and_constant_forcing_have_no_auxiliary_drift(self):
        x=np.array([.3,-.5]);b=np.array([.07,-.2]);h=.17
        dimensions=[]
        def zero(v):dimensions.append(len(v));return np.zeros_like(v)
        for _ in range(20):x,diagnostic=frozen_affine_step(zero,x,b,h,rtol=1e-8)
        np.testing.assert_allclose(x,np.array([.3,-.5])+20*h*b,atol=2e-15,rtol=0)
        self.assertEqual(set(dimensions),{2})
        second=phi_action(zero,b,h,order=2)
        np.testing.assert_allclose(second.value,h*b/2,atol=1e-16)
        json.dumps({key:value for key,value in dataclasses.asdict(diagnostic).items() if key!='value'})

    def test_nonstationary_affine_motion_matches_direct_ode(self):
        A=np.array([[-20.,5.],[0.,-.4]]);b=np.array([.03,.7]);x=np.array([.2,-.1]);h=.4
        result,diagnostic=frozen_affine_step(lambda v:A@v,x,A@x+b,h,rtol=1e-10,atol=1e-12)
        reference=solve_ivp(lambda t,y:A@y+b,(0,h),x,method='DOP853',rtol=1e-12,atol=1e-14).y[:,-1]
        np.testing.assert_allclose(result,reference,atol=2e-12,rtol=2e-11)
        self.assertTrue(diagnostic.accepted)

    def test_insufficient_space_is_reported_without_changing_method(self):
        spectrum=-np.geomspace(.01,100.,24);b=np.ones(24)
        result=phi_action(lambda x:spectrum*x,b,1.,max_dimension=2,rtol=1e-8,atol=1e-12)
        self.assertFalse(result.accepted)
        self.assertGreater(result.integrated_defect,result.tolerance)
        self.assertEqual(result.action_count,2)

    def test_invalid_and_zero_vectors(self):
        with self.assertRaises(ValueError):phi_action(lambda x:x,np.array([1j]),1.)
        with self.assertRaises(ValueError):phi_action(lambda x:x,np.ones(3),0.)
        with self.assertRaises(RuntimeError):phi_action(lambda x:np.array([np.nan]),np.ones(1),1.)
        result=phi_action(lambda x:(_ for _ in ()).throw(AssertionError()),np.zeros(3),1.)
        self.assertEqual(result.action_count,0)
        self.assertTrue(result.accepted)

    def test_etd2_has_second_order_for_nonlinear_logistic_solution(self):
        errors=[]
        exact=1/(1+4*np.exp(-1.))
        for count in (10,20,40):
            x=np.array([.2]);h=1/count
            for _ in range(count):
                x,diagnostic=etd2_trial(lambda v:-3*v,lambda y:y*(1-y),x,h,
                    rtol=1e-10,atol=1e-12)
                self.assertTrue(diagnostic['accepted'])
            errors.append(abs(x[0]-exact))
        self.assertTrue(3.5<errors[0]/errors[1]<4.5)
        self.assertTrue(3.5<errors[1]/errors[2]<4.5)

    def test_etd2_preserves_affine_rhs_without_nonlinear_correction(self):
        A=np.array([[-20.,3.],[0.,0.]]);b=np.array([.1,.2]);x=np.array([.2,.3])
        result,diagnostic=etd2_trial(lambda v:A@v,lambda y:A@y+b,x,.3,
            rtol=1e-10,atol=1e-12)
        reference=solve_ivp(lambda t,y:A@y+b,(0,.3),x,rtol=1e-12,atol=1e-14,method='DOP853').y[:,-1]
        np.testing.assert_allclose(result,reference,atol=2e-12,rtol=1e-11)
        self.assertLess(np.linalg.norm(diagnostic['correction']),1e-15)

    def test_complex_operator_is_rejected_in_real_coordinates(self):
        with self.assertRaises(RuntimeError):phi_action(lambda x:x+1j,np.ones(2),1.)


if __name__=='__main__':unittest.main()
