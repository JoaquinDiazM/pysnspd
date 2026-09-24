import unittest
import numpy as np
from scipy.optimize import brentq
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental import retarded_spatial_usadel as retarded
from pysnspd.experimental import frozen_kinetic_usadel as kinetic
from pysnspd.experimental.thermal_weak_response import SpectralTangent, ThermalKWTNormal
from pysnspd.experimental.longitudinal_reciprocal import LongitudinalReciprocalBridge


class LongitudinalReciprocalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.graph = thermal.rectangular_graph(5, 5, 3., 3.)
        cls.t = .9/8.65; cls.eps = 2*np.pi*cls.t*(np.arange(12)+.5)
        amplitude = brentq(lambda a:np.log(cls.t)+2*np.pi*cls.t*np.sum(1/cls.eps-1/np.sqrt(cls.eps**2+a*a)), .1, 3.)
        cls.d = np.full(cls.graph.n_nodes, amplitude, complex)
        cls.solutions = [thermal.solve_frequency(cls.graph, cls.d, e,
            fixed_nodes=cls.graph.boundary_nodes, fixed_u=cls.d[cls.graph.boundary_nodes]/e,
            tol=1e-11) for e in cls.eps]
        cls.tangents = [SpectralTangent(cls.graph, cls.d, s) for s in cls.solutions]
        cls.kw = ThermalKWTNormal(cls.graph, cls.d, Tc_K=8.65, T_K=.9)
        cls.energies = np.array([.2, .8, 1.8, 2.4]); cls.weights = np.array([.3, .5, .6, .5])
        cls.operators = []
        for e in cls.energies:
            a, b = retarded.uniform_fields(cls.d, .03-1j*e)
            f, ft, g = retarded.fields(a, b)
            cls.operators.append(kinetic.assemble(cls.graph, cls.d, g, f, ft))

    def bridge(self, operators=None):
        def hessian(x):
            result=2*self.graph.area_weights*np.log(self.t)*x
            for e, tangent in zip(self.eps, self.tangents):
                result += 4*np.pi*self.t*self.graph.area_weights*(x/e-tangent.apply(x).df.real)
            return result
        return LongitudinalReciprocalBridge(self.graph, self.operators if operators is None else operators,
            self.energies, self.weights, self.t, hessian, lambda g:self.kw.local_inverse(g).real,
            spectral_z=.03-1j*self.energies)

    def state(self):
        rng=np.random.default_rng(23)
        x=.002*rng.normal(size=self.graph.n_nodes)
        y=.003*rng.normal(size=(len(self.energies), self.graph.n_nodes))
        x[self.graph.boundary_nodes]=0.; y[:, self.graph.boundary_nodes]=0.
        return x,y

    def test_reference_is_equilibrium_and_charge_decouples_on_same_physical_kernels(self):
        base=thermal.evaluate_thermal(self.graph,self.d,self.t,self.solutions)
        self.assertLess(np.max(abs(base.gap_gradient)),1e-11)
        bridge=self.bridge();x,y=self.state()
        for i,op in enumerate(self.operators):
            np.testing.assert_allclose(bridge.kernel[i],-2*self.graph.area_weights*op.R[:,0,1].imag,atol=1e-14)
            np.testing.assert_array_equal(op.matrix[1::2,::2].toarray(),0.)
        self.assertEqual(bridge.response(x,y).charge_residual_max,0.)

    def test_reciprocal_work_and_dirichlet_losses_close_quadratic_availability(self):
        bridge=self.bridge();x,y=self.state();r=bridge.response(x,y)
        self.assertGreater(r.kwt_loss,0.);self.assertGreater(r.longitudinal_loss,0.)
        self.assertLess(abs(r.reciprocal_power_residual),1e-16)
        step=1e-6
        difference=(bridge.availability(x+step*r.gap_velocity,y+step*r.population_velocity)
            -bridge.availability(x-step*r.gap_velocity,y-step*r.population_velocity))/(2*step)
        self.assertAlmostEqual(difference,r.availability_rate,delta=1e-13)
        omitted_work=np.sum((bridge.weights*bridge.chi)[:,None]*y*bridge.kernel*r.gap_velocity[None,:])
        self.assertGreater(abs(omitted_work),1e-12)
        self.assertGreater(abs(r.availability_rate-omitted_work+r.kwt_loss+r.longitudinal_loss),1e-12)

    def test_normal_limit_recovers_diffusion_and_has_no_gap_work_kernel(self):
        n=self.graph.n_nodes;zero=np.zeros(n);one=np.ones(n)
        op=kinetic.assemble(self.graph,zero,one,zero,zero)
        bridge=self.bridge([op]*len(self.energies));x,y=self.state()
        r=bridge.response(np.zeros_like(x),y)
        np.testing.assert_array_equal(bridge.kernel,0.)
        tail,head=self.graph.edges.T
        expected=np.zeros_like(y)
        for row,result in zip(y,expected):
            flux=self.graph.conductance*(row[head]-row[tail])
            np.add.at(result,tail,flux);np.add.at(result,head,-flux)
        free=bridge.free
        np.testing.assert_allclose(r.population_velocity[:,free],expected[:,free]/self.graph.area_weights[free],atol=1e-15)

    def test_zero_dos_remains_an_algebraic_mass_without_a_floor(self):
        n=self.graph.n_nodes;zero=np.zeros(n);one=np.ones(n)
        # Exact causal subgap BCS values at E=0 and eta=0: g=0, f=ft=1.
        op=kinetic.assemble(self.graph,self.d,zero,one,one)
        template=self.bridge()
        bridge=LongitudinalReciprocalBridge(self.graph,[op],[0.],[1.],self.t,
            template.hessian_action,template.inverse_mobility,spectral_z=[0j])
        x,y=self.state();y=y[:1]
        np.testing.assert_array_equal(bridge.mass,0.)
        r=bridge.response(x,y,divide_mass=False)
        self.assertIsNone(r.population_velocity)
        with self.assertRaisesRegex(ValueError,'algebraic mass'):
            bridge.response(x,y)

    def test_fixed_contacts_and_phaseful_reference_are_rejected(self):
        bridge=self.bridge();x,y=self.state();x[0]=.1
        with self.assertRaisesRegex(ValueError,'contact'):
            bridge.response(x,y)
        d=self.d*np.exp(.1j)
        a,b=retarded.uniform_fields(d,.03-.8j);f,ft,g=retarded.fields(a,b)
        op=kinetic.assemble(self.graph,d,g,f,ft)
        with self.assertRaisesRegex(ValueError,'fixed-phase'):
            self.bridge([op]*len(self.energies))

    def test_energy_provenance_and_real_callback_contracts(self):
        bridge=self.bridge();x,y=self.state()
        with self.assertRaisesRegex(ValueError,'matching each'):
            LongitudinalReciprocalBridge(self.graph,self.operators,self.energies,self.weights,self.t,
                bridge.hessian_action,bridge.inverse_mobility,spectral_z=.03-2j*self.energies)
        altered=list(self.operators);altered[1]=altered[0]
        with self.assertRaisesRegex(ValueError,'declared energy'):
            self.bridge(altered)
        original=bridge.hessian_action
        for bad in (lambda x:np.ones(2),lambda x:np.full_like(x,np.nan),lambda x:original(x).astype(complex)):
            bridge.hessian_action=bad
            with self.assertRaisesRegex(ValueError,'finite real'):
                bridge.response(x,y)
        bridge.hessian_action=original;bridge.inverse_mobility=lambda x:x.astype(complex)
        with self.assertRaisesRegex(ValueError,'finite real'):
            bridge.response(x,y)


if __name__=='__main__':unittest.main()
