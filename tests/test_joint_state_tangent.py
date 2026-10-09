"""Finite-difference checks of physical state tangents, including moving bulk cuts."""
import unittest
import numpy as np
from scipy.optimize import root
from pysnspd.experimental import adiabatic_matrix_usadel as ad
from pysnspd.experimental import thermal_spatial_usadel as thermal
from pysnspd.experimental.bulk_current_reference import homogeneous_spectrum
from pysnspd.experimental.joint_state_tangent import (
    zeroth_kinetic_direction,joint_zeroth_rows_direction,
    MatsubaraStateTangent,ThermalBulkTangentBundle,
    current_carrying_thermal_cut_direction)
from pysnspd.experimental.joint_nambu_rows import JointNambuRows
from pysnspd.experimental.nambu_moments import NambuMoments


def bulk_retarded(gap,energy,eta):
    B=ad.generator(gap,energy,eta=eta).value
    scale=np.sqrt(abs(gap)**2+(eta-1j*energy)**2)
    return B/scale[:,None,None]


def bulk_retarded_direction(gap,direction,energy,eta):
    B=ad.generator(gap,energy,eta=eta).value
    dB=ad.generator(direction,0.).value
    scale=np.sqrt(abs(gap)**2+(eta-1j*energy)**2)
    ds=np.real(gap.conj()*direction)
    return dB/scale[:,None,None]-B*ds[:,None,None]/scale[:,None,None]**3


class ZerothStateTangentTests(unittest.TestCase):
    def setUp(self):
        self.graph=thermal.rectangular_graph(4,2,3.,1.)
        n=self.graph.n_nodes;x,y=self.graph.coordinates_bar.T
        self.gap=(1.1+.05*np.sin(x))*np.exp(.12j*x)
        self.dd=.08*np.cos(x)+.06j*np.sin(y+.2*x)
        self.alpha=np.linspace(-.07,.04,len(self.graph.edges))
        self.da=np.linspace(.02,-.03,len(self.graph.edges))
        self.h=np.column_stack((.35+.06*x,.04-.01*y))
        self.dh=np.column_stack((.02*np.cos(x),-.03*np.sin(y+.1)))
        self.energy=1.7;self.eta=.02

    def test_full_K_flux_and_local_residual_against_changed_state(self):
        R=bulk_retarded(self.gap,self.energy,self.eta)
        Rt=bulk_retarded_direction(self.gap,self.dd,self.energy,self.eta)
        B=ad.generator(self.gap,self.energy).value
        Bt=ad.generator(self.dd,0.).value
        tangent=zeroth_kinetic_direction(self.graph,R,Rt,self.h,self.dh,
            alpha=self.alpha,alpha_direction=self.da,B_zero=B,B_direction=Bt)
        step=2e-6;values=[]
        for sign in (-1.,1.):
            gap=self.gap+sign*step*self.dd
            h=self.h+sign*step*self.dh
            residual,flux,K=ad.kinetic_residual(self.graph,
                ad.Jet.of(bulk_retarded(gap,self.energy,self.eta)),
                ad.distribution(h[:,0],h[:,1]),
                ad.generator(gap,self.energy),alpha=self.alpha+sign*step*self.da,kappa=0.)
            values.append((K.value,flux.value,residual.value))
        for j,name in enumerate(("K_direction","flux_direction","residual_direction")):
            expected=(values[1][j]-values[0][j])/(2*step)
            np.testing.assert_allclose(getattr(tangent,name),expected,rtol=2e-8,atol=4e-10)
        self.assertIsNone(zeroth_kinetic_direction(self.graph,R,Rt,self.h,self.dh).residual_direction)
        with self.assertRaises(ValueError):
            zeroth_kinetic_direction(self.graph,R,Rt,self.h,self.dh,B_zero=B)

    def test_completed_SI_moments_against_changed_R_and_populations(self):
        energies=np.array([.7,1.7,3.])
        weights=np.array([.4,.7,1.2]);Eu=2e-22;I0=3e-6
        obs=NambuMoments(energies*Eu,weights*Eu,np.full(self.graph.n_nodes,2e-24),
            8e46,Eu,I0)
        x=self.graph.coordinates_bar[:,0]
        rows=JointNambuRows(obs,self.graph.edges,left_nodes=np.flatnonzero(x==x.min()),
            right_nodes=np.flatnonzero(x==x.max()))
        R=np.asarray([bulk_retarded(self.gap,E,self.eta) for E in energies])
        Rt=np.asarray([bulk_retarded_direction(self.gap,self.dd,E,self.eta) for E in energies])
        h=np.broadcast_to(self.h,(3,)+self.h.shape).copy()
        ht=np.broadcast_to(self.dh,h.shape).copy()
        thermal_h=np.zeros_like(h);thermal_h[:,:,0]=np.tanh(energies[:,None]/.4)
        dg=(.02+.03j)*obs.N0*Eu**2*obs.volume
        dI=np.linspace(.1e-6,-.2e-6,len(self.graph.edges))
        dIs=.3e-6
        got=joint_zeroth_rows_direction(rows,self.graph,retarded_zero=R,
            retarded_direction=Rt,h=h,h_direction=ht,thermal_h=thermal_h,
            thermal_gradient_direction_J_per_gapbar=dg,
            thermal_edge_current_direction_A=dI,
            alpha=self.alpha,alpha_direction=self.da,imposed_current_direction_A=dIs,
            reference_contract="same instantaneous spectrum and independent finite Matsubara reference")
        step=2e-6;values=[]
        for sign in (-1.,1.):
            gap=self.gap+sign*step*self.dd;links=self.alpha+sign*step*self.da
            Ks=[];Fs=[];Krefs=[];Frefs=[]
            for j,E in enumerate(energies):
                r=ad.Jet.of(bulk_retarded(gap,E,self.eta))
                for H,KK,FF in ((h[j]+sign*step*ht[j],Ks,Fs),(thermal_h[j],Krefs,Frefs)):
                    _,flux,K=ad.kinetic_residual(self.graph,r,
                        ad.distribution(H[:,0],H[:,1]),ad.generator(gap,E),
                        alpha=links,kappa=0.)
                    KK.append(K.value);FF.append(flux.value)
            z=np.zeros_like(np.asarray(Ks));fz=np.zeros_like(np.asarray(Fs))
            values.append(rows.evaluate(K_zero=np.asarray(Ks),K_first=z,
                K_thermal_zero=np.asarray(Krefs),phi_V=np.zeros(self.graph.n_nodes),
                thermal_gradient_J_per_gapbar=sign*step*dg,
                thermal_edge_current_A=np.full(len(self.graph.edges),20e-6)+sign*step*dI,
                edge_flux_zero=np.asarray(Fs),edge_flux_first=fz,
                thermal_edge_flux_zero=np.asarray(Frefs),
                imposed_current_A=20e-6+sign*step*dIs,
                reference_contract="same instantaneous spectrum and independent finite Matsubara reference"))
        for name in ("neutrality_V","gap_force_J_per_gapbar","edge_current_A",
                "port_current_residual_A","interior_current_divergence_A"):
            fd=(getattr(values[1],name)-getattr(values[0],name))/(2*step)
            np.testing.assert_allclose(getattr(got,name),fd,rtol=3e-7,atol=1e-28 if "gap" in name else 1e-13)


class ThermalStateTangentTests(unittest.TestCase):
    def setUp(self):
        self.graph=thermal.rectangular_graph(5,3,4.,2.)
        x,y=self.graph.coordinates_bar.T
        self.fixed=np.flatnonzero((x==x.min())|(x==x.max()))
        self.gap=(1.2+.03*np.cos(x))*np.exp(.13j*x)
        self.dd=.07*np.cos(.3*x)+.05j*np.sin(.4*y+.1)
        self.alpha=np.linspace(-.04,.03,len(self.graph.edges))
        self.da=np.linspace(.02,-.025,len(self.graph.edges))
        self.q=.13;self.dq=.04;self.t=.52
        self.ep=2*np.pi*self.t*(np.arange(8)+.5)
        self.energy_scale=2.7e-22;self.I0=2.1e-6
        self.solutions=self.solve(self.gap,self.alpha,self.q)
        self.bundle=ThermalBulkTangentBundle(self.graph,self.gap,self.t,self.ep,
            np.asarray([item.u for item in self.solutions]),self.alpha,
            fixed_nodes=self.fixed,bulk_q_bar=self.q,energy_scale_J=self.energy_scale,
            current_unit_A=self.I0,spectral_tolerance=1e-11)

    def solve(self,gap,alpha,q):
        solutions=[]
        for ep in self.ep:
            boundary=np.asarray([homogeneous_spectrum(abs(d),np.array([ep]),q*q)[0][0]
                *d/abs(d) for d in gap[self.fixed]])
            free=np.setdiff1d(np.arange(self.graph.n_nodes),self.fixed)
            u=gap/ep;u[self.fixed]=boundary
            def fields(v):
                out=u.copy();out[free]=v[0::2]+1j*v[1::2]
                return out
            def residual(v):
                out=thermal.spectral_energy_gradient(self.graph,gap,ep,fields(v),alpha).residual[free]
                return np.column_stack((out.real,out.imag)).ravel()
            components=np.ravel(np.column_stack((2*free,2*free+1)))
            def jacobian(v):
                _,jac=thermal.spectral_residual_jacobian(self.graph,gap,ep,fields(v),alpha)
                return jac[components][:,components].toarray()
            first=np.column_stack((u[free].real,u[free].imag)).ravel()
            solved=root(residual,first,jac=jacobian,tol=1e-12)
            error=float(np.max(abs(residual(solved.x)),initial=0))
            if error>2e-13:raise AssertionError("Independent reference root failed: "+str(error))
            field=fields(solved.x)
            ev=thermal.spectral_energy_gradient(self.graph,gap,ep,field,alpha)
            solutions.append(thermal.SpectralSolution(float(ep),field,ev.f,ev.g,ev.energy,
                error,float(np.max(abs(ev.gradient_u[free]),initial=0)),0,tuple(),
                self.fixed.copy(),thermal._signature(self.graph,gap,alpha)))
        return solutions

    def test_bulk_cut_derivative_against_actual_homogeneous_solutions(self):
        ep=self.ep[0];step=2e-6
        got=current_carrying_thermal_cut_direction(ep,self.gap[self.fixed],
            self.dd[self.fixed],self.q,q_direction=self.dq)
        values=[]
        for sign in (-1.,1.):
            gap=self.gap[self.fixed]+sign*step*self.dd[self.fixed]
            q=self.q+sign*step*self.dq
            values.append(np.asarray([homogeneous_spectrum(abs(d),np.array([ep]),q*q)[0][0]
                *d/abs(d) for d in gap]))
        np.testing.assert_allclose(got,(values[1]-values[0])/(2*step),rtol=2e-8,atol=3e-10)

    def test_finite_Matsubara_Hessian_and_current_against_new_stationary_states(self):
        got=self.bundle.apply(self.dd,alpha_direction=self.da,q_direction=self.dq)
        step=2e-5;values=[]
        for sign in (-1.,1.):
            gap=self.gap+sign*step*self.dd;alpha=self.alpha+sign*step*self.da
            solutions=self.solve(gap,alpha,self.q+sign*step*self.dq)
            values.append(thermal.evaluate_thermal(self.graph,gap,self.t,solutions,alpha))
        gradient=self.energy_scale*(values[1].gap_gradient-values[0].gap_gradient)/(2*step)
        current=self.I0*(values[1].current_bar-values[0].current_bar)/(2*step)
        np.testing.assert_allclose(got.gradient_J_per_gapbar,gradient,rtol=2e-8,atol=3e-31)
        np.testing.assert_allclose(got.edge_current_A,current,rtol=2e-8,atol=3e-15)
        self.assertLess(got.maximum_spectral_residual,1e-13)
        self.assertEqual(got.modes,8)

    def test_global_phase_preserves_current_and_rotates_thermal_force(self):
        got=self.bundle.apply(1j*self.gap)
        base=thermal.evaluate_thermal(self.graph,self.gap,self.t,self.solutions,self.alpha)
        np.testing.assert_allclose(got.gradient_J_per_gapbar,
            1j*self.energy_scale*base.gap_gradient,rtol=2e-12,atol=1e-33)
        np.testing.assert_allclose(got.edge_current_A,0.,atol=2e-20)

    def test_explicit_cut_law_is_required_and_zero_current_BCS_is_rejected(self):
        mode=MatsubaraStateTangent(self.graph,self.gap,self.ep[0],self.solutions[0].u,
            self.alpha,fixed_nodes=self.fixed)
        with self.assertRaisesRegex(ValueError,"Explicit moving"):
            mode.apply(self.dd)
        modes=np.asarray([s.u for s in self.solutions])
        modes[:,self.fixed]=self.gap[self.fixed]/self.ep[:,None]
        with self.assertRaises(ValueError):
            ThermalBulkTangentBundle(self.graph,self.gap,self.t,self.ep,modes,self.alpha,
                fixed_nodes=self.fixed,bulk_q_bar=self.q,energy_scale_J=self.energy_scale,
                current_unit_A=self.I0)


if __name__=="__main__":
    unittest.main()
