"""Variational, thermodynamic and support checks of finite-depairing closure."""
import numpy as np
import pytest
from scipy.constants import Boltzmann as kB
from test_thermal_stiffness_coupled_candidate import small_system
from pysnspd.experimental.resummed_thermal import ResummedThermalFunctional
from pysnspd.experimental.general_thermal_gradient import gradient_alpha_jacobians
from pysnspd.experimental.thermal_stiffness import entropy_compatible_carrier
from pysnspd.experimental.effective_bcs_functional import bcs_count_quadrature

def convert(old):
 D=old.ell0_m**2*2*old.Eu_J/__import__('scipy').constants.hbar
 return ResummedThermalFunctional(old.graph,energy_scale_J=old.Eu_J,N0_per_J_m3=old.N0_per_J_m3,
  ell0_m=old.ell0_m,thickness_m=old.thickness_m,count_nodes=old.xi_bar,count_weights=old.weights_bar,
  bath_temperature_K=old.bath_temperature_K,diffusion_m2_s=D)

def system_state(conductance=.1):
 m,s=small_system(conductance=conductance);m.functional=convert(m.functional)
 return m,s

def test_native_energy_chain_and_force_current_physical_fields():
 m,s=system_state();f=m.functional;d,p=s.delta_bar,s.p;alpha=np.array([.04,-.03])
 dd=np.array([.021+.012j,-.018+.022j,.011-.009j]);dp=.015*p;da=np.array([.02,-.01])
 exact=f.direction(d,p,dd,dp,alpha,da)
 errors=[]
 for h in (2e-5,1e-5):
  plus=f.evaluate(d+h*dd,p+h*dp,alpha+h*da);minus=f.evaluate(d-h*dd,p-h*dp,alpha-h*da)
  fd=(plus.energy_J-minus.energy_J)/(2*h)
  errors.append(abs(fd-exact.energy_direction_J)/max(abs(exact.energy_direction_J),1e-30))
  np.testing.assert_allclose((plus.gap_force_J_per_gapbar-minus.gap_force_J_per_gapbar)/(2*h),
    exact.gap_force_direction_J_per_gapbar,rtol=2e-6,atol=1e-31)
  np.testing.assert_allclose((plus.edge_current_A-minus.edge_current_A)/(2*h),
    exact.edge_current_direction_A,rtol=2e-6,atol=1e-18)
 assert errors[-1]<3e-6

def test_fixed_temperature_force_hessian_and_link_current_jacobians():
 m,s=system_state();f=m.functional;d,p=s.delta_bar,s.p;alpha=np.array([.04,-.03])
 T=f.equivalent_temperature_K(d,p);dd=np.array([.021+.012j,-.018+.022j,.011-.009j])
 da=np.array([.02,-.01]);cart=np.column_stack((dd.real,dd.imag)).ravel()
 value=f.evaluate(d,p,alpha,temperature_K=T);H,JI=f.hessian_and_current_jacobian(d,p,alpha,temperature_K=T)
 h=1e-5;plus=f.evaluate(d+h*dd,p,alpha,temperature_K=T);minus=f.evaluate(d-h*dd,p,alpha,temperature_K=T)
 exact=(H@cart).reshape(-1,2);exact=exact[:,0]+1j*exact[:,1]
 np.testing.assert_allclose((plus.gap_force_J_per_gapbar-minus.gap_force_J_per_gapbar)/(2*h),exact,rtol=2e-6,atol=1e-31)
 np.testing.assert_allclose((plus.edge_current_A-minus.edge_current_A)/(2*h),JI@cart,rtol=2e-6,atol=1e-18)
 np.testing.assert_allclose(H.toarray(),H.T.toarray(),rtol=1e-14,atol=1e-34)
 ge=f.gradient_evaluation(d,T,alpha);Fa,Ia=gradient_alpha_jacobians(d,ge)
 plus=f.evaluate(d,p,alpha+h*da,temperature_K=T);minus=f.evaluate(d,p,alpha-h*da,temperature_K=T)
 np.testing.assert_allclose((plus.energy_J-minus.energy_J)/(2*h),
   value.link_derivative_alpha_J@da-float(np.dot(T,ge.partials.G_L_T_per_K*(
   f.gradient_rate_geometry(d,np.zeros_like(d),alpha,da)))),rtol=2e-6,atol=1e-31)
 exact=(Fa@da).reshape(-1,2);exact=exact[:,0]+1j*exact[:,1]
 np.testing.assert_allclose((plus.gap_force_J_per_gapbar-minus.gap_force_J_per_gapbar)/(2*h),exact,rtol=2e-6,atol=1e-31)
 np.testing.assert_allclose((plus.edge_current_A-minus.edge_current_A)/(2*h),Ia@da,rtol=2e-6,atol=1e-18)
 np.testing.assert_allclose((Fa+__import__('scipy').constants.hbar/(4*__import__('scipy').constants.elementary_charge)*JI.T).toarray(),0.,atol=1e-34)

def test_continuous_independent_total_energy_balance_and_entropy():
 m,s=system_state();f=m.functional;r=m.evaluate(s);k=m.kinetics
 dp=r.dp_per_ps*1e12;dn=r.dn_per_ps*1e12;dd=r.delta_per_ps*1e12
 film=f.direction(s.delta_bar,s.p,dd,dp).energy_direction_J
 phonon=float(np.sum(k.phonon_energy_weights_J*dn))
 cmrate=r.circuit_scaled_per_ps*m.circuit_scales*1e12
 Ib,Is,vc=s.circuit_SI
 circuit=float(np.dot([m.circuit.L_bias_H*Ib,m.circuit.Lk_ext_H*Is,m.circuit.C_couple_F*vc],cmrate))
 total=film+phonon+circuit
 scale=abs(total)+abs(r.external_power_W)+1e-25
 assert abs(total-r.external_power_W)<2e-11*scale
 assert abs(r.instant_power_residual_W)<2e-11*scale
 assert abs(r.film_power_residual_W)<2e-11*scale
 td=f.thermal_data(s.delta_bar,s.p,r.temperature_K)
 eps=td['epsilon_J'];en=r.kinetics.electron_energies_J/f.Eu_J
 sigma=kB*eps/(f.Eu_J*en)*(np.log1p(-s.p)-np.log(s.p))
 P=r.kinetics.heating_power_W
 excess=np.sum(sigma*r.kinetics.heating_dp_dt_s,axis=1)-P/r.temperature_K
 assert np.all(excess>=-2e-11*(abs(P/r.temperature_K)+1e-30))
 assert np.all(td['capacity_effective_J_K']>0)
 np.testing.assert_array_equal(r.dp_per_ps[m.cut_nodes],0.)
 np.testing.assert_array_equal(r.dn_per_ps[m.cut_nodes],0.)

@pytest.mark.parametrize("endpoint",[0.,1.])
def test_closed_Pauli_endpoints_have_inward_heat_and_full_rates(endpoint):
 m,s=system_state();s.p[1,0]=endpoint;r=m.evaluate(s)
 rate=r.dp_per_ps[1,0]
 assert rate>=0 if endpoint==0 else rate<=0
 f=m.functional;td=f.thermal_data(s.delta_bar,s.p,r.temperature_K)
 en=r.kinetics.electron_energies_J/f.Eu_J;eps=td['epsilon_J']
 for power in (1e-12,-1e-12):
  carrier=entropy_compatible_carrier(s.p,eps,r.temperature_K,f.Eu_J,en,power_W=np.full(3,power))
  heat=carrier*power
  assert heat[1,0]>=0 if endpoint==0 else heat[1,0]<=0
  np.testing.assert_allclose(np.sum(eps*heat,axis=1),power,rtol=4e-15,atol=1e-28)

def test_exact_zero_gap_energy_force_corner_without_exact_curvature_claim():
 m,s=system_state();f=m.functional;T=np.array([.9,7.,.9])
 d=s.delta_bar.copy();d[1]=0.;p=f.fermi_dirac(d,T)
 value=f.evaluate(d,p,temperature_K=T)
 assert np.isfinite(value.energy_J) and np.all(np.isfinite(value.gap_force_J_per_gapbar))
 assert f.exact_zero_gap_curvature_preconditioner_count==0
 # The isolated normal node has a finite exact curvature after DC matching.
 dd=np.array([0.,.01+.02j,0.],complex);dp=np.zeros_like(p)
 exact=f.direction(d,p,dd,dp,temperature_K=T,temperature_direction_K=np.zeros(3))
 h=1e-5;plus=f.evaluate(d+h*dd,p,temperature_K=T);minus=f.evaluate(d-h*dd,p,temperature_K=T)
 np.testing.assert_allclose((plus.energy_J-minus.energy_J)/(2*h),exact.energy_direction_J,rtol=5e-5,atol=1e-29)

def test_positive_capacity_grid_with_native_full_band():
 m,s=small_system();old=m.functional
 xi,w=bcs_count_quadrature(32.,6)
 f=ResummedThermalFunctional(old.graph,energy_scale_J=old.Eu_J,N0_per_J_m3=old.N0_per_J_m3,
  ell0_m=old.ell0_m,thickness_m=old.thickness_m,count_nodes=xi,count_weights=w,
  bath_temperature_K=.9,diffusion_m2_s=5e-5)
 minimum=np.inf;count=0
 for temp in (.9,3.,6.,8.65,12.,25.):
  T=np.full(3,temp);theta=T/f.Tc_K
  for a in (.01,.1,.5,1.,2.):
   rho=np.full(3,a*a);d=np.full(3,a,dtype=complex);en=f.energies_bar(d)
   fd=f.fermi_dirac(d,T)
   C=4*f.node_energy_scale_J*np.sum(w*en**2*fd*(1-fd),axis=1)/(f.Tc_K*theta**2)
   for Q in (0.,.01,.1,.5,1.,5.,20.,100.):
    L=f.node_energy_scale_J*rho*Q
    data=f._gradient_partials(rho,L,T)
    ratio=(C-T*data.G_T_T_J_K2)/C
    assert np.all(np.isfinite(ratio)) and np.min(ratio)>0,(a,temp,Q,ratio)
    minimum=min(minimum,float(np.min(ratio)));count+=3
 assert count==720 and minimum>.8

def test_coupled_midpoint_step_refinement_and_nonzero_progress():
 m,s=system_state()
 def evolution(n):
  state=s.copy();wallpower=0.
  for _ in range(n):
   state,receipt=m.midpoint_step(state,1e-4/n,tolerance=2e-10,max_iterations=60)
   wallpower+=receipt['external_energy_increment_J']
   assert receipt['iterations']>0
  return state,wallpower
 results=[evolution(n) for n in (1,2,4)]
 derrors=[np.max(abs(results[i][0].delta_bar-results[2][0].delta_bar)) for i in (0,1)]
 assert derrors[0]/derrors[1]>3.5
 assert np.max(abs(results[0][0].p-s.p))>0
 assert np.max(abs(results[0][0].delta_bar-s.delta_bar))>0
 initial=m.evaluate(s).energy_J
 errors=[abs(m.evaluate(state).energy_J-initial-external) for state,external in results]
 assert errors[2]<errors[0]

def test_fully_normal_homogeneous_origin_uses_explicit_preconditioner_convention():
 m,s=system_state();f=m.functional;T=np.array([.9,7.,.9])
 d=np.zeros(3,complex);p=f.fermi_dirac(d,T)
 value=f.evaluate(d,p,temperature_K=T)
 assert np.isfinite(value.energy_J)
 np.testing.assert_array_equal(value.edge_current_A,0.)
 np.testing.assert_array_equal(value.gap_force_J_per_gapbar,0.)
 assert f.exact_zero_gap_curvature_preconditioner_count==3
 H,J=f.hessian_and_current_jacobian(d,p,temperature_K=T)
 assert np.all(np.isfinite(H.data)) and np.all(np.isfinite(J.data))
 # No exact Hessian claim at this direction-dependent origin.

def test_DC_matching_retains_full_positive_depairing_branch():
 from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional
 m,s=system_state();f=m.functional
 kwargs=dict(graph=f.graph,energy_scale_J=f.Eu_J,N0_per_J_m3=f.N0_per_J_m3,
     ell0_m=f.ell0_m,thickness_m=f.thickness_m,count_nodes=f.xi_bar,count_weights=f.weights_bar,
     bath_temperature_K=f.bath_temperature_K,diffusion_m2_s=5e-5,gradient_law="bcs_radial")
 base=EffectiveBCSFunctional(**kwargs)
 expected=base.bulk(temperature_K=.9,current_A=0.,width_m=80e-9)
 actual=f.bulk(temperature_K=.9,current_A=0.,width_m=80e-9)
 np.testing.assert_allclose(actual.depairing_current_A,expected.depairing_current_A,rtol=2e-12)
 np.testing.assert_allclose(actual.depairing_q_bar,expected.depairing_q_bar,rtol=5e-8)
 np.testing.assert_allclose(f.bulk_amplitude(.9,expected.depairing_q_bar),
     base.bulk_amplitude(.9,expected.depairing_q_bar),rtol=2e-12)

def test_DC_matching_preserves_reference_force_current_and_free_gradient():
 from pysnspd.experimental.effective_bcs_functional import EffectiveBCSFunctional
 m,s=system_state();f=m.functional;alpha=np.array([.07,-.04])
 base=EffectiveBCSFunctional(f.graph,energy_scale_J=f.Eu_J,N0_per_J_m3=f.N0_per_J_m3,
     ell0_m=f.ell0_m,thickness_m=f.thickness_m,count_nodes=f.xi_bar,count_weights=f.weights_bar,
     bath_temperature_K=f.bath_temperature_K,diffusion_m2_s=5e-5,gradient_law="bcs_radial")
 T=np.full(3,.9);p=f.fermi_dirac(s.delta_bar,T)
 actual=f.evaluate(s.delta_bar,p,alpha,temperature_K=T);expected=base.evaluate(s.delta_bar,p,alpha)
 np.testing.assert_allclose(actual.gap_force_J_per_gapbar,expected.gap_force_J_per_gapbar,rtol=3e-14,atol=1e-34)
 np.testing.assert_allclose(actual.edge_current_A,expected.edge_current_A,rtol=3e-14,atol=1e-24)
 np.testing.assert_allclose(actual.free_gradient_energy_J,expected.gradient_energy_J,rtol=3e-14)
