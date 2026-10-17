import numpy as np
from scipy.constants import elementary_charge as e,hbar
from pysnspd.experimental.general_thermal_gradient import *
from test_thermal_stiffness_coupled_candidate import small_system
from pysnspd.experimental.thermal_stiffness import thermal_rigidity
import pytest

def synthetic():
 rng=np.random.default_rng(113);n=6
 d=rng.normal(size=n)+1j*rng.normal(size=n)
 edges=np.array([[0,1],[0,2],[1,2],[2,3],[3,4],[4,5],[2,5]])
 c=1e-23*rng.uniform(.3,1.2,size=len(edges));alpha=rng.uniform(-1,1,size=len(edges))
 T=rng.uniform(3,12,size=n)
 return rng,d,edges,c,alpha,T

def poly(rho,L,T):
 Es=1e-23;Tc=8.65;l=L/Es;t=T/Tc
 a,b,c,d,e,f,u,v,w=.2,.3,.1,.07,.11,.09,.12,.31,.14
 return GradientPartials(Es*(a*rho*rho+b*rho*l+.5*c*l*l+d*rho*t+e*l*t+.5*f*t*t+u*rho+v*l+w*t),
   Es*(2*a*rho+b*l+d*t+u),b*rho+c*l+e*t+v,
   2*a*Es,b,c/Es,Es/Tc*(d*rho+e*l+f*t+w),Es*d/Tc,e/Tc,Es*f/Tc**2)

def ev(d,edges,c,alpha,T,provider=poly):
 geo=gradient_geometry(d,edges,c,alpha)
 return gradient_evaluate(d,T,geo,provider(geo.rho,geo.L_J,T))

def test_force_energy_and_link_current():
 rng,d,edges,c,alpha,T=synthetic();v=ev(d,edges,c,alpha,T)
 dd=rng.normal(size=len(d))+1j*rng.normal(size=len(d));da=rng.normal(size=len(c))
 h=1e-6
 fd=(ev(d+h*dd,edges,c,alpha+h*da,T).free_energy_J-ev(d-h*dd,edges,c,alpha-h*da,T).free_energy_J)/(2*h)
 analytic=2*np.real(np.vdot(v.force_J,dd))+v.link_derivative_alpha_J@da
 np.testing.assert_allclose(fd,analytic,rtol=3e-9,atol=1e-30)
 np.testing.assert_allclose(v.edge_current_A,-2*e/hbar*v.link_derivative_alpha_J,rtol=2e-15)
 for k in range(len(c)):
  ad=alpha.copy();ad[k]+=h;am=alpha.copy();am[k]-=h
  grad=(ev(d,edges,c,ad,T).free_energy_J-ev(d,edges,c,am,T).free_energy_J)/(2*h)
  np.testing.assert_allclose(grad,v.link_derivative_alpha_J[k],rtol=2e-8,atol=1e-31)

def test_hessian_current_J_and_symmetry():
 rng,d,edges,c,alpha,T=synthetic();v=ev(d,edges,c,alpha,T)
 H,J=gradient_hessian_current_jacobian(d,v)
 dd=rng.normal(size=len(d))+1j*rng.normal(size=len(d));x=np.column_stack((dd.real,dd.imag)).ravel()
 h=1e-6;vp=ev(d+h*dd,edges,c,alpha,T);vm=ev(d-h*dd,edges,c,alpha,T)
 fd=(vp.force_J-vm.force_J)/(2*h);hx=(H@x).reshape(-1,2)
 np.testing.assert_allclose(fd,hx[:,0]+1j*hx[:,1],rtol=2e-8,atol=1e-30)
 np.testing.assert_allclose((vp.edge_current_A-vm.edge_current_A)/(2*h),J@x,rtol=2e-8,atol=1e-14)
 np.testing.assert_allclose(H.toarray(),H.T.toarray(),rtol=2e-15,atol=1e-36)

def test_direction_all_fields_and_internal_entropy():
 rng,d,edges,c,alpha,T=synthetic();v=ev(d,edges,c,alpha,T)
 dd=rng.normal(size=len(d))+1j*rng.normal(size=len(d));dT=rng.normal(size=len(d));da=rng.normal(size=len(c))
 dg=gradient_direction(d,T,dd,dT,v,da)
 h=1e-6;vp=ev(d+h*dd,edges,c,alpha+h*da,T+h*dT);vm=ev(d-h*dd,edges,c,alpha-h*da,T-h*dT)
 np.testing.assert_allclose((vp.free_energy_J-vm.free_energy_J)/(2*h),dg.free_energy_direction_J,rtol=3e-9,atol=1e-30)
 np.testing.assert_allclose((vp.internal_energy_J-vm.internal_energy_J)/(2*h),dg.internal_energy_direction_J,rtol=3e-9,atol=1e-30)
 np.testing.assert_allclose((vp.force_J-vm.force_J)/(2*h),dg.force_direction_J,rtol=2e-8,atol=1e-30)
 np.testing.assert_allclose((vp.edge_current_A-vm.edge_current_A)/(2*h),dg.current_direction_A,rtol=2e-8,atol=1e-14)
 np.testing.assert_allclose(-(vp.entropy_J_K-vm.entropy_J_K)/(2*h),dg.G_T_direction_J_K,rtol=2e-8,atol=1e-31)
 np.testing.assert_allclose(dg.internal_energy_direction_J+T@dg.G_T_direction_J_K,
  2*np.real(np.vdot(v.force_J,dd))+v.link_derivative_alpha_J@da,rtol=2e-15,atol=1e-34)

def test_gauge_covariance():
 rng,d,edges,c,alpha,T=synthetic();v=ev(d,edges,c,alpha,T)
 chi=rng.normal(size=len(d));tail,head=edges.T
 w=ev(d*np.exp(1j*chi),edges,c,alpha+chi[head]-chi[tail],T)
 np.testing.assert_allclose(w.force_J,v.force_J*np.exp(1j*chi),rtol=5e-15,atol=1e-35)
 np.testing.assert_allclose(w.edge_current_A,v.edge_current_A,rtol=5e-15,atol=1e-19)
 np.testing.assert_allclose(w.internal_energy_J,v.internal_energy_J,rtol=5e-15)

def original_provider(f):
 def provider(rho,L,T):
  k=thermal_rigidity(rho,T,Tc_K=f.Tc_K,prefactor=f.radial_prefactor)
  return GradientPartials(L*k.kappa,L*k.rho,k.kappa,L*k.rho_rho,k.rho,
   np.zeros_like(rho),L*k.temperature,L*k.rho_temperature,k.temperature,L*k.temperature_temperature)
 return provider

def wrapper_from(f,provider):
 return GeneralThermalGradientFunctional(f.graph,energy_scale_J=f.Eu_J,N0_per_J_m3=f.N0_per_J_m3,
  ell0_m=f.ell0_m,thickness_m=f.thickness_m,count_nodes=f.xi_bar,count_weights=f.weights_bar,
  bath_temperature_K=f.bath_temperature_K,diffusion_m2_s=5e-5,
  gradient_partials=provider)

def test_original_G_L_kappa_complete_regression():
 m,s=small_system();f=m.functional;g=wrapper_from(f,original_provider(f))
 d,p=s.delta_bar,s.p;T=f.equivalent_temperature_K(d,p);alpha=np.array([.03,-.05])
 a=f.evaluate(d,p,alpha,temperature_K=T);b=g.evaluate(d,p,alpha,temperature_K=T)
 for key in a.__dataclass_fields__:
  if key=='force_contract':continue
  np.testing.assert_allclose(getattr(a,key),getattr(b,key),rtol=5e-13,atol=1e-32,err_msg=key)
 H,J=f.hessian_and_current_jacobian(d,p,alpha,temperature_K=T)
 h,j=g.hessian_and_current_jacobian(d,p,alpha,temperature_K=T)
 np.testing.assert_allclose(h.toarray(),H.toarray(),rtol=3e-13,atol=1e-32)
 np.testing.assert_allclose(j.toarray(),J.toarray(),rtol=3e-13,atol=1e-20)
 td=f.thermal_data(d,p,T,alpha);tg=g.thermal_data(d,p,T,alpha)
 for key in ('capacity_FD_J_K','B_J','epsilon_J','temperature_K','G_T_J_per_K',
  'G_TT_J_per_K2','G_Trho_J_per_K','capacity_effective_J_K'):
  np.testing.assert_allclose(tg[key],td[key],rtol=3e-13,atol=1e-32)
 dd=np.array([.2+.1j,-.1+.03j,.2-.03j]);dp=.03*p;da=np.array([.03,-.03])
 a=f.direction(d,p,dd,dp,alpha,da,temperature_K=T)
 b=g.direction(d,p,dd,dp,alpha,da,temperature_K=T)
 for key in a.__dataclass_fields__:
  np.testing.assert_allclose(getattr(a,key),getattr(b,key),rtol=3e-12,atol=1e-30,err_msg=key)

def test_wrapper_general_native_thermal_chain_FD():
 m,s=small_system();f=m.functional;g=wrapper_from(f,poly);d,p=s.delta_bar,s.p
 alpha=np.array([.03,-.05]);dd=np.array([.2+.1j,-.1+.03j,.2-.03j]);dp=.03*p;da=np.array([.03,-.03])
 a=g.direction(d,p,dd,dp,alpha,da)
 h=1e-6;plus=g.evaluate(d+h*dd,p+h*dp,alpha+h*da);minus=g.evaluate(d-h*dd,p-h*dp,alpha-h*da)
 np.testing.assert_allclose((plus.energy_J-minus.energy_J)/(2*h),a.energy_direction_J,rtol=2e-7,atol=1e-30)
 np.testing.assert_allclose((plus.gap_force_J_per_gapbar-minus.gap_force_J_per_gapbar)/(2*h),a.gap_force_direction_J_per_gapbar,rtol=2e-7,atol=1e-30)
 np.testing.assert_allclose((plus.edge_current_A-minus.edge_current_A)/(2*h),a.edge_current_direction_A,rtol=2e-7,atol=1e-16)

def test_alpha_jacobians_match_finite_difference_and_mixed_symmetry():
 rng,d,edges,c,alpha,T=synthetic();v=ev(d,edges,c,alpha,T)
 F,I=gradient_alpha_jacobians(d,v);da=rng.normal(size=len(c))
 h=1e-6;plus=ev(d,edges,c,alpha+h*da,T);minus=ev(d,edges,c,alpha-h*da,T)
 exact=(F@da).reshape(-1,2);exact=exact[:,0]+1j*exact[:,1]
 np.testing.assert_allclose((plus.force_J-minus.force_J)/(2*h),exact,rtol=2e-8,atol=1e-30)
 np.testing.assert_allclose((plus.edge_current_A-minus.edge_current_A)/(2*h),I@da,rtol=2e-8,atol=1e-14)
 np.testing.assert_allclose(I.toarray(),I.T.toarray(),rtol=2e-15,atol=1e-19)
 _,JI=gradient_hessian_current_jacobian(d,v)
 np.testing.assert_allclose(F.toarray(),(-hbar/(4*e)*JI.T).toarray(),rtol=3e-15,atol=1e-35)
