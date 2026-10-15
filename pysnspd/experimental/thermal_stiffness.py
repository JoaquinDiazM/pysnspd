
'''Independent mathematical prototype for the TE-stiffness candidate.
No production import or trajectory admission is implied.
G is a FREE-energy extension. Stored Ugrad=G-sum(T G_T), Sgrad=-G_T.
'''
from dataclasses import dataclass
import numpy as np
from scipy.constants import Boltzmann as kB, elementary_charge as electron_charge, hbar
from scipy.special import expit
from numpy.polynomial.polynomial import polyval,polyder

@dataclass(frozen=True)
class ThermalRigidity:
    kappa:np.ndarray
    rho:np.ndarray
    rho_rho:np.ndarray
    temperature:np.ndarray
    temperature_temperature:np.ndarray
    rho_temperature:np.ndarray

def thermal_rigidity(rho,temperature_K,*,Tc_K,prefactor=np.pi/2):
    rho,T=np.broadcast_arrays(np.asarray(rho,float),np.asarray(temperature_K,float))
    if np.any(rho<0)or np.any(T<=0)or np.any(~np.isfinite(rho+T)):
        raise ValueError("Finite rho>=0 and T>0 required")
    theta=T/Tc_K;b=1/(2*theta);a=np.sqrt(rho);x=b*a
    t=np.tanh(x);sech2=4*expit(-2*x)*expit(2*x)
    small=x<.1
    coeff=np.array([1.,-1/3,2/15,-17/315,62/2835,-1382/155925,21844/6081075,-929569/638512875])
    out=[np.empty_like(rho)for _ in range(3)]
    for order in range(3):
        out[order][small]=b[small]**(2*order+1)*polyval(x[small]**2,polyder(coeff,order))
    regular=~small;ar=a[regular];br=b[regular];tr=t[regular];sr=sech2[regular]
    out[0][regular]=tr/ar
    out[1][regular]=(br*ar*sr-tr)/(2*ar**3)
    out[2][regular]=(-2*br**2*ar**2*sr*tr-3*br*ar*sr+3*tr)/(4*ar**5)
    tanh_over_a=np.divide(t,a,out=np.array(b,copy=True),where=a>0)
    kt=-prefactor*sech2/(2*Tc_K*theta**2)
    ktt=prefactor*sech2*(1-x*t)/(Tc_K**2*theta**3)
    krt=prefactor*sech2*tanh_over_a/(4*Tc_K*theta**3)
    return ThermalRigidity(*(prefactor*v for v in out),kt,ktt,krt)

@dataclass(frozen=True)
class ElectronicTemperatureDifferentials:
    epsilon_J:np.ndarray
    sigma_J_per_K:np.ndarray
    heat_capacity_J_per_K:np.ndarray
    B_J:np.ndarray
    dT_dp_K:np.ndarray
    dT_drho_K:np.ndarray
    fd:np.ndarray

def electronic_differentials(rho,p,T,*,xi_bar,weights_bar,node_energy_scale_J,Eu_J):
    rho=np.asarray(rho,float);p=np.asarray(p,float);T=np.asarray(T,float)
    xi=np.asarray(xi_bar,float);w=np.asarray(weights_bar,float);u=np.asarray(node_energy_scale_J,float)
    if np.any(T<=0)or np.any(p<=0)or np.any(p>=1):
        raise ValueError("T>0 and strictly interior occupations required for entropy derivatives")
    Tc=Eu_J/kB;theta=T/Tc
    e=np.sqrt(rho[:,None]+xi[None,:]**2);fd=expit(-e/theta[:,None])
    epsilon=4*u[:,None]*w*e
    counts=4*u[:,None]*w/Eu_J
    sigma=kB*counts*(np.log1p(-p)-np.log(p))
    C=4*u*np.sum(w*e**2*fd*(1-fd),axis=1)/(Tc*theta**2)
    B=2*u*np.sum(w*((p-fd)/e+fd*(1-fd)/theta[:,None]),axis=1)
    if np.any(C<=0):
        raise ArithmeticError("Finite positive native FD capacity required")
    return ElectronicTemperatureDifferentials(epsilon,sigma,C,B,
        epsilon/C[:,None],B/C,fd)

@dataclass(frozen=True)
class GradientData:
    free_energy_J:float
    internal_energy_J:float
    entropy_J_per_K:np.ndarray
    L_J:np.ndarray
    G_T_J_per_K:np.ndarray
    G_TT_J_per_K2:np.ndarray
    G_Trho_J_per_K:np.ndarray
    gap_force_J:np.ndarray
    edge_current_A:np.ndarray

def gradient_data(d,T,*,edges,edge_energy_scale_J,Tc_K,prefactor=np.pi/2,alpha=None):
    d=np.asarray(d,complex);T=np.asarray(T,float);edges=np.asarray(edges,int)
    tail,head=edges.T;c=np.asarray(edge_energy_scale_J,float)
    alpha=np.zeros(len(edges))if alpha is None else np.asarray(alpha,float)
    phase=np.exp(-1j*alpha);difference=d[tail]-phase*d[head];diff2=abs(difference)**2
    kr=thermal_rigidity(abs(d)**2,T,Tc_K=Tc_K,prefactor=prefactor)
    ke=(kr.kappa[tail]+kr.kappa[head])/2
    edge=c*diff2/2
    L=np.bincount(tail,edge,minlength=len(d))+np.bincount(head,edge,minlength=len(d))
    gt=L*kr.temperature;gtt=L*kr.temperature_temperature;gtr=L*kr.rho_temperature
    force=np.zeros_like(d)
    np.add.at(force,tail,c*(ke*difference+.5*kr.rho[tail]*d[tail]*diff2))
    np.add.at(force,head,c*(-ke*np.conj(phase)*difference+.5*kr.rho[head]*d[head]*diff2))
    current=4*electron_charge/hbar*c*ke*np.imag(np.conj(d[tail])*phase*d[head])
    G=float(np.dot(L,kr.kappa))
    return GradientData(G,G-float(np.dot(T,gt)),-gt,L,gt,gtt,gtr,force,current)

def gradient_entropy_derivative_at_fixed_T(d,dd,T,*,edges,edge_energy_scale_J,Tc_K,
                                          prefactor=np.pi/2,alpha=None,alpha_rate=None):
    d=np.asarray(d,complex);dd=np.asarray(dd,complex);edges=np.asarray(edges,int)
    tail,head=edges.T;c=np.asarray(edge_energy_scale_J,float)
    alpha=np.zeros(len(edges))if alpha is None else np.asarray(alpha,float)
    da=np.zeros(len(edges))if alpha_rate is None else np.asarray(alpha_rate,float)
    phase=np.exp(-1j*alpha);diff=d[tail]-phase*d[head]
    diff_dot=dd[tail]-phase*dd[head]+1j*da*phase*d[head]
    edge_L=.5*c*abs(diff)**2;edge_Ldot=c*np.real(np.conj(diff)*diff_dot)
    L=np.bincount(tail,edge_L,minlength=len(d))+np.bincount(head,edge_L,minlength=len(d))
    Ldot=np.bincount(tail,edge_Ldot,minlength=len(d))+np.bincount(head,edge_Ldot,minlength=len(d))
    rhodot=2*np.real(np.conj(d)*dd)
    kr=thermal_rigidity(abs(d)**2,T,Tc_K=Tc_K,prefactor=prefactor)
    return kr.temperature*Ldot+L*kr.rho_temperature*rhodot

def balanced_entropy_carrier(epsilon_J,sigma_J_per_K,fd,temperature_K,*,Eu_J,e_bar):
    '''Positive heat carrier: epsilon.h=1 and sigma.h=1/TE.

    This explicit phenomenological carrier balances positive and negative
    entropy affinities of the thermal tangent. Exact FD uses its thermal
    tangent. Any signed power is permitted; endpoint support requires an
    implicit solve and step rejection, never clipping.
    '''
    eps=np.asarray(epsilon_J,float);sig=np.asarray(sigma_J_per_K,float)
    T=np.asarray(temperature_K,float);fd=np.asarray(fd,float);eb=np.asarray(e_bar,float)
    r=T[:,None]*sig/eps-1
    v=Eu_J*eb*fd*(1-fd)/(kB*T[:,None]**2)
    rp=np.maximum(r,0.);rm=np.maximum(-r,0.)
    Ap=np.sum(eps*v*rp,axis=1);Am=np.sum(eps*v*rm,axis=1)
    close=np.max(abs(r),axis=1)<256*np.finfo(float).eps
    bad=(~close)&((Ap<=0)|(Am<=0))
    if np.any(bad):
        raise ArithmeticError("TE must match the native energy moment; no balanced carrier bracket")
    factor=np.where(r>0,Am[:,None],np.where(r<0,Ap[:,None],((Ap+Am)/2)[:,None]))
    raw=np.where(close[:,None],v,v*factor)
    normalizer=np.sum(eps*raw,axis=1)
    if np.any(normalizer<=0):
        raise ArithmeticError("No finite entropy-compatible heat carrier")
    h=raw/normalizer[:,None]
    return h,dict(energy_moment=np.sum(eps*h,axis=1),
        entropy_moment_T=np.sum(sig*h,axis=1)*T,
        exact_FD_branch=close)

@dataclass(frozen=True)
class ReversibleLedger:
    temperature_rate_K_per_s:np.ndarray
    source_power_W:np.ndarray
    reversible_power_W:np.ndarray
    heat_capacity_effective_J_per_K:np.ndarray
    G_T_rate_W_per_K:np.ndarray

def reversible_ledger(T,C,B,rho_rate,G_TT,G_T_rate_fixed_T,base_native_electron_power_W,
                      dissipative_power_W):
    '''Algebraic continuous ledger; powers are analytic constitutive terms.
    base_native_electron_power includes collisions, transport, BGK, QDelta+PJ.
    dissipative_power is QDelta+PJ only. No energy residual enters this call.
    '''
    T=np.asarray(T,float);C=np.asarray(C,float);B=np.asarray(B,float)
    Ceff=C-T*G_TT
    if np.any(Ceff<=0)or np.any(~np.isfinite(Ceff)):
        bad=np.flatnonzero((Ceff<=0)|(~np.isfinite(Ceff)))
        raise ArithmeticError("Candidate left Ceff>0 domain at nodes "+str(bad[:12]))
    Tdot=(base_native_electron_power_W+B*rho_rate+T*G_T_rate_fixed_T)/Ceff
    gtdot=G_T_rate_fixed_T+G_TT*Tdot
    rev=T*gtdot
    return ReversibleLedger(Tdot,np.asarray(dissipative_power_W)+rev,rev,Ceff,gtdot)

from .effective_bcs_functional import EffectiveBCSFunctional,EffectiveBCSEvaluation,EffectiveBCSDirection

@dataclass(frozen=True)
class ThermalStiffnessEvaluation(EffectiveBCSEvaluation):
    free_gradient_energy_J:float
    gradient_entropy_J_per_K:np.ndarray
    gradient_G_T_J_per_K:np.ndarray
    gradient_G_TT_J_per_K2:np.ndarray
    gradient_L_J:np.ndarray
    temperature_K:np.ndarray
    force_contract:str

def entropy_compatible_carrier(p,epsilon_J,temperature_K,Eu_J,enbar,weights_bar=None):
    '''Positive carrier: epsilon.h=1, sigma.h=1/TE. No residual argument.'''
    p=np.asarray(p,float);epsilon=np.asarray(epsilon_J,float)
    T=np.asarray(temperature_K,float);eb=np.asarray(enbar,float)
    if np.any(p<=0)or np.any(p>=1)or np.any(T<=0):
        raise ValueError("Carrier requires TE>0 and interior Pauli support")
    counts=epsilon/(Eu_J*eb)
    sigma=kB*counts*(np.log1p(-p)-np.log(p))
    fd=expit(-Eu_J*eb/(kB*T[:,None]))
    return balanced_entropy_carrier(epsilon,sigma,fd,T,Eu_J=Eu_J,e_bar=eb)[0]

class ThermalStiffnessFunctional(EffectiveBCSFunctional):
    '''Candidate K(|Delta|,TE) FREE gradient with complete internal U/S.

    energy_J and node_energy_J are U_BCS+G-sum(T G_T).
    gap_force=U_d-sum(T S_d)=d(U_BCS+G)/dd* at fixed T, not U_d at fixed p.
    Current derives from that same G. Hessian/current Jacobian fix T and serve
    only as preconditioners of the fully coupled RHS. No global admission.
    '''
    closure_id="TE_free_gradient_internal_entropy_candidate_v1"
    def __init__(self,*args,**kwargs):
        law=kwargs.pop("gradient_law","bcs_radial")
        if law not in ("bcs_radial","bcs_radial_TE_free_energy"):
            raise ValueError("Radial dirty-limit rigidity required")
        super().__init__(*args,gradient_law="bcs_radial",**kwargs)
        self.gradient_law="bcs_radial_TE_free_energy"

    def rigidity(self,rho,temperature_K=None):
        T=self.bath_temperature_K if temperature_K is None else temperature_K
        r=thermal_rigidity(rho,T,Tc_K=self.Tc_K,prefactor=self.radial_prefactor)
        return r.kappa,r.rho,r.rho_rho

    def equivalent_temperature_K(self,gapbar,p):
        d,p,_=self._state(gapbar,p,None);en=self.energies_bar(d);w=self.weights_bar
        target=np.sum(w*en*p,axis=1);upper=.5*np.sum(w*en,axis=1)
        if np.any(target<=0)or np.any(target>=upper):
            raise ValueError("Finite positive energy-equivalent temperature required")
        lo=np.zeros(len(d));hi=np.ones(len(d))
        for _ in range(128):
            bad=np.sum(w*en*expit(-en/hi[:,None]),axis=1)<target
            if not np.any(bad):break
            hi=np.where(bad,2*hi,hi)
        else:raise ArithmeticError("No finite TE bracket")
        for _ in range(64):
            mid=(lo+hi)/2
            bad=np.sum(w*en*expit(-en/mid[:,None]),axis=1)<target
            lo=np.where(bad,mid,lo);hi=np.where(bad,hi,mid)
        return (lo+hi)/2*self.Tc_K

    def _temperature(self,d,p,temperature_K):
        if temperature_K is None:return self.equivalent_temperature_K(d,p)
        T=np.broadcast_to(np.asarray(temperature_K,float),(self.graph.n_nodes,)).copy()
        if np.any(~np.isfinite(T))or np.any(T<=0):
            raise ValueError("Finite positive TE required; no bath floor for stiffness")
        return T

    def _fixed_T(self,T):
        # A separate base-class object prevents temperature mutation or recursion.
        obj=object.__new__(EffectiveBCSFunctional);obj.__dict__.update(self.__dict__)
        def fixed_rigidity(rho):
            kr=thermal_rigidity(rho,T,Tc_K=self.Tc_K,prefactor=self.radial_prefactor)
            return kr.kappa,kr.rho,kr.rho_rho
        obj.rigidity=fixed_rigidity
        return obj

    def _L(self,d,alpha=None):
        tail,head=self.graph.edges.T
        phase=1 if alpha is None else np.exp(-1j*np.asarray(alpha,float))
        edge=.5*self.edge_energy_scale_J*abs(d[tail]-phase*d[head])**2
        return np.bincount(tail,edge,minlength=len(d))+np.bincount(head,edge,minlength=len(d))

    def gradient_rate_geometry(self,d,dz_per_s,alpha=None,alpha_rate=None):
        d=np.asarray(d,complex);dd=np.asarray(dz_per_s,complex);tail,head=self.graph.edges.T
        alpha=np.zeros(len(tail))if alpha is None else np.asarray(alpha,float)
        da=np.zeros(len(tail))if alpha_rate is None else np.asarray(alpha_rate,float)
        phase=np.exp(-1j*alpha);diff=d[tail]-phase*d[head]
        ddiff=dd[tail]-phase*dd[head]+1j*da*phase*d[head]
        edge=self.edge_energy_scale_J*np.real(np.conj(diff)*ddiff)
        return np.bincount(tail,edge,minlength=len(d))+np.bincount(head,edge,minlength=len(d))

    def thermal_data(self,d,p,temp=None,alpha=None):
        d,p,alpha=self._state(d,p,alpha);T=self._temperature(d,p,temp)
        rho=abs(d)**2;en=self.energies_bar(d);theta=T/self.Tc_K
        fd=expit(-en/theta[:,None]);u=self.node_energy_scale_J;w=self.weights_bar
        C=4*u*np.sum(w*en**2*fd*(1-fd),axis=1)/(self.Tc_K*theta**2)
        B=2*u*np.sum(w*((p-fd)/en+fd*(1-fd)/theta[:,None]),axis=1)
        L=self._L(d,alpha)
        r=thermal_rigidity(rho,T,Tc_K=self.Tc_K,prefactor=self.radial_prefactor)
        return dict(L_J=L,kappa=r.kappa,kappa_rho=r.rho,kappa_rho_rho=r.rho_rho,
            kappa_T_per_K=r.temperature,kappa_TT_per_K2=r.temperature_temperature,
            kappa_rho_T_per_K=r.rho_temperature,capacity_FD_J_K=C,B_J=B,
            capacity_effective_J_K=C-T*L*r.temperature_temperature,
            epsilon_J=4*u[:,None]*w*en,temperature_K=T,
            G_T_J_per_K=L*r.temperature,G_TT_J_per_K2=L*r.temperature_temperature,
            G_Trho_J_per_K=L*r.rho_temperature,fd=fd)

    def evaluate(self,gapbar,p,alpha=None,*,temperature_K=None):
        d,p,alpha=self._state(gapbar,p,alpha);T=self._temperature(d,p,temperature_K)
        raw=self._fixed_T(T).evaluate(d,p,alpha);td=self.thermal_data(d,p,T,alpha)
        rho=abs(d)**2;en=raw.electron_energies_J/self.Eu_J
        local=self.node_energy_scale_J*(rho*self.inverse_pairing_coupling-
            2*np.sum(self.weights_bar*rho[:,None]/(en+self.xi_bar),axis=1)+
            4*np.sum(self.weights_bar*en*p,axis=1))
        gradient_node=td["L_J"]*(td["kappa"]-T*td["kappa_T_per_K"])
        node_energy=local+gradient_node
        return ThermalStiffnessEvaluation(float(node_energy.sum()),raw.local_energy_J,
            float(gradient_node.sum()),node_energy,raw.electron_energies_J,
            raw.gap_force_J_per_gapbar,raw.gap_gradient_J_per_gapbar,
            raw.edge_current_A,raw.current_divergence_A,raw.link_derivative_alpha_J,
            raw.local_coefficient_bar,raw.gradient_kappa,raw.gradient_energy_J,
            -td["G_T_J_per_K"],td["G_T_J_per_K"],td["G_TT_J_per_K2"],
            td["L_J"],T,"Fthermo=U_d-sum(T*S_d)=d(U_BCS+G)/dd* at fixed TE; energy is U")

    def hessian_and_current_jacobian(self,gapbar,p,alpha=None,*,temperature_K=None):
        d,p,alpha=self._state(gapbar,p,alpha);T=self._temperature(d,p,temperature_K)
        return self._fixed_T(T).hessian_and_current_jacobian(d,p,alpha)

    def direction(self,gapbar,p,gap_direction,p_direction=None,alpha=None,
                  alpha_direction=None,*,temperature_K=None,temperature_direction_K=None):
        d,p,alpha=self._state(gapbar,p,alpha);T=self._temperature(d,p,temperature_K)
        dd=np.asarray(gap_direction,complex);dp=np.zeros_like(p)if p_direction is None else np.asarray(p_direction,float)
        da=np.zeros_like(alpha)if alpha_direction is None else np.asarray(alpha_direction,float)
        td=self.thermal_data(d,p,T,alpha);drho=2*np.real(np.conj(d)*dd)
        if temperature_direction_K is None:
            dT=(td["B_J"]*drho+np.sum(td["epsilon_J"]*dp,axis=1))/td["capacity_FD_J_K"]
        else:dT=np.broadcast_to(np.asarray(temperature_direction_K,float),T.shape)
        fixed=self._fixed_T(T).direction(d,p,dd,dp,alpha,da)
        value=self.evaluate(d,p,alpha,temperature_K=T)
        Ldot=self.gradient_rate_geometry(d,dd,alpha,da)
        GTdot=td["kappa_T_per_K"]*Ldot+td["G_Trho_J_per_K"]*drho+td["G_TT_J_per_K2"]*dT
        du=(2*np.real(np.vdot(value.gap_force_J_per_gapbar,dd))+
            np.dot(value.link_derivative_alpha_J,da)+np.sum(td["epsilon_J"]*dp)-np.dot(T,GTdot))
        tail,head=self.graph.edges.T;phase=np.exp(-1j*alpha)
        diff=d[tail]-phase*d[head];diff2=abs(diff)**2;c=self.edge_energy_scale_J
        dk=td["kappa_T_per_K"]*dT;dkp=td["kappa_rho_T_per_K"]*dT;dke=(dk[tail]+dk[head])/2
        df=fixed.gap_force_direction_J_per_gapbar.copy()
        np.add.at(df,tail,c*(dke*diff+.5*dkp[tail]*d[tail]*diff2))
        np.add.at(df,head,c*(-dke*np.conj(phase)*diff+.5*dkp[head]*d[head]*diff2))
        di=fixed.edge_current_direction_A+4*electron_charge/hbar*c*dke*np.imag(np.conj(d[tail])*phase*d[head])
        return EffectiveBCSDirection(float(du),df,di,fixed.electron_energy_direction_J)

    def bulk_coefficient(self,amplitude_bar,temperature_K,q_bar=0.,phase_step_m=None):
        return self._fixed_T(temperature_K).bulk_coefficient(amplitude_bar,temperature_K,q_bar,phase_step_m)

    def bulk_current(self,temperature_K,q_bar,width_m,phase_step_m=None):
        return self._fixed_T(temperature_K).bulk_current(temperature_K,q_bar,width_m,phase_step_m)

    def bulk(self,*,temperature_K,current_A,width_m,phase_step_m=None):
        return self._fixed_T(temperature_K).bulk(temperature_K=temperature_K,current_A=current_A,
                                                width_m=width_m,phase_step_m=phase_step_m)
