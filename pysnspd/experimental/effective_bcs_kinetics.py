"""Neutral native-count BCS kinetics in physical SI units.

This is the explicit effective alternative to the rejected ultrafast Moyal
closure. Fixed xi states have E=sqrt(xi**2+|Delta|**2). Shared scattering and
recombination events conserve native excitation+phonon energy exactly and
obey FD/BE balance. Transport is at equal physical E, with barycentric
stoichiometry and geometric Pauli activities; it does not diffuse xi labels.
No Usadel root, population clipping, DOS floor or energy repair is used.
"""
from __future__ import annotations
from dataclasses import dataclass
import os
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.constants import Boltzmann as kB, hbar
from scipy.special import expit
from .kinetic_events import PhononGrid, _signed_exponential_difference
try:
    from . import effective_bcs_kinetics_numba as _fast
except ImportError:
    _fast=None
_fused=_fast


@dataclass(frozen=True)
class QuadraticAcousticCoupling:
    """Exact declared acoustic alpha2F(omega_bar)=coefficient*omega_bar**2."""
    coefficient: float

    def __call__(self,omega):
        return self.coefficient*np.asarray(omega)**2



def _array(value, name, shape=None, positive=False):
    raw=np.asarray(value)
    if np.iscomplexobj(raw): raise ValueError(name+' must be real')
    a=np.array(raw,float,copy=True)
    if (shape is not None and a.shape!=shape) or np.any(~np.isfinite(a)):
        raise ValueError(name+' has invalid shape or nonfinite entries')
    if positive and np.any(a<=0): raise ValueError(name+' must be positive')
    return a


def _log(a):
    out=np.full_like(a,-np.inf,float)
    np.log(a,out=out,where=a>0)
    return out


def _weighted_pair_logs(logs, lower, upper, beta):
    a=np.take_along_axis(logs,lower,axis=1)
    b=np.take_along_axis(logs,upper,axis=1)
    out=np.zeros_like(a)
    np.multiply(1-beta,a,out=out,where=beta<1)
    term=np.zeros_like(b)
    np.multiply(beta,b,out=term,where=beta>0)
    return out+term


@dataclass(frozen=True)
class EffectiveBCSKineticResult:
    dp_dt_s: np.ndarray
    dn_dt_s: np.ndarray
    temperature_K: np.ndarray
    electron_energies_J: np.ndarray
    node_phonon_energy_J: np.ndarray
    phonon_energy_J: float
    electron_population_power_W: np.ndarray
    phonon_power_W: np.ndarray
    collision_dp_dt_s: np.ndarray
    collision_dn_dt_s: np.ndarray
    transport_dp_dt_s: np.ndarray
    bgk_dp_dt_s: np.ndarray
    heating_dp_dt_s: np.ndarray
    escape_dn_dt_s: np.ndarray
    collision_electron_power_W: np.ndarray
    collision_phonon_power_W: np.ndarray
    ep_energy_residual_W: np.ndarray
    transport_face_power_W: np.ndarray
    transport_node_power_W: np.ndarray
    transport_energy_residual_W: float
    transport_number_residual_per_s: float
    bgk_energy_residual_W: np.ndarray
    heating_power_W: np.ndarray
    heating_residual_W: np.ndarray
    escape_power_W: np.ndarray
    balance_residual_W: float
    forward_euler_time_limit_s: float
    limited_event_energy_J: float=0.
    limited_event_fraction: float=0.
    totalabs_event_energy_J: float=0.
    max_flux_reduction: float=0.
    equilibrium_collision_nodes: int=0


class EffectiveBCSKinetics:
    """Vectorized neutral p_xi/n_Omega material and transport RHS.

    xi_bar and weights_bar use Eu=kB*Tc; capacities in supplied PhononGrid
    are (gph/N0)*dOmega_J/Eu. Collision prefactor is 8*pi*Eu/hbar per second.
    tau_ee_Tc_s is a free calibrator: tau_BGK(T)=tau_ee_Tc*Tc/max(TE,Tb).
    infinity disables BGK/escape. No material relaxation time is inferred.
    """
    def __init__(self, *, xi_bar, weights_bar, phonons, Eu_J,
                 N0_per_J_m3, D_m2_s, volumes_m3, edges, face_areas_m2,
                 edge_lengths_m, bath_temperature_K, tau_ee_Tc_s=5e-12,
                 tau_escape_s=15e-12, transport_order=48, accelerated=True,numba_threads=None,recognize_equilibrium=True):
        self.accelerated=bool(accelerated and _fast is not None)
        self.recognize_equilibrium=bool(recognize_equilibrium)
        self.numba_threads=int(os.environ.get('NUMBA_NUM_THREADS','1') if numba_threads is None else numba_threads)
        if self.numba_threads<1:raise ValueError('Positive explicit numba thread count required')
        if self.accelerated:_fast.set_num_threads(self.numba_threads)
        self.xi_bar=_array(xi_bar,'xi_bar',positive=True)
        self.weights_bar=_array(weights_bar,'weights_bar',self.xi_bar.shape,True)
        if self.xi_bar.ndim!=1 or len(self.xi_bar)<2 or np.any(np.diff(self.xi_bar)<=0):
            raise ValueError('Increasing positive native xi quadrature required')
        self.phonons=phonons
        if not isinstance(phonons,PhononGrid):raise TypeError('Explicit PhononGrid required')
        for key,val in [('Eu_J',Eu_J),('N0_per_J_m3',N0_per_J_m3),
                        ('D_m2_s',D_m2_s),('bath_temperature_K',bath_temperature_K)]:
            if not np.isfinite(val) or val<=0:raise ValueError(key+' must be positive finite')
            setattr(self,key,float(val))
        self.Tc_K=self.Eu_J/kB
        self.volumes_m3=_array(volumes_m3,'volumes_m3',positive=True)
        if self.volumes_m3.ndim!=1:raise ValueError('One volume per node required')
        self.node_count=len(self.volumes_m3);self.electron_size=len(self.xi_bar)
        self.phonon_size=len(phonons.energies)
        ed=np.asarray(edges)
        if ed.size==0:ed=np.empty((0,2),int)
        if (ed.ndim!=2 or ed.shape[1]!=2 or ed.dtype.kind not in 'iu'
            or np.any(ed<0) or np.any(ed>=self.node_count) or np.any(ed[:,0]==ed[:,1])):
            raise ValueError('Integer oriented graph edges required')
        if len(np.unique(np.sort(ed,axis=1),axis=0))!=len(ed):
            raise ValueError('Duplicate transport faces')
        self.edges=ed.astype(int,copy=True)
        self.face_areas_m2=_array(face_areas_m2,'face_areas_m2',(len(ed),))
        self.edge_lengths_m=_array(edge_lengths_m,'edge_lengths_m',(len(ed),),True)
        if np.any(self.face_areas_m2<0):raise ValueError('Face areas must be nonnegative')
        self.face_rates_m3_s=self.D_m2_s*self.face_areas_m2/self.edge_lengths_m
        for name,value in [('tau_ee_Tc_s',tau_ee_Tc_s),('tau_escape_s',tau_escape_s)]:
            a=np.broadcast_to(np.asarray(value,float),(self.node_count,)).copy()
            if np.any(np.isnan(a)) or np.any(a<=0):raise ValueError(name+' must be positive seconds or infinity')
            setattr(self,name,a)
        if isinstance(transport_order,bool) or int(transport_order)!=transport_order or transport_order<2:
            raise ValueError('transport_order must be integer >=2')
        x,w=leggauss(int(transport_order));self.face_rule=(x,w)
        si,sj=np.triu_indices(self.electron_size,1)
        ri,rj=np.triu_indices(self.electron_size,0)
        self.event_i=np.r_[si,ri];self.event_j=np.r_[sj,rj]
        self.recombination=np.r_[np.zeros(len(si),bool),np.ones(len(ri),bool)]
        pairweights=self.weights_bar[self.event_i]*self.weights_bar[self.event_j]
        pairweights[self.recombination & (self.event_i==self.event_j)]*=.5
        self.pairweights=pairweights
        self._deposit_i=self._deposit_j=None
        if not self.accelerated:
            nodes=np.arange(self.node_count)[:,None]
            self._deposit_i=nodes*self.electron_size+self.event_i
            self._deposit_j=nodes*self.electron_size+self.event_j
        self.node_energy_scale_J=self.N0_per_J_m3*self.Eu_J**2*self.volumes_m3
        self.phonon_energy_weights_J=(self.node_energy_scale_J[:,None]*
                phonons.capacities*phonons.energies)
        self.bath_phonons=self.bose(self.bath_temperature_K)

    @classmethod
    def from_graph(cls,graph,*,ell0_m,thickness_m,**kwargs):
        edge_lengths=np.linalg.norm(graph.coordinates_bar[graph.edges[:,1]]-
            graph.coordinates_bar[graph.edges[:,0]],axis=1)*ell0_m
        return cls(volumes_m3=graph.area_weights*ell0_m**2*thickness_m,
            edges=graph.edges,edge_lengths_m=edge_lengths,
            face_areas_m2=graph.conductance*edge_lengths*thickness_m,**kwargs)

    def metadata(self):
        return dict(model='neutral native-count local BCS effective kinetics',accelerated=self.accelerated,numba_threads=self.numba_threads,
            equilibrium_recognition_relative_epsilon=32 if self.recognize_equilibrium else None,
            Eu_J=self.Eu_J,Tc_K=self.Tc_K,electron_states=self.electron_size,
            phonon_modes=self.phonon_size,transport_order=len(self.face_rule[0]),
            rate_prefactor_per_s=8*np.pi*self.Eu_J/hbar,
            tau_ee_Tc_s=self.tau_ee_Tc_s.tolist(),
            bgk_temperature_law='tau_ee(T)=tau_ee(Tc)*Tc/max(TE,Tb)',
            tau_escape_s=self.tau_escape_s.tolist(),
            energy_convention='N0*Eu^2*V*(4 sum w Ebar p + sum cph Omegabar n)',
            phonon_capacities='(gph/N0)*dOmega_J/Eu',
            spectral_approximation='BCS local neutral; D_L=1 above local gap',
            no_population_clipping=True,no_Usadel_queries=True)

    def energies_bar(self,delta_bar):
        d=np.asarray(delta_bar,complex)
        if d.shape!=(self.node_count,) or np.any(~np.isfinite(d)):
            raise ValueError('One finite gap per node required')
        return np.sqrt(self.xi_bar[None,:]**2+abs(d[:,None])**2)

    def fermi_dirac(self,delta_bar,temperature_K):
        t=np.broadcast_to(np.asarray(temperature_K,float),(self.node_count,))
        if np.any(~np.isfinite(t)) or np.any(t<0):raise ValueError('Nonnegative finite temperature required')
        e=self.energies_bar(delta_bar);out=np.zeros_like(e)
        np.divide(-e,kB*t[:,None]/self.Eu_J,out=out,where=t[:,None]>0)
        return np.where(t[:,None]>0,expit(out),0.)

    def bose(self,temperature_K):
        t=float(temperature_K)
        if not np.isfinite(t) or t<0:raise ValueError('Nonnegative finite temperature required')
        if t==0:return np.zeros(self.phonon_size)
        x=self.phonons.energies*self.Eu_J/(kB*t)
        return np.exp(-x)/(-np.expm1(-x))

    def _validate(self,delta_bar,p,n):
        e=self.energies_bar(delta_bar)
        p=_array(p,'p',(self.node_count,self.electron_size))
        n=_array(n,'n',(self.node_count,self.phonon_size))
        if np.any((p<0)|(p>1)) or np.any(n<0):
            raise ValueError('Native occupations leave Pauli/phonon support; no clipping')
        return e,p,n

    def _temperature(self,e,p):
        target=np.sum(self.weights_bar*e*p,axis=1)
        upperlimit=.5*np.sum(self.weights_bar*e,axis=1)
        if np.any(target>=upperlimit):raise ValueError('Population energy has no finite positive FD temperature')
        lo=np.zeros(self.node_count);hi=np.ones(self.node_count)
        def moments(t):
            f=expit(-e/t[:,None])
            value=np.sum(self.weights_bar*e*f,axis=1)
            derivative=np.sum(self.weights_bar*e*e*f*(1-f),axis=1)/t**2
            return value,derivative
        for _ in range(128):
            need=moments(hi)[0]<target
            if not np.any(need):break
            hi[need]*=2
        else:raise FloatingPointError('No finite energy-equivalent temperature bracket')
        previous=getattr(self,'_temperature_guess',np.full(self.node_count,kB*self.bath_temperature_K/self.Eu_J))
        theta=np.minimum(previous,hi)
        theta=np.where(theta>0,theta,hi/2)
        tolerance=64*np.finfo(float).eps*np.maximum(target,np.finfo(float).tiny)
        for _ in range(64):
            value,derivative=moments(theta)
            converged=(abs(value-target)<=tolerance)|(target==0)
            if np.all(converged):break
            below=value<target
            lo=np.where(below&~converged,theta,lo)
            hi=np.where((~below)&~converged,theta,hi)
            correction=np.divide(value-target,derivative,out=np.full_like(target,np.inf),where=derivative>0)
            candidate=theta-correction
            safe=np.isfinite(candidate)&(candidate>lo)&(candidate<hi)
            candidate=np.where(safe,candidate,(lo+hi)/2)
            theta=np.where(converged,theta,candidate)
        else:raise FloatingPointError('Energy-equivalent FD inversion did not converge')
        theta=np.where(target==0,0.,theta)
        self._temperature_guess=np.where(theta>0,theta,kB*self.bath_temperature_K/self.Eu_J)
        return theta*self.Eu_J/kB

    def equivalent_temperature_K(self,delta_bar,p):
        e=self.energies_bar(delta_bar)
        p=_array(p,'p',e.shape)
        if np.any((p<0)|(p>1)):raise ValueError('Pauli support required')
        return self._temperature(e,p)

    def energy_J(self,delta_bar,p,n):
        e,p,n=self._validate(delta_bar,p,n)
        return float(np.sum(self.node_energy_scale_J[:,None]*4*self.weights_bar*e*p)
            +np.sum(self.phonon_energy_weights_J*n))

    def _collisions(self,delta,e,p,n,return_events=False,equilibrium_mask=None):
        i,j,rec=self.event_i,self.event_j,self.recombination
        if equilibrium_mask is None:equilibrium_mask=np.zeros(self.node_count,bool)
        if self.accelerated and isinstance(self.phonons.alpha2F,QuadraticAcousticCoupling):
            dp,dn,rate,omega,lower,upper,beta=_fast.acoustic_collisions(
                e,abs(delta),p,n,_log(p),_log(1-p),_log(n),np.log1p(n),
                i,j,rec,self.pairweights,self.phonons.energies,self.phonons.capacities,
                self.weights_bar,self.phonons.alpha2F.coefficient,*self.phonons.coupling_support,
                8*np.pi*self.Eu_J/hbar,equilibrium_mask)
            return (dp,dn,dict(rate=rate,omega=omega,lower=lower,upper=upper,beta=beta)) if return_events else (dp,dn)
        ei,ej=e[:,i],e[:,j]
        omega=np.where(rec,ei+ej,ej-ei)
        ratio=abs(delta[:,None])**2/(ei*ej)
        coherence=np.where(rec,1+ratio,1-ratio)
        coeff=(8*np.pi*self.Eu_J/hbar)*self.pairweights*self.phonons.coupling(omega)*coherence
        coeff[equilibrium_mask]=0.
        active=coeff>0
        if np.any(active & ((omega<self.phonons.energies[0])|
                            (omega>self.phonons.energies[-1]))):
            raise ValueError('Coupled event lacks phonon bracket; extend physical quadrature')
        upper=np.searchsorted(self.phonons.energies,omega,side='left')
        # Inactive compact-support events have zero physical rate; these indices
        # only make array gathers legal and do not modify coupled events.
        upper=np.where(active,np.maximum(upper,1),1)
        lower=np.maximum(upper-1,0)
        beta=(omega-self.phonons.energies[lower])/(self.phonons.energies[upper]-self.phonons.energies[lower])
        beta=np.where(active,beta,0.)
        emission=_weighted_pair_logs(np.log1p(n),lower,upper,beta)
        absorption=_weighted_pair_logs(_log(n),lower,upper,beta)
        lp,lh=_log(p),_log(1-p)
        forward=np.where(rec,lp[:,i]+lp[:,j],lh[:,i]+lp[:,j])+emission
        reverse=np.where(rec,lh[:,i]+lh[:,j],lp[:,i]+lh[:,j])+absorption
        rate=coeff*_signed_exponential_difference(forward,reverse)
        result=self._collision_deposit(rate,lower,upper,beta,e.shape,n.shape)
        return (*result,dict(rate=rate,lower=lower,upper=upper,beta=beta,omega=omega)) if return_events else result

    def _collision_deposit(self,rate,lower,upper,beta,eshape,nshape):
        if self._deposit_i is None:
            nodes=np.arange(self.node_count)[:,None]
            self._deposit_i=nodes*self.electron_size+self.event_i
            self._deposit_j=nodes*self.electron_size+self.event_j
        rec=self.recombination
        ne=self.node_count*self.electron_size
        dp=(np.bincount(self._deposit_i.ravel(),np.where(rec,-rate,rate).ravel(),minlength=ne)
            +np.bincount(self._deposit_j.ravel(),(-rate).ravel(),minlength=ne)).reshape(eshape)
        dp/=4*self.weights_bar
        offsets=np.arange(self.node_count)[:,None]*self.phonon_size
        dn=(np.bincount((offsets+lower).ravel(),((1-beta)*rate).ravel(),
                minlength=self.node_count*self.phonon_size)
            +np.bincount((offsets+upper).ravel(),(beta*rate).ravel(),
                minlength=self.node_count*self.phonon_size)).reshape(nshape)
        dn/=self.phonons.capacities
        return dp,dn

    def _transport(self,delta,e,p,return_events=False):
        nf=len(self.edges)
        if not nf:
            result=(np.zeros_like(p),np.empty(0),np.empty(0))
            return (*result,None) if return_events else result
        if self.accelerated:
            dp,power,number,flux,q,lower,upper,beta,lp,lh=_fast.transport(
                e,abs(delta),_log(p),_log(1-p),self.xi_bar,self.weights_bar,self.edges,
                self.volumes_m3,self.face_rates_m3_s,*self.face_rule,self.N0_per_J_m3,self.Eu_J)
            data=[(lower[k],upper[k],beta[k],lp[k],lh[k])for k in range(2)]
            return (dp,power,number,dict(flux=flux,q=q,data=data)) if return_events else (dp,power,number)
        left,right=self.edges.T
        low=np.maximum(e[left,0],e[right,0])
        high=np.minimum(e[left,-1],e[right,-1])
        if np.any(high<=low):raise ValueError('No shared electronic energy support at a face')
        x,w=self.face_rule
        q=(low[:,None]+high[:,None])/2+(high-low)[:,None]*x/2
        measure=(high-low)[:,None]*w/2
        data=[]
        for nodes in (left,right):
            a=abs(delta[nodes,None])
            xi=np.sqrt((q-a)*(q+a))
            upper=np.searchsorted(self.xi_bar,xi,side='right')
            if np.any((upper<1)|(upper>=self.electron_size)):
                raise FloatingPointError('Transport quadrature lacks native brackets')
            lower=upper-1
            elo=np.sqrt(self.xi_bar[lower]**2+a*a)
            ehi=np.sqrt(self.xi_bar[upper]**2+a*a)
            beta=(q-elo)/(ehi-elo)
            if np.any((beta<0)|(beta>1)):raise FloatingPointError('Invalid transport barycentric activity')
            data.append((lower,upper,beta,
                _weighted_pair_logs(_log(p[nodes]),lower,upper,beta),
                _weighted_pair_logs(_log(1-p[nodes]),lower,upper,beta)))
        l,r=data
        flux=4*self.face_rates_m3_s[:,None]*measure*_signed_exponential_difference(l[3]+r[4],r[3]+l[4])
        result=self._transport_deposit(flux,q,data,e.shape)
        return (*result,dict(flux=flux,q=q,data=data)) if return_events else result

    def _transport_deposit(self,flux,q,data,eshape):
        left,right=self.edges.T
        l,r=data
        count=np.zeros(self.node_count*self.electron_size)
        for nodes,(lower,upper,beta,_,_),sign in ((left,l,-1),(right,r,1)):
            offsets=nodes[:,None]*self.electron_size
            count+=np.bincount((offsets+lower).ravel(),(sign*(1-beta)*flux).ravel(),minlength=len(count))
            count+=np.bincount((offsets+upper).ravel(),(sign*beta*flux).ravel(),minlength=len(count))
        dp=count.reshape(eshape)/(4*self.weights_bar*self.volumes_m3[:,None])
        power=self.N0_per_J_m3*self.Eu_J**2*np.sum(q*flux,axis=1)
        number=self.N0_per_J_m3*self.Eu_J*np.sum(flux,axis=1)
        return dp,power,number

    def evaluate(self,delta_bar,p_nodes,n_nodes,heating_power_W=None,temperature_K=None):
        return self._evaluate(delta_bar,p_nodes,n_nodes,heating_power_W,
            temperature_K=temperature_K)

    def evaluate_limited(self,delta_bar,p_nodes,n_nodes,*,budget_p,budget_n,
                         dt_s,heating_power_W=None,temperature_K=None):
        """Finite-step production/destruction limiter on SHARED event extents.

        Availability of old-state particles, holes and phonons bounds outgoing
        event fluxes. Every event retains exactly its common stoichiometry and
        energy. Populations are never clipped or repaired. This is a first-order
        numerical time treatment, not a rate floor or modified continuous RHS;
        limited energy and flux fractions must be checked by step refinement.
        BGK/escape use exponential-fitted midpoint coefficients 2*tanh(dt/(2*tau))/dt. Their full energy
        moments and the complete supplied heat are retained before event limits.
        """
        if not np.isfinite(dt_s) or dt_s<=0:raise ValueError('Positive physical dt_s required')
        _,oldp,oldn=self._validate(delta_bar,budget_p,budget_n)
        return self._evaluate(delta_bar,p_nodes,n_nodes,heating_power_W,
            budget=(oldp,oldn,float(dt_s)),temperature_K=temperature_K)

    def _evaluate(self,delta_bar,p_nodes,n_nodes,heating_power_W=None,budget=None,temperature_K=None):
        delta=np.asarray(delta_bar,complex)
        e,p,n=self._validate(delta,p_nodes,n_nodes)
        temperature=(self._temperature(e,p) if temperature_K is None else
            _array(temperature_K,'temperature_K',(self.node_count,)))
        if np.any(temperature<0):raise ValueError('Nonnegative energy-equivalent temperatures required')
        equilibrium_mask=self._equilibrium_mask(e,p,n,temperature)
        native_node_path=self.accelerated and isinstance(self.phonons.alpha2F,QuadraticAcousticCoupling)
        if native_node_path:
            collision_dp=collision_dn=collision_events=None
        else:
            collision_dp,collision_dn,collision_events=self._collisions(delta,e,p,n,True,equilibrium_mask)
        transport_dp,face_power,face_number,transport_events=self._transport(delta,e,p,True)
        theta=kB*np.maximum(temperature,self.bath_temperature_K)/self.Eu_J
        target=expit(-e/theta[:,None])
        tau=self.tau_ee_Tc_s*self.Tc_K/np.maximum(temperature,self.bath_temperature_K)
        # BGK target uses TE (not Tmob); bath floor affects the rate only.
        physical_theta=kB*temperature/self.Eu_J
        bgk_target=np.zeros_like(e)
        np.divide(-e,physical_theta[:,None],out=bgk_target,where=physical_theta[:,None]>0)
        bgk_target=np.where(physical_theta[:,None]>0,expit(bgk_target),0.)
        bgk=(bgk_target-p)/tau[:,None]
        if budget is not None:
            dt=budget[2]
            bgk=(bgk_target-p)*(2*np.tanh(dt/(2*tau[:,None])))/dt
        source=np.zeros(self.node_count) if heating_power_W is None else _array(heating_power_W,'heating_power_W',(self.node_count,))
        if np.any(source<0):raise ValueError('D18 supplied heat must be nonnegative')
        shape=e*target*(1-p)
        normalizer=4*np.sum(self.weights_bar*e*shape,axis=1)*self.node_energy_scale_J
        if np.any((normalizer<=0)&(source>0)):raise FloatingPointError('No available states for supplied heating')
        factor=np.divide(source,normalizer,out=np.zeros_like(source),where=normalizer>0)
        heating=shape*factor[:,None]
        escape=-(n-self.bath_phonons)/self.tau_escape_s[:,None]
        limiter_metrics=dict(limited_event_energy_J=0.,limited_event_fraction=0.,
            totalabs_event_energy_J=0.,max_flux_reduction=0.)
        if budget is None and native_node_path:
            collision_dp,collision_dn=_fused.acoustic_unlimited(*self._node_arguments(delta,e,p,n,equilibrium_mask))
        if budget is not None:
            oldp,oldn,dt=budget
            escape=(self.bath_phonons-n)*(2*np.tanh(dt/(2*self.tau_escape_s[:,None])))/dt
            ap=np.tanh(dt/(2*tau[:,None]))
            an=np.tanh(dt/(2*self.tau_escape_s[:,None]))
            particle_budget=4*self.weights_bar*((1-ap)*oldp+2*ap*bgk_target+dt*heating)/dt
            hole_budget=4*self.weights_bar*((1-ap)*(1-oldp)+2*ap*(1-bgk_target)-dt*heating)/dt
            phonon_budget=self.phonons.capacities*((1-an)*oldn+2*an*self.bath_phonons)/dt
            if native_node_path:
                collision_dp,collision_dn,transport_dp,face_power,face_number,limiter_metrics=(
                    self._node_limited(delta,e,p,n,equilibrium_mask,transport_events,
                        dt,particle_budget,hole_budget,phonon_budget))
            else:
                collision_dp,collision_dn,transport_dp,face_power,face_number,limiter_metrics=(
                    self._limit_events(collision_events,transport_events,e.shape,n.shape,
                        oldp,oldn,dt,particle_budget,hole_budget,phonon_budget))
        dp=collision_dp+transport_dp+bgk+heating;dn=collision_dn+escape
        electron_weights=4*self.node_energy_scale_J[:,None]*self.weights_bar*e
        emoment=lambda a:np.sum(electron_weights*a,axis=1)
        pmoment=lambda a:np.sum(self.phonon_energy_weights_J*a,axis=1)
        ce,cp=emoment(collision_dp),pmoment(collision_dn)
        transport_power=emoment(transport_dp)
        transport_number=float(np.sum(4*self.N0_per_J_m3*self.Eu_J*
            self.volumes_m3[:,None]*self.weights_bar*transport_dp))
        bgk_residual=emoment(bgk);heat_residual=emoment(heating)-source
        escape_power=-pmoment(escape)
        erate,prate=emoment(dp),pmoment(dn)
        if np.any(~np.isfinite(dp)) or np.any(~np.isfinite(dn)):
            raise FloatingPointError('Nonfinite physical BCS kinetic RHS')
        phenergy=pmoment(n)
        return EffectiveBCSKineticResult(dp,dn,temperature,e*self.Eu_J,phenergy,float(np.sum(phenergy)),
            erate,prate,collision_dp,collision_dn,transport_dp,bgk,heating,escape,ce,cp,ce+cp,
            face_power,transport_power,float(np.sum(transport_power)),transport_number,
            bgk_residual,source.copy(),heat_residual,escape_power,
            float(np.sum(erate+prate+escape_power-source)),
            forward_euler_time_limit_s(p if budget is None else budget[0],
                n if budget is None else budget[1],dp,dn),equilibrium_collision_nodes=int(np.count_nonzero(equilibrium_mask)),**limiter_metrics)


    def _node_arguments(self,delta,e,p,n,equilibrium_mask):
        return (e,abs(delta),p,n,_log(p),_log(1-p),_log(n),np.log1p(n),
            self.event_i,self.event_j,self.recombination,self.pairweights,
            self.phonons.energies,self.phonons.capacities,self.weights_bar,
            self.phonons.alpha2F.coefficient,*self.phonons.coupling_support,
            8*np.pi*self.Eu_J/hbar,equilibrium_mask)

    def _node_limited(self,delta,e,p,n,equilibrium_mask,te,dt,bp,bh,bn):
        if te is None:
            flux=np.empty((0,0));q=np.empty((0,0))
            tlo=np.empty((2,0,0),np.int64);tup=tlo.copy();tbeta=np.empty((2,0,0))
        else:
            flux,q=te['flux'],te['q']
            tlo=np.array([r[0]for r in te['data']]);tup=np.array([r[1]for r in te['data']])
            tbeta=np.array([r[2]for r in te['data']])
        if not hasattr(self,'_node_incidence'):
            self._node_incidence=_fused.ordered_incidence(self.edges,self.node_count)
        cdp,cdn,fp,fh,changed,total,count,active,reduction=_fused.acoustic_limited_nodes(
            *self._node_arguments(delta,e,p,n,equilibrium_mask),
            flux,tlo,tup,tbeta,self.volumes_m3,self.node_energy_scale_J,
            dt,bp,bh,bn,*self._node_incidence)
        values=_fused.apply_transport(flux,q,tlo,tup,tbeta,self.edges,self.volumes_m3,
            self.weights_bar,self.N0_per_J_m3,self.Eu_J,dt,
            fp,fh,changed,total,count,active,reduction)
        return cdp,cdn,*values[:3],dict(limited_event_energy_J=values[3],
            totalabs_event_energy_J=values[4],limited_event_fraction=values[5],
            max_flux_reduction=values[6])

    def _equilibrium_mask(self,e,p,n,temperature):
        """FD/BE detailed-balance identity recognized only at IEEE accuracy.

        Relative tests have NO absolute occupation floor. Transport remains
        active at every node, and any resolved nonthermal deviation activates
        the complete collision operator immediately.
        """
        if not self.recognize_equilibrium:return np.zeros(self.node_count,bool)
        epsilon=32*np.finfo(float).eps
        def equal(a,b):
            return np.all(abs(a-b)<=epsilon*np.maximum(abs(a),abs(b)),axis=1)
        theta=kB*temperature/self.Eu_J
        fd=np.zeros_like(e)
        np.divide(-e,theta[:,None],out=fd,where=theta[:,None]>0)
        fd=np.where(theta[:,None]>0,expit(fd),0.)
        x=np.zeros_like(n)
        np.divide(self.phonons.energies,theta[:,None],out=x,where=theta[:,None]>0)
        be=np.zeros_like(n)
        positive=theta[:,None]>0
        np.divide(np.exp(-x),-np.expm1(-x),out=be,where=positive)
        common=equal(p,fd)&equal(n,be)
        bathfd=expit(-e/(kB*self.bath_temperature_K/self.Eu_J))
        return common|(equal(p,bathfd)&equal(n,self.bath_phonons))

    def _limit_events(self,ce,te,eshape,nshape,oldp,oldn,dt,particle_budget,hole_budget,phonon_budget):
        """Conservative outgoing-flux correction using native support budgets."""
        if self.accelerated:
            if te is None:
                flux=np.empty((0,0));q=np.empty((0,0))
                tlo=np.empty((2,0,0),np.int64);tup=tlo.copy();tbeta=np.empty((2,0,0))
            else:
                flux,q=te['flux'],te['q']
                tlo=np.array([r[0]for r in te['data']]);tup=np.array([r[1]for r in te['data']])
                tbeta=np.array([r[2]for r in te['data']])
            values=_fast.limit_events(ce['rate'],ce['omega'],ce['lower'],ce['upper'],ce['beta'],
                self.event_i,self.event_j,self.recombination,flux,q,tlo,tup,tbeta,
                self.edges,self.volumes_m3,self.weights_bar,self.phonons.capacities,
                self.node_energy_scale_J,self.N0_per_J_m3,self.Eu_J,oldp,oldn,dt,particle_budget,hole_budget,phonon_budget)
            return (*values[:5],dict(limited_event_energy_J=values[5],
                totalabs_event_energy_J=values[6],limited_event_fraction=values[7],
                max_flux_reduction=values[8]))
        ne=self.node_count*self.electron_size
        npop=self.node_count*self.phonon_size
        loss_p=np.zeros(ne);loss_h=np.zeros(ne);loss_n=np.zeros(npop)
        rate=ce['rate'];rec=self.recombination
        def electron_loss(indices,component):
            nonlocal loss_p,loss_h
            loss_p+=np.bincount(indices.ravel(),np.maximum(-component,0).ravel(),minlength=ne)
            loss_h+=np.bincount(indices.ravel(),np.maximum(component,0).ravel(),minlength=ne)
        electron_loss(self._deposit_i,np.where(rec,-rate,rate))
        electron_loss(self._deposit_j,-rate)
        offsets=np.arange(self.node_count)[:,None]*self.phonon_size
        for ix,beta in ((ce['lower'],1-ce['beta']),(ce['upper'],ce['beta'])):
            loss_n+=np.bincount((offsets+ix).ravel(),np.maximum(-beta*rate,0).ravel(),minlength=npop)
        if te is not None:
            left,right=self.edges.T
            flux=te['flux']
            for nodes,(lower,upper,beta,_,_),sign in ((left,te['data'][0],-1),(right,te['data'][1],1)):
                for ix,weight in ((lower,1-beta),(upper,beta)):
                    electron_loss(nodes[:,None]*self.electron_size+ix,
                        sign*weight*flux/self.volumes_m3[nodes,None])
        bp,bh,bn=particle_budget,hole_budget,phonon_budget
        # Exact midpoint elimination supplies these positive production budgets:
        # p_new=[(1-a)p_old+2a fFD+h H+h C]/(1+a), a=tanh(h/(2 tau)).
        # The analogous phonon numerator includes the bath's 2a nBE. Newly
        # produced phonons can escape within the SAME step without requiring
        # their old exponentially tiny population to pay the entire escape.
        if np.any(bp<0) or np.any(bh<0) or np.any(bn<0):
            raise FloatingPointError('Supplied heating exhausts midpoint native-hole budget; reduce dt')
        def fractions(available,loss):
            # No outgoing loss consumes no budget.  Otherwise the safety
            # fraction is continuous through available/loss=1.
            q=np.full_like(available,np.inf)
            np.divide(available,loss,out=q,where=loss>0)
            return np.minimum(1.,.95*q)
        fp=fractions(bp.ravel(),loss_p).reshape(eshape)
        fh=fractions(bh.ravel(),loss_h).reshape(eshape)
        fn=fractions(bn.ravel(),loss_n).reshape(nshape)
        def efactor(nodes,ix,component):
            return np.where(component<0,fp[nodes,ix],np.where(component>0,fh[nodes,ix],1.))
        nodes=np.arange(self.node_count)[:,None]
        cf=np.minimum(efactor(nodes,self.event_i,np.where(rec,-rate,rate)),
                      efactor(nodes,self.event_j,-rate))
        for ix,weight in ((ce['lower'],1-ce['beta']),(ce['upper'],ce['beta'])):
            cf=np.minimum(cf,np.where((rate<0)&(weight>0),np.take_along_axis(fn,ix,axis=1),1.))
        limited=rate*cf
        cdp,cdn=self._collision_deposit(limited,ce['lower'],ce['upper'],ce['beta'],eshape,nshape)
        energy_factor=self.node_energy_scale_J[:,None]*dt*ce['omega']
        changed=float(np.sum(energy_factor*abs(rate-limited)))
        total=float(np.sum(energy_factor*abs(rate)))
        count=int(np.count_nonzero((cf<1)&(rate!=0)))
        active=int(np.count_nonzero(rate))
        reduction=float(np.max(1-cf,initial=0.))
        if te is not None:
            tf=np.ones_like(te['flux'])
            left,right=self.edges.T
            for nodes,(lower,upper,beta,_,_),sign in ((left,te['data'][0],-1),(right,te['data'][1],1)):
                for ix,weight in ((lower,1-beta),(upper,beta)):
                    tf=np.minimum(tf,efactor(nodes[:,None],ix,sign*weight*te['flux']))
            newflux=te['flux']*tf
            tdp,power,number=self._transport_deposit(newflux,te['q'],te['data'],eshape)
            unit=self.N0_per_J_m3*self.Eu_J**2*dt*te['q']
            changed+=float(np.sum(unit*abs(te['flux']-newflux)))
            total+=float(np.sum(unit*abs(te['flux'])))
            count+=int(np.count_nonzero((tf<1)&(te['flux']!=0)))
            active+=int(np.count_nonzero(te['flux']))
            reduction=max(reduction,float(np.max(1-tf,initial=0.)))
        else:tdp,power,number=np.zeros(eshape),np.empty(0),np.empty(0)
        return cdp,cdn,tdp,power,number,dict(limited_event_energy_J=changed,
            limited_event_fraction=count/max(active,1),totalabs_event_energy_J=total,
            max_flux_reduction=reduction)


def acoustic_phonon_grid(material,Eu_J,*,nodes=96,minimum_bar=1e-8):
    """Explicit nonuniform native phonon grid covering the entire material band.

    Low-energy nodes resolve acoustic heat capacity near the bath; a logarithmic
    high-energy section includes all coupled Debye events. Trapezoid capacities
    define this finite model and are refined as a numerical quadrature.
    """
    if nodes<24:raise ValueError('At least24 phonon nodes required')
    cutoff=material.energy_J[-1]/Eu_J
    nlow=max(5,nodes//12);nmid=2*nodes//3
    split=min(4.,cutoff/3)
    energy=np.r_[np.geomspace(minimum_bar,.01,nlow,endpoint=False),
        np.linspace(.01,split,nmid,endpoint=False),
        np.geomspace(split,cutoff,nodes-nlow-nmid)]
    weights=np.empty_like(energy);weights[1:-1]=(energy[2:]-energy[:-2])/2
    weights[0]=(energy[1]-energy[0])/2;weights[-1]=(energy[-1]-energy[-2])/2
    grid=material.phonon_grid(energy*Eu_J,weights*Eu_J,Eu_J)
    if material.manifest.get('analytic_acoustic_Debye'):
        from dataclasses import replace
        grid=replace(grid,alpha2F=QuadraticAcousticCoupling(
            material.manifest['alpha2F_coefficient_per_J2']*Eu_J**2))
    return grid


def forward_euler_time_limit_s(p,n,dp_dt_s,dn_dt_s):
    """Global native-support bound; zero identifies outward exhausted support.

    The caller applies its safety fraction and may replace explicitly fixed
    reservoir derivatives before requesting this bound. No vector is repaired.
    """
    p=np.asarray(p,float);n=np.asarray(n,float)
    dp=np.asarray(dp_dt_s,float);dn=np.asarray(dn_dt_s,float)
    if p.shape!=dp.shape or n.shape!=dn.shape:
        raise ValueError('Matching native populations and RHS required')
    bounds=[np.inf]
    if np.any(dp<0):bounds.append(float(np.min(p[dp<0]/(-dp[dp<0]))))
    if np.any(dp>0):bounds.append(float(np.min((1-p[dp>0])/dp[dp>0])))
    if np.any(dn<0):bounds.append(float(np.min(n[dn<0]/(-dn[dn<0]))))
    return min(bounds)
