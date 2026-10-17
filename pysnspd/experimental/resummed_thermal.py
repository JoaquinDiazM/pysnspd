"""DC-matched finite-depairing thermal gradient; effective, not microscopic.

G = kappa(rho,Tb)*L + u*(phi(rho,T,Q)-phi(rho,Tb,Q)),
Q=L/(u*rho). The temperature-independent matching term preserves the
validated DC elastic reference; only the finite-Q thermal correction is
resummed. H=kappa(rho,Tb)*L-u*phi(rho,Tb,Q) persists at fixed rho,L as
Te grows and can dominate hot strong-gradient states. This is a strong
mesoscopic hypothesis, not microscopic high-temperature stiffness.
It introduces no fitted latency parameter or numerical capacity floor.
Native collision spectra remain the explicitly declared q=0 approximation.
"""
import numpy as np
from scipy.constants import elementary_charge as e,hbar
from scipy.optimize import brentq,minimize_scalar
from .general_thermal_gradient import GeneralThermalGradientFunctional,GradientPartials
from .finite_depairing_sums import finite_depairing_partials,_partials_one
from .thermal_stiffness import thermal_rigidity
from scipy.special import zeta
from .effective_bcs_functional import EffectiveBCSBulk,EffectiveBCSFunctional

class ResummedThermalFunctional(GeneralThermalGradientFunctional):
    closure_id="TE_DC_matched_finite_depairing_thermal_gradient_v1"
    def __init__(self,*args,**kwargs):
        kwargs.pop("gradient_partials",None)
        super().__init__(*args,gradient_partials=self._gradient_partials,
                         closure_id=self.closure_id,**kwargs)
        self.gradient_law="DC_matched_finite_depairing_thermal_free_energy"
        if not np.isclose(self.radial_prefactor,np.pi/2,rtol=1e-12,atol=0.):
            raise ValueError("Finite-depairing mapping requires ell0=sqrt(hbar D/(2 Eu))")
    def _gradient_partials(self,rho,L_J,temperature_K):
        rho=np.asarray(rho,float);L=np.asarray(L_J,float);T=np.asarray(temperature_K,float)
        for oldrho,oldL,oldT,result in getattr(self,"_gradient_cache",[]):
            if np.array_equal(rho,oldrho) and np.array_equal(L,oldL) and np.array_equal(T,oldT):
                self._record_gradient_diagnostics(rho,L,T,result)
                return result
        u=self.node_energy_scale_J;fields=np.zeros((10,len(rho)))
        kb=thermal_rigidity(rho,self.bath_temperature_K,Tc_K=self.Tc_K)
        active=rho>0
        if np.any(active):
            r0=rho[active];l=L[active];uu=u[active];tt=T[active];Q=l/(uu*r0)
            vals=finite_depairing_partials(r0,tt/self.Tc_K,Q)
            ref=finite_depairing_partials(r0,np.full(len(r0),self.bath_temperature_K/self.Tc_K),Q)
            (v,r,q,t1,rr,rq,qq,rt,qt,t2)=vals
            dv=v-ref[0];dr=r-ref[1];dq=q-ref[2]
            drr=rr-ref[4];drq=rq-ref[5];dqq=qq-ref[6]
            fields[:,active]=np.array((
                l*kb.kappa[active]+uu*dv,
                l*kb.rho[active]+uu*(dr-Q*dq/r0),
                kb.kappa[active]+dq/r0,
                l*kb.rho_rho[active]+uu*(drr-2*Q/r0*drq+(Q*Q*dqq+2*Q*dq)/r0**2),
                kb.rho[active]+drq/r0-(dq+Q*dqq)/r0**2,
                dqq/(uu*r0**2),
                uu*t1/self.Tc_K,uu*(rt-Q*qt/r0)/self.Tc_K,
                qt/(r0*self.Tc_K),uu*t2/self.Tc_K**2))
        zero=~active;homogeneous=zero&(L==0);theta=T/self.Tc_K
        # At the fully normal homogeneous origin the Cartesian curvature
        # depends on the perturbation direction. The weak-gradient curvature
        # below is a preconditioner convention, not an exact Hessian claim.
        kt=thermal_rigidity(rho,T,Tc_K=self.Tc_K)
        fields[2,homogeneous]=kt.kappa[homogeneous]
        fields[4,homogeneous]=kt.rho[homogeneous]
        fields[8,homogeneous]=kt.temperature[homogeneous]
        isolated=zero&(L>0);b=self.bath_temperature_K/self.Tc_K
        # Exact rho->0 at fixed nonzero L: the rho*log(rho) terms cancel.
        fields[0,isolated]=L[isolated]*kb.kappa[isolated]
        fields[1,isolated]=L[isolated]*kb.rho[isolated]+u[isolated]*np.log(self.bath_temperature_K/T[isolated])
        fields[2,isolated]=kb.kappa[isolated]
        fields[3,isolated]=(L[isolated]*kb.rho_rho[isolated]-
            7*zeta(3)*u[isolated]/(8*np.pi**2)*(theta[isolated]**-2-b**-2))
        fields[4,isolated]=kb.rho[isolated]
        fields[7,isolated]=-u[isolated]/T[isolated]
        self.exact_zero_gap_curvature_preconditioner_count=int(np.count_nonzero(homogeneous))
        principal=fields[2]+2*L*fields[5]
        if (np.any(~np.isfinite(fields)) or np.any(fields[2]<=0) or np.any(principal<=0)):
            node=int(np.argmin(principal))
            raise ArithmeticError("Thermal gradient left elliptic constitutive domain: node="+str(node)+
                ", Te_K="+str(float(T[node]))+", G_L="+str(float(fields[2,node]))+
                ", G_L_plus_2L_G_LL="+str(float(principal[node]))+"; no stiffness floor")
        result=GradientPartials(*fields)
        self._record_gradient_diagnostics(rho,L,T,result)
        cache=getattr(self,"_gradient_cache",[])
        cache.insert(0,(rho.copy(),L.copy(),T.copy(),result));self._gradient_cache=cache[:4]
        return result
    def _record_gradient_diagnostics(self,rho,L,T,result):
        principal=result.G_L+2*L*result.G_L_L_per_J
        count=int(np.count_nonzero((rho==0)&(L==0)))
        self.exact_zero_gap_curvature_preconditioner_count=count
        self.last_gradient_diagnostics=dict(
            minimum_gradient_stiffness=float(np.min(result.G_L)),
            minimum_principal_stiffness=float(np.min(principal)),
            minimum_electron_temperature_K=float(np.min(T)),
            homogeneous_zero_gap_preconditioner_nodes=count,
            no_temperature_floor=True,no_stiffness_floor=True,
            matching_reference_temperature_K=float(self.bath_temperature_K),
            persistent_temperature_independent_matching_background=True,
            matching_background_is_effective_hypothesis=True)
    def _phi(self,rho,temperature_K,Q):
        t=temperature_K/self.Tc_K;b=self.bath_temperature_K/self.Tc_K
        actual=_partials_one(float(rho),float(t),float(Q),64,32.)
        ref=_partials_one(float(rho),float(b),float(Q),64,32.)
        k=thermal_rigidity(np.asarray(rho),self.bath_temperature_K,Tc_K=self.Tc_K)
        delta=actual-ref
        delta[0]+=rho*Q*k.kappa
        delta[1]+=Q*(k.kappa+rho*k.rho)
        delta[2]+=rho*k.kappa
        delta[4]+=Q*(2*k.rho+rho*k.rho_rho)
        delta[5]+=k.kappa+rho*k.rho
        # The matching subtraction is independent of the evolved T.
        delta[3]=actual[3];delta[7]=actual[7];delta[8]=actual[8];delta[9]=actual[9]
        return tuple(float(v) for v in delta)
    def _bulk_partials(self,rho,temperature_K,Q):
        if float(temperature_K)==self.bath_temperature_K:
            k=thermal_rigidity(np.asarray(rho),self.bath_temperature_K,Tc_K=self.Tc_K)
            return (float(Q*(k.kappa+rho*k.rho)),float(rho*k.kappa),
                float(Q*(2*k.rho+rho*k.rho_rho)),float(k.kappa+rho*k.rho),0.)
        values=self._phi(rho,temperature_K,Q)
        return values[1],values[2],values[4],values[5],values[6]
    def bulk_coefficient(self,amplitude_bar,temperature_K,q_bar=0.,phase_step_m=None):
        a=float(amplitude_bar);t=float(temperature_K)/self.Tc_K
        if a<0 or not np.isfinite(a+q_bar+t) or t<=0:
            raise ValueError("Finite nonnegative amplitude/q and positive T required")
        en=np.sqrt(self.xi_bar**2+a*a)
        native=self.inverse_pairing_coupling-np.dot(self.weights_bar,np.tanh(en/(2*t))/en)
        Q,_,_=self._bulk_q_terms(q_bar,phase_step_m)
        return float(native+self._bulk_partials(a*a,temperature_K,Q)[0])
    def bulk_amplitude(self,temperature_K,q_bar=0.,phase_step_m=None):
        # The matched elastic reference can have a metastable positive branch
        # even when the zero-amplitude coefficient is positive. Preserve the
        # inherited largest stable-root scan; do not discard that branch.
        return EffectiveBCSFunctional.bulk_amplitude(self,temperature_K,q_bar,phase_step_m)
    def bulk_current(self,temperature_K,q_bar,width_m,phase_step_m=None):
        a=self.bulk_amplitude(temperature_K,q_bar,phase_step_m)
        Q,currentq,_=self._bulk_q_terms(q_bar,phase_step_m)
        phiQ=self._bulk_partials(a*a,temperature_K,Q)[1]
        pref=4*e/hbar*self.N0_per_J_m3*self.Eu_J**2*self.ell0_m*float(width_m)*self.thickness_m
        return float(pref*phiQ*currentq),a
    def bulk(self,*,temperature_K,current_A,width_m,phase_step_m=None):
        temp=float(temperature_K);current=float(current_A);width=float(width_m)
        if not 0<temp<self.Tc_K or current<0 or width<=0:
            raise ValueError("Superconducting reference requires positive width/T and nonnegative current")
        limit=3. if phase_step_m is None else min(3.,np.pi*self.ell0_m/float(phase_step_m))
        axis=np.linspace(0.,limit,151)
        vals=np.array([self.bulk_current(temp,q,width,phase_step_m)[0] for q in axis])
        peak=int(np.argmax(vals))
        if peak in (0,len(axis)-1):raise RuntimeError("Finite-depairing maximum not bracketed")
        opt=minimize_scalar(lambda q:-self.bulk_current(temp,q,width,phase_step_m)[0],
            bounds=(axis[peak-1],axis[peak+1]),method="bounded",options={"xatol":2e-10})
        qdep=float(opt.x);idep=float(-opt.fun)
        if current>=idep:raise ValueError("Requested current exceeds finite-depairing current")
        q=0. if current==0 else brentq(lambda q:self.bulk_current(temp,q,width,phase_step_m)[0]-current,
            0.,qdep,xtol=2e-13,rtol=2e-13)
        measured,a=self.bulk_current(temp,q,width,phase_step_m)
        en=np.sqrt(self.xi_bar**2+a*a);t=temp/self.Tc_K;th=np.tanh(en/(2*t))
        local_da=np.dot(self.weights_bar,a*(th/en**3-(1-th*th)/(2*t*en**2)))
        Q,currentq,dcurrentq=self._bulk_q_terms(q,phase_step_m)
        _,phiQ,phirr,phirq,phiqq=self._bulk_partials(a*a,temp,Q)
        ca=float(local_da+2*a*phirr);cq=float(2*currentq*phirq)
        da_dq=-cq/ca
        pref=4*e/hbar*self.N0_per_J_m3*self.Eu_J**2*self.ell0_m*width*self.thickness_m
        di_dq=pref*(phiQ*dcurrentq+currentq*(2*a*phirq*da_dq+2*currentq*phiqq))
        return EffectiveBCSBulk(a,self.Eu_J*a,q,q/self.ell0_m,measured,idep,qdep,
            self.Eu_J*da_dq/di_dq,1/(self.ell0_m*di_dq),
            abs(self.bulk_coefficient(a,temp,q,phase_step_m)),bool(di_dq>0),temp,
            "Effective native BCS plus finite-depairing equilibrium free-gradient completion; Qeff includes radial gradients",phase_step_m)
