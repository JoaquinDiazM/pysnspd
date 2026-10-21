"""Gated late-time cold-wire manifold of the SAME finite-link BCS functional.

New explicit approximation: once the entire film is superconducting and its
native excess energy is negligible, retain only the three actual circuit
states. The film follows its isothermal finite-current bulk branch. This is
not applied to a hot belt, a latched wire or detection latency.
"""
from dataclasses import replace
import numpy as np
from scipy.constants import hbar,elementary_charge
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
from .effective_coupled import EffectiveState
from .dimensional_transition import film_energy,population_distance

class ColdRecovery:
    minimum_time_ps=100.
    qualification_ps=10.
    gap_margin=.01
    native_margin=.01
    current_margin=.005
    voltage_margin_V=1e-6
    normal_fraction_margin=.001
    maximum_step_ps=50.
    def __init__(self,system,geometry,dc,reference_J):
        self.s,self.g,self.dc,self.reference=system,geometry,dc,float(reference_J)
        self.ready=False
        self.receipt=dict(approximation="late cold isothermal BCS branch; negligible native excess discarded explicitly",
            minimum_time_ps=self.minimum_time_ps,qualification_ps=self.qualification_ps,
            gap_margin=self.gap_margin,native_energy_over_retained_photon=self.native_margin,
            supercurrent_margin=self.current_margin,voltage_margin_V=self.voltage_margin_V,
            normal_current_fraction_margin=self.normal_fraction_margin,
            maximum_step_ps=self.maximum_step_ps,
            limitation="not admitted for residual hot/resistive states; no late quasiparticle spectroscopy or retrigger claim")
    def state(self,cm):
        c=self.s.lead_curve;I=float(cm[1]);a,_=c.amplitude(I)
        q=float(c.wave(I));q0=float(c.wave(self.s.current_scale_A))
        x=self.g.coordinates_m[:,0]-self.g.length_m/2
        d=self.dc.delta_bar*(a/c.reference_amplitude)*np.exp(1j*(q-q0)*x)
        p=self.s.functional.fermi_dirac(d,self.s.kinetics.bath_temperature_K)
        n=np.broadcast_to(self.s.kinetics.bath_phonons,self.dc.n.shape).copy()
        return EffectiveState(d,p,n,np.array(cm,copy=True))
    def inductance(self,I):
        L=hbar/(2*elementary_charge)*self.g.length_m*float(self.s.lead_curve.wave.derivative()(I))
        if not np.isfinite(L) or L<=0:raise ArithmeticError("Nonpositive cold differential wire inductance")
        return L
    def rate(self,cm):
        c=self.s.circuit;Ib,Is,vc=cm
        if min(Ib,Is)<=0:raise ArithmeticError("Cold recovery branch reversal")
        L=self.inductance(Is)
        dIs=(vc+c.R_load_ohm*(Ib-Is))/(c.Lk_ext_H+L)
        dIb=(c.V_bias_V-c.R_bias_ohm*Ib-vc-c.R_load_ohm*(Ib-Is))/c.L_bias_H
        dvc=(Ib-Is)/c.C_couple_F
        return np.array([dIb,dIs,dvc]),L*dIs
    def eligibility(self,state,rhs,t_ps):
        if t_ps<self.minimum_time_ps:return None
        cold=self.state(state.circuit_SI);dcamp=abs(cold.delta_bar)
        amp=float(np.max(abs(abs(state.delta_bar)/dcamp-1)))
        if amp>self.gap_margin:return None
        native=population_distance(self.s,state,cold,self.reference)
        if native>self.native_margin:return None
        current=float(np.max(abs(rhs.evaluation.edge_current_A-state.circuit_SI[1]))/self.s.current_scale_A)
        if current>self.current_margin:return None
        rate,voltage=self.rate(state.circuit_SI)
        if abs(rhs.electrical.Vdev_V-voltage)>self.voltage_margin_V:return None
        R=self.s.ports.sheet_resistance_ohm*self.g.length_m/self.g.width_m
        normal=abs(voltage)/(R*state.circuit_SI[1])
        if normal>self.normal_fraction_margin:return None
        exchange=film_energy(self.s,state)-film_energy(self.s,cold)
        if abs(exchange)>self.native_margin*self.reference:return None
        return dict(gap_amplitude_relative=amp,native_energy_relative=native,
            supercurrent_relative=current,device_voltage_error_V=abs(rhs.electrical.Vdev_V-voltage),
            omitted_normal_current_fraction=normal,projection_bath_release_J=exchange)
    def prepare_table(self):
        if self.ready:return
        # E(I) uses the actual native functional; F_q gives electrical L.
        # Their difference is the isothermal bath entropy exchange, retained.
        I=self.s.lead_curve.current
        energies=np.array([film_energy(self.s,self.state([v,v,0.])) for v in I])
        z=.5*I*I
        self.energy=PchipInterpolator(z,energies,extrapolate=False)
        self.ready=True
    def step(self,state,rhs,h_ps,*,energy_tolerance,reference,absolute_work,**kwargs):
        self.prepare_table();s=self.s;c=s.circuit;scale=s.circuit_scales
        start=state.circuit_SI
        def ode(t,y):
            cm=y[:3]*scale
            rate,V=self.rate(cm);Is=cm[1]
            LE=float(self.energy.derivative()(.5*Is*Is))
            # Lq*I*dI is circuit work. LE*I*dI is actual resolved-film energy.
            bath=(self.inductance(Is)-LE)*Is*rate[1]
            external=(c.V_bias_V*cm[0]-c.R_bias_ohm*cm[0]**2-
                      c.R_load_ohm*(cm[0]-cm[1])**2-bath)
            resistance=self.s.ports.sheet_resistance_ohm*self.g.length_m/self.g.width_m
            return np.r_[rate/scale*1e-12,external/self.reference*1e-12,
                         V*V/resistance/self.reference*1e-12]
        sol=solve_ivp(ode,(0.,h_ps),np.r_[start/scale,0.,0.],method="DOP853",
            rtol=2e-9,atol=2e-11,max_step=min(10.,h_ps))
        if not sol.success:raise ArithmeticError("Cold circuit integration failed: "+sol.message)
        result=self.state(sol.y[:3,-1]*scale)
        rr=s.evaluate(result)
        rate,V=self.rate(result.circuit_SI)
        R=self.s.ports.sheet_resistance_ohm*self.g.length_m/self.g.width_m
        if abs(V)/(R*result.circuit_SI[1])>self.normal_fraction_margin:
            raise ArithmeticError("Cold manifold lost its normal-current scale separation")
        el=replace(rr.electrical,Vdev_V=float(V))
        rr=replace(rr,electrical=el,circuit_scaled_per_ps=rate/scale*1e-12)
        work=float(sol.y[3,-1]*self.reference)
        omitted_joule=float(sol.y[4,-1]*self.reference)
        defect=abs(rr.energy_J-rhs.energy_J-work)/(reference+absolute_work+abs(work))
        if defect>energy_tolerance:raise ArithmeticError("Cold independent energy acceptance failed: "+str(defect))
        rr.cold_reconstruction=True
        return result,rr,dict(temporal_error=0.,accepted_local_energy_fraction=defect,iterations=int(sol.nfev),
            external_energy_increment_J=work,absolute_external_energy_increment_J=abs(work),
            absolute_cut_energy_increment_J=0.,absolute_device_work_increment_J=abs(
                .5*(start[1]+result.circuit_SI[1])*.5*(rhs.electrical.Vdev_V+V)*h_ps*1e-12),
            regime="cold equilibrium reconstruction; 3-state actual circuit",
            omitted_normal_joule_J=omitted_joule,
            cold_solver_rtol=2e-9,cold_solver_atol=2e-11)


def advance_qualification(cold,state,rhs,t_ps,m,args):
    """Independent 10ps native/cold overlap before an official reduction."""
    info=m.setdefault("cold_recovery",dict(active=False))
    eligible=cold.eligibility(state,rhs,t_ps)
    if eligible is None:
        if not info.get("active"):
            info.clear();info.update(active=False)
        return None
    if "since_ps" not in info:
        initial=cold.state(state.circuit_SI)
        info.update(since_ps=t_ps,shadow_time_ps=t_ps,shadow_cm=initial.circuit_SI.tolist(),
            max_voltage_error_V=0.,max_current_relative=0.,checks=0,admission=cold.receipt)
        return None
    shadow=cold.state(np.array(info["shadow_cm"]))
    srhs=cold.s.evaluate(shadow)
    stime=info["shadow_time_ps"]
    while stime<t_ps-1e-10:
        h=min(cold.maximum_step_ps,t_ps-stime)
        try:
            shadow,srhs,_=cold.step(shadow,srhs,h,energy_tolerance=args.local_energy_tolerance,
                reference=cold.reference,absolute_work=m["absolute_external_J"])
        except (ArithmeticError,ValueError) as exc:
            info.clear();info.update(active=False,last_qualification_rejection=str(exc))
            return None
        stime+=h
    verror=abs(cold.s.circuit.R_load_ohm*((state.circuit_SI[0]-state.circuit_SI[1])-
                                            (shadow.circuit_SI[0]-shadow.circuit_SI[1])))
    ierror=float(np.max(abs(state.circuit_SI[:2]-shadow.circuit_SI[:2]))/cold.s.current_scale_A)
    if verror>cold.voltage_margin_V or ierror>cold.current_margin:
        info.clear();info.update(active=False,last_rejected_voltage_error_V=verror,
                                last_rejected_current_relative=ierror)
        return None
    info.update(shadow_cm=shadow.circuit_SI.tolist(),shadow_time_ps=t_ps,
        max_voltage_error_V=max(info["max_voltage_error_V"],verror),
        max_current_relative=max(info["max_current_relative"],ierror),checks=info["checks"]+1,
        latest_eligibility=eligible)
    if t_ps-info["since_ps"]<cold.qualification_ps:return None
    # Preserve the actual circuit at official activation; never replace it
    # with the independently compared shadow circuit.
    official=cold.state(state.circuit_SI)
    cold_rhs=cold.s.evaluate(official)
    rate,voltage=cold.rate(official.circuit_SI)
    cold_rhs=replace(cold_rhs,electrical=replace(cold_rhs.electrical,Vdev_V=float(voltage)),
                     circuit_scaled_per_ps=rate/cold.s.circuit_scales*1e-12)
    exchange=rhs.energy_J-cold_rhs.energy_J
    if abs(exchange)>cold.native_margin*cold.reference:return None
    info.update(active=True,activation_time_ps=t_ps,projection_bath_release_J=exchange,
                approximation=cold.receipt)
    m["external_J"]-=exchange
    m["absolute_external_J"]+=abs(exchange)
    print(__import__("json").dumps(dict(status="COLD_RECOVERY_ADMITTED",time_ps=t_ps,
        projection_bath_release_J=exchange,overlap=info)),flush=True)
    return official,cold_rhs
