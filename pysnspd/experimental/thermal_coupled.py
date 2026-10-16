"""Conserving, admissibility-guarded electron-temperature stiffness candidate.

The v1.0 solver is unchanged. This subclass uses its implicit midpoint
algorithm, circuit, ports and shared event limiter with a new full RHS.
No residual-derived energy correction or thermal-capacity floor is used.
"""
from __future__ import annotations
from dataclasses import replace
import numpy as np
from scipy.constants import elementary_charge as e, hbar
from .effective_coupled import EffectiveCoupledSystem, EffectiveRHS
from .thermal_stiffness import entropy_compatible_carrier

class ThermalCapacityError(ArithmeticError):
    """The gradient free-energy closure left its positive-capacity domain."""

class ThermalCoupledSystem(EffectiveCoupledSystem):
    def evaluate(self, state, *, kinetic_budget=None, kinetic_dt_ps=None):
        f,k=self.functional,self.kinetics
        d,p,n,cm=state.delta_bar,state.p,state.n,state.circuit_SI
        if np.any(~np.isfinite(cm)) or min(cm[:2])<=0:
            raise ValueError("Global current reversal is outside the permitted window")
        temp=k.equivalent_temperature_K(d,p)
        ev=f.evaluate(d,p,temperature_K=temp)
        thermal=f.thermal_data(d,p,temp)
        ce=thermal["capacity_effective_J_K"]
        cfd=thermal["capacity_FD_J_K"]
        if np.any(~np.isfinite(ce)) or np.any(ce<=0):
            bad=int(np.argmin(ce/cfd))
            raise ThermalCapacityError("Nonpositive thermal capacity: node="+str(bad)+
                ", Ceff/C_FD="+str(float(ce[bad]/cfd[bad]))+
                ", Te_K="+str(float(temp[bad]))+"; no capacity floor is permitted")
        self.minimum_evaluated_capacity_ratio=min(getattr(self,"minimum_evaluated_capacity_ratio",np.inf),float(np.min(ce/cfd)))
        el=self.ports.solve(ev.edge_current_A,float(cm[1]))
        ratio=np.maximum(temp,k.bath_temperature_K)/f.Tc_K
        taupsi=1/(ratio/self.tau_ee_Tc_ps+ratio**3/self.tau_ep_Tc_ps)
        R=np.sqrt(1+(2*f.Eu_J*taupsi*1e-12/hbar)**2*abs(d)**2)
        A=np.sqrt((1+ratio)/2)
        radial=1/(A*self.tau0_ps*R);tangential=R/(A*self.tau0_ps)
        direction=np.ones(len(d),complex);nz=abs(d)>0
        direction[nz]=d[nz]/abs(d[nz])
        projected=np.conj(direction)*ev.gap_force_J_per_gapbar/f.node_energy_scale_J
        material=-direction*(radial*projected.real+1j*tangential*projected.imag)
        material[self.cut_nodes]=0.
        dz=material-1j*self.kappa_per_ps_V*el.phi_V*d
        coord=np.empty(self.ngap)
        coord[:2*len(self.free)]=np.column_stack((dz[self.free].real,dz[self.free].imag)).ravel()
        for j,nodes in enumerate(self.cuts):
            coord[2*len(self.free)+j]=-self.kappa_per_ps_V*el.phi_V[nodes[0]]
        q=2*f.node_energy_scale_J*(radial*projected.real**2+
            tangential*projected.imag**2)*1e12
        q[self.cut_nodes]=0.
        if np.any(q<0) or np.any(~np.isfinite(q)):
            raise ArithmeticError("Invalid KWT dissipation")
        # Bare conserving events/transport/BGK. The derived signed heat is
        # inserted below; it is not the old nonnegative heat-shape closure.
        zero=np.zeros(len(d))
        if kinetic_budget is None:
            bare=k.evaluate(d,p,n,heating_power_W=zero,temperature_K=temp)
        else:
            bare=k.evaluate_limited(d,p,n,budget_p=kinetic_budget.p,
                budget_n=kinetic_budget.n,dt_s=float(kinetic_dt_ps)*1e-12,
                heating_power_W=zero,temperature_K=temp)
        en=bare.electron_energies_J/f.Eu_J
        eps=4*f.node_energy_scale_J[:,None]*f.weights_bar*en
        rhodot=2*np.real(np.conj(d)*dz)*1e12
        Ldot=f.gradient_rate_geometry(d,dz*1e12)
        GTdot_fixedT=(thermal["kappa_T_per_K"]*Ldot+
            thermal["L_J"]*thermal["kappa_rho_T_per_K"]*rhodot)
        barepower=np.sum(eps*bare.dp_dt_s,axis=1)
        Tdot=(barepower+q+el.nodal_joule_W+thermal["B_J"]*rhodot+
            temp*GTdot_fixedT)/ce
        # Thermal reservoirs hold both native populations and |Delta| fixed.
        Tdot[self.cut_nodes]=0.
        GTdot=GTdot_fixedT+thermal["L_J"]*thermal["kappa_TT_per_K2"]*Tdot
        source=q+el.nodal_joule_W+temp*GTdot
        carrier=entropy_compatible_carrier(p,eps,temp,f.Eu_J,en,weights_bar=f.weights_bar,power_W=source)
        heat=carrier*source[:,None]
        full_dp=bare.dp_dt_s+heat
        actual_heat=np.sum(eps*heat,axis=1)
        heat_defect=actual_heat-source
        kin=replace(bare,dp_dt_s=full_dp,
            electron_population_power_W=np.sum(eps*full_dp,axis=1),
            heating_dp_dt_s=heat,heating_power_W=actual_heat,
            heating_residual_W=heat_defect,
            balance_residual_W=float(bare.balance_residual_W+np.sum(heat_defect)))
        dp=kin.dp_dt_s.copy();dn=kin.dn_dt_s.copy()
        rawcutpower=np.sum(eps*kin.dp_dt_s,axis=1)+kin.phonon_power_W
        reservoir=float(np.sum(rawcutpower[self.cut_nodes]))
        dp[self.cut_nodes]=0.;dn[self.cut_nodes]=0.
        cmrate=self.circuit.rhs(cm,el.Vdev_V)
        cb=self.circuit.power_balance(cm,el.Vdev_V,derivative=cmrate)
        fieldpower=2*np.real(np.vdot(ev.gap_force_J_per_gapbar,dz))*1e12
        epower=float(np.sum(eps*dp))
        phpower=float(np.sum(kin.phonon_power_W)-np.sum(kin.phonon_power_W[self.cut_nodes]))
        reversible=float(np.sum(temp*GTdot))
        escape=float(np.sum(kin.escape_power_W))
        domainrate=fieldpower+epower+phpower-reversible
        expectedfilm=el.port_power_W-escape-reservoir
        external=cb.source_power_W-cb.bias_joule_W-cb.load_joule_W-escape-reservoir
        residual=domainrate+cb.stored_energy_rate_W-external
        energy=ev.energy_J+kin.phonon_energy_J
        limits=[]
        for val,rate,upper in ((p,dp,1.),(n,dn,None)):
            negative=rate<0
            if np.any(negative):limits.append(float(np.min(-val[negative]/rate[negative])))
            if upper is not None:
                positive=rate>0
                if np.any(positive):limits.append(float(np.min((upper-val[positive])/rate[positive])))
        limit=min(limits,default=np.inf)*1e12
        kin=replace(kin,forward_euler_time_limit_s=limit*1e-12)
        self.last_thermal_diagnostics=dict(
            minimum_capacity_ratio=float(np.min(ce/cfd)),
            maximum_abs_reversible_gradient_power_W=float(np.max(abs(temp*GTdot))),
            maximum_abs_temperature_rate_K_ps=float(np.max(abs(Tdot))*1e-12),
            maximum_abs_heat_moment_defect_W=float(np.max(abs(heat_defect))),
            signed_population_source_allowed=True,
            heat_carrier_mobility="FD f*(1-p) heating; p*(1-f) cooling; nonnegative excess entropy",
            closed_Pauli_limits_supported=True,
            occupation_floor_used=False,
            capacity_floor_used=False)
        return EffectiveRHS(coord,dp*1e-12,dn*1e-12,cmrate/self.circuit_scales*1e-12,
            dz,material,q,temp,energy+self.circuit.energy(cm),energy,external,residual,
            domainrate-expectedfilm,reservoir,escape,limit,ev,el,kin)

    def coupled_iteration_damping(self,state,du,dp,dn,dc,cu,cp,cn,cc):
        """Keep internal nonlinear iterates in the positive-capacity domain."""
        f=self.functional;fraction=1.
        for _ in range(40):
            valid=True
            for scale in (.5,1.):
                trial=self.shifted(state,du+fraction*cu,dp+fraction*cp,
                    dn+fraction*cn,dc+fraction*cc,scale)
                if (np.any(~np.isfinite(trial.delta_bar)) or
                        np.any(trial.p<0) or np.any(trial.p>1) or
                        np.any(trial.n<0) or min(trial.circuit_SI[:2])<=0):
                    valid=False;break
                try:
                    td=f.thermal_data(trial.delta_bar,trial.p)
                except (ValueError,ArithmeticError):
                    valid=False;break
                ce=td["capacity_effective_J_K"]
                if np.any(~np.isfinite(ce)) or np.any(ce<=0):
                    valid=False;break
            if valid:return fraction
            fraction*=.5
        raise ArithmeticError("No admissible internal thermal Newton correction")

    def population_gap_response(self,state,rhs,p_direction):
        '''Analytic RHS response to a signed population residual direction.

        The direction need not define a supported trial population. We never
        evaluate p+direction or invert an FD temperature at such a trial.
        Force/current include the full TE chain; neutrality response is exact.
        Mobility uses its analytic temperature derivative, with the explicit
        directional increment of max(TE,Tb) at the bath corner.
        '''
        f=self.functional;d=state.delta_bar;p=state.p;temp=rhs.temperature_K
        dp=np.asarray(p_direction,float).copy()
        if dp.shape!=p.shape or np.any(~np.isfinite(dp)):
            raise ValueError("Finite matching population direction required")
        dp[self.cut_nodes]=0.
        thermal=f.thermal_data(d,p,temp)
        dtemp=np.sum(thermal["epsilon_J"]*dp,axis=1)/thermal["capacity_FD_J_K"]
        if np.any(~np.isfinite(dtemp)):
            raise ArithmeticError("Nonfinite thermal population direction")
        derivative=f.direction(d,p,np.zeros_like(d),dp,temperature_K=temp,
            temperature_direction_K=dtemp)
        dphi=self.ports.solve(derivative.edge_current_direction_A,0.).phi_V
        ratio=np.maximum(temp,self.kinetics.bath_temperature_K)/f.Tc_K
        dratio=(np.maximum(temp+dtemp,self.kinetics.bath_temperature_K)-
            np.maximum(temp,self.kinetics.bath_temperature_K))/f.Tc_K
        inverse_tau=ratio/self.tau_ee_Tc_ps+ratio**3/self.tau_ep_Tc_ps
        taupsi=1/inverse_tau
        R=np.sqrt(1+(2*f.Eu_J*taupsi*1e-12/hbar)**2*abs(d)**2)
        A=np.sqrt((1+ratio)/2)
        radial=1/(A*self.tau0_ps*R);tangential=R/(A*self.tau0_ps)
        log_tau_derivative=-(1/self.tau_ee_Tc_ps+3*ratio**2/self.tau_ep_Tc_ps)/inverse_tau
        log_R_derivative=(1-1/R**2)*log_tau_derivative
        dradial=radial*(-.5/(1+ratio)-log_R_derivative)*dratio
        dtangential=tangential*(-.5/(1+ratio)+log_R_derivative)*dratio
        z=np.ones(len(d),complex);nz=abs(d)>0;z[nz]=d[nz]/abs(d[nz])
        projected=np.conj(z)*rhs.evaluation.gap_force_J_per_gapbar/f.node_energy_scale_J
        dprojected=np.conj(z)*derivative.gap_force_direction_J_per_gapbar/f.node_energy_scale_J
        material=-z*(dradial*projected.real+radial*dprojected.real+
            1j*(dtangential*projected.imag+tangential*dprojected.imag))
        material[self.cut_nodes]=0.
        response=material-1j*self.kappa_per_ps_V*dphi*d
        if np.any(~np.isfinite(response)) or np.any(~np.isfinite(dphi)):
            raise ArithmeticError("Nonfinite full thermal gap-response indicator")
        return response,dphi,dtemp

    def thermal_population_gap_defect(self,state,rhs,p_residual,h_ps):
        response,dphi,_=self.population_gap_response(state,rhs,.5*p_residual)
        induced=float(h_ps*np.max(abs(response[self.free]),initial=0.))
        for nodes in self.cuts:
            induced=max(induced,float(h_ps*self.kappa_per_ps_V*abs(dphi[nodes[0]])))
        return induced

    def midpoint_step(self,state,h_ps,*,tolerance=2e-9,max_iterations=40):
        """Also check the full thermal gap response to population residual."""
        # Bound accumulated nonlinear error per unit physical time. A tiny
        # adaptive step must not receive the whole nominal-step error budget.
        solve_tolerance=tolerance*min(1.,float(h_ps)/.02)
        result,receipt=super().midpoint_step(state,h_ps,tolerance=solve_tolerance,
            max_iterations=max_iterations)
        mid=state.copy()
        mid.delta_bar=(state.delta_bar+result.delta_bar)/2
        for nodes in self.cuts:
            phase=np.angle(result.delta_bar[nodes]/state.delta_bar[nodes])
            mid.delta_bar[nodes]=state.delta_bar[nodes]*np.exp(.5j*phase)
        mid.p=(state.p+result.p)/2;mid.n=(state.n+result.n)/2
        mid.circuit_SI=(state.circuit_SI+result.circuit_SI)/2
        rhs=self.evaluate(mid,kinetic_budget=state,kinetic_dt_ps=h_ps)
        gp=h_ps*rhs.dp_per_ps-(result.p-state.p)
        response,dphi,dtemp=self.population_gap_response(mid,rhs,.5*gp)
        induced=float(h_ps*np.max(abs(response[self.free]),initial=0.))
        for nodes in self.cuts:
            induced=max(induced,float(h_ps*self.kappa_per_ps_V*abs(dphi[nodes[0]])))
        if induced>=solve_tolerance:
            raise ArithmeticError("Full thermal gap response to population residual exceeds tolerance: "+str(induced))
        receipt["nonlinear_nominal_tolerance"]=float(tolerance)
        receipt["nonlinear_actual_tolerance"]=float(solve_tolerance)
        receipt["nonlinear_reference_step_ps"]=.02
        receipt["thermal_induced_gap_defect"]=induced
        receipt["thermal_response_method"]="analytic signed population direction; no unsupported trial"
        receipt["thermal_population_direction_max_temperature_K"]=float(np.max(abs(dtemp),initial=0.))
        receipt["thermal_capacity_min_ratio"]=self.last_thermal_diagnostics["minimum_capacity_ratio"]
        return result,receipt
