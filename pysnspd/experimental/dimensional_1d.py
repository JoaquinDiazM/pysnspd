"""v1.2.2 1D cold-continuation cuts and polar midpoint coordinates.

The native BCS/KWT/p/n/circuit RHS is retained. The declared new exterior
assumption is an isothermal stable bulk branch at instantaneous Is. Its
prescribed gap/population work is included in the reservoir power.
Polar coordinates remove the midpoint phase-rotation error; Cartesian
coordinates remain for nearly vanishing condensate. No occupation/capacity
floor or relaxed energy tolerance is introduced.
"""
from dataclasses import replace
import numpy as np
from scipy.constants import elementary_charge as e,hbar
from scipy.interpolate import PchipInterpolator
from scipy.sparse import coo_matrix,csr_matrix,eye,bmat
from scipy.sparse.linalg import splu
from .thermal_coupled import ThermalCoupledSystem,ThermalCapacityError
from .effective_coupled import EffectiveRHS,EffectiveState,EffectiveCoupledSystem
from .thermal_stiffness import entropy_compatible_carrier
from .dimensional_transition import ObservedThermalSystem

class ColdLeadCurve:
    def __init__(self,functional,temperature_K,width_m,spacing_m,current_A):
        f=functional
        ref=f.bulk(temperature_K=temperature_K,current_A=current_A,width_m=width_m,
                   phase_step_m=spacing_m)
        q=np.unique(np.r_[np.linspace(0.,.98*ref.depairing_q_bar,257),
                          ref.q_m_inv*f.ell0_m])
        rows=np.array([f.bulk_current(temperature_K,v,width_m,spacing_m) for v in q])
        self.current=rows[:,0];self.gap=rows[:,1];self.q=q/f.ell0_m
        if np.any(np.diff(self.current)<=0):raise ArithmeticError("Nonmonotone stable bulk table")
        self.amp=PchipInterpolator(self.current,self.gap,extrapolate=False)
        self.wave=PchipInterpolator(self.current,self.q,extrapolate=False)
        self.reference_amplitude=float(self.amp(current_A))
        self.receipt=dict(assumption="instantaneous-current cold stable exterior; prescribed boundary work counted",
            temperature_K=temperature_K,knots=len(q),max_current_A=float(self.current[-1]),
            finite_link_spacing_m=spacing_m)
    def amplitude(self,current_A):
        a=float(self.amp(current_A));der=float(self.amp.derivative()(current_A))
        if not np.isfinite(a) or a<=0:raise ArithmeticError("Cold exterior left its stable branch")
        return a,der

class ColdPolar1DSystem(ObservedThermalSystem):
    physical_temporal_norm=True
    def configure(self,geometry,dc):
        self.lead_curve=ColdLeadCurve(self.functional,self.kinetics.bath_temperature_K,
            geometry.width_m,float(np.median(geometry.edge_lengths_m)),self.current_scale_A)
        self.cut_amplitude_scale=abs(dc.delta_bar[self.cut_nodes])/self.lead_curve.reference_amplitude
        self._polar=np.ones(len(self.free),bool)
        self.gauge_weights=self.functional.volumes_m3/self.functional.volumes_m3.sum()
        # The same gradient is requested repeatedly by force, heat capacity
        # and the independent population-response test. Reuse EXACT inputs;
        # never round states or interpolate this nonlinear thermal closure.
        f=self.functional
        original=getattr(f,"_uncached_gradient_evaluation",f.gradient_evaluation)
        f._uncached_gradient_evaluation=original
        cache={}
        def exact_gradient(d,T,alpha=None):
            key=(np.asarray(d).tobytes(),np.asarray(T).tobytes(),
                 None if alpha is None else np.asarray(alpha).tobytes())
            if key in cache:
                result,diagnostics=cache[key]
                if diagnostics is not None:f.last_gradient_diagnostics=dict(diagnostics)
                return result
            result=original(d,T,alpha)
            cache[key]=(result,dict(f.last_gradient_diagnostics) if hasattr(f,"last_gradient_diagnostics") else None)
            if len(cache)>8:cache.pop(next(iter(cache)))
            return result
        f.gradient_evaluation=exact_gradient
        self.solver_receipt=dict(coordinates="polar above 1e-3 gapbar; Cartesian below; full same midpoint RHS",
            gradient_cache="exact full d,T,alpha bytes; eight entries; no rounding",
            gap_gauge="zero volume-weighted mean scalar potential; electrostatic output remains right-cut referenced",
            temporal_norm="gap amplitude, physical supercurrent, device voltage, native energy, circuit; complex phase diagnostic",
            boundary=self.lead_curve.receipt)
        return self
    def impose_cold_cuts(self,state):
        st=state.copy();a,_=self.lead_curve.amplitude(st.circuit_SI[1])
        st.delta_bar[self.cut_nodes]=a*self.cut_amplitude_scale*np.exp(1j*np.angle(st.delta_bar[self.cut_nodes]))
        fd=self.functional.fermi_dirac(st.delta_bar,self.kinetics.bath_temperature_K)
        st.p[self.cut_nodes]=fd[self.cut_nodes]
        st.n[self.cut_nodes]=self.kinetics.bath_phonons
        return st
    def shifted(self,state,du,dp,dn,dc,fraction=1.):
        st=super().shifted(state,du,dp,dn,dc,fraction)
        change=du[:2*len(self.free)].reshape(-1,2)
        ix=self.free[self._polar]
        amplitude=abs(state.delta_bar[ix])+fraction*change[self._polar,0]
        if np.any(amplitude<=0):raise ArithmeticError("Internal polar amplitude crossed zero")
        st.delta_bar[ix]=amplitude*np.exp(1j*(np.angle(state.delta_bar[ix])+fraction*change[self._polar,1]))
        return self.impose_cold_cuts(st)
    def coordinates_rate(self,d,dz,coord):
        # Bootstrap RHS evaluations may precede midpoint_step's coordinate
        # selection. An exactly normal node always has Cartesian rates.
        polar=self._polar & (abs(d[self.free])>0)
        ix=self.free[polar]
        direction=d[ix]/abs(d[ix])
        v=np.conj(direction)*dz[ix]
        c=coord[:2*len(self.free)].reshape(-1,2)
        c[polar,0]=v.real;c[polar,1]=v.imag/abs(d[ix])
        return coord
    def coordinate_tangent(self,d):
        T=super().coordinate_tangent(d).tolil()
        for j,i in enumerate(self.free):
            if not self._polar[j]:continue
            z=d[i]/abs(d[i])
            T[2*i,2*j]=z.real;T[2*i+1,2*j]=z.imag
            T[2*i,2*j+1]=-d[i].imag;T[2*i+1,2*j+1]=d[i].real
        return T.tocsr()
    def midpoint_preconditioner(self,state,rhs,h_ps):
        d=state.delta_bar;nf=len(self.free);ng=self.ngap
        H,JI=self.functional.hessian_and_current_jacobian(d,state.p)
        T=self.coordinate_tangent(d);mob=self.mobility_blocks(d,rhs.temperature_K)
        selector=coo_matrix((np.ones(2*nf),(np.arange(2*nf),self._free_cart)),
            shape=(ng,2*len(d))).tolil()
        for j,i in enumerate(self.free):
            if not self._polar[j]:continue
            z=d[i]/abs(d[i])
            selector[2*j,2*i]=z.real;selector[2*j,2*i+1]=z.imag
            selector[2*j+1,2*i]=-z.imag/abs(d[i]);selector[2*j+1,2*i+1]=z.real/abs(d[i])
        forceJ=selector.tocsr()@mob@H@T
        rr=[];cc=[];gg=[]
        for j,i in enumerate(self.free):
            if self._polar[j]:continue
            rr.extend((2*j,2*j+1));cc.extend((2*j+1,2*j))
            phi=rhs.electrical.phi_V-float(np.dot(self.gauge_weights,rhs.electrical.phi_V))
            gg.extend((self.kappa_per_ps_V*phi[i],-self.kappa_per_ps_V*phi[i]))
        gaugeJ=coo_matrix((gg,(rr,cc)),shape=(ng,ng)).tocsr()
        A=eye(ng,format="csr")+h_ps/2*(forceJ-gaugeJ)
        unknown=self.ports.unknown
        group=self.ports.node_to_group
        gu=np.full(self.ports.group_count,-1,int);gu[np.flatnonzero(unknown)]=np.arange(np.count_nonzero(unknown))
        cr=[];cl=[];cv=[]
        gauge_vector=np.zeros(ng)
        for j,i in enumerate(self.free):
            k=gu[group[i]]
            if k<0:continue
            if self._polar[j]:
                cr.append(2*j+1);cl.append(k);cv.append(self.kappa_per_ps_V)
                gauge_vector[2*j+1]=self.kappa_per_ps_V
            else:
                cr.extend((2*j,2*j+1));cl.extend((k,k))
                cv.extend((-self.kappa_per_ps_V*d[i].imag,self.kappa_per_ps_V*d[i].real))
                gauge_vector[2*j:2*j+2]=[-self.kappa_per_ps_V*d[i].imag,self.kappa_per_ps_V*d[i].real]
        for j,nodes in enumerate(self.cuts):
            k=gu[group[nodes[0]]]
            if k>=0:cr.append(2*nf+j);cl.append(k);cv.append(self.kappa_per_ps_V)
            gauge_vector[2*nf+j]=self.kappa_per_ps_V
        nphi=np.count_nonzero(unknown)
        C=coo_matrix((cv,(cr,cl)),shape=(ng,nphi)).tocsr()*h_ps
        mean=np.asarray(self.gauge_weights@self.ports.prolongation[:,unknown]).ravel()
        C-=csr_matrix(gauge_vector[:,None]*mean[None,:])*h_ps
        B=self.ports.reduced_incidence[unknown];V=(B@JI@T)/2
        lap=self.ports.normal_laplacian[unknown][:,unknown]
        inj=self.ports.group_injection_per_A[unknown]
        W=csr_matrix(np.column_stack((np.zeros(len(inj)),-inj*self.current_scale_A/2,np.zeros(len(inj)))))
        c=self.circuit;dt=h_ps*1e-12;Is=self.current_scale_A;Vs=self.voltage_scale_V
        CM=np.array([[1+dt*(c.R_bias_ohm+c.R_load_ohm)/(2*c.L_bias_H),-dt*c.R_load_ohm/(2*c.L_bias_H),dt*Vs/(2*c.L_bias_H*Is)],
                     [-dt*c.R_load_ohm/(2*c.Lk_ext_H),1+dt*c.R_load_ohm/(2*c.Lk_ext_H),-dt*Vs/(2*c.Lk_ext_H*Is)],
                     [-dt*Is/(2*c.C_couple_F*Vs),dt*Is/(2*c.C_couple_F*Vs),1.]])
        port=np.zeros(nphi)
        for sign,g in ((1,self.ports.left_group),(-1,self.ports.right_group)):
            k=gu[g]
            if k>=0:port[k]+=sign
        D=csr_matrix(np.vstack((np.zeros_like(port),dt*port/(c.Lk_ext_H*Is),np.zeros_like(port))))
        return splu(bmat([[A,C,None],[V,lap,W],[None,D,csr_matrix(CM)]],format="csc")),nphi
    def population_gap_response(self,state,rhs,p_direction):
        response,dphi,dtemp=super().population_gap_response(state,rhs,p_direction)
        mean=float(np.dot(self.gauge_weights,dphi))
        response=response+1j*self.kappa_per_ps_V*mean*state.delta_bar
        return response,dphi-mean,dtemp

    def coupled_iteration_damping(self,state,du,dp,dn,dc,cu,cp,cn,cc):
        fraction=1.
        for _ in range(40):
            try:
                value=super().coupled_iteration_damping(state,du,dp,dn,dc,
                    fraction*cu,fraction*cp,fraction*cn,fraction*cc)
                return fraction*value
            except ArithmeticError:fraction*=.5
        raise ArithmeticError("No positive polar thermal correction")
    def midpoint_step(self,state,h_ps,*,tolerance=2e-7,max_iterations=60):
        self._polar=abs(state.delta_bar[self.free])>1e-3
        solve_tolerance=tolerance*min(1.,float(h_ps)/.02)
        result,receipt=EffectiveCoupledSystem.midpoint_step(self,state,h_ps,
            tolerance=solve_tolerance,max_iterations=max_iterations)
        # The converged polar midpoint is NOT the Cartesian endpoint average.
        # Reuse its exact accepted RHS for the same independent thermal test.
        rhs=self.last_observed_rhs;mid=self.last_observed_state
        gp=h_ps*rhs.dp_per_ps-(result.p-state.p)
        response,dphi,dtemp=self.population_gap_response(mid,rhs,.5*gp)
        induced=float(h_ps*np.max(abs(response[self.free]),initial=0.))
        for nodes in self.cuts:
            induced=max(induced,float(h_ps*self.kappa_per_ps_V*abs(dphi[nodes[0]])))
        if induced>=solve_tolerance:
            raise ArithmeticError("Full thermal gap response exceeds nonlinear tolerance")
        f=self.functional
        eps=4*f.node_energy_scale_J[:,None]*f.weights_bar*f.energies_bar(mid.delta_bar)
        raw=np.sum(eps*rhs.kinetics.dp_dt_s,axis=1)+rhs.kinetics.phonon_power_W
        flux=np.array([np.sum(raw[n]) for n in self.cuts])-self.last_boundary_work_W
        receipt.update(absolute_cut_energy_increment_J=h_ps*1e-12*float(np.sum(abs(flux))),
            cut_left_power_W=float(flux[0]),cut_right_power_W=float(flux[1]),
            cold_lead_work_increment_J=h_ps*1e-12*float(self.last_boundary_work_W.sum()),
            absolute_device_work_increment_J=h_ps*1e-12*abs(mid.circuit_SI[1]*rhs.electrical.Vdev_V),
            nonlinear_actual_tolerance=solve_tolerance,thermal_induced_gap_defect=induced,
            coordinates="mixed polar/Cartesian; exact phase rotation")
        return result,receipt

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
            phi_gap=el.phi_V-float(np.dot(self.gauge_weights,el.phi_V))
            cmrate=self.circuit.rhs(cm,el.Vdev_V)
            amp,derivative=self.lead_curve.amplitude(cm[1])
            adot=derivative*cmrate[1]*1e-12
            boundary_dd=np.zeros(len(d),complex)
            boundary_dd[self.cut_nodes]=(d[self.cut_nodes]/abs(d[self.cut_nodes]))*adot*self.cut_amplitude_scale
            boundary_dp=np.zeros_like(p)
            en0=f.energies_bar(d)
            fd=f.fermi_dirac(d,k.bath_temperature_K)
            drho0=2*np.real(np.conj(d)*boundary_dd)
            boundary_dp[self.cut_nodes]=-fd[self.cut_nodes]*(1-fd[self.cut_nodes])*drho0[self.cut_nodes,None]/(
                2*en0[self.cut_nodes]*(k.bath_temperature_K/f.Tc_K))
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
            dz=material-1j*self.kappa_per_ps_V*phi_gap*d+boundary_dd
            coord=np.empty(self.ngap)
            coord[:2*len(self.free)]=np.column_stack((dz[self.free].real,dz[self.free].imag)).ravel()
            for j,nodes in enumerate(self.cuts):
                coord[2*len(self.free)+j]=-self.kappa_per_ps_V*phi_gap[nodes[0]]
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
            GTdot_fixedT=(thermal["G_TL_per_K"]*Ldot+
                thermal["G_Trho_J_per_K"]*rhodot)
            barepower=np.sum(eps*bare.dp_dt_s,axis=1)
            Tdot=(barepower+q+el.nodal_joule_W+thermal["B_J"]*rhodot+
                temp*GTdot_fixedT)/ce
            # Thermal reservoirs hold both native populations and |Delta| fixed.
            Tdot[self.cut_nodes]=0.
            GTdot=GTdot_fixedT+thermal["G_TT_J_per_K2"]*Tdot
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
            dp[self.cut_nodes]=boundary_dp[self.cut_nodes]*1e12;dn[self.cut_nodes]=0.
            prescribed_power=(2*np.real(np.conj(ev.gap_force_J_per_gapbar)*boundary_dd)*1e12+
                              np.sum(eps*boundary_dp,axis=1)*1e12)
            boundary_work=np.array([float(np.sum(prescribed_power[nodes])) for nodes in self.cuts])
            reservoir-=float(boundary_work.sum())
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
            if hasattr(f,"last_gradient_diagnostics"):
                self.last_thermal_diagnostics.update(f.last_gradient_diagnostics)
            coord=self.coordinates_rate(d,dz,coord)
            result=EffectiveRHS(coord,dp*1e-12,dn*1e-12,cmrate/self.circuit_scales*1e-12,
                dz,material,q,temp,energy+self.circuit.energy(cm),energy,external,residual,
                domainrate-expectedfilm,reservoir,escape,limit,ev,el,kin)
            self.last_observed_rhs=result;self.last_observed_delta=d;self.last_observed_state=state
            self.last_boundary_work_W=boundary_work
            return result
