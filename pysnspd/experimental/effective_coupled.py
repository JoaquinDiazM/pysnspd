"""Unified effective neutral BCS--KWT--kinetic--CM evolution, in SI.

All physical rates come from one state evaluation. Time stepping may use a
sparse block preconditioner; it does not define separate thermal time windows.
Rigid superconducting cuts have D_t Delta=0, fixed DC magnitude and independent
Josephson phases. Their population reservoirs exchange explicitly counted heat.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
import os,json
from scipy.constants import elementary_charge as e, hbar, Boltzmann as kB
from scipy.sparse import coo_matrix, csr_matrix, diags, eye, bmat
from scipy.sparse.linalg import splu


@dataclass
class EffectiveState:
    delta_bar: np.ndarray
    p: np.ndarray
    n: np.ndarray
    circuit_SI: np.ndarray

    def copy(self):
        return EffectiveState(*(np.array(x, copy=True) for x in
            (self.delta_bar, self.p, self.n, self.circuit_SI)))


@dataclass
class EffectiveRHS:
    gap_coordinate_per_ps: np.ndarray
    dp_per_ps: np.ndarray
    dn_per_ps: np.ndarray
    circuit_scaled_per_ps: np.ndarray
    delta_per_ps: np.ndarray
    material_delta_per_ps: np.ndarray
    kwt_heat_W: np.ndarray
    temperature_K: np.ndarray
    energy_J: float
    domain_energy_J: float
    external_power_W: float
    instant_power_residual_W: float
    film_power_residual_W: float
    reservoir_power_out_W: float
    escape_power_W: float
    forward_euler_limit_ps: float
    evaluation: object
    electrical: object
    kinetics: object


class EffectiveCoupledSystem:
    """One evolving gap/p/n/circuit state and its independent power observer."""
    def __init__(self, functional, kinetics, ports, circuit, *,
                 tau_ee_Tc_ps=5., tau_ep_Tc_ps=24.7,
                 current_scale_A=21.29720614489798e-6, voltage_scale_V=1e-3,
                 kinetic_residual_energy_scale_J=None):
        self.functional, self.kinetics = functional, kinetics
        self.ports, self.circuit = ports, circuit
        self.tau_ee_Tc_ps, self.tau_ep_Tc_ps = float(tau_ee_Tc_ps), float(tau_ep_Tc_ps)
        if min(self.tau_ee_Tc_ps, self.tau_ep_Tc_ps) <= 0:
            raise ValueError("Positive free KWT calibrators required")
        self.current_scale_A, self.voltage_scale_V = float(current_scale_A), float(voltage_scale_V)
        self.circuit_scales = np.array([self.current_scale_A]*2+[self.voltage_scale_V])
        self.free = np.asarray(ports.free, int)
        self.cuts = (np.asarray(ports.left, int), np.asarray(ports.right, int))
        self.cut_nodes = np.concatenate(self.cuts)
        self.ngap = 2*len(self.free)+2
        self.kappa_per_ps_V = 2*e/hbar*1e-12
        self._free_cart = np.ravel(np.column_stack((2*self.free,2*self.free+1)))
        self.tau0_ps = np.pi*hbar/(8*functional.Eu_J)*1e12
        reference=float(np.sum(functional.node_energy_scale_J))
        if kinetic_residual_energy_scale_J is not None:
            supplied=float(kinetic_residual_energy_scale_J)
            if not np.isfinite(supplied) or supplied<=0:
                raise ValueError('Positive physical kinetic residual energy scale required')
            reference=supplied
        self.kinetic_residual_energy_scale_J=reference

    def coordinate_tangent(self, d):
        n, nf = len(d), len(self.free)
        rows = list(self._free_cart); cols = list(range(2*nf)); values = [1.]*(2*nf)
        for j, nodes in enumerate(self.cuts):
            for i in nodes:
                rows.extend((2*i,2*i+1)); cols.extend((2*nf+j,2*nf+j))
                values.extend((-d[i].imag,d[i].real))
        return coo_matrix((values,(rows,cols)),shape=(2*n,self.ngap)).tocsr()

    def shifted(self, state, du, dp, dn, dc, fraction=1.):
        d = state.delta_bar.copy()
        interior = du[:2*len(self.free)].reshape(-1,2)
        d[self.free] += fraction*(interior[:,0]+1j*interior[:,1])
        for j, nodes in enumerate(self.cuts):
            d[nodes] *= np.exp(1j*fraction*du[2*len(self.free)+j])
        return EffectiveState(d, state.p+fraction*dp, state.n+fraction*dn,
            state.circuit_SI+fraction*dc*self.circuit_scales)

    def mobility_blocks(self, d, temperature_K):
        f = self.functional
        ratio = np.maximum(temperature_K, self.kinetics.bath_temperature_K)/f.Tc_K
        taupsi = 1/(ratio/self.tau_ee_Tc_ps+ratio**3/self.tau_ep_Tc_ps)
        stretch = np.sqrt(1+(2*f.Eu_J*taupsi*1e-12/hbar)**2*abs(d)**2)
        A = np.sqrt((1+ratio)/2)
        radial = 1/(A*self.tau0_ps*stretch)
        tangential = stretch/(A*self.tau0_ps)
        direction = np.ones(len(d), complex)
        nz = abs(d)>0
        direction[nz] = d[nz]/abs(d[nz])
        x,y=direction.real,direction.imag
        aa=tangential+(radial-tangential)*x*x
        bb=(radial-tangential)*x*y
        cc=tangential+(radial-tangential)*y*y
        scale=f.node_energy_scale_J
        rows=np.ravel(np.column_stack((2*self.free,2*self.free,2*self.free+1,2*self.free+1)))
        cols=np.ravel(np.column_stack((2*self.free,2*self.free+1,2*self.free,2*self.free+1)))
        vals=np.column_stack((aa[self.free]/scale[self.free],bb[self.free]/scale[self.free],
                              bb[self.free]/scale[self.free],cc[self.free]/scale[self.free])).ravel()
        return coo_matrix((vals,(rows,cols)),shape=(2*len(d),2*len(d))).tocsr()

    def evaluate(self, state, *, kinetic_budget=None, kinetic_dt_ps=None):
        d,p,n,cm=state.delta_bar,state.p,state.n,state.circuit_SI
        if min(cm[:2])<=0:
            raise ValueError("Global current reversal is outside the permitted window")
        ev=self.functional.evaluate(d,p)
        el=self.ports.solve(ev.edge_current_A,float(cm[1]))
        temp=self.kinetics.equivalent_temperature_K(d,p)
        mob=self.mobility_blocks(d,temp)
        force=np.column_stack((ev.gap_force_J_per_gapbar.real,ev.gap_force_J_per_gapbar.imag)).ravel()
        ratio=np.maximum(temp,self.kinetics.bath_temperature_K)/self.functional.Tc_K
        taupsi=1/(ratio/self.tau_ee_Tc_ps+ratio**3/self.tau_ep_Tc_ps)
        R=np.sqrt(1+(2*self.functional.Eu_J*taupsi*1e-12/hbar)**2*abs(d)**2)
        A=np.sqrt((1+ratio)/2)
        radial=1/(A*self.tau0_ps*R);tangential=R/(A*self.tau0_ps)
        direction=np.ones(len(d),complex);nz=abs(d)>0;direction[nz]=d[nz]/abs(d[nz])
        projected=np.conj(direction)*ev.gap_force_J_per_gapbar/self.functional.node_energy_scale_J
        material=-direction*(radial*projected.real+1j*tangential*projected.imag)
        material[self.cut_nodes]=0.
        dz=material-1j*self.kappa_per_ps_V*el.phi_V*d
        coord=np.empty(self.ngap)
        coord[:2*len(self.free)]=np.column_stack((dz[self.free].real,dz[self.free].imag)).ravel()
        for j,nodes in enumerate(self.cuts):
            coord[2*len(self.free)+j]=-self.kappa_per_ps_V*el.phi_V[nodes[0]]
        q=2*self.functional.node_energy_scale_J*(radial*projected.real**2+
            tangential*projected.imag**2)*1e12
        q[self.cut_nodes]=0.
        if np.any(q < -1e-26):
            raise ArithmeticError("KWT dissipation is negative")
        if kinetic_budget is None:
            kin=self.kinetics.evaluate(d,p,n,heating_power_W=q+el.nodal_joule_W,temperature_K=temp)
        else:
            kin=self.kinetics.evaluate_limited(d,p,n,budget_p=kinetic_budget.p,
                budget_n=kinetic_budget.n,dt_s=float(kinetic_dt_ps)*1e-12,
                heating_power_W=q+el.nodal_joule_W,temperature_K=temp)
        dp=kin.dp_dt_s.copy();dn=kin.dn_dt_s.copy()
        popscale=4*self.functional.node_energy_scale_J[:,None]*self.functional.weights_bar
        rawcutpower=(np.sum(popscale*kin.electron_energies_J/self.functional.Eu_J
                              *kin.dp_dt_s,axis=1)+kin.phonon_power_W)
        reservoir=float(np.sum(rawcutpower[self.cut_nodes]))
        dp[self.cut_nodes]=0.;dn[self.cut_nodes]=0.
        cmrate=self.circuit.rhs(cm,el.Vdev_V)
        cb=self.circuit.power_balance(cm,el.Vdev_V,derivative=cmrate)
        fieldpower=2*np.real(np.vdot(ev.gap_force_J_per_gapbar,dz))*1e12
        epower=float(np.sum(popscale*kin.electron_energies_J/self.functional.Eu_J*dp))
        phpower=float(np.sum(kin.phonon_power_W)-np.sum(kin.phonon_power_W[self.cut_nodes]))
        escape=float(np.sum(kin.escape_power_W))
        domainrate=fieldpower+epower+phpower
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
        return EffectiveRHS(coord,dp*1e-12,dn*1e-12,cmrate/self.circuit_scales*1e-12,
            dz,material,q,temp,energy+self.circuit.energy(cm),energy,external,residual,
            domainrate-expectedfilm,reservoir,escape,limit,ev,el,kin)

    def midpoint_preconditioner(self,state,rhs,h_ps):
        """Sparse Newton preconditioner; the converged midpoint uses full RHS."""
        d=state.delta_bar; nf=len(self.free); ng=self.ngap
        H,JI=self.functional.hessian_and_current_jacobian(d,state.p)
        T=self.coordinate_tangent(d)
        mob=self.mobility_blocks(d,rhs.temperature_K)
        selector=coo_matrix((np.ones(2*nf),(np.arange(2*nf),self._free_cart)),
            shape=(ng,2*len(d))).tocsr()
        forceJ=selector@mob@H@T
        # Ordinary-time gauge term derivative at fixed phi, on interior Cartesian DOFs.
        rr=np.ravel(np.column_stack((np.arange(0,2*nf,2),np.arange(1,2*nf,2))))
        cc=np.ravel(np.column_stack((np.arange(1,2*nf,2),np.arange(0,2*nf,2))))
        gg=np.column_stack((self.kappa_per_ps_V*rhs.electrical.phi_V[self.free],
                           -self.kappa_per_ps_V*rhs.electrical.phi_V[self.free])).ravel()
        gaugeJ=coo_matrix((gg,(rr,cc)),shape=(ng,ng)).tocsr()
        A=eye(ng,format="csr")+h_ps/2*(forceJ-gaugeJ)
        unknown=self.ports.unknown
        P=self.ports.prolongation[:,unknown]
        group=self.ports.node_to_group
        group_unknown=np.full(self.ports.group_count,-1,int)
        group_unknown[np.flatnonzero(unknown)]=np.arange(np.count_nonzero(unknown))
        cr=[];cl=[];cv=[]
        for j,i in enumerate(self.free):
            g=group_unknown[group[i]]
            if g>=0:
                cr.extend((2*j,2*j+1));cl.extend((g,g))
                cv.extend((-self.kappa_per_ps_V*d[i].imag,self.kappa_per_ps_V*d[i].real))
        for j,nodes in enumerate(self.cuts):
            g=group_unknown[group[nodes[0]]]
            if g>=0: cr.append(2*nf+j);cl.append(g);cv.append(self.kappa_per_ps_V)
        C=coo_matrix((cv,(cr,cl)),shape=(ng,P.shape[1])).tocsr()*h_ps
        B=self.ports.reduced_incidence[unknown]
        V=(B@JI@T)/2
        lap=self.ports.normal_laplacian[unknown][:,unknown]
        inj=self.ports.group_injection_per_A[unknown]
        W=csr_matrix(np.column_stack((np.zeros(len(inj)),-inj*self.current_scale_A/2,np.zeros(len(inj)))))
        c=self.circuit;dt=h_ps*1e-12;Is=self.current_scale_A;Vs=self.voltage_scale_V
        CM=np.array([[1+dt*(c.R_bias_ohm+c.R_load_ohm)/(2*c.L_bias_H),
                      -dt*c.R_load_ohm/(2*c.L_bias_H),dt*Vs/(2*c.L_bias_H*Is)],
                     [-dt*c.R_load_ohm/(2*c.Lk_ext_H),
                      1+dt*c.R_load_ohm/(2*c.Lk_ext_H),-dt*Vs/(2*c.Lk_ext_H*Is)],
                     [-dt*Is/(2*c.C_couple_F*Vs),dt*Is/(2*c.C_couple_F*Vs),1.]])
        vport=np.zeros(P.shape[1])
        for sign,g in ((1,self.ports.left_group),(-1,self.ports.right_group)):
            gi=group_unknown[g]
            if gi>=0:vport[gi]+=sign
        D=csr_matrix(np.vstack((np.zeros_like(vport),dt*vport/(c.Lk_ext_H*Is),np.zeros_like(vport))))
        matrix=bmat([[A,C,None],[V,lap,W],[None,D,csr_matrix(CM)]],format="csc")
        return splu(matrix),P.shape[1]

    def midpoint_step(self,state,h_ps,*,tolerance=2e-9,max_iterations=40):
        """Solve ALL midpoint equations, including kinetics and CM, to tolerance.

        Gap/circuit/potential corrections are preconditioned together. Kinetic
        corrections are fixed point iterations of the SAME full midpoint RHS.
        Shared-event availability limits are a finite-step Patankar treatment.
        They approach the continuous RHS under refinement and are audited by
        limited event energy. Population convergence uses unsigned physical
        energy plus an independent induced-gap-force bound at the same tolerance.
        Pauli/Bose support and the complete independent energy audit remain strict. Internal positivity damping changes only iteration.
        This is neither population clipping nor an energy-defect correction.
        """
        h=float(h_ps);ng=self.ngap
        du=np.zeros(ng);dp=np.zeros_like(state.p);dn=np.zeros_like(state.n);dc=np.zeros(3)
        factor=None;previous=np.inf
        # Legacy models accelerate populations; thermal models accelerate the full
        # preconditioned gap/population/circuit map. Acceptance checks the same
        # unpreconditioned physical defect and constrained internal iterates.
        history_x=[];history_r=[]
        coupled_anderson=bool(getattr(self,"coupled_anderson",False))
        nscale=1+state.n
        for iteration in range(max_iterations):
            mid=self.shifted(state,du,dp,dn,dc,.5)
            rhs=self.evaluate(mid,kinetic_budget=state,kinetic_dt_ps=h)
            gu=h*rhs.gap_coordinate_per_ps-du
            gp=h*rhs.dp_per_ps-dp
            gn=h*rhs.dn_per_ps-dn
            gc=h*rhs.circuit_scaled_per_ps-dc
            # Population defects use unsigned physical energy, not equal
            # accuracy for every nearly empty cold occupation. The induced
            # gap-force correction is guarded independently at the gap tolerance.
            f=self.functional
            ep_weights=(4*f.node_energy_scale_J[:,None]*f.weights_bar
                        *rhs.kinetics.electron_energies_J/f.Eu_J)
            p_energy_defect=float(np.sum(ep_weights*abs(gp)))
            n_energy_defect=float(np.sum(self.kinetics.phonon_energy_weights_J*abs(gn)))
            delta_force=(f.node_energy_scale_J*mid.delta_bar
                *np.sum(f.weights_bar*gp/(rhs.kinetics.electron_energies_J/f.Eu_J),axis=1))
            force_cart=np.column_stack((delta_force.real,delta_force.imag)).ravel()
            induced_gap=h*(self.mobility_blocks(mid.delta_bar,rhs.temperature_K)@force_cart)
            induced_gap_defect=float(np.max(abs(induced_gap[self._free_cart]),initial=0.))
            if hasattr(self,"thermal_population_gap_defect"):
                induced_gap_defect=max(induced_gap_defect,
                    self.thermal_population_gap_defect(mid,rhs,gp,h))
            defect=max(float(np.max(abs(gu))),float(np.max(abs(gc))),
                p_energy_defect/self.kinetic_residual_energy_scale_J,
                n_energy_defect/self.kinetic_residual_energy_scale_J,induced_gap_defect)
            if os.environ.get('PYSNSPD_MIDPOINT_DEBUG')=='1':
                print(json.dumps(dict(midpoint_iteration=iteration,h_ps=h,defect=defect,
                    gap_defect=float(np.max(abs(gu))),p_defect=float(np.max(abs(gp))),
                    p_energy_defect_J=p_energy_defect,n_energy_defect_J=n_energy_defect,
                    kinetic_reference_energy_J=self.kinetic_residual_energy_scale_J,
                    induced_gap_defect=induced_gap_defect,
                    n_defect=float(np.max(abs(gn)/(1+mid.n))),cm_defect=float(np.max(abs(gc))),
                    p_min=float(mid.p.min()),n_min=float(mid.n.min()))),flush=True)
            # A small h scales every residual below an absolute tolerance.
            # Accepting the untouched predictor then advances time with zero
            # dynamics. Always perform a coupled correction before acceptance.
            if iteration > 0 and defect<tolerance:
                result=self.shifted(state,du,dp,dn,dc)
                if (np.any(result.p<0) or np.any(result.p>1) or np.any(result.n<0)
                        or np.any(~np.isfinite(result.delta_bar)) or min(result.circuit_SI[:2])<=0):
                    raise ArithmeticError("Midpoint endpoint outside physical support")
                return result,dict(iterations=iteration,defect=defect,
                    anderson_scope="gap,p,n,circuit" if coupled_anderson else "p,n",
                    electron_residual_L1_energy_J=p_energy_defect,
                    phonon_residual_L1_energy_J=n_energy_defect,
                    kinetic_reference_energy_J=self.kinetic_residual_energy_scale_J,
                    induced_gap_defect=induced_gap_defect,
                    maximum_electron_occupation_residual=float(np.max(abs(gp))),
                    external_energy_increment_J=h*1e-12*rhs.external_power_W,
                    midpoint_power_residual_W=rhs.instant_power_residual_W,
                    limited_event_energy_J=getattr(rhs.kinetics,'limited_event_energy_J',0.),
                    totalabs_event_energy_J=getattr(rhs.kinetics,'totalabs_event_energy_J',0.))
            if factor is None or defect>previous*.9:
                factor,nphi=self.midpoint_preconditioner(mid,rhs,h)
            correction=factor.solve(np.r_[gu,np.zeros(nphi),gc])
            # Bound internal iterative corrections coordinate by coordinate.
            # A microscopic tail must not damp the entire gap/current Newton step.
            # This changes only the nonlinear solver path: the accepted equations
            # and the event stoichiometry remain those in the residual above.
            def feasible_fraction(value, increment, upper=None):
                alpha=np.ones_like(value)
                negative=increment<0
                alpha[negative]=np.minimum(1.,.95*value[negative]/(-increment[negative]))
                if upper is not None:
                    positive=increment>0
                    alpha[positive]=np.minimum(1.,.95*(upper-value[positive])/increment[positive])
                if np.any(alpha<0):
                    raise ArithmeticError("Internal iterate already outside physical support")
                return alpha
            if coupled_anderson:
                kinetic_x=np.r_[du,dp.ravel(),(dn/nscale).ravel(),dc]
                kinetic_r=np.r_[correction[:ng],gp.ravel(),(gn/nscale).ravel(),correction[ng+nphi:]]
            else:
                kinetic_x=np.r_[dp.ravel(),(dn/nscale).ravel()]
                kinetic_r=np.r_[gp.ravel(),(gn/nscale).ravel()]
            kinetic_correction=kinetic_r.copy()
            if defect>previous*2:
                history_x.clear();history_r.clear()
            if history_x and iteration>=2:
                xs=history_x+[kinetic_x];rs=history_r+[kinetic_r]
                DX=np.column_stack([xs[j+1]-xs[j] for j in range(len(xs)-1)])
                DR=np.column_stack([rs[j+1]-rs[j] for j in range(len(rs)-1)])
                gamma=np.linalg.lstsq(DR,kinetic_r,rcond=1e-10)[0]
                accelerated=kinetic_r-(DX+DR)@gamma
                if np.all(np.isfinite(accelerated)):
                    kinetic_correction=accelerated
            history_x.append(kinetic_x);history_r.append(kinetic_r)
            if len(history_x)>(10 if coupled_anderson else 5):history_x.pop(0);history_r.pop(0)
            split=state.p.size
            if coupled_anderson:
                cu=kinetic_correction[:ng];cc=kinetic_correction[-3:]
                cp=kinetic_correction[ng:ng+split].reshape(state.p.shape)
                cn=kinetic_correction[ng+split:-3].reshape(state.n.shape)*nscale
            else:
                cu=correction[:ng];cc=correction[ng+nphi:]
                cp=kinetic_correction[:split].reshape(state.p.shape)
                cn=kinetic_correction[split:].reshape(state.n.shape)*nscale
            ap=feasible_fraction(state.p+dp,cp,1.)
            an=feasible_fraction(state.n+dn,cn)
            cp=ap*cp;cn=an*cn
            damping=(self.coupled_iteration_damping(state,du,dp,dn,dc,cu,cp,cn,cc)
                if hasattr(self,"coupled_iteration_damping") else 1.)
            du+=damping*cu;dp+=damping*cp;dn+=damping*cn
            dc+=damping*cc
            if min((state.circuit_SI+dc*self.circuit_scales)[:2])<=0:
                raise ArithmeticError("Newton current step would reverse global current")
            previous=defect
        raise ArithmeticError("Coupled midpoint did not converge")
