"""First physical Rosenbrock-Euler candidate of the coupled photon DAE.

A first-step IMEX convention treats C(y_old)=0 at the actual prepared common
local FD/BE state explicitly. It is forbidden for subsequent states. Born
spectral broadening is still omitted in the inherited weak-Born spectral
ordering; no assumption of general material/time/energy admission is made.

All zeroth spatial kinetic derivatives enter M+dt*Lh. Gap/charge/current
moments include their physical state derivative; CM uses its exact linear
Jacobian. The only allowed numerical reduction is the verified voltage gauge.
"""
from __future__ import annotations
from dataclasses import replace
from copy import copy
import numpy as np
from scipy.constants import elementary_charge as e, hbar, Boltzmann as kB
from scipy.sparse import coo_matrix
from . import adiabatic_matrix_usadel as ad
from . import causal_kinetic_regulator as reg
from .joint_evolution import JointSnapshotSystem, JointRateValue


def assemble_zeroth_kinetic_operator(graph, R, B, *, numerical_eta, kappa_s, alpha=None):
    """Exact sparse L_h = d Pi(F0) / dh at fixed physical spectrum.

    Its coefficient has no seconds; M+dt L_h retains the sign of the
    original kinetic residual, including its normal diffusion limit.
    """
    n=graph.n_nodes
    R=ad.Jet.of(R.value); B=ad.Jet.of(B.value)
    A=ad.advanced(R)
    tail,head=graph.edges.T
    alpha=np.zeros(len(tail)) if alpha is None else np.asarray(alpha,float)
    if alpha.shape!=tail.shape or np.any(~np.isfinite(alpha)):
        raise ValueError("One finite alpha per edge required")
    P=ad.links(alpha)
    Ri,Ai=R.take(tail),A.take(tail)
    Rj=ad.transport(P,R.take(head),kappa_s)
    Aj=ad.transport(P,A.take(head),kappa_s)
    local=np.zeros((n,2,2)); tail_blocks=np.zeros((len(tail),2,2)); head_blocks=np.zeros_like(tail_blocks)
    for c,basis in enumerate((ad.I2,ad.TAU3)):
        H=ad.Jet.of(np.broadcast_to(basis,R.value.shape))
        K=ad.keldysh(R,H,kappa_s)
        term=ad.commutator(B,K,kappa_s).scale(-graph.area_weights/2)
        term=term-reg.regulator_term(R,H,numerical_eta,kappa_s=kappa_s).scale(graph.area_weights/2)
        local[:,:,c]=ad.projections(term.value)
        Ki=K.take(tail); Kj=ad.transport(P,K.take(head),kappa_s)
        tf=(ad.star(Ki,Aj,kappa_s)-ad.star(Rj,Ki,kappa_s)).scale(graph.conductance/2)
        hf=(ad.star(Ri,Kj,kappa_s)-ad.star(Kj,Ai,kappa_s)).scale(graph.conductance/2)
        tail_blocks[:,:,c]=ad.projections(tf.value)
        head_blocks[:,:,c]=ad.projections(hf.value)
    nodes=np.arange(n); rows=[];cols=[];values=[]
    for out in range(2):
        for c in range(2):
            rows += [2*nodes+out,2*tail+out,2*tail+out,2*head+out,2*head+out]
            cols += [2*nodes+c,2*tail+c,2*head+c,2*tail+c,2*head+c]
            values += [local[:,out,c],tail_blocks[:,out,c],head_blocks[:,out,c],-tail_blocks[:,out,c],-head_blocks[:,out,c]]
    result=coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(cols))),shape=(2*n,2*n)).tocsr()
    result.eliminate_zeros()
    return result


class JointRosenbrockSystem(JointSnapshotSystem):
    """Coupled Schur elimination for ONE linearized first photon step.

    Reconstructing a finite candidate and re-solving its static spectral
    equations is mandatory. Neither the affine solve nor a small dt admits
    the material model or a nonlinear trajectory.
    """
    def __init__(self, snapshot, dt_s, thermal_tangent, *, tau_escape_s=np.inf):
        if not np.isfinite(dt_s) or dt_s<=0:
            raise ValueError("Positive finite physical step required")
        self.__dict__.update(snapshot.__dict__)
        self.slices=tuple(copy(sl) for sl in snapshot.slices)
        self.dt_s=float(dt_s); self.thermal_tangent=thermal_tangent
        if self.phonon_rate is not None and np.any(self.phonon_rate!=0):
            raise ValueError("First FD/BE detailed-balance step expects zero phonon rate")
        if not hasattr(self,"initial_data") or self.initial_data["receipt"].get("source_must_not_be_applied_again") is not True:
            raise ValueError("Actual single-absorption FD/BE initial preparation required")
        for sl in self.slices:
            if np.any(sl.source!=0) or np.any(sl.SigmaB!=0):
                raise ValueError("This first IMEX adapter uses declared initial zero collision, not general material")
            sl.configure_linearized_step(self.dt_s)
        c=self.circuit
        J=np.array([[-(c.R_bias_ohm+c.R_load_ohm)/c.L_bias_H,c.R_load_ohm/c.L_bias_H,-1/c.L_bias_H],
            [c.R_load_ohm/c.Lk_ext_H,-c.R_load_ohm/c.Lk_ext_H,1/c.Lk_ext_H],
            [1/c.C_couple_F,-1/c.C_couple_F,0.]])
        self.circuit_step_inverse=np.linalg.inv(np.eye(3)-self.dt_s*J)
        self.phonon_energy_J,self.phonon_weights_J,self.phonon_DOS_per_J_m3=self.initial_data["effective_material"]._quadrature()
        self.phonon_occupation_initial=np.zeros((len(self.phonon_energy_J),self.graph.n_nodes))
        xx=self.phonon_energy_J[:,None]/(kB*self.initial_data["common_temperature_K"][None,:])
        np.divide(1.,np.expm1(np.minimum(xx,700)),out=self.phonon_occupation_initial,where=xx<700)
        if np.isnan(tau_escape_s) or tau_escape_s<=0:
            raise ValueError("Positive escape time or infinity for an explicit closed control required")
        self.tau_escape_s=float(tau_escape_s)
        Tb=self.initial_data["inputs"]["summary"]["plan"]["T_K"]
        xb=self.phonon_energy_J[:,None]/(kB*Tb)
        bath=np.zeros_like(xb)
        np.divide(1.,np.expm1(np.minimum(xb,700)),out=bath,where=xb<700)
        self.phonon_bath=np.broadcast_to(bath,self.phonon_occupation_initial.shape).copy()
        # Closed-form Schur solution of the same linear implicit escape row.
        # Explicit e-ph exchange is zero only at the actual local FD/BE input.
        ratio=dt_s/self.tau_escape_s
        self.phonon_occupation=self.phonon_bath+(self.phonon_occupation_initial-self.phonon_bath)/(1+ratio)
        self.phonon_rate=-(self.phonon_occupation-self.phonon_bath)/self.tau_escape_s
        self.phonon_escape_power_W=float(np.sum(
            (self.phonon_energy_J*self.phonon_weights_J*self.phonon_DOS_per_J_m3)[:,None]*
            (self.phonon_occupation-self.phonon_bath)/self.tau_escape_s*
            self.rows.moments.volume[None,:]))

    def unpack(self,x):
        dd,phi,_=super().unpack(x)
        cy=self.circuit_step_inverse@self.circuit.rhs(self.y,phi[self.rows.left[0]]-phi[self.rows.right[0]])
        dd[self.fixed]=self.gap[self.fixed]*(self.bulk.gap_bar_derivative_per_A*cy[1]/self.bulk.gap_bar)
        return dd,phi,cy

    def evaluate(self,x,*,save_fields=False):
        dd,phi,cy=self.unpack(x); dt=self.dt_s
        alpha_t=2*e/hbar*(phi[self.graph.edges[:,1]]-phi[self.graph.edges[:,0]])
        reference=self.thermal_tangent.apply(dd,alpha_direction=alpha_t,
            q_direction=self.bulk.q_bar_derivative_per_A*cy[1])
        K0=[];K1=[];Kt=[];Ktht=[];F0=[];F1=[];Ft=[];Ftht=[];ht=[]
        metrics=[];Rs=[];RE=[];RT=[];R1=[];kinetic=0.;full=0.
        for sl in self.slices:
            R,K,F,dH,m,ke,fu=sl.evaluate(dd,phi,cy[1])
            Hth=ad.distribution(np.full(self.graph.n_nodes,
                np.tanh(sl.energy*self.Eu/(2*self.initial_data["inputs"]["summary"]["plan"]["T_K"]*kB))),
                np.zeros(self.graph.n_nodes))
            _,Fth,Kth=ad.kinetic_residual(self.graph,R,Hth,
                ad.generator(self.gap,sl.energy,d_time=dd),alpha=self.alpha,
                alpha_time=alpha_t,kappa=sl.kappa)
            K0.append(K.value);K1.append(K.first);Kt.append(K.time);Ktht.append(Kth.time)
            F0.append(F.value);F1.append(F.first);Ft.append(F.time);Ftht.append(Fth.time)
            ht.append(dH);metrics.append(m);kinetic=max(kinetic,ke);full=max(full,fu)
            if save_fields:
                Rs.append(R.value);RE.append(R.energy);RT.append(R.time);R1.append(R.first)
        K0,K1,Kt,Ktht,F0,F1,Ft,Ftht,ht=map(np.asarray,(K0,K1,Kt,Ktht,F0,F1,Ft,Ftht,ht))
        value=self.rows.evaluate(K_zero=K0+dt*Kt,K_first=K1,
            K_thermal_zero=self.Kth+dt*Ktht,phi_V=np.zeros_like(phi),
            thermal_gradient_J_per_gapbar=self.Gth+dt*reference.gradient_J_per_gapbar,
            thermal_edge_current_A=self.Ith+dt*reference.edge_current_A,
            edge_flux_zero=F0+dt*Ft,edge_flux_first=F1,thermal_edge_flux_zero=self.Fth+dt*Ftht,
            imposed_current_A=self.y[1]+dt*cy[1],reference_contract=self.contract)
        value=replace(value,port_voltage_V=float(phi[self.rows.left[0]]-phi[self.rows.right[0]]))
        packed=value.pack(neutrality_nodes=self.free,force_nodes=self.free,scales=self.rows.scales)
        fields=None if not save_fields else dict(R0=np.asarray(Rs),RE=np.asarray(RE),Rt=np.asarray(RT),R1=np.asarray(R1),
            K0=K0,K1=K1,K0_time_per_s=Kt,edge_flux_zero=F0,edge_flux_first=F1,
            covariant_gap_rate_per_s=dd,alpha_time_per_s=alpha_t)
        return JointRateValue(packed,dd-2j*e*phi*self.gap/hbar,phi,ht,cy,value,tuple(metrics),kinetic,full,
            float(np.max(abs(value.neutrality_V[self.fixed]),initial=0)),fields,
            self.phonon_rate,
            "First step only: C_e-ph(y0)=0 and dn_e-ph(y0)=0 explicit by local FD/BE; linear phonon escape implicit at supplied tau; no general material admission",
            "Single linearly implicit candidate; not converged backward Euler or admitted trajectory")
