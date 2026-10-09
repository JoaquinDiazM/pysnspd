"""One affine first-gradient gap/kinetic/neutrality/CM system at a state.

Full hL/hT fields are retained at every energy/node. Spectral and kinetic
linear solves are exact Schur eliminations, not temperature/profile closures.
The unknown vector contains two covariant gap rates per interior node and one potential
per interior node plus two equipotential superconducting cuts. Circuit rates
and moving cut amplitudes/q are functions of those same potentials.
This is a snapshot DAE assembler; stage5.2 energy and stage5.3 time admission
are not implied. Sources must be explicitly identified, never fitted to a
residual. Numerical causal eta follows H and is not an inelastic bath.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
from types import SimpleNamespace
import numpy as np
from scipy.constants import elementary_charge as e, hbar
from scipy.sparse.linalg import splu
from . import adiabatic_matrix_usadel as ad
from . import retarded_spatial_usadel as ret
from . import causal_kinetic_regulator as reg
from .first_spectral_response import ChartFirstSpectralResponse
from .joint_nambu_rows import JointNambuRows, checked_gauge_reduction
from .electrical_ports import ThesisCircuitParameters


@dataclass(frozen=True)
class BulkCutLaw:
    gap_bar: float
    q_bar: float
    gap_bar_derivative_per_A: float
    q_bar_derivative_per_A: float
    contract: str

    def __post_init__(self):
        if (not np.all(np.isfinite([self.gap_bar,self.q_bar,
                self.gap_bar_derivative_per_A,self.q_bar_derivative_per_A]))
                or self.gap_bar <= 0 or not self.contract.strip()):
            raise ValueError("Finite superconducting bulk law and provenance required")


@dataclass(frozen=True)
class JointRateValue:
    packed_rows: np.ndarray
    gap_rate_per_s: np.ndarray
    phi_V: np.ndarray
    distribution_rate_per_s: np.ndarray
    circuit_rate_SI_per_s: np.ndarray
    moments: object
    spectral_metrics: tuple
    kinetic_residual_max: float
    full_nambu_residual_max_per_mass: float
    cut_neutrality_max_V: float
    fields: dict | None
    phonon_rate_per_s: np.ndarray | None = None
    material_rhs_contract: str = "No general phonon RHS supplied"
    scope: str = "unified instantaneous candidate; no temporal/energy admission"


class JointSpectralSlice:
    """Reusable spectral/kinetic eliminations at one physical energy."""
    def __init__(self, graph, gap, alpha, R0, energy_J, eta_J, Eu_J,
                 h, h_energy, fixed_nodes, bulk_law, *,
                 projected_source, source_contract, retarded_generator_first=None,
                 full_nambu_source=None,
                 first_tolerance=1e-10, first_maxiter=3000):
        self.graph=graph; self.gap=np.asarray(gap,complex).copy()
        self.alpha=np.asarray(alpha,float).copy(); self.fixed=np.asarray(fixed_nodes,int)
        self.Eu=float(Eu_J); self.energy=float(energy_J)/self.Eu
        self.eta=float(eta_J)/self.Eu; self.kappa=hbar/(2*self.Eu)
        self.h=np.asarray(h,float).copy();self.hE=np.asarray(h_energy,float).copy()
        self.source=np.asarray(projected_source,float).copy()
        n=graph.n_nodes
        if any(v.shape!=(n,2) or np.any(~np.isfinite(v)) for v in (self.h,self.hE,self.source)):
            raise ValueError("Finite two-mode fields/source at every node required")
        if np.max(abs(self.h[:,0]+self.h[:,1]))>1+32*np.finfo(float).eps or np.max(abs(self.h[:,0]-self.h[:,1]))>1+32*np.finfo(float).eps:
            raise ValueError("Initial fields leave Pauli support; no clipping")
        if not isinstance(source_contract,str) or not source_contract.strip():
            raise ValueError("Explicit source contract required")
        self.source_contract=source_contract;self.bulk=bulk_law
        self.full_source=np.zeros((n,2,2),complex) if full_nambu_source is None else np.asarray(full_nambu_source,complex)
        if self.full_source.shape!=(n,2,2) or np.any(~np.isfinite(self.full_source)):
            raise ValueError("Finite complete physical Nambu source required")
        if not np.allclose(ad.projections(self.full_source),self.source,rtol=1e-12,atol=1e-30):
            raise ValueError("Physical trace source must project the supplied full Nambu collision; no inferred anomalous entries")
        R0=np.asarray(R0,complex)
        a,b=R0[:,0,1]/(1+R0[:,0,0]),R0[:,1,0]/(1+R0[:,0,0])
        root=SimpleNamespace(z=self.eta-1j*self.energy,a=a,b=b,g=R0[:,0,0],f=R0[:,0,1],f_tilde=R0[:,1,0])
        self.tangent=ret.RetardedTangent(graph,self.gap,root,self.alpha,fixed_nodes=self.fixed,spectral_tolerance=1e-7)
        self.phase=np.angle(self.gap[self.fixed])
        self.theta=2*np.arctan(a[self.fixed]*np.exp(-1j*self.phase))
        da,db=ret.current_carrying_cut_tangent(self.theta,self.eta-1j*self.energy,
            bulk_law.gap_bar,bulk_law.q_bar,self.phase,z_direction=-1j)
        self.RE=self.tangent.apply(z_direction=-1j,fixed_a_direction=da,
            fixed_b_direction=db).matrix
        self.Rbase=ad.Jet.of(R0,self.RE)
        self.Bbase=ad.generator(self.gap,self.energy,eta=self.eta)
        self.spectral=ChartFirstSpectralResponse(graph,self.Rbase,self.Bbase,
            self.alpha,self.fixed,kappa_s=self.kappa,retarded_tangent=self.tangent)
        self.SigmaB=np.zeros_like(R0) if retarded_generator_first is None else np.asarray(retarded_generator_first,complex)
        if self.SigmaB.shape!=R0.shape or np.any(~np.isfinite(self.SigmaB)):
            raise ValueError("Finite full first-order retarded generator required")
        self.first_tolerance=first_tolerance;self.first_maxiter=first_maxiter
        M=reg.assemble_storage(graph,self.Rbase,ad.generator(self.gap,self.energy),
            numerical_eta=self.eta,kappa_s=self.kappa,alpha=self.alpha)
        free=np.ones(n,bool);free[self.fixed]=False
        self.free=free;self.components=np.flatnonzero(np.repeat(free,2))
        self.storage=M
        # A singular storage is exposed; eta is not increased to make it invertible.
        self.storage_factor=splu(M[self.components][:,self.components].tocsc()) if len(self.components) else None

    def configure_linearized_step(self, dt_s):
        """Prepare M + dt*Lh, with the physical zeroth kinetic Jacobian.

        This is Rosenbrock-Euler about this state, not a converged nonlinear
        backward-Euler trajectory. Initial detailed-balance collisions can
        be explicit for its first step only. No physical damping is added.
        """
        if not np.isfinite(dt_s) or dt_s < 0:
            raise ValueError("Nonnegative finite linearized step required")
        self.linearized_step_s = float(dt_s)
        if dt_s == 0:
            self.step_storage = self.storage
            self.step_factor = self.storage_factor
        else:
            from .joint_rosenbrock import assemble_zeroth_kinetic_operator
            if not hasattr(self, "zeroth_operator"):
                self.zeroth_operator = assemble_zeroth_kinetic_operator(
                    self.graph, self.Rbase, ad.generator(self.gap, self.energy),
                    numerical_eta=self.eta, kappa_s=self.kappa, alpha=self.alpha)
            self.step_storage = self.storage + dt_s*self.zeroth_operator
            self.step_factor = splu(self.step_storage[self.components][:,self.components].tocsc()) if len(self.components) else None

    def evaluate(self, covariant_gap_rate_per_s, phi_V, current_rate_A_per_s):
        g=self.graph;fixed=self.fixed;b=self.bulk
        dd=b.gap_bar_derivative_per_A*current_rate_A_per_s
        dq=b.q_bar_derivative_per_A*current_rate_A_per_s
        # Instantaneous local temporal gauge (Wilson time links):
        # phi'=0, Dt d is the gap derivative and the electric field appears
        # in alpha_t'=2e(phi_head-phi_tail)/hbar. This defines hL/hT relative
        # to the local electrochemical coordinate, not a bath temperature.
        phase_rate=np.zeros(len(fixed))
        alpha_time=2*e/hbar*(np.asarray(phi_V)[g.edges[:,1]]-
                            np.asarray(phi_V)[g.edges[:,0]])
        da,db=ret.current_carrying_cut_tangent(self.theta,self.eta-1j*self.energy,
            b.gap_bar,b.q_bar,self.phase,gap_direction=dd,q_direction=dq,
            phase_direction=phase_rate)
        Rt=self.tangent.apply(gap_direction=covariant_gap_rate_per_s,alpha_direction=alpha_time,
            fixed_a_direction=da,fixed_b_direction=db).matrix
        raw=ad.Jet.of(self.Rbase.value,self.RE,Rt)
        Bs=ad.generator(self.gap,self.energy,d_time=covariant_gap_rate_per_s,eta=self.eta)
        Bs=ad.Jet.of(Bs.value,Bs.energy,Bs.time,self.SigmaB)
        # Exterior continuation: the same current-carrying bulk equation,
        # with moving q, amplitude and phase, supplies all four R1 entries.
        cut_R=raw.take(fixed);cut_B=Bs.take(fixed)
        completed_cut,cut_metrics=ad.bulk_first_spectral(cut_R,cut_B,b.q_bar,
            q_time_per_s=dq,kappa_s=self.kappa,tolerance=self.first_tolerance)
        R,metrics=self.spectral.apply(raw,Bs,alpha_time_per_s=alpha_time,fixed_first=completed_cut.first,
            tolerance=self.first_tolerance,maxiter=self.first_maxiter)
        B=ad.generator(self.gap,self.energy,d_time=covariant_gap_rate_per_s)
        H0=ad.distribution(self.h[:,0],self.h[:,1],energy_L=self.hE[:,0],energy_T=self.hE[:,1])
        base,_,_,_=reg.kinetic_residual(g,R,H0,B,numerical_eta=self.eta,
            kappa_s=self.kappa,alpha=self.alpha,alpha_time_per_s=alpha_time)
        dt = getattr(self, "linearized_step_s", 0.)
        rhs=(self.source-ad.projections(base.value+base.first+dt*base.time)).ravel()
        ht=np.zeros_like(self.h)
        # The cuts continue the thermal external superconducting bank;
        # their distribution function is prescribed, not inferred from Te.
        factor = getattr(self, "step_factor", self.storage_factor)
        if factor is not None:
            ht.ravel()[self.components]=factor.solve(rhs[self.components])
        if np.any(~np.isfinite(ht)):
            raise FloatingPointError("Nonfinite eliminated population derivative")
        H=ad.distribution(self.h[:,0],self.h[:,1],energy_L=self.hE[:,0],energy_T=self.hE[:,1],
            time_L=ht[:,0],time_T=ht[:,1])
        residual,flux,K,_=reg.kinetic_residual(g,R,H,B,numerical_eta=self.eta,
            kappa_s=self.kappa,alpha=self.alpha,alpha_time_per_s=alpha_time)
        pr=ad.projections(residual.value+residual.first+dt*residual.time)-self.source
        err=float(np.max(abs(pr[self.free]),initial=0))
        # Full anomalous residual is retained for a Ward test, not discarded.
        full=float(np.max(abs((residual.value+residual.first+dt*residual.time-self.full_source)[self.free]/
            g.area_weights[self.free,None,None]),initial=0))
        metrics=dict(metrics,cut_metrics=cut_metrics)
        return R,K,flux,ht,metrics,err,full


class JointSnapshotSystem:
    """Square Schur system for gap rates, neutrality and both CM ports.

    Fields use the local electrochemical/Wilson energy coordinate. The
    common voltage is an exact gauge mode, checked before eliminating it. Full-cut neutrality and interior current
    divergence are independent diagnostics of continuation/Ward/UV accuracy.
    """
    def __init__(self, graph, gap_bar, alpha, slices, rows:JointNambuRows, *,
            K_thermal_zero, thermal_gradient_J_per_gapbar,thermal_edge_current_A,
            thermal_edge_flux_zero,circuit:ThesisCircuitParameters,circuit_state,
            bulk_law:BulkCutLaw,reference_contract,phonon_rate_per_s=None,material_rhs_contract=None):
        self.graph=graph;self.gap=np.asarray(gap_bar,complex).copy()
        self.alpha=np.asarray(alpha,float).copy();self.slices=tuple(slices);self.rows=rows
        self.fixed=np.concatenate((rows.left,rows.right))
        free=np.ones(graph.n_nodes,bool);free[self.fixed]=False
        self.free=np.flatnonzero(free);self.nf=len(self.free);self.size=3*self.nf+2
        if len(self.slices)!=len(rows.moments.energy):
            raise ValueError("One kinetic slice per physical quadrature energy required")
        self.Kth=np.asarray(K_thermal_zero,complex)
        self.Gth=np.asarray(thermal_gradient_J_per_gapbar,complex)
        self.Ith=np.asarray(thermal_edge_current_A,float)
        self.Fth=np.asarray(thermal_edge_flux_zero,complex)
        self.circuit=circuit;self.y=np.asarray(circuit_state,float).copy();self.bulk=bulk_law
        if self.y.shape!=(3,) or not np.all(np.isfinite(self.y)):
            raise ValueError("All three finite CM states required")
        self.contract=reference_contract
        self.phonon_rate=None if phonon_rate_per_s is None else np.asarray(phonon_rate_per_s,float).copy()
        self.material_contract="No general phonon RHS supplied" if material_rhs_contract is None else str(material_rhs_contract)
        if self.phonon_rate is not None:
            if (self.phonon_rate.ndim!=2 or self.phonon_rate.shape[1]!=graph.n_nodes
                    or np.any(~np.isfinite(self.phonon_rate)) or not self.material_contract.strip()
                    or material_rhs_contract is None):
                raise ValueError("Finite full phonon rate and shared material collision contract required")
        self.Eu=rows.moments.Eu;self.kappa=hbar/(2*self.Eu)
        self.Vunit=self.Eu/e

    def unpack(self,x):
        x=np.asarray(x,float)
        if x.shape!=(self.size,) or np.any(~np.isfinite(x)):
            raise ValueError("One finite full joint unknown vector required")
        dd=np.zeros(self.graph.n_nodes,complex)
        dd[self.free]=(x[:2*self.nf:2]+1j*x[1:2*self.nf:2])/self.kappa
        phi=np.zeros(self.graph.n_nodes)
        phi[self.free]=x[2*self.nf:3*self.nf]*self.Vunit
        phi[self.rows.left]=x[-2]*self.Vunit;phi[self.rows.right]=x[-1]*self.Vunit
        c_rate=self.circuit.rhs(self.y,(x[-2]-x[-1])*self.Vunit)
        amp=self.bulk.gap_bar_derivative_per_A*c_rate[1]/self.bulk.gap_bar
        dd[self.fixed]=self.gap[self.fixed]*amp
        return dd,phi,c_rate

    def evaluate(self,x,*,save_fields=False):
        dd,phi,cy=self.unpack(x)
        K0=[];K1=[];F0=[];F1=[];ht=[];metrics=[];Rs=[];RE=[];RT=[];R1=[]
        kinetic=0.;full=0.
        for sl in self.slices:
            R,K,F,dH,m,ke,fu=sl.evaluate(dd,phi,cy[1])
            K0.append(K.value);K1.append(K.first);F0.append(F.value);F1.append(F.first)
            ht.append(dH);metrics.append(m);kinetic=max(kinetic,ke);full=max(full,fu)
            if save_fields:
                Rs.append(R.value);RE.append(R.energy);RT.append(R.time);R1.append(R.first)
        K0,K1,F0,F1,ht=map(np.asarray,(K0,K1,F0,F1,ht))
        value=self.rows.evaluate(K_zero=K0,K_first=K1,K_thermal_zero=self.Kth,
            phi_V=np.zeros_like(phi),thermal_gradient_J_per_gapbar=self.Gth,
            thermal_edge_current_A=self.Ith,edge_flux_zero=F0,edge_flux_first=F1,
            thermal_edge_flux_zero=self.Fth,imposed_current_A=self.y[1],
            reference_contract=self.contract)
        value=replace(value,port_voltage_V=float(np.mean(phi[self.rows.left])-np.mean(phi[self.rows.right])),
            left_equipotential_residual_V=phi[self.rows.left]-np.mean(phi[self.rows.left]),
            right_equipotential_residual_V=phi[self.rows.right]-np.mean(phi[self.rows.right]))
        packed=value.pack(neutrality_nodes=self.free,force_nodes=self.free,
            scales=self.rows.scales)
        fields=None if not save_fields else dict(R0=np.asarray(Rs),RE=np.asarray(RE),
            Rt=np.asarray(RT),R1=np.asarray(R1),K0=K0,K1=K1,
            edge_flux_zero=F0,edge_flux_first=F1,covariant_gap_rate_per_s=dd)
        actual_dd=dd-2j*e*phi*self.gap/hbar
        return JointRateValue(packed,actual_dd,phi,ht,cy,value,tuple(metrics),kinetic,full,
            float(np.max(abs(value.neutrality_V[self.fixed]),initial=0)),fields,
            self.phonon_rate,self.material_contract)

    @property
    def gauge_direction(self):
        result=np.zeros(self.size)
        result[2*self.nf:]=1.
        return result

    def affine_matrix(self,*,progress=None):
        base=self.evaluate(np.zeros(self.size)).packed_rows
        matrix=np.empty((self.size,self.size))
        for j in range(self.size):
            v=np.zeros(self.size);v[j]=1
            matrix[:,j]=self.evaluate(v).packed_rows-base
            if progress is not None:progress(j+1,self.size)
        return matrix,base


def solve_affine_candidate(matrix,base,*,relative_rank_tolerance=1e-7,
                           residual_tolerance=1e-7,gauge_direction=None,gauge_reference_index=None):
    """Expose all singular values; never add friction or discard a physical row.

    Unit-scaled inputs are mandatory. A minimum-norm candidate is returned for
    diagnostic inspection when deficient, but it cannot be marked admitted.
    Re-evaluation through JointSnapshotSystem is required after this solve.
    """
    A=np.asarray(matrix,float);b=np.asarray(base,float)
    if A.ndim!=2 or A.shape[0]!=A.shape[1] or b.shape!=(len(A),) or np.any(~np.isfinite(A)) or np.any(~np.isfinite(b)):
        raise ValueError("Finite square SI-scaled affine system required")
    if (not 0 < relative_rank_tolerance < 1 or not np.isfinite(relative_rank_tolerance)
            or residual_tolerance <= 0 or not np.isfinite(residual_tolerance)):
        raise ValueError("Positive explicit numerical accuracy/rank limits required")
    original=A
    reduction=None
    if gauge_direction is not None:
        reduction=checked_gauge_reduction(A,np.asarray(gauge_direction,float),
            reference_index=gauge_reference_index,row_scales=np.ones(len(A)),
            column_scales=np.ones(A.shape[1]),tolerance=1e-10)
        A=np.asarray(reduction.matrix)
    U,s,Vh=np.linalg.svd(A,full_matrices=False)
    limit=relative_rank_tolerance*s[0] if len(s) else 0
    keep=s>limit;coeff=np.zeros_like(s)
    np.divide(-(U.T@b),s,out=coeff,where=keep)
    reduced_x=Vh.T@coeff
    x=reduced_x if reduction is None else reduction.expand(reduced_x)
    residual=original@x+b
    rank=int(np.count_nonzero(keep))
    ok=rank==A.shape[1] and np.max(abs(residual),initial=0)<=residual_tolerance
    return x,dict(rank=rank,size=original.shape[1],physical_columns=A.shape[1],
        singular_values=s,scaled_residual_max=float(np.max(abs(residual),initial=0)),
        condition_number=float(s[0]/s[-1]) if s[-1]>0 and np.isfinite(s[0]/s[-1]) else None,
        relative_rank_tolerance=relative_rank_tolerance,
        algebraically_solved_full_rank=bool(ok),temporal_admission=False,
        energy_admission=False,no_regularization=True,
        gauge_reduced=bool(reduction is not None),
        gauge_ward_max=None if reduction is None else reduction.scaled_ward_residual,
        gauge_reference_index=gauge_reference_index,
        no_physical_rows_removed=True,
        rank_scope="Numerically resolved at declared spectral/affine accuracy; no claim of exact physical rank")
