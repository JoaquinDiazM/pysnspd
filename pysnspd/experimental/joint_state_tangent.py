"""Physical zeroth-state tangents for a linearly implicit joint detector step.

The derivative here is an ordinary state differential, not a Moyal correction:
dK0 = dR0 H0 + R0 dH0 - dH0 A0 - H0 dA0.
The corresponding link-flux and local-commutator derivatives retain moving
electromagnetic links. No gap friction, Ohmic current, DOS floor, collision
source, or time-step admission is supplied by this module.

The Matsubara bundle differentiates the SAME finite thermal reference, with
the actual current-carrying superconducting continuation at both cuts. It
does not replace populations by temperature, prescribe an interior profile,
or use zero-current BCS cuts. All directions can be rates per SI second or
ordinary variations, provided the caller uses one consistent convention.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy.sparse.linalg import splu
from . import adiabatic_matrix_usadel as ad
from . import thermal_spatial_usadel as thermal
from .bulk_current_reference import homogeneous_spectrum


def _real(values,shape,name):
    raw=np.asarray(values)
    if (raw.shape!=shape or np.iscomplexobj(raw) or raw.dtype.kind=="b"
            or np.any(~np.isfinite(raw))):
        raise ValueError(name+" requires finite real fields of shape "+str(shape))
    return np.asarray(raw,float)


def _complex(values,shape,name):
    result=np.asarray(values,complex)
    if result.shape!=shape or np.any(~np.isfinite(result)):
        raise ValueError(name+" requires finite complex fields of shape "+str(shape))
    return result


@dataclass(frozen=True)
class ZerothKineticDirection:
    K_zero: np.ndarray
    K_direction: np.ndarray
    flux_zero: np.ndarray
    flux_direction: np.ndarray
    residual_zero: np.ndarray | None
    residual_direction: np.ndarray | None
    scope: str="ordinary zeroth state differential; no independent first-gradient correction"


def zeroth_kinetic_direction(graph,retarded_zero,retarded_direction,h,h_direction,*,
        alpha=None,alpha_direction=None,B_zero=None,B_direction=None):
    """Differentiate the full zeroth Nambu flux, with optional physical local B.

    Inputs R,dR have shape(nodes,2,2), h,dh shape(nodes,2). A supplied B must
    exclude numerical causal eta, as in the physical kinetic equation.
    Derivatives stored in Jet.time are used ONLY as ordinary differentials;
    kappa=0 removes every Moyal product. R1 and energy derivatives do not enter.
    """
    n=graph.n_nodes;shape=(n,2,2)
    R=_complex(retarded_zero,shape,"retarded_zero")
    dR=_complex(retarded_direction,shape,"retarded_direction")
    H=_real(h,(n,2),"h");dH=_real(h_direction,(n,2),"h_direction")
    links=np.zeros(len(graph.edges)) if alpha is None else _real(
        alpha,(len(graph.edges),),"alpha")
    dl=np.zeros(len(graph.edges)) if alpha_direction is None else _real(
        alpha_direction,(len(graph.edges),),"alpha_direction")
    if (B_zero is None)!=(B_direction is None):
        raise ValueError("Supply both physical B_zero and B_direction or neither")
    supplied=B_zero is not None
    B=np.zeros(shape,complex) if not supplied else _complex(B_zero,shape,"B_zero")
    dB=np.zeros(shape,complex) if not supplied else _complex(B_direction,shape,"B_direction")
    Rjet=ad.Jet.of(R,time=dR)
    Hjet=ad.distribution(H[:,0],H[:,1],time_L=dH[:,0],time_T=dH[:,1])
    residual,flux,K=ad.kinetic_residual(graph,Rjet,Hjet,ad.Jet.of(B,time=dB),
        alpha=links,alpha_time=dl,kappa=0.)
    return ZerothKineticDirection(K.value,K.time,flux.value,flux.time,
        residual.value if supplied else None,residual.time if supplied else None)


def joint_zeroth_rows_direction(rows,graph,*,retarded_zero,retarded_direction,
        h,h_direction,thermal_h,thermal_gradient_direction_J_per_gapbar,
        thermal_edge_current_direction_A,reference_contract,thermal_h_direction=None,
        alpha=None,alpha_direction=None,phi_direction_V=None,
        imposed_current_direction_A=0.):
    """Integrate a full state differential of G0, Q0 and completed I0 in SI.

    The thermal spectrum equals the actual R0 at the instantaneous same gap,
    links and gauge. The thermal H derivative defaults to zero only when its
    bath reference is fixed in that coordinate. No diagonal-voltage term is
    silently inserted for a caller using the local Wilson energy coordinate.
    """
    ne,n=len(rows.moments.energy),graph.n_nodes
    R=_complex(retarded_zero,(ne,n,2,2),"retarded_zero")
    dR=_complex(retarded_direction,R.shape,"retarded_direction")
    H=_real(h,(ne,n,2),"h");dH=_real(h_direction,H.shape,"h_direction")
    Hth=_real(thermal_h,H.shape,"thermal_h")
    dHth=np.zeros_like(Hth) if thermal_h_direction is None else _real(
        thermal_h_direction,H.shape,"thermal_h_direction")
    phi=np.zeros(n) if phi_direction_V is None else _real(
        phi_direction_V,(n,),"phi_direction_V")
    actual=[];reference=[]
    for j in range(ne):
        actual.append(zeroth_kinetic_direction(graph,R[j],dR[j],H[j],dH[j],
            alpha=alpha,alpha_direction=alpha_direction))
        reference.append(zeroth_kinetic_direction(graph,R[j],dR[j],Hth[j],dHth[j],
            alpha=alpha,alpha_direction=alpha_direction))
    Kdot=np.asarray([item.K_direction for item in actual])
    Fdot=np.asarray([item.flux_direction for item in actual])
    Kthdot=np.asarray([item.K_direction for item in reference])
    Fthdot=np.asarray([item.flux_direction for item in reference])
    return rows.linear_direction(K_first_direction=np.zeros_like(Kdot),
        edge_flux_first_direction=np.zeros_like(Fdot),phi_direction_V=phi,
        K_zero_direction=Kdot,K_thermal_zero_direction=Kthdot,
        edge_flux_zero_direction=Fdot,thermal_edge_flux_zero_direction=Fthdot,
        thermal_gradient_direction_J_per_gapbar=thermal_gradient_direction_J_per_gapbar,
        thermal_edge_current_direction_A=thermal_edge_current_direction_A,
        imposed_current_direction_A=imposed_current_direction_A,
        reference_contract=reference_contract)


@dataclass(frozen=True)
class MatsubaraStateDirection:
    du: np.ndarray
    df: np.ndarray
    dg: np.ndarray
    current_bar_direction: np.ndarray
    residual_direction: np.ndarray
    free_equation_residual_max: float


class MatsubaraStateTangent:
    """Exact Cartesian sparse tangent with arbitrary prescribed spectral cuts."""
    def __init__(self,graph,gap_bar,epsilon_bar,u,alpha=None,*,fixed_nodes=None,
            spectral_tolerance=1e-7):
        self.graph=graph;n=graph.n_nodes
        self.gap=_complex(gap_bar,(n,),"gap_bar").copy()
        self.u=_complex(u,(n,),"u").copy()
        self.epsilon=float(epsilon_bar)
        if not np.isfinite(self.epsilon) or self.epsilon<=0:
            raise ValueError("Positive Matsubara frequency required")
        self.alpha=np.zeros(len(graph.edges)) if alpha is None else _real(
            alpha,(len(graph.edges),),"alpha").copy()
        fixed=graph.boundary_nodes if fixed_nodes is None else fixed_nodes
        if fixed is None:raise ValueError("Explicit superconducting cut nodes required")
        self.fixed=thermal._fixed_indices(fixed,n)
        self.free=np.ones(n,bool);self.free[self.fixed]=False
        self.components=np.flatnonzero(np.repeat(self.free,2))
        self.fixed_components=np.flatnonzero(np.repeat(~self.free,2))
        self.value,self.jacobian=thermal.spectral_residual_jacobian(
            graph,self.gap,self.epsilon,self.u,self.alpha)
        tol=float(spectral_tolerance)
        if not np.isfinite(tol) or tol<=0:raise ValueError("Positive spectral_tolerance required")
        error=float(np.max(abs(self.value.residual[self.free])/
            (graph.area_weights[self.free]*np.maximum(1.,abs(self.gap[self.free]))),initial=0))
        if error>tol:raise ValueError("Thermal reference is not stationary on its free nodes")
        self.reference_spectral_residual=error
        self.factor=splu(self.jacobian[self.components][:,self.components].tocsc()) if len(self.components) else None

    def apply(self,gap_direction,*,alpha_direction=None,fixed_u_direction=None):
        n=self.graph.n_nodes
        dd=_complex(gap_direction,(n,),"gap_direction")
        da=np.zeros(len(self.graph.edges)) if alpha_direction is None else _real(
            alpha_direction,(len(self.graph.edges),),"alpha_direction")
        if len(self.fixed) and fixed_u_direction is None:
            raise ValueError("Explicit moving spectral cut direction required, including a declared zero")
        fixed=np.empty(0,complex) if not len(self.fixed) else _complex(
            fixed_u_direction,(len(self.fixed),),"fixed_u_direction")
        source=self.graph.area_weights*dd
        tail,head=self.graph.edges.T;phase=np.exp(-1j*self.alpha)
        # R_alpha = +i*c*P*f_j at tail, -i*c*P^dag*f_i at head.
        np.add.at(source,tail,-1j*self.graph.conductance*phase*self.value.f[head]*da)
        np.add.at(source,head,+1j*self.graph.conductance*phase.conj()*self.value.f[tail]*da)
        rhs=np.column_stack((source.real,source.imag)).ravel()
        answer=np.zeros(2*n)
        answer[self.fixed_components]=np.column_stack((fixed.real,fixed.imag)).ravel()
        if self.factor is not None:
            right=rhs[self.components]-self.jacobian[self.components][:,self.fixed_components]@answer[self.fixed_components]
            answer[self.components]=self.factor.solve(right)
        du=answer[0::2]+1j*answer[1::2]
        dg=-self.value.g**3*np.real(self.u.conj()*du)
        df=self.value.g*du+self.u*dg
        current=2*self.graph.conductance*np.imag(
            df[tail].conj()*phase*self.value.f[head]+
            self.value.f[tail].conj()*phase*df[head]-
            1j*da*self.value.f[tail].conj()*phase*self.value.f[head])
        real_residual=self.jacobian@answer-rhs
        residual=real_residual[0::2]+1j*real_residual[1::2]
        return MatsubaraStateDirection(du,df,dg,current,residual,
            float(np.max(abs(residual[self.free]),initial=0)))


def current_carrying_thermal_cut_direction(epsilon_bar,gap_boundary,
        gap_boundary_direction,q_bar,*,q_direction=0.):
    """du_cut from epsilon*u+q^2*u/sqrt(1+u^2)=|d_cut| and its phase.

    Each cut can have its own phase direction. Its finite amplitude is needed;
    no arg(d) is evaluated at an interior zero, where this boundary law is unused.
    """
    d=np.asarray(gap_boundary,complex)
    dd=_complex(gap_boundary_direction,d.shape,"gap_boundary_direction")
    q,dq=float(q_bar),float(q_direction)
    ep=float(epsilon_bar)
    if (np.any(~np.isfinite(d)) or np.any(abs(d)<=0) or not np.isfinite(q)
            or not np.isfinite(dq) or not np.isfinite(ep) or ep<=0):
        raise ValueError("Finite nonzero superconducting boundary gaps and frequency required")
    out=np.empty_like(d)
    for index in np.ndindex(d.shape):
        amplitude=abs(d[index]);phase=d[index]/amplitude
        u,f,g=(item[0] for item in homogeneous_spectrum(amplitude,np.array([ep]),q*q))
        radial=float(np.real(d[index].conj()*dd[index])/amplitude)
        angular=float(np.imag(d[index].conj()*dd[index])/amplitude**2)
        up=(radial-2*q*f*dq)/(ep+q*q*g**3)
        out[index]=phase*(up+1j*u*angular)
    return out


@dataclass(frozen=True)
class ThermalReferenceDirection:
    gradient_J_per_gapbar: np.ndarray
    edge_current_A: np.ndarray
    maximum_spectral_residual: float
    modes: int
    scope: str="same finite Matsubara reference derivative; no analytic tail or KWT term"

    @property
    def gradient_direction_J_per_gapbar(self):
        return self.gradient_J_per_gapbar

    @property
    def current_direction_A(self):
        return self.edge_current_A


class ThermalBulkTangentBundle:
    """Resident finite Matsubara Hessian/current tangent with moving bulk cuts.

    energy_scale_J converts the graph free energy to J. current_unit_A is I0,
    the same Eu/(2eRsheet) factor used by the independent Matsubara baseline.
    The caller controls memory/CPU for constructing this resident factor bundle.
    """
    def __init__(self,graph,gap_bar,temperature_ratio,epsilon_bar,u_modes,
            alpha=None,*,fixed_nodes=None,bulk_q_bar,energy_scale_J,current_unit_A,
            spectral_tolerance=1e-7,cut_tolerance=1e-10):
        self.graph=graph
        self.gap=_complex(gap_bar,(graph.n_nodes,),"gap_bar").copy()
        self.t=float(temperature_ratio);self.q=float(bulk_q_bar)
        self.scale=float(energy_scale_J);self.I0=float(current_unit_A)
        ep=np.asarray(epsilon_bar,float);u=np.asarray(u_modes,complex)
        if (not np.isfinite(self.t) or self.t<=0 or not np.isfinite(self.q)
                or not np.isfinite(self.scale) or self.scale<=0 or not np.isfinite(self.I0)
                or self.I0<=0 or ep.ndim!=1 or not len(ep)
                or u.shape!=(len(ep),graph.n_nodes) or np.any(~np.isfinite(u))):
            raise ValueError("Finite complete Matsubara modes and positive SI scales required")
        expected=2*np.pi*self.t*(np.arange(len(ep))+.5)
        if not np.allclose(ep,expected,rtol=1e-13,atol=0):
            raise ValueError("Consecutive Matsubara frequencies of the same thermal reference required")
        self.epsilon=ep.copy()
        self.modes=tuple(MatsubaraStateTangent(graph,self.gap,epsilon,field,alpha,
            fixed_nodes=fixed_nodes,spectral_tolerance=spectral_tolerance)
            for epsilon,field in zip(ep,u))
        self.fixed=self.modes[0].fixed
        if not np.isfinite(cut_tolerance) or cut_tolerance<=0:
            raise ValueError("Positive cut_tolerance required")
        amplitudes=abs(self.gap[self.fixed])
        if np.any(amplitudes<=0):
            raise ValueError("Thermal continuation requires finite nonzero gap on its cuts")
        self.cut_phase=self.gap[self.fixed]/amplitudes
        self.cut_amplitudes=amplitudes
        self.cut_coefficients=[]
        for mode in self.modes:
            # The monotone homogeneous equation has a unique nonnegative
            # root. Check that equation directly on the supplied fixed field,
            # then cache its exact linear coefficients once. No scalar roots
            # are repeated for every later state/Jacobian direction.
            radial=mode.u[self.fixed]/self.cut_phase
            ur=radial.real
            gc=1/np.hypot(1.,ur)
            error=mode.epsilon*radial+self.q*self.q*radial*gc-amplitudes
            if (np.any(ur<0) or np.max(abs(radial.imag),initial=0)>cut_tolerance
                    or np.max(abs(error)/np.maximum(1.,amplitudes),initial=0)>cut_tolerance):
                raise ValueError("Thermal cuts do not match the declared current-carrying continuation")
            response=1/(mode.epsilon+self.q*self.q*gc**3)
            self.cut_coefficients.append((ur,ur*gc,response))
        self.cut_coefficients=tuple(self.cut_coefficients)

    def apply(self,gap_direction,*,alpha_direction=None,q_direction=0.):
        dd=_complex(gap_direction,(self.graph.n_nodes,),"gap_direction")
        gradient=2*self.graph.area_weights*dd*np.log(self.t)
        current=np.zeros(len(self.graph.edges))
        error=0.
        dq=float(q_direction)
        if not np.isfinite(dq):raise ValueError("Finite bulk q direction required")
        radial=np.real(self.gap[self.fixed].conj()*dd[self.fixed])/self.cut_amplitudes
        angular=np.imag(self.gap[self.fixed].conj()*dd[self.fixed])/self.cut_amplitudes**2
        for mode,(uc,fc,response) in zip(self.modes,self.cut_coefficients):
            boundary=self.cut_phase*((radial-2*self.q*fc*dq)*response+1j*uc*angular)
            item=mode.apply(dd,alpha_direction=alpha_direction,fixed_u_direction=boundary)
            gradient += 4*np.pi*self.t*self.graph.area_weights*(dd/mode.epsilon-item.df)
            current += 2*np.pi*self.t*item.current_bar_direction
            error=max(error,item.free_equation_residual_max)
        return ThermalReferenceDirection(self.scale*gradient,self.I0*current,error,len(self.modes))
