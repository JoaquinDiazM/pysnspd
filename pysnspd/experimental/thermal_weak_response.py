"""Full-node tangents of the admitted thermal Usadel/KWT/normal-potential system.

No nonthermal energy law is supplied. Spectral factors are reused for arbitrary
Cartesian gap directions; probe profiles never restrict the physical DOFs.
The gap reference need not be stationary: the baseline RHS and all derivatives
of KWT mobility and covariant phase velocity are retained.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from .energy_catalog import HBAR_J_S as hbar, K_B_J_K as k
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
from . import thermal_spatial_usadel as thermal


def _vector(value, size, name):
    result = np.asarray(value, complex)
    if result.shape != (size,) or np.any(~np.isfinite(result)):
        raise ValueError(name+' requires one finite complex value per node')
    return result


def divergence(graph, current):
    result = np.zeros(graph.n_nodes, dtype=np.asarray(current).dtype)
    np.add.at(result, graph.edges[:, 0], current)
    np.add.at(result, graph.edges[:, 1], -np.asarray(current))
    return result


@dataclass(frozen=True)
class SpectralDirection:
    du: np.ndarray
    df: np.ndarray
    dg: np.ndarray
    current_derivative: np.ndarray
    equation_residual: float


class SpectralTangent:
    """Factor the exact real residual Jacobian once, then apply J^-1 m dd."""
    def __init__(self, graph, d, solution, alpha=None):
        self.graph, self.d, self.solution = graph, _vector(d, graph.n_nodes, 'd'), solution
        self.alpha = thermal._alpha(graph, alpha)
        if not solution.converged or solution.source_signature != thermal._signature(graph, self.d, self.alpha):
            raise ValueError('A converged spectral solution of these fields is required')
        self.free = np.ones(graph.n_nodes, bool); self.free[solution.fixed_nodes] = False
        self.components = np.flatnonzero(np.repeat(self.free, 2))
        value, jac = thermal.spectral_residual_jacobian(graph, self.d, solution.epsilon, solution.u, self.alpha)
        self.f, self.g, self.u = value.f, value.g, solution.u
        self.jacobian = jac[self.components][:, self.components].tocsc()
        self.factor = splu(self.jacobian) if len(self.components) else None

    def apply(self, direction):
        direction = _vector(direction, self.graph.n_nodes, 'direction')
        if np.any(direction[~self.free] != 0):
            raise ValueError('This fixed-contact tangent requires zero gap directions on contacts')
        rhs = self.graph.area_weights*direction
        real_rhs = np.column_stack((rhs.real, rhs.imag)).ravel()[self.components]
        real_du = np.zeros(2*self.graph.n_nodes)
        if self.factor is not None:
            real_du[self.components] = self.factor.solve(real_rhs)
        du = real_du[0::2]+1j*real_du[1::2]
        projection = np.real(np.conj(self.u)*du)
        dg = -self.g**3*projection
        df = self.g*du+self.u*dg
        tail, head = self.graph.edges.T
        phase = np.exp(-1j*self.alpha)
        # Minus derivative of the spectral action with respect to the link.
        current = 2*self.graph.conductance*np.imag(
            np.conj(df[tail])*phase*self.f[head]+np.conj(self.f[tail])*phase*df[head])
        error = self.jacobian@real_du[self.components]-real_rhs
        return SpectralDirection(du, df, dg, current, float(np.max(abs(error), initial=0.)))


def thermal_hessian_action(graph, temperature_ratio, directions, spectral_directions, epsilons):
    """Assemble any number of full-node H*direction actions in complex notation."""
    values = np.asarray(directions, complex)
    if values.ndim == 1: values = values[None, :]
    if values.ndim != 2 or values.shape[1] != graph.n_nodes or np.any(~np.isfinite(values)):
        raise ValueError('Directions must have shape (directions,nodes)')
    t = float(temperature_ratio)
    if not np.isfinite(t) or t <= 0 or len(spectral_directions) != len(epsilons):
        raise ValueError('Positive temperature and one spectral bundle per frequency required')
    out = 2*graph.area_weights*values*np.log(t)
    currents = np.zeros((len(values),len(graph.edges)))
    for n,(bundle,epsilon) in enumerate(zip(spectral_directions,epsilons)):
        if len(bundle) != len(values) or not np.isclose(epsilon,2*np.pi*t*(n+.5),rtol=1e-13,atol=0):
            raise ValueError('Consecutive frequencies and identical direction counts required')
        for index,item in enumerate(bundle):
            out[index] += 4*np.pi*t*graph.area_weights*(values[index]/epsilon-item.df)
            currents[index] += 2*np.pi*t*item.current_derivative
    return out,currents


class ThermalKWTNormal:
    """Full-node inherited KWT law with the existing normal Ohmic graph.

    d=Delta/(kBTc), tau=t/[hbar/(2kBTc)], v=2e*phi/(kBTc).
    Fixed gap and zero potential contacts are constraints of this core control.
    The same physical KWT times at Tc are retained; there is no fitted rate.
    """
    def __init__(self, graph, d, *, Tc_K, T_K, tau_ee_Tc_ps=.50, tau_ep_Tc_ps=2.47):
        self.graph, self.d = graph, _vector(d,graph.n_nodes,'d').copy()
        values = np.array([Tc_K,T_K,tau_ee_Tc_ps,tau_ep_Tc_ps],float)
        if np.any(~np.isfinite(values)) or np.any(values<=0):
            raise ValueError('Finite positive temperatures and KWT times required')
        if graph.boundary_nodes is None or not len(graph.boundary_nodes):
            raise ValueError('Explicit fixed gap/potential contacts required')
        self.free=np.ones(graph.n_nodes,bool);self.free[graph.boundary_nodes]=False
        ratio=T_K/Tc_K
        self.tD_ps=hbar/(2*k*Tc_K)*1e12
        self.taupsi_ps=1/(ratio/tau_ee_Tc_ps+ratio**3/tau_ep_Tc_ps)
        self.kappa=4*(k*Tc_K*self.taupsi_ps*1e-12/hbar)**2
        self.denominator=graph.area_weights*2*np.sqrt((1+ratio)/2)*(np.pi/4)
        self.stretch=np.sqrt(1+self.kappa*abs(self.d)**2)
        tail,head=graph.edges.T;c=graph.conductance
        self.laplacian=coo_matrix((np.r_[c,c,-c,-c],
            (np.r_[tail,head,tail,head],np.r_[tail,head,head,tail])),
            shape=(graph.n_nodes,graph.n_nodes)).tocsr()
        self.potential_factor=splu(self.laplacian[self.free][:,self.free].tocsc())

    def local_inverse(self, force):
        f=_vector(force,self.graph.n_nodes,'force')
        projection=np.real(np.conj(self.d)*f)
        out=(self.stretch*f-self.kappa/self.stretch*self.d*projection)/self.denominator
        out[~self.free]=0
        return out

    def local_inverse_derivative(self, direction, force):
        v=_vector(direction,self.graph.n_nodes,'direction');f=_vector(force,self.graph.n_nodes,'force')
        if np.any(v[~self.free]!=0):raise ValueError('Fixed gap contacts cannot move')
        dot=np.real(np.conj(self.d)*v)
        dR=self.kappa*dot/self.stretch
        df_projection=np.real(np.conj(self.d)*f)
        vf_projection=np.real(np.conj(v)*f)
        out=(dR*f+self.kappa*dR/self.stretch**2*self.d*df_projection
             -self.kappa/self.stretch*(v*df_projection+self.d*vf_projection))/self.denominator
        out[~self.free]=0
        return out

    def potential(self,current):
        current=np.asarray(current,float)
        if current.shape!=(len(self.graph.edges),) or np.any(~np.isfinite(current)):
            raise ValueError('One finite real superconducting current per edge required')
        result=np.zeros(self.graph.n_nodes)
        rhs=-divergence(self.graph,current)
        result[self.free]=self.potential_factor.solve(rhs[self.free])
        return result

    def response(self,gradient,current):
        force=_vector(gradient,self.graph.n_nodes,'gradient');v=self.potential(current)
        material=-self.local_inverse(force)
        velocity=material-.5j*v*self.d
        tail,head=self.graph.edges.T
        normal=self.graph.conductance*(v[tail]-v[head])
        continuity=divergence(self.graph,np.asarray(current)+normal)
        kwt=float(-np.real(np.vdot(force,material)))
        joule=float(.5*np.dot(self.graph.conductance,(v[tail]-v[head])**2))
        rate=float(np.real(np.vdot(force,velocity)))
        noether=np.imag(np.conj(self.d)*force)+divergence(self.graph,current)
        return dict(velocity=velocity,material_velocity=material,potential_v=v,normal_current=normal,
            kwt_loss=kwt,normal_loss=joule,free_energy_rate=rate,
            dissipation_residual=rate+kwt+joule,
            continuity_max=float(np.max(abs(continuity[self.free]),initial=0.)),
            noether_max=float(np.max(abs(noether[self.free]),initial=0.)))

    def rhs_tangent(self,direction,hessian_direction,gradient,current,current_direction):
        direction=_vector(direction,self.graph.n_nodes,'direction')
        hd=_vector(hessian_direction,self.graph.n_nodes,'hessian_direction')
        if np.any(direction[~self.free]!=0):raise ValueError('Fixed gap contacts cannot move')
        v=self.potential(current);dv=self.potential(current_direction)
        fixed_metric=-self.local_inverse(hd)-.5j*dv*self.d
        baseline_correction=-self.local_inverse_derivative(direction,gradient)-.5j*v*direction
        result=fixed_metric+baseline_correction
        result[~self.free]=0
        return dict(velocity_direction=result,potential_direction=dv,
            baseline_mobility_gauge_correction=baseline_correction,
            continuity_direction_max=float(np.max(abs((divergence(self.graph,current_direction)+
                self.laplacian@dv)[self.free]),initial=0.)))
