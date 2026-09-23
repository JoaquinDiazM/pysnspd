"""Open 1D common-energy functional with compatible nodal spectral elements.

Degree-four Legendre--Gauss--Lobatto elements share their endpoint nodes.
M is the assembled positive quadrature, MD the assembled weak derivative,
and K the assembled positive stiffness. Jensen's inequality gives
K-D^*MD >= 0, so the regularized gradient remainder is nonnegative without
giving an alternating mode zero stiffness. Each global node owns one spectrum
and population. D is the mass-weighted derivative at a shared endpoint.

This removes the periodic diagnostic closure. Endpoint amplitudes and phase
fluxes are supplied by the caller; no Dirichlet phase is silently imposed.
No 2D interface or time evolution is certified by this longitudinal operator.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from numpy.polynomial.legendre import Legendre

from .energy_catalog import E_CHARGE_C, HBAR_J_S
from .spatial_functional import PeriodicSpatialFunctional, SpatialAdmissibilityError


def lobatto_operators(elements, length_bar, degree=4):
    if type(elements) is not int or elements < 1 or type(degree) is not int or degree < 2:
        raise ValueError('positive element count and degree at least two required')
    if not np.isfinite(length_bar) or length_bar <= 0:
        raise ValueError('length must be positive')
    polynomial=Legendre.basis(degree)
    roots=np.real_if_close(polynomial.deriv().roots())
    if np.iscomplexobj(roots):raise ValueError('Lobatto roots are not numerically real')
    reference=np.r_[-1.,np.sort(roots.astype(float)),1.]
    weights=2/(degree*(degree+1)*polynomial(reference)**2)
    barycentric=np.array([1/np.prod(reference[i]-np.delete(reference,i)) for i in range(degree+1)])
    derivative=np.zeros((degree+1,degree+1))
    for i in range(degree+1):
        for j in range(degree+1):
            if i!=j:derivative[i,j]=barycentric[j]/barycentric[i]/(reference[i]-reference[j])
        derivative[i,i]=-np.sum(derivative[i])
    size=elements*degree+1
    x=np.empty(size);mass=np.zeros(size);weak=np.zeros((size,size));stiffness=np.zeros_like(weak)
    width=length_bar/elements
    local_mass=weights*width/2;local_d=derivative*2/width
    for element in range(elements):
        indices=element*degree+np.arange(degree+1)
        x[indices]=element*width+(reference+1)*width/2
        mass[indices]+=local_mass
        weak[np.ix_(indices,indices)]+=local_mass[:,None]*local_d
        stiffness[np.ix_(indices,indices)]+=local_d.T@(local_mass[:,None]*local_d)
    return x,mass,weak/mass[:,None],stiffness


@dataclass(frozen=True)
class OpenSpatialEvaluation:
    energy_bar: float
    free_energy_bar: float | None
    energy_J: float
    gradient_cartesian_bar: np.ndarray
    link_current_bar: np.ndarray
    current_A: np.ndarray
    amplitude_nodes_bar: np.ndarray
    gamma_nodes_bar: np.ndarray
    q_delta_nodes_bar: np.ndarray
    electronic_derivatives: np.ndarray
    p_nodes: np.ndarray
    principal_symbols: tuple | None
    gradient_remainder_bar: float
    global_phase_residual_bar: float


class OpenSpatialFunctional(PeriodicSpatialFunctional):
    """Open longitudinal energy; populations are fixed in field variations."""

    def __init__(self,catalog,length_m,cross_section_m2,elements,*,degree=4,Tc_K=None):
        super().__init__(catalog,length_m,cross_section_m2,elements*degree+1,Tc_K=Tc_K)
        self.elements=elements;self.degree=degree
        self.x_bar,self.mass_bar,self.derivative_matrix,self.stiffness_matrix=lobatto_operators(
            elements,self.length_bar,degree)

    def _operators(self,delta_bar,link_phases):
        z=np.asarray(delta_bar,dtype=complex)
        if z.shape!=(self.cells,) or np.any(~np.isfinite(z)):
            raise ValueError('one finite complex field per node required')
        phases=np.zeros(self.cells-1) if link_phases is None else np.asarray(link_phases)
        if np.iscomplexobj(phases) or phases.shape!=(self.cells-1,) or np.any(~np.isfinite(phases)):
            raise ValueError('one finite real phase per open-chain link required')
        transport=np.exp(1j*np.r_[0.,np.cumsum(phases)])
        conjugation=np.conj(transport)[:,None]*transport[None,:]
        return z,self.derivative_matrix*conjugation,self.stiffness_matrix*conjugation

    def covariant_derivative(self,delta_bar,*,link_phases=None):
        z,d,_=self._operators(delta_bar,link_phases)
        return d@z

    def evaluate(self,delta_bar,p_nodes,*,link_phases=None,require_stability=True,on_node=None):
        return self._evaluate(delta_bar,p_nodes,None,link_phases,require_stability,on_node)

    def evaluate_thermal(self,delta_bar,bath_theta,*,link_phases=None,
                         require_stability=True,on_node=None):
        """Thermal free energy and its field gradient, including entropy.

        The thermal variational identity equals the internal-energy field
        derivative at fixed instantaneous FD occupations. Differentiating the
        thermal internal energy alone would incorrectly include T*dS.
        """
        if not np.isfinite(bath_theta) or bath_theta < 0:
            raise ValueError('bath_theta must be finite and nonnegative')
        return self._evaluate(delta_bar,None,float(bath_theta),link_phases,require_stability,on_node)

    def _evaluate(self,delta_bar,p_nodes,theta,link_phases,require_stability,on_node):
        z,dmat,kmat=self._operators(delta_bar,link_phases)
        d=dmat@z;mass=self.mass_bar;cache={};rows=[];free=[]
        if theta is None:
            raw=np.asarray(p_nodes)
            if np.iscomplexobj(raw):raise ValueError('populations must be real')
            p=np.asarray(raw,dtype=float)
            if p.shape!=(self.cells,len(self.catalog.count_nodes)) or np.any(~np.isfinite(p)) or np.any((p<0)|(p>1)):
                raise ValueError('physical nodal populations required')
        else:p=np.zeros((self.cells,len(self.catalog.count_nodes)))
        symbols=[] if require_stability else None
        for i in range(self.cells):
            if theta is not None and theta > 0:
                a=float(abs(z[i]));q=float(np.imag(np.conj(z[i])*d[i])/(a*a+.01));g=q*q/self.gap_ratio
                key=(a,g)
                if key not in cache:
                    cache[key]=[np.asarray(self.catalog.vacuum.evaluate(a,g)),
                                tuple(np.asarray(k) for k in self.catalog.energy_kernel(a,g))]
                vacuum,kernels=cache[key]
                exponent=np.exp(-kernels[0]/theta)
                p[i]=exponent/(1+exponent)
                free.append(float(vacuum[0]-4*theta*np.dot(self.catalog.count_weights,
                                                          np.logaddexp(0.,-kernels[0]/theta))))
            row=self._local(z[i],d[i],p[i],cache);rows.append(row)
            if theta==0:free.append(row.electronic_derivatives[0])
            if require_stability:
                symbol=self._symbol(z[i],d[i],p[i],cache);symbols.append(symbol)
                if not symbol.stable:raise SpatialAdmissibilityError(i,symbol)
            if on_node is not None:on_node(i,self.cells)
        electronic=np.asarray([r.electronic_derivatives for r in rows])
        q=np.asarray([r.q_delta_bar for r in rows]);gamma=np.asarray([r.gamma_bar for r in rows])
        rho=abs(z)**2
        v=(2*q*electronic[:,2]/self.gap_ratio-2*self.kappa*rho*q)/(rho+.01)
        gz=np.asarray([r.gradient_delta_bar for r in rows]);gd=1j*v*z
        gradient=mass*gz+dmat.conj().T@(mass*gd)+2*self.kappa*(kmat@z)
        # Every gauge link on a tree is generated by a nodal phase variation.
        # Noether gives G_theta_i + J_i - J_(i-1)=0, with no exterior link.
        phase_gradient=np.imag(np.conj(z)*gradient)
        current=-np.cumsum(phase_gradient)[:-1]
        remainder=float(self.kappa*(np.real(np.vdot(z,kmat@z))-np.dot(mass,rho*q*q)))
        energy=float(np.dot(mass,electronic[:,0])+remainder)
        free_energy=None if theta is None else float(np.dot(mass,free)+remainder)
        scale=(2*E_CHARGE_C/HBAR_J_S)*self.energy_scale_J
        return OpenSpatialEvaluation(energy,free_energy,self.energy_scale_J*energy,
            np.column_stack((gradient.real,gradient.imag)),current,scale*current,abs(z),gamma,q,
            electronic,p,tuple(symbols) if symbols is not None else None,remainder,
            float(np.sum(phase_gradient)))
