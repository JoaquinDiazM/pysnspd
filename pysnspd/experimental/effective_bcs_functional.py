"""Effective neutral BCS energy and dirty-limit gradient rigidity.

This is an explicit mesoscopic closure, not the rejected first-Moyal backend.
Quasiparticle labels xi are fixed native states and E=sqrt(xi**2+|d|**2).
The finite pairing counterterm uses exactly the same quadrature and imposes Tc
without a singular vacuum logarithm. Pair breaking enters the gauge gradient
only; the local DOS has Gamma=0.

The optional radial rigidity extends Vodolazov (2017), Eq. (33), to ONE energy:
K(|Delta|;Tb)=pi*N0*hbar*D*tanh(|Delta|/(2*kB*Tb))/(4*|Delta|).
Tb is fixed and the full derivative of K is retained. This extension is our
effective hypothesis, not a derivation of nonequilibrium Usadel. It retains
finite positive Cartesian gradient stiffness at Delta=0 without a core cutoff.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from numpy.polynomial.legendre import leggauss
from numpy.polynomial.polynomial import polyval, polyder
from scipy.constants import Boltzmann as kB, hbar, elementary_charge as e
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expit
from .thermal_spatial_usadel import ThermalGraph

VODOLAZOV_SOURCE = "https://arxiv.org/pdf/1611.06060"
BCS_SOURCE = "https://doi.org/10.1103/PhysRev.108.1175"


def bcs_count_quadrature(cutoff_bar=32., order=6):
    """Positive Gauss nodes/weights in xi/Eu, Eu=kB*Tc, with no state at xi=0."""
    cutoff = float(cutoff_bar)
    if not np.isfinite(cutoff) or cutoff <= 0 or type(order) is not int or order < 2:
        raise ValueError("positive cutoff and integer Gauss order >=2 required")
    knots = np.array([0., .125, .25, .5, 1., 2., 4., 8., 16., 32., 64., 128.])
    knots = np.r_[knots[knots < cutoff], cutoff]
    if knots[0] != 0.:
        knots = np.r_[0., knots]
    x, w = leggauss(order)
    lo, hi = knots[:-1, None], knots[1:, None]
    nodes = ((lo+hi)/2+(hi-lo)*x/2).ravel()
    weights = ((hi-lo)*w/2).ravel()
    return nodes, weights


def _positive(value, name):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError(name+" must be finite and positive")
    return value


def _axis(values, name):
    raw = np.asarray(values)
    if np.iscomplexobj(raw):
        raise ValueError(name+" must be real")
    arr = np.asarray(raw, float).copy()
    if (arr.ndim != 1 or len(arr) < 2 or np.any(~np.isfinite(arr))
            or np.any(arr <= 0)):
        raise ValueError(name+" must be a positive finite vector")
    arr.setflags(write=False)
    return arr


@dataclass(frozen=True)
class EffectiveBCSEvaluation:
    energy_J: float
    local_energy_J: float
    gradient_energy_J: float
    node_energy_J: np.ndarray
    electron_energies_J: np.ndarray
    gap_force_J_per_gapbar: np.ndarray
    gap_gradient_J_per_gapbar: np.ndarray
    edge_current_A: np.ndarray
    current_divergence_A: np.ndarray
    link_derivative_alpha_J: np.ndarray
    local_coefficient_bar: np.ndarray
    gradient_kappa: np.ndarray


@dataclass(frozen=True)
class EffectiveBCSDirection:
    energy_direction_J: float
    gap_force_direction_J_per_gapbar: np.ndarray
    edge_current_direction_A: np.ndarray
    electron_energy_direction_J: np.ndarray


@dataclass(frozen=True)
class EffectiveBCSBulk:
    amplitude_bar: float
    delta_J: float
    q_bar: float
    q_m_inv: float
    current_A: float
    depairing_current_A: float
    depairing_q_bar: float
    gap_derivative_J_per_A: float
    q_derivative_m_inv_per_A: float
    stationarity_residual_bar: float
    longitudinal_stable: bool
    temperature_K: float
    contract: str
    phase_step_m: float | None = None


class EffectiveBCSFunctional:
    """One SI energy, its Wirtinger gap force and oriented physical current.

    gapbar=Delta/Eu, xi and weights are dimensionless in Eu=kB*Tc.
    p has shape (node, xi). dU=2*Re(vdot(gap_force,dd))+
    dot(link_derivative_alpha,da)+sum(4*N0*Eu**2*V*w*E*dp).
    Graph conductance is dual-length/primal-length, without material factors.
    """
    def __init__(self, graph, *, energy_scale_J, N0_per_J_m3, ell0_m,
                 thickness_m, count_nodes, count_weights,
                 gradient_stiffness_SI=None, gradient_law="constant",
                 bath_temperature_K=None, diffusion_m2_s=None):
        if not isinstance(graph, ThermalGraph):
            raise ValueError("ThermalGraph geometry required")
        self.graph = graph
        self.Eu_J = _positive(energy_scale_J, "energy_scale_J")
        self.N0_per_J_m3 = _positive(N0_per_J_m3, "N0_per_J_m3")
        self.ell0_m = _positive(ell0_m, "ell0_m")
        self.thickness_m = _positive(thickness_m, "thickness_m")
        self.Tc_K = self.Eu_J/kB
        self.xi_bar = _axis(count_nodes, "count_nodes")
        self.weights_bar = _axis(count_weights, "count_weights")
        if self.xi_bar.shape != self.weights_bar.shape or np.any(np.diff(self.xi_bar) <= 0):
            raise ValueError("ordered native nodes and matching positive weights required")
        self.count_nodes, self.count_weights = self.xi_bar, self.weights_bar
        self.inverse_pairing_coupling = float(np.dot(
            self.weights_bar, np.tanh(self.xi_bar/2)/self.xi_bar))
        self.volumes_m3 = graph.area_weights*self.ell0_m**2*self.thickness_m
        self.node_energy_scale_J = self.N0_per_J_m3*self.Eu_J**2*self.volumes_m3
        self.edge_energy_scale_J = (self.N0_per_J_m3*self.Eu_J**2
                                   *self.ell0_m**2*self.thickness_m*graph.conductance)
        if gradient_law not in ("constant", "bcs_radial"):
            raise ValueError("gradient_law must be constant or bcs_radial")
        self.gradient_law = gradient_law
        self.bath_temperature_K = None
        if gradient_law == "constant":
            stiffness = _positive(gradient_stiffness_SI, "gradient_stiffness_SI")
            self.constant_kappa = stiffness/(self.N0_per_J_m3*self.ell0_m**2)
            self.radial_prefactor = None
        else:
            self.bath_temperature_K = _positive(bath_temperature_K, "bath_temperature_K")
            diffusion = _positive(diffusion_m2_s, "diffusion_m2_s")
            if gradient_stiffness_SI is not None:
                raise ValueError("bcs_radial stiffness follows D; no extra fitted stiffness allowed")
            self.radial_prefactor = np.pi*hbar*diffusion/(4*self.Eu_J*self.ell0_m**2)
            self.constant_kappa = None

    def rigidity(self, rho):
        """Return dimensionless kappa and its first two rho derivatives."""
        raw = np.asarray(rho)
        if np.iscomplexobj(raw):
            raise ValueError("rho must be real")
        rho = np.asarray(raw, float)
        if np.any(~np.isfinite(rho)) or np.any(rho < 0):
            raise ValueError("rho must be finite and nonnegative")
        if self.gradient_law == "constant":
            return (np.full_like(rho, self.constant_kappa),
                    np.zeros_like(rho), np.zeros_like(rho))
        b = self.Eu_J/(2*kB*self.bath_temperature_K)
        z = b*np.sqrt(rho)
        result = [np.empty_like(rho) for _ in range(3)]
        small = z < .1
        # tanh(z)/z, differentiated as a polynomial in z**2 near the origin.
        coeff = np.array([1., -1/3, 2/15, -17/315, 62/2835,
                          -1382/155925, 21844/6081075, -929569/638512875])
        for order in range(3):
            result[order][small] = b**(2*order+1)*polyval(
                z[small]**2, polyder(coeff, order))
        regular = ~small
        a = np.sqrt(rho[regular])
        t = np.tanh(z[regular])
        s = 1-t*t
        result[0][regular] = t/a
        result[1][regular] = (b*a*s-t)/(2*a**3)
        result[2][regular] = (-2*b*b*a*a*s*t-3*b*a*s+3*t)/(4*a**5)
        return tuple(self.radial_prefactor*r for r in result)

    def _state(self, gapbar, p, alpha):
        d = np.asarray(gapbar, complex)
        rawp = np.asarray(p)
        if (d.shape != (self.graph.n_nodes,) or np.any(~np.isfinite(d))
                or np.iscomplexobj(rawp)):
            raise ValueError("finite node gap and real native occupations required")
        p = np.asarray(rawp, float)
        if (p.shape != (self.graph.n_nodes, len(self.xi_bar))
                or np.any(~np.isfinite(p)) or np.any(p < 0) or np.any(p > 1)):
            raise ValueError("native occupations must have shape (node,xi) and lie in [0,1]")
        alpha = np.zeros(len(self.graph.edges)) if alpha is None else np.asarray(alpha)
        if (np.iscomplexobj(alpha) or alpha.shape != (len(self.graph.edges),)
                or np.any(~np.isfinite(alpha))):
            raise ValueError("one finite real alpha per edge required")
        return d, p, np.asarray(alpha, float)

    def energies_bar(self, gapbar):
        d = np.asarray(gapbar, complex)
        if np.any(~np.isfinite(d)):
            raise ValueError("finite gap required")
        return np.sqrt(abs(d)[..., None]**2+self.xi_bar**2)

    def fermi_dirac(self, gapbar, temperature_K):
        temp = np.asarray(temperature_K, float)
        if np.any(~np.isfinite(temp)) or np.any(temp <= 0):
            raise ValueError("positive temperature required")
        energies = self.energies_bar(gapbar)
        if temp.ndim:
            temp = temp[..., None]
        return expit(-self.Eu_J*energies/(kB*temp))

    def evaluate(self, gapbar, p, alpha=None):
        d, p, alpha = self._state(gapbar, p, alpha)
        rho = abs(d)**2
        energy = self.energies_bar(d)
        # Stable E-xi avoids cancellation at small amplitude.
        energy_increment = rho[:, None]/(energy+self.xi_bar)
        local = (rho*self.inverse_pairing_coupling
                 -2*np.sum(self.weights_bar*energy_increment, axis=1)
                 +4*np.sum(self.weights_bar*energy*p, axis=1))
        coefficient = (self.inverse_pairing_coupling
                       -np.sum(self.weights_bar*(1-2*p)/energy, axis=1))
        force = self.node_energy_scale_J*coefficient*d
        node_energy = self.node_energy_scale_J*local
        tail, head = self.graph.edges.T
        phase = np.exp(-1j*alpha)
        difference = d[tail]-phase*d[head]
        difference2 = abs(difference)**2
        kap, kp, _ = self.rigidity(rho)
        ke = (kap[tail]+kap[head])/2
        c = self.edge_energy_scale_J
        edge_energy = c*ke*difference2
        np.add.at(force, tail, c*(ke*difference+.5*kp[tail]*d[tail]*difference2))
        np.add.at(force, head, c*(-ke*np.conj(phase)*difference+.5*kp[head]*d[head]*difference2))
        np.add.at(node_energy, tail, edge_energy/2)
        np.add.at(node_energy, head, edge_energy/2)
        derivative = -2*c*ke*np.imag(np.conj(d[tail])*phase*d[head])
        current = -2*e/hbar*derivative
        divergence = np.zeros(self.graph.n_nodes)
        np.add.at(divergence, tail, current)
        np.add.at(divergence, head, -current)
        return EffectiveBCSEvaluation(
            float(node_energy.sum()), float(np.dot(self.node_energy_scale_J, local)),
            float(edge_energy.sum()), node_energy, self.Eu_J*energy, force, 2*force,
            current, divergence, derivative, coefficient, kap)

    def direction(self, gapbar, p, gap_direction, p_direction=None, alpha=None,
                  alpha_direction=None):
        """Analytic directional derivative, including the complete radial K law."""
        d, p, alpha = self._state(gapbar, p, alpha)
        dd = np.asarray(gap_direction, complex)
        dp = np.zeros_like(p) if p_direction is None else np.asarray(p_direction)
        da = np.zeros_like(alpha) if alpha_direction is None else np.asarray(alpha_direction)
        if (dd.shape != d.shape or np.any(~np.isfinite(dd)) or np.iscomplexobj(dp)
                or dp.shape != p.shape or np.any(~np.isfinite(dp))
                or np.iscomplexobj(da) or da.shape != alpha.shape or np.any(~np.isfinite(da))):
            raise ValueError("finite matching directional fields required")
        dp, da = np.asarray(dp, float), np.asarray(da, float)
        value = self.evaluate(d, p, alpha)
        rho = abs(d)**2
        drho = 2*np.real(np.conj(d)*dd)
        energy = value.electron_energies_J/self.Eu_J
        de = drho[:, None]/(2*energy)
        dc = (np.sum(self.weights_bar*(1-2*p)*de/energy**2, axis=1)
              +2*np.sum(self.weights_bar*dp/energy, axis=1))
        df = self.node_energy_scale_J*(value.local_coefficient_bar*dd+dc*d)
        kap, kp, kpp = self.rigidity(rho)
        dk, dkp = kp*drho, kpp*drho
        tail, head = self.graph.edges.T
        phase = np.exp(-1j*alpha)
        diff = d[tail]-phase*d[head]
        d_diff = dd[tail]-phase*dd[head]+1j*da*phase*d[head]
        diff2 = abs(diff)**2
        d_diff2 = 2*np.real(np.conj(diff)*d_diff)
        ke, dke = (kap[tail]+kap[head])/2, (dk[tail]+dk[head])/2
        c = self.edge_energy_scale_J
        common = dke*diff+ke*d_diff
        np.add.at(df, tail, c*(common+.5*(dkp[tail]*d[tail]*diff2
                  +kp[tail]*dd[tail]*diff2+kp[tail]*d[tail]*d_diff2)))
        np.add.at(df, head, c*(-np.conj(phase)*(common+1j*da*ke*diff)
                  +.5*(dkp[head]*d[head]*diff2+kp[head]*dd[head]*diff2
                       +kp[head]*d[head]*d_diff2)))
        product = np.conj(d[tail])*phase*d[head]
        dproduct = np.conj(dd[tail])*phase*d[head]+np.conj(d[tail])*phase*dd[head]-1j*da*product
        di = 4*e/hbar*c*(dke*np.imag(product)+ke*np.imag(dproduct))
        du = (2*np.real(np.vdot(value.gap_force_J_per_gapbar, dd))
              +np.dot(value.link_derivative_alpha_J, da)
              +4*np.sum(self.node_energy_scale_J[:, None]*self.weights_bar*energy*dp))
        return EffectiveBCSDirection(float(du), df, di, self.Eu_J*de)


    def hessian_and_current_jacobian(self, gapbar, p, alpha=None):
        """Sparse real Jacobians of the Wirtinger force and edge current.

        Coordinates are interleaved (Re d0,Im d0,Re d1,Im d1,...).
        H is dF/dd, i.e. one half of the real Hessian of U, symmetric at fixed
        p and alpha. No finite differences, dense assembly or row removal.
        """
        from scipy.sparse import coo_matrix
        d,p,alpha=self._state(gapbar,p,alpha)
        n=self.graph.n_nodes
        z=np.column_stack((d.real,d.imag))
        rho=abs(d)**2
        en=self.energies_bar(d)
        local_c=self.inverse_pairing_coupling-np.sum(self.weights_bar*(1-2*p)/en,axis=1)
        local_s=np.sum(self.weights_bar*(1-2*p)/en**3,axis=1)
        identity=np.eye(2)
        local=(self.node_energy_scale_J[:,None,None]*
               (local_c[:,None,None]*identity+
                local_s[:,None,None]*np.einsum("ei,ej->eij",z,z)))
        t,h=self.graph.edges.T
        phase=np.exp(-1j*alpha)
        rot=np.stack((np.stack((phase.real,-phase.imag),axis=1),
                      np.stack((phase.imag,phase.real),axis=1)),axis=1)
        zi,zj=z[t],z[h]
        rzj=np.einsum("eij,ej->ei",rot,zj)
        w=zi-rzj
        rtw=np.einsum("eji,ej->ei",rot,w)
        w2=np.sum(w*w,axis=1)
        k,kp,kpp=self.rigidity(rho)
        ke=(k[t]+k[h])/2
        c=self.edge_energy_scale_J
        outer=lambda a,b:np.einsum("ei,ej->eij",a,b)
        hii=(ke+.5*kp[t]*w2)[:,None,None]*identity
        hii+=kp[t,None,None]*(outer(w,zi)+outer(zi,w))
        hii+=kpp[t,None,None]*w2[:,None,None]*outer(zi,zi)
        hjj=(ke+.5*kp[h]*w2)[:,None,None]*identity
        hjj-=kp[h,None,None]*(outer(rtw,zj)+outer(zj,rtw))
        hjj+=kpp[h,None,None]*w2[:,None,None]*outer(zj,zj)
        hij=-ke[:,None,None]*rot+kp[h,None,None]*outer(w,zj)-kp[t,None,None]*outer(zi,rtw)
        blocks=[local,c[:,None,None]*hii,c[:,None,None]*hij,
                c[:,None,None]*np.transpose(hij,(0,2,1)),c[:,None,None]*hjj]
        row_nodes=[np.arange(n),t,t,h,h]
        col_nodes=[np.arange(n),t,h,t,h]
        rows,cols,vals=[],[],[]
        for rr,cc,block in zip(row_nodes,col_nodes,blocks):
            rindex=2*rr[:,None]+np.arange(2)
            cindex=2*cc[:,None]+np.arange(2)
            rows.append(np.broadcast_to(rindex[:,:,None],block.shape).ravel())
            cols.append(np.broadcast_to(cindex[:,None,:],block.shape).ravel())
            vals.append(block.ravel())
        hessian=coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),
                          shape=(2*n,2*n)).tocsr()
        # Im(conj(d_i)*U*d_j)=z_i^T J R z_j, J=[[0,1],[-1,0]].
        imag_product=zi[:,0]*rzj[:,1]-zi[:,1]*rzj[:,0]
        j_rzj=np.column_stack((rzj[:,1],-rzj[:,0]))
        j_zi=np.column_stack((zi[:,1],-zi[:,0]))
        rt_jzi=np.einsum("eji,ej->ei",rot,j_zi)
        scale=4*e/hbar*c
        it=scale[:,None]*(kp[t,None]*zi*imag_product[:,None]+ke[:,None]*j_rzj)
        ih=scale[:,None]*(kp[h,None]*zj*imag_product[:,None]-ke[:,None]*rt_jzi)
        edge_rows=np.repeat(np.arange(len(t)),4)
        edge_cols=np.column_stack((2*t,2*t+1,2*h,2*h+1)).ravel()
        current_jac=coo_matrix((np.column_stack((it,ih)).ravel(),(edge_rows,edge_cols)),
                              shape=(len(t),2*n)).tocsr()
        return hessian,current_jac

    def _bulk_q_terms(self, q_bar, phase_step_m):
        if phase_step_m is None:
            return q_bar*q_bar, q_bar, 1.
        step = _positive(phase_step_m, "phase_step_m")/self.ell0_m
        angle = q_bar*step
        return (4*np.sin(angle/2)**2/step**2,
                np.sin(angle)/step, np.cos(angle))

    def bulk_coefficient(self, amplitude_bar, temperature_K, q_bar=0.,
                         phase_step_m=None):
        """Uniform FD stationarity divided by Delta, at the same finite quadrature."""
        a = float(amplitude_bar)
        t = _positive(temperature_K, "temperature_K")/self.Tc_K
        q = float(q_bar)
        if not np.isfinite(a) or a < 0 or not np.isfinite(q):
            raise ValueError("finite nonnegative amplitude and finite q required")
        en = np.sqrt(self.xi_bar**2+a*a)
        coefficient = self.inverse_pairing_coupling-np.dot(
            self.weights_bar, np.tanh(en/(2*t))/en)
        kap, kp, _ = self.rigidity(np.asarray(a*a))
        q_squared, _, _ = self._bulk_q_terms(q, phase_step_m)
        return float(coefficient+q_squared*(kap+a*a*kp))

    def bulk_amplitude(self, temperature_K, q_bar=0., phase_step_m=None):
        """Largest locally stable positive FD amplitude, or zero if no branch."""
        if self.bulk_coefficient(0., temperature_K, q_bar, phase_step_m) == 0:
            return 0.
        upper = 2.
        while self.bulk_coefficient(upper, temperature_K, q_bar, phase_step_m) <= 0 and upper < 16:
            upper *= 2
        axis = np.r_[0., np.geomspace(1e-7, upper, 220)]
        values = np.array([self.bulk_coefficient(a, temperature_K, q_bar, phase_step_m) for a in axis])
        crossings = np.flatnonzero((values[:-1] < 0) & (values[1:] >= 0))
        if not len(crossings):
            return 0.
        k = crossings[-1]
        return float(brentq(lambda a:self.bulk_coefficient(a, temperature_K, q_bar, phase_step_m),
                            axis[k], axis[k+1], xtol=5e-14, rtol=5e-14))

    def bulk_current(self, temperature_K, q_bar, width_m, phase_step_m=None):
        a = self.bulk_amplitude(temperature_K, q_bar, phase_step_m)
        kap = self.rigidity(np.asarray(a*a))[0]
        prefactor = (4*e/hbar*self.N0_per_J_m3*self.Eu_J**2*self.ell0_m
                     *_positive(width_m, "width_m")*self.thickness_m)
        _, current_q, _ = self._bulk_q_terms(q_bar, phase_step_m)
        return float(prefactor*kap*a*a*current_q), a

    def bulk(self, *, temperature_K, current_A, width_m, phase_step_m=None):
        """Fresh stable homogeneous reference, with its own depairing admission."""
        temp = _positive(temperature_K, "temperature_K")
        current = float(current_A)
        if not np.isfinite(current) or current < 0:
            raise ValueError("nonnegative current required")
        if temp >= self.Tc_K:
            raise ValueError("superconducting reference requires temperature below Tc")
        width = _positive(width_m, "width_m")
        # A declared scan brackets the physical low-q maximum; no inherited DC admission.
        q_limit = 3. if phase_step_m is None else min(
            3., np.pi*self.ell0_m/_positive(phase_step_m, "phase_step_m"))
        axis = np.linspace(0., q_limit, 151)
        currents = np.array([self.bulk_current(temp, q, width, phase_step_m)[0] for q in axis])
        peak = int(np.argmax(currents))
        if peak in (0, len(axis)-1):
            raise RuntimeError("depairing maximum was not bracketed")
        opt = minimize_scalar(lambda q:-self.bulk_current(temp, q, width, phase_step_m)[0],
                              bounds=(axis[peak-1], axis[peak+1]), method="bounded",
                              options={"xatol":2e-10})
        qdep, idep = float(opt.x), float(-opt.fun)
        if current >= idep:
            raise ValueError("requested current exceeds this effective model depairing current")
        q = (0. if current == 0 else brentq(
             lambda q:self.bulk_current(temp,q,width,phase_step_m)[0]-current, 0.,qdep,
             xtol=2e-13, rtol=2e-13))
        measured, a = self.bulk_current(temp,q,width,phase_step_m)
        en = np.sqrt(self.xi_bar**2+a*a)
        t = temp/self.Tc_K
        th = np.tanh(en/(2*t))
        local_da = float(np.dot(self.weights_bar,
            a*(th/en**3-(1-th*th)/(2*t*en**2))))
        kap, kp, kpp = [float(v) for v in self.rigidity(np.asarray(a*a))]
        q_squared,current_q,d_current_q = self._bulk_q_terms(q,phase_step_m)
        ca = local_da+q_squared*(4*a*kp+2*a**3*kpp)
        cq = 2*current_q*(kap+a*a*kp)
        da_dq = -cq/ca
        prefactor = 4*e/hbar*self.N0_per_J_m3*self.Eu_J**2*self.ell0_m*width*self.thickness_m
        di_dq = prefactor*(kap*a*a*d_current_q+current_q*2*a*(kap+a*a*kp)*da_dq)
        return EffectiveBCSBulk(
            a,self.Eu_J*a,q,q/self.ell0_m,measured,idep,qdep,
            self.Eu_J*da_dq/di_dq,1/(self.ell0_m*di_dq),
            abs(self.bulk_coefficient(a,temp,q,phase_step_m)),bool(di_dq>0),temp,
            "effective neutral native BCS; current pair breaking only in shared gradient; fresh DC",phase_step_m)
