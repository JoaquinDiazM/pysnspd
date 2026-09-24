"""Weak longitudinal gap/population coupling at a real, fixed-phase reference.

This is a reciprocal mesoscopic linear prototype. Its quadratic availability
is not total internal energy. KWT dissipation is reported once; it is not also
injected into first-order populations. No collision bath, photon, charge mode,
electrostatic closure or circuit is introduced by this module.
"""
from dataclasses import dataclass
import numpy as np


def thermal_susceptibility(energy, temperature):
    energy = np.asarray(energy, float)
    if (energy.ndim != 1 or np.any(~np.isfinite(energy)) or np.any(energy < 0)
            or not np.isfinite(temperature) or temperature <= 0):
        raise ValueError('Finite nonnegative energies and positive temperature required')
    q = np.exp(-energy/temperature)
    chi = 2*q/(temperature*(1+q)**2)
    if np.any(chi == 0):
        raise ValueError('Thermal susceptibility underflow; no artificial population capacity added')
    return chi


@dataclass(frozen=True)
class LongitudinalResponse:
    gap_velocity: np.ndarray
    population_mass_rhs: np.ndarray
    population_velocity: object
    effective_gap_force: np.ndarray
    availability_rate: float
    kwt_loss: float
    longitudinal_loss: float
    reciprocal_power_residual: float
    charge_residual_max: float


class LongitudinalReciprocalBridge:
    """All-node radial-in-gap response, with y(E)=delta hL(E)/chi(E).

    hessian_action must be the thermal Cartesian Schur Hessian and
    inverse_mobility the inherited radial KWT map on integrated forces.
    Both contacts and distribution increments are fixed to zero. The affine
    drift of a nonequilibrium reference is outside this homogeneous tangent;
    supplying this prototype never redefines an admitted nonzero G0 as zero.
    """
    def __init__(self, graph, operators, energies, weights, temperature,
                 hessian_action, inverse_mobility, *, spectral_z,
                 spectral_residual_tolerance=1e-7):
        self.graph = graph; self.operators = tuple(operators)
        self.energies = np.asarray(energies, float)
        self.weights = np.asarray(weights, float)
        if (not self.operators or self.weights.shape != self.energies.shape
                or len(self.operators) != len(self.energies)
                or np.any(~np.isfinite(self.weights)) or np.any(self.weights <= 0)):
            raise ValueError('One positive finite quadrature weight per spectral operator required')
        self.chi = thermal_susceptibility(self.energies, temperature)
        self.spectral_z = np.asarray(spectral_z, complex)
        if (self.spectral_z.shape != self.energies.shape or np.any(~np.isfinite(self.spectral_z))
                or np.any(self.spectral_z.real < 0)
                or np.any(abs(self.spectral_z.imag+self.energies) > 64*np.finfo(float).eps*np.maximum(1.,self.energies))
                or not np.isfinite(spectral_residual_tolerance) or spectral_residual_tolerance <= 0):
            raise ValueError('One causal z=eta-iE matching each quadrature energy is required')
        self.hessian_action = hessian_action; self.inverse_mobility = inverse_mobility
        self.free = np.ones(graph.n_nodes, bool); self.free[graph.boundary_nodes] = False
        self.reference_gap = self.operators[0].d.copy()
        self.mass = []; self.kernel = []; self.longitudinal = []; self.edge_rates = []
        for operator, z in zip(self.operators,self.spectral_z):
            if (operator.graph is not graph or np.any(operator.d.imag != 0)
                    or not np.array_equal(operator.d,self.reference_gap)):
                raise ValueError('One common graph and strictly real fixed-phase gap required')
            f, ft = operator.R[:, 0, 1], operator.R[:, 1, 0]
            g=operator.R[:,0,0];a=f/(1+g);b=ft/(1+g)
            tail,head=graph.edges.T;c=graph.conductance;m=graph.area_weights
            hz=m*z;hp=m*operator.d;hm=m*np.conj(operator.d)
            np.add.at(hz,tail,c*g[head]);np.add.at(hz,head,c*g[tail])
            np.add.at(hp,tail,c*f[head]);np.add.at(hp,head,c*f[tail])
            np.add.at(hm,tail,c*ft[head]);np.add.at(hm,head,c*ft[tail])
            residual=np.column_stack((hm*a*a+2*hz*a-hp,hp*b*b+2*hz*b-hm))
            if np.max(abs(residual[self.free])/(m*np.maximum(1.,abs(operator.d)))[self.free,None],initial=0.) > spectral_residual_tolerance:
                raise ValueError('Spectrum is not stationary for its declared energy/reference with alpha=0')
            magnitude = max(1., float(np.max(abs(operator.matrix.data), initial=0.)))
            tolerance = 256*np.finfo(float).eps*magnitude
            if np.max(abs(f-ft), initial=0.) > tolerance:
                raise ValueError('The spectral reference is not in the real fixed-phase sector')
            lt = operator.matrix[::2, 1::2]; tl = operator.matrix[1::2, ::2]
            if max(np.max(abs(lt.data), initial=0.), np.max(abs(tl.data), initial=0.)) > tolerance:
                raise ValueError('Charge and longitudinal modes do not decouple at this reference')
            rate = operator.edge_head_blocks[:, 0, 0]
            if (np.any(rate < 0) or np.max(abs(operator.edge_tail_blocks[:, 0, 0]+rate), initial=0.) > tolerance
                    or np.max(abs(operator.pair_blocks[:, 0, 0]), initial=0.) > tolerance):
                raise ValueError('Longitudinal operator lacks an admitted nonnegative Dirichlet form; no clipping')
            probe = np.column_stack((np.ones(graph.n_nodes), np.zeros(graph.n_nodes)))
            force = operator.evaluate(probe).gap_force_increment
            if np.max(abs(force.imag), initial=0.) > tolerance:
                raise ValueError('Longitudinal forcing is not purely amplitude-like')
            rho = operator.R[:, 0, 0].real
            if np.any(rho < 0):
                raise ValueError('Negative DOS cannot define population storage; no floor')
            self.mass.append(graph.area_weights*rho)
            self.kernel.append(force.real)
            self.longitudinal.append(operator.matrix[::2, ::2].tocsr())
            self.edge_rates.append(rate.copy())
        self.mass = np.asarray(self.mass); self.kernel = np.asarray(self.kernel)
        self.edge_rates = np.asarray(self.edge_rates)

    def _callback_vector(self, callback, argument, name):
        raw=np.asarray(callback(argument))
        if (np.iscomplexobj(raw) or raw.shape != (self.graph.n_nodes,)
                or np.any(~np.isfinite(raw))):
            raise ValueError(name+' must return one finite real value per node; no complex part is discarded')
        return raw.astype(float)

    def _state(self, x, y):
        xraw, yraw = np.asarray(x), np.asarray(y)
        if (np.iscomplexobj(xraw) or np.iscomplexobj(yraw)
                or xraw.shape != (self.graph.n_nodes,) or yraw.shape != self.mass.shape
                or np.any(~np.isfinite(xraw)) or np.any(~np.isfinite(yraw))):
            raise ValueError('Finite real full-node amplitude and spectral population directions required')
        x, y = xraw.astype(float), yraw.astype(float)
        if np.any(x[~self.free] != 0) or np.any(y[:, ~self.free] != 0):
            raise ValueError('Fixed contact increments must remain zero')
        return x, y

    def availability(self, x, y):
        x, y = self._state(x, y)
        Hx = self._callback_vector(self.hessian_action,x,'Thermal Hessian')
        return float(.5*np.dot(x, Hx)+np.sum(self.weights[:, None]*self.chi[:, None]*self.mass*y*y))

    def response(self, x, y, *, divide_mass=True):
        x, y = self._state(x, y)
        Hx = self._callback_vector(self.hessian_action,x,'Thermal Hessian')
        force = Hx+np.sum((self.weights*self.chi)[:, None]*self.kernel*y, axis=0)
        velocity = -self._callback_vector(self.inverse_mobility,force,'Inverse mobility')
        if velocity.shape != x.shape or np.any(~np.isfinite(velocity)) or np.any(velocity[~self.free] != 0):
            raise ValueError('The inherited mobility must return finite velocities and preserve fixed contacts')
        transport = np.asarray([operator@row for operator, row in zip(self.longitudinal, y)])
        rhs = transport+.5*self.kernel*velocity[None, :]
        rhs[:, ~self.free] = 0.
        ydot = None
        if divide_mass:
            if np.any(self.mass[:, self.free] <= 0):
                raise ValueError('Zero DOS is an algebraic mass constraint; use the mass form, never a DOS floor')
            ydot = np.zeros_like(y); ydot[:, self.free] = rhs[:, self.free]/self.mass[:, self.free]
        kwt_loss = -float(np.dot(force, velocity))
        tail, head = self.graph.edges.T
        longitudinal_loss = float(2*np.sum((self.weights*self.chi)[:, None]*self.edge_rates*(y[:, head]-y[:, tail])**2))
        # This contraction uses the mass RHS, so zero DOS needs no division.
        rate = float(np.dot(Hx, velocity)+2*np.sum((self.weights*self.chi)[:, None]*y*rhs))
        charge = 0.
        for operator, chi, row in zip(self.operators, self.chi, y):
            residual = operator.evaluate(np.column_stack((chi*row, np.zeros_like(row)))).residual
            charge = max(charge, float(np.max(abs(residual[:, 1]), initial=0.)))
        return LongitudinalResponse(velocity, rhs, ydot, force, rate, kwt_loss,
            longitudinal_loss, rate+kwt_loss+longitudinal_loss, charge)
