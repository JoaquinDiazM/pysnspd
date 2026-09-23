"""Registered, non-photon static controls; no physical time integration.

Small-amplitude support is a NEW direct-spectrum reference, not extrapolation
of the archived occupation table. The same causal solver and implicit fixed-
count derivatives are reused. A cache uses exact field coordinates only.
"""
from dataclasses import dataclass
from functools import lru_cache
import time
import numpy as np
from scipy.optimize import brentq

from pysnspd.experimental.energy_catalog import (
    build_vacuum_catalog, E_CHARGE_C, HBAR_J_S, K_B_J_K)
from pysnspd.experimental.refined_cells import refined_count_catalog
from pysnspd.experimental.cell_closures import CellScales, KWTMobility
from pysnspd.experimental.rectangular_spatial import RectangularSpatialFunctional
from pysnspd.experimental.spatial_dynamics import kwt_spatial_response


@dataclass(frozen=True)
class DirectReference:
    vacuum: object
    eta: float


class CachedDirect:
    def __init__(self, source):
        self.source = source
        self.vacuum, self.eta = source.vacuum, source.eta
        self.count_nodes, self.count_weights = source.count_nodes, source.count_weights
        self.calls = 0
        self.seconds = 0.
        self.energy_kernel = lru_cache(maxsize=4096)(self._kernel)

    def _kernel(self, amplitude, gamma):
        start = time.monotonic()
        result = self.source.energy_kernel(amplitude, gamma)
        self.calls += 1
        self.seconds += time.monotonic()-start
        for value in result:
            value.setflags(write=False)
        return result

    def diagnostics(self):
        return dict(count_states=len(self.count_nodes), direct_queries=self.calls,
                    direct_seconds=self.seconds, exact_cache=self.energy_kernel.cache_info()._asdict())


def material_reference(plan):
    m = plan['material']
    sigma = 1/(m['R_sheet_ohm']*m['thickness_m'])
    n0 = sigma/(2*E_CHARGE_C**2*m['D_m2_s'])
    return dict(**m, conductivity_S_m=sigma, N0_per_J_m3=n0,
                ell0_m=float(np.sqrt(HBAR_J_S*m['D_m2_s']/(2*K_B_J_K*m['Tc_K']))),
                N0_convention='single-spin; sigma=2e^2 N0 D')


def catalogue(plan, resolution):
    m = material_reference(plan)
    # Only closed-form vacuum values use this tiny array. No excitation-energy
    # interpolant is built; direct retarded roots are checked on every miss.
    vac = build_vacuum_catalog(np.array([1e-4, .01, .1, 1.2]),
        np.array([0., .05, .5, 2.]), Tc_K=m['Tc_K'], N0_per_J_m3=m['N0_per_J_m3'],
        D_m2_s=m['D_m2_s'], analytic=True)
    vac.metadata['scope'] = 'stage4 direct diagnostic support; no physical core or trajectory admission'
    spec = plan['spectral_resolutions'][resolution]
    base = DirectReference(vac, spec['eta'])
    return CachedDirect(refined_count_catalog(base, refinement=spec['refinement']))


def scales(cat, plan):
    m = plan['material']
    return CellScales(cat.vacuum.delta0_J, cat.vacuum.N0_per_J_m3,
                      m['Tc_K'], m['Tb_K'], 1.)


def fd_population(energy, theta):
    if theta == 0:
        return np.zeros_like(energy)
    e = np.exp(-energy/theta)
    return e/(1+e)


def population(cat, kind):
    x = cat.count_nodes
    if kind == 'vacuum':
        return np.zeros_like(x)
    if kind == 'synthetic_fixed_count':
        return .04*np.exp(-x/.35)+.035*np.exp(-((x-.95)/.24)**2)
    raise ValueError('Unregistered population '+kind)


def equivalent_temperature(cat, amplitude, gamma, p):
    if not np.any(p):
        return 0.
    energy = cat.energy_kernel(float(amplitude), float(gamma))[0]
    target = np.dot(cat.count_weights*energy, p)
    def defect(theta):
        return np.dot(cat.count_weights*energy, fd_population(energy, theta))-target
    upper = 1.
    while defect(upper) < 0:
        upper *= 2
        if upper > 1e6:
            raise ValueError('Equivalent-temperature bracket not found')
    return float(brentq(defect, 0., upper, xtol=1e-12))


def response_records(plan, cell_scales, z, force, temperature):
    values = {}
    for name, pair in plan['mobility_pairs_ps'].items():
        mobility = KWTMobility(cell_scales, *pair)
        result = mobility.tensor_response(np.array([z.real, z.imag]), force, temperature)
        values[name] = dict(metadata=mobility.metadata(), velocity=result.velocity.tolist(),
            heat=result.heat, matrix=result.mobility.tolist(), taupsi_ps=result.taupsi_ps,
            Tmob_K=result.Tmob_K,
            power_identity_residual=float(np.dot(force, result.velocity)+result.heat))
    return values


def local_case(plan, case, cat, progress):
    # Tensor2D geometry is only a carrier for local constitutive methods here.
    model = RectangularSpatialFunctional(cat, 160e-9, 80e-9, 7e-9,
        elements_x=1, elements_y=1, degree=2, delta_regularizer_bar=case['delta_reg'])
    state = plan['local_controls'][case['control']]
    z = complex(*state['z'])
    derivatives = np.array([complex(*value) for value in state['derivatives']])
    p = population(cat, case['population'])
    q = np.imag(np.conj(z)*derivatives)/(abs(z)**2+case['delta_reg']**2)
    gamma = float(np.dot(q, q)/model.gap_ratio)
    # A 2D local energy is u(a,sum Gamma)+kappa*(sum|grad|²-a²sum q²).
    def density(zvalue, dvalue):
        a = abs(zvalue)
        qv = np.imag(np.conj(zvalue)*dvalue)/(a*a+case['delta_reg']**2)
        gv = float(np.dot(qv, qv)/model.gap_ratio)
        u = model._potential(a, gv, p, {})[0]
        return float(u+model.kappa*(np.vdot(dvalue, dvalue).real-a*a*np.dot(qv, qv)))
    symbol = model.principal_symbol(z, derivatives, p, spatial_dimensions=2)
    u, fa, fg = model._potential(abs(z), gamma, p, {})
    rho = abs(z)**2
    v = (2*q*fg/model.gap_ratio-2*model.kappa*rho*q)/(rho+case['delta_reg']**2)
    gz = fa*z/abs(z) if abs(z) else 0j
    for derivative, qi, vi in zip(derivatives, q, v):
        gz += -2*model.kappa*qi*qi*z-1j*vi*derivative-2*vi*qi*z
    gd = 2*model.kappa*derivatives+1j*v*z
    # Representative independent directional check at fixed p. Skip amplitude
    # FD at the exact normal point: its neighbourhood must remain supported.
    direction_z = .6+.8j if abs(z) else 0j
    direction_d = np.array([.3+.4j, -.2+.1j])
    expected = float(np.real(np.conj(gz)*direction_z+np.vdot(gd, direction_d)))
    estimates = []
    for h in (2e-5, 1e-5):
        estimates.append((density(z+h*direction_z, derivatives+h*direction_d)
                          -density(z-h*direction_z, derivatives-h*direction_d))/(2*h))
    theta = equivalent_temperature(cat, abs(z), gamma, p)
    force = np.array([gz.real, gz.imag])
    result = dict(kind='local_constitutive_2D', state=state, delta_reg=case['delta_reg'],
        population=case['population'], fixed_count_in_all_variations=True,
        energy_density_bar=density(z, derivatives), force_cartesian=force.tolist(),
        gradient_conjugate=[[v.real, v.imag] for v in gd], gamma_bar=gamma, q_delta=q.tolist(),
        symbol_matrix=symbol.matrix.tolist(), symbol_eigenvalues=symbol.eigenvalues.tolist(),
        symbol_uncertainty=symbol.uncertainty, symbol_stable=symbol.stable,
        verdict=('POSITIVE_LOCAL_SYMBOL' if symbol.stable else
                 'NEGATIVE_LOCAL_SYMBOL' if symbol.eigenvalues[0] < -symbol.uncertainty else 'UNRESOLVED_LOCAL_SYMBOL'),
        fd_identity=dict(expected=expected, estimates=estimates,
                         absolute_errors=[abs(v-expected) for v in estimates]),
        equivalent_temperature_K=cell_temperature(cat, theta),
        mobility=response_records(plan, scales(cat, plan), z, force, theta),
        mobility_scope='Local algebraic response to partial density derivative at fixed gradient; not full spatial Euler-Lagrange RHS',
        normal_point_fd_scope='At exactly zero only gradient variations are tested, not unsupported amplitude neighbourhood',
        no_trajectory=True, physical_core_admitted=False)
    progress(1, 1, 'local result')
    return result


def cell_temperature(cat, theta):
    return float(theta*cat.vacuum.delta0_J/K_B_J_K)


def spatial_field(model, profile):
    x, y = model.dof_coordinates_bar.T
    x = x-model.length_bar/2
    if profile == 'smooth':
        envelope = np.exp(-((x/8)**2+(y/6)**2))
        amplitude = .95-.08*envelope
    elif profile == 'suppressed':
        envelope = np.exp(-((x/2.5)**2+(y/2.5)**2))
        amplitude = .95-.91*envelope
    else:
        raise ValueError('Unregistered field profile')
    phase = .07*x+.025*np.sin(y/5)*np.exp(-(x/10)**2)
    return amplitude*np.exp(1j*phase)


def spatial_case(plan, case, cat, progress, output):
    model = RectangularSpatialFunctional(cat, plan['geometry']['length_m'],
        plan['material']['width_m'], plan['material']['thickness_m'],
        elements_x=case['elements_x'], elements_y=case['elements_y'],
        degree=case['degree'], delta_regularizer_bar=case['delta_reg'])
    z = spatial_field(model, case['profile'])
    fields = model.sample_fields(z)
    p = np.broadcast_to(population(cat, case['population']), (model.cells, len(cat.count_nodes))).copy()
    # Diagnostic fields may have a negative symbol. Record that sign without
    # evolving, clipping or turning it into an execution exception.
    base = model.evaluate(z, p, require_stability=False,
        on_quadrature=lambda i, n: progress(i+1, n, 'energy and force'))
    eigen, uncertainty, stable = [], [], []
    for i in range(model.cells):
        s = model.principal_symbol(z[i], fields.derivative_quadrature_bar[i], p[i])
        eigen.append(s.eigenvalues); uncertainty.append(s.uncertainty); stable.append(s.stable)
        progress(i+1, model.cells, 'local stability')
    temperatures = []
    for i, (v, g, pp) in enumerate(zip(z, fields.gamma_quadrature_bar, p)):
        temperatures.append(equivalent_temperature(cat, float(abs(v)), float(g), pp))
        progress(i+1, model.cells, 'equivalent temperatures')
    theta = np.asarray(temperatures)
    heat_fields = {}
    mobility_values = {}
    for name, pair in plan['mobility_pairs_ps'].items():
        mob = KWTMobility(scales(cat, plan), *pair)
        response = kwt_spatial_response(z, base.gradient_cartesian_bar,
            model.mass_bar, mob, theta, np.zeros(model.cells))
        mobility_values[name] = dict(metadata=mob.metadata(),
            total_heat_bar=response.condensate_heat_rate_bar,
            field_work_bar=response.field_energy_rate_bar,
            identity_residual_bar=response.identity_residual_bar,
            maximum_velocity_per_ps=float(np.max(abs(response.material_velocity_bar))),
            interpretation='Unconstrained instantaneous descent including boundary DOFs; no reservoir/circuit/trajectory')
        heat_fields['heat_'+name] = response.heat_density_bar
    # One field and one gauge direction, tested at two steps. p is literally
    # fixed; this is internal energy, not thermal internal-energy differentiation.
    xx, yy = model.dof_coordinates_bar.T
    dz = np.sin(np.pi*xx/model.length_bar)*(1+.2*np.cos(yy))*complex(.6,.8)
    links = .1*np.sin(np.arange(len(model.graph_edges))*.7)
    expected_z = float(np.sum(base.gradient_cartesian_bar[:,0]*dz.real+base.gradient_cartesian_bar[:,1]*dz.imag))
    expected_link = float(np.dot(base.link_current_bar, links))
    identities = []
    for h in (2e-5, 1e-5):
        def callback(label):
            return lambda i,n:progress(i+1,n,f'{label} h={h:g}')
        plus = model.evaluate(z+h*dz, p, require_stability=False,
            on_quadrature=callback('field+')).energy_bar
        minus = model.evaluate(z-h*dz, p, require_stability=False,
            on_quadrature=callback('field-')).energy_bar
        plus_link = model.evaluate(z, p, link_phases=h*links, require_stability=False,
            on_quadrature=callback('gauge+')).energy_bar
        minus_link = model.evaluate(z, p, link_phases=-h*links, require_stability=False,
            on_quadrature=callback('gauge-')).energy_bar
        identities.append(dict(step=h, field_derivative=(plus-minus)/(2*h),
            expected_field=expected_z, gauge_derivative=(plus_link-minus_link)/(2*h),
            expected_gauge=expected_link))
        progress(len(identities), 2, 'field/current identities')
    eigen = np.asarray(eigen);uncertainty = np.asarray(uncertainty)
    np.savez_compressed(output/'fields.npz', coordinates_m=model.dof_coordinates_bar*model.ell0_m,
        shape=np.array(model.rectangle_shape), delta=z, gamma=fields.gamma_quadrature_bar,
        cartesian_force=base.gradient_cartesian_bar, symbol_eigenvalues=eigen,
        symbol_uncertainty=uncertainty, temperature_equivalent_K=theta*cat.vacuum.delta0_J/K_B_J_K,
        quadrature_volume_m3=model.node_volumes_m3, **heat_fields)
    return dict(kind='full2D_static_prescribed_field', profile=case['profile'], nodes=model.cells,
        elements=[case['elements_x'], case['elements_y']], degree=case['degree'],
        delta_reg=case['delta_reg'], population=case['population'],
        minimum_amplitude=float(min(abs(z))), energy_bar=base.energy_bar, energy_J=base.energy_J,
        integrated_force_norm=float(np.linalg.norm(base.gradient_cartesian_bar)),
        maximum_current_A=float(np.max(abs(base.current_A))),
        noether_max=float(np.max(abs(base.noether_residual_bar))),
        symbol_min_eigenvalue=float(eigen[:,0].min()), maximum_symbol_uncertainty=float(uncertainty.max()),
        positive_symbols=int(sum(stable)), negative_symbols=int(np.sum(eigen[:,0]<-uncertainty)),
        unresolved_symbols=int(np.sum(abs(eigen[:,0])<=uncertainty)),
        identity_checks=identities, mobility=mobility_values,
        fields_file='fields.npz', no_trajectory=True, physical_core_admitted=False,
        photon=False, boundaries='Independent prescribed input values; no current-driven reservoir or dynamic boundary admission')
