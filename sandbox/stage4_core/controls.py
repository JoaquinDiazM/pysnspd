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
        self.prefetched = {}
        self.prefetch_hits = 0
        self.strict_prefetch = False
        self.energy_kernel = lru_cache(maxsize=4096)(self._kernel)

    def _kernel(self, amplitude, gamma):
        key = (float(amplitude), float(gamma))
        if key in self.prefetched:
            self.prefetch_hits += 1
            return self.prefetched[key]
        if self.strict_prefetch:
            raise RuntimeError('Unscheduled spectral key; no implicit serial fallback: '+repr(key))
        start = time.monotonic()
        result = self.source.energy_kernel(amplitude, gamma)
        self.calls += 1
        self.seconds += time.monotonic()-start
        for value in result:
            value.setflags(write=False)
        return result

    def install_prefetched(self, values):
        """Exact evaluated keys only: no rounding, interpolation or projection."""
        checked = {}
        for key, kernels in values.items():
            if len(key) != 2 or len(kernels) != 3:
                raise ValueError('Expected exact (amplitude,Gamma) and three kernels')
            arrays = tuple(np.asarray(value, dtype=float) for value in kernels)
            if any(value.shape != self.count_nodes.shape or np.any(~np.isfinite(value)) for value in arrays):
                raise ValueError('Invalid prefetched spectral kernel')
            for value in arrays:
                value.setflags(write=False)
            checked[(float(key[0]), float(key[1]))] = arrays
        self.energy_kernel.cache_clear()
        self.prefetched = checked
        self.strict_prefetch = True

    def diagnostics(self):
        return dict(count_states=len(self.count_nodes), direct_queries=self.calls,
                    direct_seconds=self.seconds, exact_cache=self.energy_kernel.cache_info()._asdict(),
                    prefetched_keys=len(self.prefetched), prefetch_hits=self.prefetch_hits,
                    strict_prefetch=self.strict_prefetch)


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


def spatial_profile(coordinates_bar, length_bar, profile):
    """Prescribed complex field and its exact Cartesian derivatives in X,Y.

    These are diagnostic input functions, not a stationary solution or a
    photon deposition profile. Their exact derivatives provide a reference
    independent of the GLL differentiation of their sampled values.
    """
    x, y = np.asarray(coordinates_bar, float).T
    x = x-length_bar/2
    if profile == 'smooth':
        sx, sy, suppression = 8., 6., .08
    elif profile == 'suppressed':
        sx, sy, suppression = 2.5, 2.5, .91
    else:
        raise ValueError('Unregistered field profile')
    envelope = np.exp(-((x/sx)**2+(y/sy)**2))
    amplitude = .95-suppression*envelope
    phase_envelope = np.exp(-(x/10)**2)
    phase = .07*x+.025*np.sin(y/5)*phase_envelope
    da = np.column_stack((2*suppression*x/sx**2*envelope,
                          2*suppression*y/sy**2*envelope))
    dp = np.column_stack((.07-.025*np.sin(y/5)*phase_envelope*2*x/100,
                          .005*np.cos(y/5)*phase_envelope))
    rotation = np.exp(1j*phase)
    return amplitude*rotation, rotation[:, None]*(da+1j*amplitude[:, None]*dp)


def spatial_field(model, profile):
    return spatial_profile(model.dof_coordinates_bar, model.length_bar, profile)[0]


def boundary_mask(model):
    """Four prescribed sides, with each corner counted exactly once."""
    mask = np.zeros(model.cells, bool)
    mask[np.concatenate(tuple(model.boundary_dofs.values()))] = True
    return mask


def volume_force_metrics(model, cartesian_gradient, constrained):
    """Norms of the volume derivative G/m, distinguishing boundary reactions."""
    gradient = np.asarray(cartesian_gradient)
    density = gradient/model.mass_bar[:, None]
    square = np.sum(density*density, axis=1)
    result = {}
    for name, mask in (('all_dofs_including_boundary_reactions', np.ones(model.cells, bool)),
                       ('interior', ~constrained)):
        measure = float(np.sum(model.mass_bar[mask]))
        integral = float(np.dot(model.mass_bar[mask], square[mask]))
        result[name] = dict(mass_bar=measure, l2_volume_bar=float(np.sqrt(integral)),
            rms_volume_bar=float(np.sqrt(integral/measure)) if measure else None,
            maximum_volume_force_bar=float(np.sqrt(square[mask].max())) if np.any(mask) else None)
    return result, density


def section_currents(model, current_A):
    """Oriented graph flux through every x cut; positive means left to right.

    The incidence convention is B_tail=+1, B_head=-1. Edge currents already
    include transverse quadrature weights, so they are summed without another
    area factor. A maximum individual edge current is not a section current.
    """
    x = model.dof_coordinates_bar[:, 0]
    unique_x = np.unique(x)
    cuts = (unique_x[:-1]+unique_x[1:])/2
    tail, head = x[model.graph_edges[:, 0]], x[model.graph_edges[:, 1]]
    current = np.asarray(current_A)
    flux = np.array([np.sum(current[(tail < cut)&(head > cut)])
                     -np.sum(current[(head < cut)&(tail > cut)]) for cut in cuts])
    return cuts*model.ell0_m, flux


def derivative_diagnostics(model, profile, sampled, constrained):
    _, exact = spatial_profile(model.dof_coordinates_bar, model.length_bar, profile)
    error = sampled.derivative_quadrature_bar-exact
    result = {}
    for name, mask in (('all_nodes', np.ones(model.cells, bool)), ('interior', ~constrained)):
        integral = float(np.dot(model.mass_bar[mask], np.sum(abs(error[mask])**2, axis=1)))
        norm = float(np.dot(model.mass_bar[mask], np.sum(abs(exact[mask])**2, axis=1)))
        result[name] = dict(l2_error_bar=float(np.sqrt(integral)),
            relative_l2_error=float(np.sqrt(integral/norm)) if norm else None,
            maximum_component_error=float(np.max(abs(error[mask]))) if np.any(mask) else None)
    z = sampled.delta_quadrature_bar
    exact_q = np.imag(np.conj(z)[:, None]*exact)/(abs(z[:, None])**2+model.delta_regularizer_bar**2)
    exact_gamma = np.sum(exact_q*exact_q, axis=1)/model.gap_ratio
    center = int(np.argmin(np.sum((model.dof_coordinates_bar-[model.length_bar/2, 0.])**2, axis=1)))
    result['center'] = dict(node=center, coordinates_bar=model.dof_coordinates_bar[center].tolist(),
        derivative_discrete=[[v.real, v.imag] for v in sampled.derivative_quadrature_bar[center]],
        derivative_analytic=[[v.real, v.imag] for v in exact[center]],
        derivative_error_norm=float(np.linalg.norm(error[center])),
        gamma_discrete=float(sampled.gamma_quadrature_bar[center]),
        gamma_analytic=float(exact_gamma[center]),
        gamma_ratio_discrete_to_analytic=float(sampled.gamma_quadrature_bar[center]/exact_gamma[center]) if exact_gamma[center] else None)
    return result, exact, exact_gamma


def mobility_record(response, constrained):
    return dict(total_heat_bar=response.condensate_heat_rate_bar,
        interior_heat_bar=float(np.dot(response.quadrature_mass_bar[~constrained], response.heat_density_bar[~constrained])),
        boundary_heat_bar=float(np.dot(response.quadrature_mass_bar[constrained], response.heat_density_bar[constrained])),
        field_work_bar=response.field_energy_rate_bar,
        boundary_work_bar=response.boundary_work_rate_bar,
        identity_residual_bar=response.identity_residual_bar,
        maximum_velocity_per_ps=float(np.max(abs(response.material_velocity_bar))),
        maximum_interior_velocity_per_ps=float(np.max(abs(response.material_velocity_bar[~constrained]))) if np.any(~constrained) else None,
        maximum_boundary_velocity_per_ps=float(np.max(abs(response.material_velocity_bar[constrained]))))


def spatial_case(plan, case, cat, progress, output):
    boundary_policy = case.get('boundary_response', 'unconstrained_instantaneous')
    if boundary_policy not in ('unconstrained_instantaneous', 'fixed_prescribed_all_sides'):
        raise ValueError('Unregistered boundary response '+boundary_policy)
    model = RectangularSpatialFunctional(cat, plan['geometry']['length_m'],
        plan['material']['width_m'], plan['material']['thickness_m'],
        elements_x=case['elements_x'], elements_y=case['elements_y'],
        degree=case['degree'], delta_regularizer_bar=case['delta_reg'])
    z = spatial_field(model, case['profile'])
    fields = model.sample_fields(z)
    constrained = boundary_mask(model)
    derivative_check, exact_derivative, exact_gamma = derivative_diagnostics(
        model, case['profile'], fields, constrained)
    p = np.broadcast_to(population(cat, case['population']), (model.cells, len(cat.count_nodes))).copy()
    # Diagnostic fields may have a negative symbol. Record that sign without
    # evolving, clipping or turning it into an execution exception.
    base = model.evaluate(z, p, require_stability=False,
        on_quadrature=lambda i, n: progress(i+1, n, 'energy and force'))
    force_metrics, force_density = volume_force_metrics(model, base.gradient_cartesian_bar, constrained)
    cuts_m, cuts_A = section_currents(model, base.current_A)
    constraint_load = np.zeros_like(base.gradient_cartesian_bar)
    if boundary_policy == 'fixed_prescribed_all_sides':
        # External conjugate load +G cancels the force -G at constrained DOFs.
        # With phi=0 this yields zero material/field velocity without overwriting
        # a computed velocity, and the constraint work is identically zero.
        constraint_load[constrained] = base.gradient_cartesian_bar[constrained]
    eigen, uncertainty, stable = [], [], []
    for i in range(model.cells):
        s = model.principal_symbol(z[i], fields.derivative_quadrature_bar[i], p[i])
        eigen.append(s.eigenvalues); uncertainty.append(s.uncertainty); stable.append(s.stable)
        progress(i+1, model.cells, 'local stability')
    center = derivative_check['center']['node']
    center_analytic = model.principal_symbol(z[center], exact_derivative[center], p[center])
    derivative_check['center'].update(
        discrete_symbol_eigenvalues=np.asarray(eigen[center]).tolist(),
        discrete_symbol_uncertainty=float(uncertainty[center]),
        analytic_symbol_eigenvalues=center_analytic.eigenvalues.tolist(),
        analytic_symbol_uncertainty=float(center_analytic.uncertainty),
        interpretation='Signs refer to their specified gradients; a coarse discrete negative sign alone does not establish failure of this analytic profile')
    temperatures = []
    for i, (v, g, pp) in enumerate(zip(z, fields.gamma_quadrature_bar, p)):
        temperatures.append(equivalent_temperature(cat, float(abs(v)), float(g), pp))
        progress(i+1, model.cells, 'equivalent temperatures')
    theta = np.asarray(temperatures)
    heat_fields = {}
    mobility_values = {}
    fixed_mobility_values = {}
    for name, pair in plan['mobility_pairs_ps'].items():
        mob = KWTMobility(scales(cat, plan), *pair)
        response = kwt_spatial_response(z, base.gradient_cartesian_bar,
            model.mass_bar, mob, theta, np.zeros(model.cells))
        mobility_values[name] = dict(metadata=mob.metadata(), **mobility_record(response, constrained),
            interpretation='Unconstrained instantaneous descent including boundary DOFs; no reservoir/circuit/trajectory')
        heat_fields['heat_'+name] = response.heat_density_bar
        heat_fields['velocity_unconstrained_'+name] = response.material_velocity_bar
        if boundary_policy == 'fixed_prescribed_all_sides':
            fixed = kwt_spatial_response(z, base.gradient_cartesian_bar,
                model.mass_bar, mob, theta, np.zeros(model.cells), boundary_load_bar=constraint_load)
            fixed_mobility_values[name] = dict(metadata=mob.metadata(), **mobility_record(fixed, constrained),
                interpretation='Instantaneous interior descent with all four sides fixed by conjugate reaction loads; diagnostic constraint, not device boundary conditions')
            heat_fields['heat_fixed_boundary_'+name] = fixed.heat_density_bar
            heat_fields['velocity_fixed_boundary_'+name] = fixed.material_velocity_bar
    # One field and one gauge direction, tested at two steps. p is literally
    # fixed; this is internal energy, not thermal internal-energy differentiation.
    xx, yy = model.dof_coordinates_bar.T
    dz = np.sin(np.pi*xx/model.length_bar)*(1+.2*np.cos(yy))*complex(.6,.8)
    links = .1*np.sin(np.arange(len(model.graph_edges))*.7)
    expected_z = float(np.sum(base.gradient_cartesian_bar[:,0]*dz.real+base.gradient_cartesian_bar[:,1]*dz.imag))
    expected_link = float(np.dot(base.link_current_bar, links))
    identities = []
    for h in ((2e-5, 1e-5) if case.get('finite_differences', True) else ()):
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
        quadrature_volume_m3=model.node_volumes_m3, quadrature_mass_bar=model.mass_bar,
        boundary_mask=constrained, interior_mask=~constrained, constraint_load_cartesian_bar=constraint_load,
        volume_force_cartesian_bar=force_density, derivative_discrete=fields.derivative_quadrature_bar,
        derivative_analytic=exact_derivative, gamma_analytic=exact_gamma,
        graph_edges=model.graph_edges, edge_current_A=base.current_A,
        section_x_m=cuts_m, section_current_A=cuts_A, **heat_fields)
    return dict(kind='full2D_static_prescribed_field', profile=case['profile'], nodes=model.cells,
        elements=[case['elements_x'], case['elements_y']], degree=case['degree'],
        delta_reg=case['delta_reg'], population=case['population'],
        minimum_amplitude=float(min(abs(z))), energy_bar=base.energy_bar, energy_J=base.energy_J,
        integrated_force_norm=float(np.linalg.norm(base.gradient_cartesian_bar)),
        maximum_current_A=float(np.max(abs(base.current_A))),
        historical_nodal_metrics=dict(integrated_force_norm='Unweighted Euclidean norm of integrated DOF gradients, not a volume-force norm',
            maximum_current_A='Maximum single weighted graph-edge current, not device section current',
            mobility='Unconstrained historical diagnostic including boundary reactions'),
        volume_force_norms=force_metrics,
        section_currents=dict(x_m=cuts_m.tolist(), current_A=cuts_A.tolist(),
            maximum_absolute_A=float(np.max(abs(cuts_A))), orientation='Positive left-to-right; sum of oriented graph edges crossing each x cut'),
        derivative_diagnostics=derivative_check,
        response_boundary_policy=boundary_policy,
        constraint_reaction_ledger=dict(active=boundary_policy == 'fixed_prescribed_all_sides',
            unique_boundary_nodes=int(np.sum(constrained)),
            external_conjugate_load_sum_cartesian_bar=np.sum(constraint_load, axis=0).tolist(),
            integrated_load_l2_bar=float(np.linalg.norm(constraint_load)),
            phase_torque_bar=float(np.sum(np.imag(np.conj(z)*(constraint_load[:, 0]+1j*constraint_load[:, 1])))),
            sign='Applied conjugate load +G at fixed nodes cancels intrinsic force -G; each corner counted once',
            array_file='fields.npz:constraint_load_cartesian_bar',
            scope='Prescribed fixed control values, not a superconducting reservoir or device boundary admission'),
        noether_max=float(np.max(abs(base.noether_residual_bar))),
        symbol_min_eigenvalue=float(eigen[:,0].min()), maximum_symbol_uncertainty=float(uncertainty.max()),
        positive_symbols=int(sum(stable)), negative_symbols=int(np.sum(eigen[:,0]<-uncertainty)),
        unresolved_symbols=int(np.sum(abs(eigen[:,0])<=uncertainty)),
        identity_checks=identities, mobility=mobility_values,
        mobility_fixed_boundary=fixed_mobility_values,
        identity_scope=('Direct two-step field/current checks' if identities else
                        'Uses the same tested assembly; finite differences only at registered coarse anchors'),
        fields_file='fields.npz', no_trajectory=True, physical_core_admitted=False,
        photon=False, boundaries='Independent prescribed input values; no current-driven reservoir or dynamic boundary admission')
