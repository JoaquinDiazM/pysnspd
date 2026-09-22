"""Experimental conservative event-limited SSP steps; NOT an admitted solver.

Physical kernels remain unchanged. Each forward-Euler map reserves the regular
BGK/heating/escape update, then limits a COMMON extent for each kinetic event
against global electron/hole/phonon inventories. All face and reaction demands
share those inventories. The same extents update populations and ledgers.

The limiter changes the finite-step approximation, not the continuous kernel.
Accuracy/order and the nonlinear moving-condensate energy ledger MUST be
revalidated. Frozen-field kinetic energy/count invariants are algebraic.
"""
from __future__ import annotations
import numpy as np
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.cell_transport import NativeTransportEvents

SAFETY = .9


def _add(stats, key, value):
    stats[key] = stats.get(key, 0) + value


def _ratio(available, demand, h):
    result = np.ones_like(available)
    active = demand > 0
    raw = np.divide(available, h*demand, out=np.ones_like(available), where=active)
    constrained = active & (raw <= 1.)
    result[constrained] = SAFETY*raw[constrained]
    if np.any(~np.isfinite(result)) or np.any(result < 0):
        raise FloatingPointError("invalid availability ratio; no population repair")
    return result


def _demands(indices, coefficients, rates, size):
    particle, hole = np.zeros(size), np.zeros(size)
    for ix, coefficient in zip(indices, coefficients):
        signed = coefficient*rates
        particle += np.bincount(ix, weights=np.maximum(-signed, 0.), minlength=size)
        hole += np.bincount(ix, weights=np.maximum(signed, 0.), minlength=size)
    return particle, hole


def _electron_stoichiometry(event):
    """At most four slots/event; merge repeated nodes before donor budgets."""
    if hasattr(event, "electron_lower_i"):
        indices = [event.electron_lower_i, event.electron_upper_i,
                   event.electron_lower_j, event.electron_upper_j]
        sign = np.where(event.recombination, -1., 1.)
        coefficients = [sign*event.electron_beta_lower_i,
                        sign*event.electron_beta_upper_i,
                        -event.electron_beta_lower_j, -event.electron_beta_upper_j]
    else:
        indices = [event.index_i, event.index_j]
        coefficients = [np.where(event.recombination, -1., 1.), -np.ones_like(event.omega)]
    coefficients = [np.array(c, copy=True) for c in coefficients]
    for k in range(1, len(indices)):
        for j in range(k):
            same = indices[k] == indices[j]
            coefficients[j] += np.where(same, coefficients[k], 0.)
            coefficients[k][same] = 0.
    return indices, coefficients


def _limit(theta, indices, coefficients, rates, particle_ratio, hole_ratio=None):
    for ix, coefficient in zip(indices, coefficients):
        signed = coefficient*rates
        theta = np.minimum(theta, np.where(signed < 0, particle_ratio[ix], 1.))
        if hole_ratio is not None:
            theta = np.minimum(theta, np.where(signed > 0, hole_ratio[ix], 1.))
    return theta


def _scatter(indices, coefficients, rates, size):
    result = np.zeros(size)
    for ix, coefficient in zip(indices, coefficients):
        result += np.bincount(ix, weights=coefficient*rates, minlength=size)
    return result


def _record_flux(stats, rates, theta, energies, h, kind):
    weight = stats.get("_stage_weight", 1.)
    limited = theta < 1.
    _add(stats, "event_evaluations", int(len(rates)))
    _add(stats, "limited_event_evaluations", int(np.count_nonzero(limited)))
    _add(stats, kind+"_limited_evaluations", int(np.count_nonzero(limited)))
    stats["minimum_event_factor"] = min(stats.get("minimum_event_factor", 1.), float(theta.min()) if len(theta) else 1.)
    omitted = abs(rates)*(1-theta)
    _add(stats, "integrated_absolute_extent_defect", float(weight*h*np.sum(omitted)))
    _add(stats, "integrated_energy_weighted_flux_defect", float(weight*h*np.dot(energies, omitted)))
    _add(stats, "integrated_absolute_raw_extent", float(weight*h*np.sum(abs(rates))))


def step(system, t, z, h, stats=None):
    """One physical limited FE map. Fail rather than repair a bad baseline."""
    if not np.isfinite(h) or h <= 0:
        raise ValueError("positive finite step required")
    stats = {} if stats is None else stats
    amplitudes, p, n = system.unpack(z)
    _add(stats, "forward_euler_stages", 1)
    system.rhs_calls += 1
    cells = [ElectronicCell(system.catalog, amplitudes[i], system.gammas[i]) for i in range(system.cell_count)]
    capacity = np.asarray([4*cell.weights for cell in cells])
    phcap = system.phonons.capacities
    base = np.array(z, copy=True)
    base_ledger = base[system.population_size:-1].reshape(system.cell_count, 4)
    stage_weight = stats.get("_stage_weight", 1.)
    for i, cell in enumerate(cells):
        temperature = cell.equivalent_temperature(p[i])
        force = cell.moments(p[i])[1]
        motion = system.mobility.amplitude_response(amplitudes[i], force, temperature)
        heating = cell.heating(p[i], system.external_powers[i]+motion.heat, system.bath_temperature)
        # Positive affine form avoids subtracting a nearly depleted BGK tail.
        fraction = h/system.tau_kin
        if fraction > 1:
            raise ValueError("regular BGK forward-Euler step exceeds its convex interval")
        pbase = (1-fraction)*p[i]+fraction*cell.fermi_dirac(temperature)+h*heating
        escape_fraction = h/system.tau_escape
        if escape_fraction > 1:
            raise ValueError("regular escape step exceeds its convex interval")
        nbase = (1-escape_fraction)*n[i]+escape_fraction*system.bath_phonons
        escaping = (n[i]-system.bath_phonons)/system.tau_escape
        start = i*system.block_size
        base[start] = amplitudes[i]+h*motion.velocity
        base[start+1:start+1+system.electron_size] = pbase
        base[start+1+system.electron_size:start+system.block_size] = nbase
        escape_power = float(np.dot(phcap*system.phonons.energies, escaping))
        base_ledger[i] += h*np.array([escape_power, system.external_powers[i], 0., motion.heat])
        _add(stats, "integrated_external_input", stage_weight*h*system.external_powers[i])
        _add(stats, "integrated_condensate_heat", stage_weight*h*motion.heat)
    try:
        abase, pbase, nbase = system.unpack(base)
        if np.any(abase < 0):
            raise ValueError("regular amplitude left its physical domain")
    except ValueError as exc:
        _add(stats, "regular_baseline_rejections", 1)
        raise ValueError("regular BGK/heat/escape baseline is not physical; reduce the whole step") from exc
    face = None
    face_rates = None
    face_particle, face_hole = np.zeros_like(p), np.zeros_like(p)
    if system.cell_count == 2:
        face = NativeTransportEvents(*cells, system.diffusion_over_length_squared, system.face_order)
        face_rates = face.rates(p)
        for side in (0, 1):
            ix = [face.indices[side, :, k] for k in (0, 1)]
            coeff = [(2*side-1)*face.barycentric[side, :, k] for k in (0, 1)]
            face_particle[side], face_hole[side] = _demands(ix, coeff, face_rates, system.electron_size)
    particle_ratios, hole_ratios = [], []
    result = base.copy()
    ledger = result[system.population_size:-1].reshape(system.cell_count, 4)
    for i, cell in enumerate(cells):
        # Keep exactly one large reaction network live at a time.
        event = system.reaction_events(cell)
        rates = event.rates(p[i], n[i])
        ix, coeff = _electron_stoichiometry(event)
        demand_p, demand_h = _demands(ix, coeff, rates, system.electron_size)
        pix = [event.phonon_lower, event.phonon_upper]
        pcoeff = [event.beta_lower, event.beta_upper]
        demand_n, _ = _demands(pix, pcoeff, rates, system.phonon_size)
        rp = _ratio(capacity[i]*pbase[i], demand_p+face_particle[i], h)
        rh = _ratio(capacity[i]*(1-pbase[i]), demand_h+face_hole[i], h)
        rn = _ratio(phcap*nbase[i], demand_n, h)
        particle_ratios.append(rp)
        hole_ratios.append(rh)
        theta = _limit(np.ones_like(rates), ix, coeff, rates, rp, rh)
        theta = _limit(theta, pix, pcoeff, rates, rn)
        actual = theta*rates
        count_change = _scatter(ix, coeff, actual, system.electron_size)
        phonon_change = _scatter(pix, pcoeff, actual, system.phonon_size)
        start = i*system.block_size
        result[start+1:start+1+system.electron_size] += h*count_change/capacity[i]
        result[start+1+system.electron_size:start+system.block_size] += h*phonon_change/phcap
        ledger[i, 2] += h*float(np.dot(event.omega, actual))
        _record_flux(stats, rates, theta, event.omega, h, "reaction")
        del event, rates, ix, coeff, pix, pcoeff, theta, actual
    if face is not None:
        theta = np.ones_like(face_rates)
        for side in (0, 1):
            ix = [face.indices[side, :, k] for k in (0, 1)]
            coeff = [(2*side-1)*face.barycentric[side, :, k] for k in (0, 1)]
            theta = _limit(theta, ix, coeff, face_rates, particle_ratios[side], hole_ratios[side])
        actual = theta*face_rates
        for side in (0, 1):
            ix = [face.indices[side, :, k] for k in (0, 1)]
            coeff = [(2*side-1)*face.barycentric[side, :, k] for k in (0, 1)]
            start = side*system.block_size+1
            result[start:start+system.electron_size] += h*_scatter(ix, coeff, actual, system.electron_size)/capacity[side]
        result[-1] += h*float(np.dot(face.energies, actual))
        _record_flux(stats, face_rates, theta, face.energies, h, "transport")
    system.unpack(result)
    if np.any(~np.isfinite(result)):
        raise FloatingPointError("nonfinite limited step; no repair applied")
    return result


def integrate(system, initial, duration, steps, *, callback=None, stats=None):
    """SSPRK3 convex composition of the limited FE map; no adaptive retry."""
    if not np.isfinite(duration) or duration <= 0 or type(steps) is not int or steps < 1:
        raise ValueError("positive duration and integer step count required")
    stats = {} if stats is None else stats
    stats.update({"method": "experimental conservative-event-limited SSPRK3", "safety": SAFETY,
        "status": "IN_PROGRESS_NOT_NUMERICALLY_ADMITTED",
        "notes": ["Physical kernels unchanged; common event extents preserve frozen-field invariants.",
                  "Formal SSPRK3 order holds only where limiter is inactive and the RHS is sufficiently smooth.",
                  "Moving-condensate ledger and temporal convergence require independent measurement.",
                  "Event/flux defects are numerical modifications, not hidden energy corrections."]})
    times = np.linspace(0., duration, steps+1)
    states = np.empty((steps+1, len(initial)))
    states[0] = initial
    system.unpack(initial)
    h = duration/steps
    for k in range(steps):
        t, z = times[k], states[k]
        stats["_stage_weight"] = 1/6
        first = step(system, t, z, h, stats)
        stats["_stage_weight"] = 1/6
        second = .75*z+.25*step(system, t+h, first, h, stats)
        system.unpack(second)
        stats["_stage_weight"] = 2/3
        states[k+1] = (1/3)*z+(2/3)*step(system, t+h/2, second, h, stats)
        system.unpack(states[k+1])
        stats.pop("_stage_weight", None)
        stats["completed_steps"] = k+1
        if callback is not None and (k == 0 or (k+1) % 10 == 0 or k+1 == steps):
            callback({"event": "limited_ssp_progress", "macrostep": k+1,
                "time": float(times[k+1]), "steps": steps,
                "limited_event_evaluations": stats.get("limited_event_evaluations", 0),
                "integrated_absolute_extent_defect": stats.get("integrated_absolute_extent_defect", 0.),
                "integrated_energy_weighted_flux_defect": stats.get("integrated_energy_weighted_flux_defect", 0.),
                "minimum_event_factor": stats.get("minimum_event_factor", 1.)})
    stats["status"] = "COMPLETED_REQUIRES_CONVERGENCE_ASSESSMENT"
    return times, states, stats
