"""Preregistered fixed-field boundary experiments, using the frozen SSP map.

No production source is changed. A directional view selects one exact positive
forward/reverse activity from an existing projected event network. This is an
isolated suboperator test, not the equilibrium of the full reversible physics.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(ROOT/'sandbox/stage2_cells'),
               str(ROOT/'sandbox/stage2_cells/recovery_20260921')]
import numpy as np
import scipy
from run_coupled import setup
from limited_ssp import integrate
from pysnspd.experimental.cell_validation import ElectronicCell
from pysnspd.experimental.coupled_cells import CoupledCellSystem
from pysnspd.experimental.kinetic_events import ProjectedElectronPhononEvents

OUT = ROOT/'docs/implementation/stage2/closure_prep_20260922'
CRITERIA = ROOT/'docs/implementation/stage2/acceptance_criteria.json'
PLAN = OUT/'isolated_registration.json'
FIELDS = ('target_energy_i', 'target_energy_j', 'omega', 'coefficients',
          'recombination', 'phonon_lower', 'phonon_upper', 'beta_lower', 'beta_upper',
          'electron_lower_i', 'electron_upper_i', 'electron_beta_lower_i',
          'electron_beta_upper_i', 'electron_lower_j', 'electron_upper_j',
          'electron_beta_lower_j', 'electron_beta_upper_j')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', encoding='utf-8') as f:
        json.dump(value, f, indent=2, allow_nan=False)
        f.write('\n')


def sources():
    paths = list((ROOT/'pysnspd/experimental').glob('*.py'))
    paths += [Path(__file__), ROOT/'sandbox/stage2_cells/run_coupled.py',
              ROOT/'sandbox/stage2_cells/recovery_20260921/limited_ssp.py', CRITERIA,
              ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz']
    return {p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(paths)}


def registration():
    return dict(schema='pysnspd.stage2.isolated-registration.v1',
        status='PREREGISTERED_NOT_EXECUTED', source_hashes=sources(),
        candidate=dict(electron_refinement=1, electron_states=630, phonon_nodes=1025,
                       infrared=.005, face_order=2, reaction_order=2,
                       reaction_layout='resolved_panels'),
        isolated=dict(cases=['emission', 'absorption', 'recombination', 'creation', 'transport'],
                      amplitudes=[.5], transport_amplitudes=[.5, 1.], gammas=0.,
                      duration=.02, steps=20, saved_times=21,
                      filled_count_band=[1.1, 1.5],
                      population='p=1 strictly inside the band, p=0 outside; creation starts at p=0',
                      phonons='n=0 for emission/recombination/transport; n=.05 for absorption/creation',
                      disabled=['BGK', 'condensate motion and heating', 'external heating', 'escape',
                                'every reaction channel outside the named directional suboperator'],
                      transport='full bidirectional native-energy transport; all reactions disabled'),
        escape=dict(duration=.2, steps=[10, 20, 40], tau=.7, amplitude=.5,
                    electron_temperature=.2, phonons='bath+.02*exp(-Omega)',
                    reference='bath+(n_initial-bath)*exp(-t/tau)',
                    expected_order=3, minimum_error_reduction=6.,
                    reduction_reason='8 times 0.75: frozen 25% smooth-order margin'),
        gates=dict(energy_ledger_scaled_max=1e-7, count_invariant_scaled_max=1e-12,
                   instantaneous_energy_scaled_max=1e-11, event_energy_relative_max=1e-12,
                   strict_populations=True, minimum_samples=10,
                   require_nonzero_evolution=True, clipping=False, energy_repair=False),
        execution=dict(threads=1, wall_timeout_seconds=240,
                       expected_seconds='30-90; uncertain', max_attempts=1),
        limitations=['Directional tests disable the inverse reaction intentionally; no detailed-balance claim.',
                     'Fields are fixed and networks are cached only because fields never change.',
                     'These are boundary/invariant tests, not full moving-condensate time or mesh admission.',
                     'Mixed exact Pauli faces are tested; globally saturated electrons are outside the positive-temperature inversion used by this SSP adapter.'])


class DirectionalEvents(ProjectedElectronPhononEvents):
    """Keep exact geometry; select a family and one original log-activity."""
    def __init__(self, parent, mask, direction):
        self.__dict__.update(parent.__dict__)
        for name in FIELDS:
            data = np.array(getattr(parent, name)[mask], copy=True)
            data.setflags(write=False)
            setattr(self, name, data)
        self.direction = direction

    def rates(self, p, n):
        forward, backward = super().log_activities(p, n)
        if self.direction == 'forward':
            return self.coefficients*np.exp(forward)
        if self.direction == 'backward':
            return -self.coefficients*np.exp(backward)
        if self.direction == 'disabled' and len(self.omega) == 0:
            return np.zeros(0)
        raise ValueError('unregistered directional component')


class FixedSystem(CoupledCellSystem):
    def reaction_events(self, cell):
        if cell.gamma != 0 or cell.amplitude not in self.fixed_events:
            raise ValueError('a fixed-field experiment changed its field')
        return self.fixed_events[cell.amplitude]


def fixed_system(template, amplitudes, events, *, transport=False, escape=np.inf):
    mobility = SimpleNamespace(scales=template.mobility.scales,
        amplitude_response=lambda amplitude, force, temperature: SimpleNamespace(velocity=0., heat=0.))
    system = FixedSystem(template.catalog, template.phonons, mobility,
        gammas=(0.,)*len(amplitudes), tau_kin=np.inf, tau_escape=escape,
        diffusion_over_length_squared=.2 if transport else 0., face_order=2,
        external_powers=(0.,)*len(amplitudes), rate_prefactor=template.rate_prefactor,
        reaction_method='projected', reaction_layout='resolved_panels')
    system.fixed_events = dict(zip(amplitudes, events))
    return system


def measures(system, state):
    a, p, n = system.unpack(state)
    qe, ee, full_e = [], [], []
    for i in range(system.cell_count):
        cell = ElectronicCell(system.catalog, a[i], 0.)
        qe.append(float(np.dot(4*cell.weights, p[i])))
        ee.append(float(np.dot(4*cell.weights*cell.energies, p[i])))
        full_e.append(float(cell.moments(p[i])[0]))
    phonon_counts = n @ system.phonons.capacities
    phonon_energy = n @ (system.phonons.capacities*system.phonons.energies)
    ledger, transport = system.ledgers(state)
    return dict(electron_count=float(sum(qe)), phonon_count=float(sum(phonon_counts)),
                excitation_energy=float(sum(ee)), phonon_energy=float(sum(phonon_energy)),
                total=float(sum(full_e)+sum(phonon_energy)), escape=float(sum(ledger[:, 0])),
                eph=float(sum(ledger[:, 2])), transport=float(transport),
                external=float(sum(ledger[:, 1])), heat=float(sum(ledger[:, 3])),
                amplitudes=a.tolist())


def event_geometry(event):
    if not len(event.omega):
        return 0.
    e = event.electron_energies
    ei = e[event.electron_lower_i]*event.electron_beta_lower_i+e[event.electron_upper_i]*event.electron_beta_upper_i
    ej = e[event.electron_lower_j]*event.electron_beta_lower_j+e[event.electron_upper_j]*event.electron_beta_upper_j
    ph = event.phonons.energies[event.phonon_lower]*event.beta_lower+event.phonons.energies[event.phonon_upper]*event.beta_upper
    residual = np.where(event.recombination, -ei-ej, ei-ej)+ph
    return float(np.max(abs(residual)/(abs(ei)+abs(ej)+abs(ph))))


def save_case(name, system, initial, duration, steps):
    started = time.perf_counter()
    instantaneous = abs(system.instantaneous_energy_residual(initial))
    geometry = max(event_geometry(e) for e in system.fixed_events.values())
    times, states, stats = integrate(system, initial, duration, steps)
    rows = [measures(system, state) for state in states]
    a0, p0, n0 = system.unpack(initial)
    energy_scale = max(1., abs(rows[0]['total']))
    energy_error = max(abs(r['total']+r['escape']-rows[0]['total'])/energy_scale for r in rows)
    number = [(r['electron_count']+2*r['phonon_count']) if name in ('recombination','creation')
              else r['electron_count'] for r in rows]
    count_error = max(abs(v-number[0]) for v in number)/max(1., abs(number[0]))
    population_change = max(float(np.sum(abs(system.unpack(s)[1]-p0)*(4*system.catalog.count_weights)))
                            +float(np.sum(abs(system.unpack(s)[2]-n0)*system.phonons.capacities)) for s in states)
    eph_error = max(abs(r['phonon_energy']-rows[0]['phonon_energy']+r['escape']-r['eph']) for r in rows)/energy_scale
    disabled_error = max(abs(r['external'])+abs(r['heat']) for r in rows)
    fields_unchanged = all(np.array_equal(system.unpack(s)[0], a0) for s in states)
    ph_change=rows[-1]['phonon_energy']-rows[0]['phonon_energy']
    q_change=rows[-1]['electron_count']-rows[0]['electron_count']
    direction_ok = (ph_change>0 if name in ('emission','recombination') else
                    ph_change<0 if name in ('absorption','creation') else
                    rows[-1]['transport']>0 if name=='transport' else rows[-1]['escape']>0)
    if name=='recombination':
        direction_ok = direction_ok and q_change<0
    if name=='creation':
        direction_ok = direction_ok and q_change>0
    target = OUT/f'isolated_{name}.npz'
    if target.exists():
        raise FileExistsError(target)
    np.savez_compressed(target, times=times, states=states, electron_count=system.catalog.count_nodes,
                        electron_weights=system.catalog.count_weights, phonon_energies=system.phonons.energies,
                        phonon_capacities=system.phonons.capacities)
    gates = dict(physical_populations=all(np.all((system.unpack(s)[1]>=0)&(system.unpack(s)[1]<=1))
                       and np.all(system.unpack(s)[2]>=0) for s in states),
                 energy=energy_error<=1e-7, count=count_error<=1e-12,
                 instantaneous=instantaneous<=1e-11, geometry=geometry<=1e-12,
                 actual_evolution=population_change>100*np.finfo(float).eps,
                 fields_fixed=fields_unchanged, sample_count=len(times)>=10,
                 shared_energy_ledger=eph_error<=1e-12, disabled_channels=disabled_error==0.,
                 expected_transfer_direction=bool(direction_ok))
    record = dict(case=name, status='PASS' if all(gates.values()) else 'FAIL', gates=gates,
                  energy_error=energy_error, count_error=count_error, instantaneous_energy_absolute=instantaneous,
                  event_energy_relative=geometry, eph_ledger_error=eph_error,
                  population_change_capacity_L1=population_change, initial=rows[0], final=rows[-1],
                  initial_exact_faces=dict(electrons_zero=int(np.count_nonzero(p0==0)),
                                          electrons_full=int(np.count_nonzero(p0==1)),
                                          phonons_zero=int(np.count_nonzero(n0==0))),
                  samples=rows, times=times.tolist(), integrator_stats=stats,
                  trajectory_sha256=sha(target), registration_sha256=sha(PLAN),
                  runtime_seconds=time.perf_counter()-started)
    write_new(OUT/f'isolated_{name}.json', record)
    print(json.dumps({k:record[k] for k in ('case','status','runtime_seconds','energy_error','count_error')}), flush=True)
    if not all(gates.values()):
        raise ValueError('isolated case failed: '+name)
    return record, states


def run():
    plan = json.loads(PLAN.read_text(encoding='utf-8'))
    if plan != registration():
        raise ValueError('registration/source contract changed')
    names = plan['isolated']['cases']+[f'escape_{n}' for n in plan['escape']['steps']]
    paths = [OUT/'isolated_results.json']
    paths += [OUT/f'isolated_{name}{suffix}' for name in names for suffix in ('.json','.npz')]
    if any(path.exists() for path in paths):
        raise FileExistsError('experiment outputs already exist; no reexecution or overwrite')
    started = time.perf_counter()
    output = dict(schema='pysnspd.stage2.isolated-results.v1', status='INCOMPLETE',
                  source_hashes=sources(), registration_sha256=sha(PLAN), cases=[],
                  python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                  limitations=plan['limitations'])
    try:
        template, _, _ = setup('one', heating=0., escape=np.inf)
        if template.electron_size != 630 or template.phonon_size != 1025:
            raise ValueError('selected candidate grid changed')
        cell = ElectronicCell(template.catalog, .5, 0.)
        parent = template.reaction_events(cell)
        band = ((template.catalog.count_nodes>1.1)&(template.catalog.count_nodes<1.5)).astype(float)
        for name in plan['isolated']['cases']:
            if name == 'transport':
                empty = DirectionalEvents(parent, np.zeros(len(parent.omega), bool), 'disabled')
                system = fixed_system(template, [.5, 1.], [empty, empty], transport=True)
                initial = system.pack([.5, 1.], [band, np.zeros_like(band)], np.zeros((2,1025)))
            else:
                rec = name in ('recombination','creation')
                forward = name in ('emission','recombination')
                event = DirectionalEvents(parent, parent.recombination==rec,
                                          'forward' if forward else 'backward')
                system = fixed_system(template, [.5], [event])
                p = np.zeros_like(band) if name=='creation' else band
                n = np.zeros(1025) if forward else np.full(1025,.05)
                initial = system.pack([.5], [p], [n])
                # Independent check: the chosen rate is precisely the original
                # activity on this subset, not a clipped signed net rate.
                lf, lb = parent.log_activities(p,n)
                mask = parent.recombination==rec
                expected = parent.coefficients[mask]*np.exp((lf if forward else lb)[mask])
                if not np.array_equal(event.rates(p,n), expected if forward else -expected):
                    raise ValueError('directional activity differs from physical event law')
            record, _ = save_case(name, system, initial, .02, 20)
            output['cases'].append({k:record[k] for k in ('case','status','energy_error','count_error','runtime_seconds')})
            del system
            if name != 'transport':
                del event
        empty = DirectionalEvents(parent, np.zeros(len(parent.omega), bool), 'disabled')
        system = fixed_system(template, [.5], [empty], escape=.7)
        p = cell.fermi_dirac(.2)
        n0 = system.bath_phonons+.02*np.exp(-system.phonons.energies)
        expected = system.bath_phonons+(n0-system.bath_phonons)*np.exp(-.2/.7)
        errors = []
        for steps in plan['escape']['steps']:
            initial = system.pack([.5], [p], [n0])
            # Escape is stationary in electron count, but an evolving phonon state.
            record, states = save_case(f'escape_{steps}', system, initial, .2, steps)
            errors.append(float(np.dot(system.phonons.capacities,abs(system.unpack(states[-1])[2][0]-expected))))
        ratios = [errors[i]/errors[i+1] for i in (0,1)]
        output['escape_order'] = dict(errors_capacity_L1=errors, reductions=ratios,
                                     status='PASS' if min(ratios)>=6. else 'FAIL')
        if min(ratios)<6.:
            raise ValueError('SSP escape smooth-order gate failed')
        output['status']='PASS_ISOLATED_BOUNDARIES_AND_SSP_ESCAPE_ONLY'
    except Exception as exc:
        output['status']='FAIL_OR_INCOMPLETE'
        output['error']=repr(exc)
        raise
    finally:
        output['runtime_seconds']=time.perf_counter()-started
        write_new(OUT/'isolated_results.json', output)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--register', action='store_true')
    args=parser.parse_args()
    if args.register:
        write_new(PLAN,registration())
        print('PREREGISTERED', PLAN)
    else:
        run()


if __name__=='__main__':
    main()
