"""One fresh coupled-cell trajectory, with explicit parameters and failure output."""
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog, K_B_J_K
from pysnspd.experimental.cell_closures import CellScales, KWTMobility, DebyePhonons, bose_occupation
from pysnspd.experimental.cell_validation import ElectronicCell, rk4_trajectory
from pysnspd.experimental.kinetic_events import PhononGrid
from pysnspd.experimental.coupled_cells import CoupledCellSystem
from pysnspd.experimental.refined_cells import refined_count_catalog

CRITERIA = ROOT/'docs/implementation/stage2/acceptance_criteria.json'
CATALOG = ROOT/'docs/implementation/stage1_r2/catalogs/occupation_catalog.npz'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def jsonable(value):
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {k: jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    return value


def setup(case, phonon_nodes=513, infrared=.005, face_order=2, escape=15., heating=.01,
          scenario='driven', reaction_order=2, reaction_method='projected',
          electron_refinement=1, reaction_layout='resolved_panels', reaction_outer_order=2,
          reaction_max_panel=.5, reaction_max_energy_panel=.125):
    catalog = OccupationEnergyCatalog.load(CATALOG)
    if electron_refinement:
        catalog = refined_count_catalog(catalog, refinement=electron_refinement)
    vacuum = catalog.vacuum
    scales = CellScales(vacuum.delta0_J, vacuum.N0_per_J_m3, 8.65,
                        .12*vacuum.delta0_J/K_B_J_K, 1.)
    debye = DebyePhonons(scales, cutoff_energy_bar=4.,
                        atom_density_m3=10*scales.N0_per_J_m3*scales.delta0_J,
                        lambda_eph=.1, synthetic_label='Stage2 synthetic Debye, 10 ions/(N0 Delta0)',
                        infrared_cutoff_bar=infrared)
    quadrature = debye.quadrature(phonon_nodes)
    phonons = PhononGrid(quadrature.energies, quadrature.capacities, debye.alpha2F,
                        debye.coupling_support, 'Stage2 Debye with declared infrared cutoff')
    amplitudes = [.6] if case == 'one' else [.55, .95]
    cells = [ElectronicCell(catalog, a, 0.) for a in amplitudes]
    if case == 'one':
        x = catalog.count_nodes
        p = [cells[0].fermi_dirac(.12)+.2*np.exp(-((x-.8)/.3)**2)]
    else:
        p = [c.fermi_dirac(t) for c, t in zip(cells, (.35, .12))]
    n = [bose_occupation(phonons.energies, .12) for _ in cells]
    # B.39: bubble proportional to alpha^2=alpha^2F/g, normalized to energy.
    shape = debye.alpha2F(phonons.energies)/quadrature.dos
    n[0] = n[0]+.15*shape/np.dot(phonons.capacities*phonons.energies, shape)
    gammas = tuple(0. for _ in cells)
    if scenario == 'equilibrium':
        if heating != 0:
            raise ValueError('equilibrium scenario requires explicitly disabled external heating')
        gammas = (0.,) if case == 'one' else (0., .1)
        amplitudes = [brentq(lambda a: ElectronicCell(catalog, a, g).moments(
            ElectronicCell(catalog, a, g).fermi_dirac(.12))[1], .08, 1.5,
            xtol=1e-14) for g in gammas]
        cells = [ElectronicCell(catalog, a, g) for a, g in zip(amplitudes, gammas)]
        p = [c.fermi_dirac(.12) for c in cells]
        n = [bose_occupation(phonons.energies, .12) for _ in cells]
    elif scenario == 'phonon_vacuum':
        n = [np.zeros_like(phonons.energies) for _ in cells]
    elif scenario == 'sparse_electrons':
        p = [np.where((catalog.count_nodes > .65) & (catalog.count_nodes < .85), 1., 0.)
             for _ in cells]
    elif scenario != 'driven':
        raise ValueError('unknown scenario')
    powers = tuple([heating]+[0.]*(len(cells)-1))
    system = CoupledCellSystem(catalog, phonons, KWTMobility(scales),
                              gammas=gammas, tau_kin=.7,
                              tau_escape=escape, diffusion_over_length_squared=.2,
                              face_order=face_order, external_powers=powers,
                              rate_prefactor=8*np.pi*scales.quantum_rate_scale,
                              reaction_quadrature_order=reaction_order, reaction_method=reaction_method,
                              reaction_layout=reaction_layout, reaction_outer_order=reaction_outer_order,
                              reaction_max_panel=reaction_max_panel,
                              reaction_max_energy_panel=reaction_max_energy_panel)
    initial = system.pack(amplitudes, p, n)
    return system, initial, debye


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('one', 'two'), required=True)
    parser.add_argument('--steps', type=int, default=100)
    parser.add_argument('--duration', type=float, default=2.)
    parser.add_argument('--phonon-nodes', type=int, default=513)
    parser.add_argument('--infrared', type=float, default=.005)
    parser.add_argument('--face-order', type=int, default=2)
    parser.add_argument('--reaction-order', type=int, default=2)
    parser.add_argument('--reaction-layout', choices=('global','native_intervals','energy_panels','resolved_panels'), default='resolved_panels')
    parser.add_argument('--reaction-max-energy-panel', type=float, default=.125)
    parser.add_argument('--reaction-outer-order', type=int, default=2)
    parser.add_argument('--reaction-max-panel', type=float, default=.5)
    parser.add_argument('--electron-refinement', type=int, default=1,
                        help='0 selects archived R2; positive levels select complementary meshes')
    parser.add_argument('--reaction-method', choices=('projected','native_pairs'), default='projected')
    parser.add_argument('--escape', type=float, default=15.)
    parser.add_argument('--heating', type=float, default=.01)
    parser.add_argument('--method', choices=('rk4', 'dop853'), default='rk4')
    parser.add_argument('--scenario', choices=('driven', 'equilibrium', 'phonon_vacuum', 'sparse_electrons'), default='driven')
    parser.add_argument('--rtol', type=float, default=1e-9)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not CRITERIA.is_file():
        raise RuntimeError('Freeze acceptance criteria before running any new diagnostic.')
    started = time.perf_counter()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    sources = {p.relative_to(ROOT).as_posix(): sha(p) for p in
               [*sorted((ROOT/'pysnspd/experimental').glob('*.py')), Path(__file__)]}
    record = dict(schema='pysnspd.stage2.coupled_run.v1', status='INCOMPLETE',
                  parameters={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
                  source_hashes=sources, criteria_sha256=sha(CRITERIA),
                  catalog_sha256=sha(CATALOG), host=platform.node(),
                  python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)
    try:
        system, initial, material = setup(args.case, args.phonon_nodes, args.infrared,
                                           args.face_order, args.escape, args.heating, args.scenario,
                                           args.reaction_order, args.reaction_method, args.electron_refinement,
                                           args.reaction_layout, args.reaction_outer_order, args.reaction_max_panel,
                                           args.reaction_max_energy_panel)
        record['occupation_mesh'] = dict(electron_states=system.electron_size,
            evaluation='causal direct complementary' if args.electron_refinement else 'archived R2',
            refinement=args.electron_refinement)
        record['material'] = material.metadata()
        # Persist the actual initial populations and source contract BEFORE the
        # first numerical evaluation; an interrupted run remains reviewable.
        initial_path = args.output.with_name(args.output.stem+'_initial.npz')
        np.savez_compressed(initial_path, state=initial,
                            electron_count=system.catalog.count_nodes,
                            electron_weights=system.catalog.count_weights,
                            phonon_energies=system.phonons.energies,
                            phonon_capacities=system.phonons.capacities)
        record['initial_conditions_sha256'] = sha(initial_path)
        record['initial'] = system.snapshot(initial)
        initial_a, initial_p, _ = system.unpack(initial)
        record['initial_excitation_energy'] = float(sum(
            ElectronicCell(system.catalog, a, g).excitation_energy(p)
            for a, g, p in zip(initial_a, system.gammas, initial_p)))
        args.output.write_text(json.dumps(jsonable(record), indent=2, allow_nan=False)+'\n', encoding='utf-8')
        record['initial_instantaneous_residual'] = system.instantaneous_energy_residual(initial)
        if args.method == 'rk4':
            times, states = rk4_trajectory(system.rhs, initial, args.duration, args.steps)
        else:
            times = np.linspace(0, args.duration, args.steps+1)
            solution = solve_ivp(system.rhs, (0, args.duration), initial, method='DOP853',
                                 t_eval=times, rtol=args.rtol, atol=1e-13,
                                 max_step=args.duration/args.steps)
            if not solution.success:
                raise RuntimeError(solution.message)
            states = solution.y.T
        snapshots = [system.snapshot(state) for state in states]
        ledger = np.array([s['conserved'] for s in snapshots])
        scale = max(1., abs(snapshots[0]['total']), float(np.sum(snapshots[-1]['input'])))
        checkpoint_indices = np.unique(np.linspace(0, len(states)-1, 11, dtype=int))
        instantaneous = [system.instantaneous_energy_residual(states[i]) for i in checkpoint_indices]
        record.update(status='COMPLETED_NOT_YET_ADJUDICATED',
                      rhs_calls=system.rhs_calls,
                      energy_ledger_scaled_max=float(np.max(abs(ledger-ledger[0]))/scale),
                      instantaneous_residual_max=float(max(abs(np.asarray(instantaneous)))),
                      minimum_electron=min(s['minimum_electron'] for s in snapshots),
                      maximum_electron=max(s['maximum_electron'] for s in snapshots),
                      minimum_phonon=min(s['minimum_phonon'] for s in snapshots),
                      initial=snapshots[0], final=snapshots[-1],
                      times=times, snapshots=snapshots)
        np.savez_compressed(args.output.with_suffix('.npz'), times=times, states=states,
                            electron_count=system.catalog.count_nodes,
                            electron_weights=system.catalog.count_weights,
                            phonon_energies=system.phonons.energies,
                            phonon_capacities=system.phonons.capacities)
        record['trajectory_sha256'] = sha(args.output.with_suffix('.npz'))
    except Exception as exc:
        record.update(status='FAILED', exception=type(exc).__name__, reason=str(exc),
                      traceback=traceback.format_exc())
    record['runtime_seconds'] = time.perf_counter()-started
    args.output.write_text(json.dumps(jsonable(record), indent=2, allow_nan=False)+'\n', encoding='utf-8')
    summary = {k: record[k] for k in ('status', 'runtime_seconds')}
    for key in ('reason', 'energy_ledger_scaled_max', 'instantaneous_residual_max', 'rhs_calls', 'final'):
        if key in record:
            summary[key] = record[key]
    print(json.dumps(jsonable(summary), indent=2))
    if record['status'] == 'FAILED':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
