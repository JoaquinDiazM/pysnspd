"""Fresh, bounded one/two-cell experiments against the immutable R2 catalogue."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import platform
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
import scipy
from scipy.special import expit
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pysnspd.experimental.energy_catalog import OccupationEnergyCatalog
from pysnspd.experimental.cell_validation import ElectronicCell, EnergyFacePair, rk4_trajectory

OUT = ROOT / "docs/implementation/stage1_closure"
CATALOG = ROOT / "docs/implementation/stage1_r2/catalogs/occupation_catalog.npz"
CRITERIA = OUT / "acceptance_criteria.json"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_csv(path, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def nonthermal(name, x):
    if name == "low_energy":
        return .2*np.exp(-x/.25)
    if name == "finite_energy_band":
        return .16*np.exp(-((x-.95)/.3)**2)
    raise ValueError(name)


def run_bgk(catalog):
    cell = ElectronicCell(catalog, .72, .2)
    tau, duration = .7, 3.5
    result = []; trajectories = []
    for profile in ("low_energy", "finite_energy_band"):
        initial = nonthermal(profile, catalog.count_nodes)
        target = cell.fermi_dirac(cell.equivalent_temperature(initial))
        exact = target+(initial-target)*np.exp(-duration/tau)
        for steps in (80, 160, 320):
            times, states = rk4_trajectory(lambda t, p: cell.bgk(p, tau)[0], initial, duration, steps)
            energy = np.asarray([cell.excitation_energy(p) for p in states])
            entropy = np.asarray([cell.entropy(p) for p in states])
            moments = [abs(4*np.dot(cell.weights*cell.energies, cell.bgk(p, tau)[0])) for p in states]
            result.append(dict(profile=profile, steps=steps,
                               energy_scaled_error=float(np.max(abs(energy-energy[0]))/max(1., abs(energy[0]))),
                               instantaneous_energy_moment_max=float(max(moments)),
                               final_population_error=float(np.max(abs(states[-1]-exact))),
                               entropy_decrease_max=float(max(0., -np.min(np.diff(entropy)))),
                               initial_count=cell.quasiparticle_count(initial), final_count=cell.quasiparticle_count(states[-1]),
                               minimum_population=float(states.min()), maximum_population=float(states.max())))
            if steps == 320:
                for t, p, u, s in zip(times, states, energy, entropy):
                    trajectories.append(dict(profile=profile, time=float(t), excitation_energy=float(u), entropy=float(s),
                                             distance_to_FD=float(np.max(abs(p-target))), temperature=cell.equivalent_temperature(p)))
    zero = np.zeros_like(cell.energies)
    thermal = cell.fermi_dirac(.25)
    stationary = dict(vacuum_rhs_max=float(np.max(abs(cell.bgk(zero, tau)[0]))),
                      thermal_rhs_max=float(np.max(abs(cell.bgk(thermal, tau)[0]))))
    return dict(amplitude=.72, gamma=.2, tau=tau, duration=duration, resolutions=result, stationarity=stationary), trajectories


def imposed_fields(t):
    return .72+.04*t, .1+.03*np.sin(np.pi*t/2), .04, .03*np.pi/2*np.cos(np.pi*t/2)


def run_moving(catalog):
    duration, power, tau, bath = 2., .04, .7, .08
    initial = np.zeros(len(catalog.count_nodes)+2)
    result = []; trajectory = []; largest_source_moment = 0.

    def rhs(t, state):
        nonlocal largest_source_moment
        a, g, adot, gdot = imposed_fields(t)
        cell = ElectronicCell(catalog, a, g)
        p = state[:-2]
        bgk, _ = cell.bgk(p, tau)
        heat = cell.heating(p, power, bath)
        _, fa, fg = cell.moments(p)
        largest_source_moment = max(largest_source_moment, abs(4*np.dot(cell.weights*cell.energies, heat)-power))
        return np.r_[bgk+heat, fa*adot+fg*gdot, power]

    for steps in (80, 160, 320):
        times, states = rk4_trajectory(rhs, initial, duration, steps)
        rows = []
        u0 = ElectronicCell(catalog, *imposed_fields(0)[:2]).moments(initial[:-2])[0]
        for t, state in zip(times, states):
            a, g = imposed_fields(t)[:2]
            cell = ElectronicCell(catalog, a, g)
            u, fa, fg = cell.moments(state[:-2])
            rows.append(dict(time=float(t), amplitude=a, gamma=g, temperature=cell.equivalent_temperature(state[:-2]),
                             total_energy=u, field_work=float(state[-2]), heat_work=float(state[-1]),
                             energy_residual=u-u0-state[-2]-state[-1],
                             incorrect_residual_without_field_work=u-u0-state[-1]))
        errors = np.asarray([row['energy_residual'] for row in rows])
        result.append(dict(steps=steps, energy_ledger_scaled_max=float(np.max(abs(errors))/max(1., abs(u0))),
                           energy_ledger_final_abs=float(abs(errors[-1])),
                           negative_control_final_abs=float(abs(rows[-1]['incorrect_residual_without_field_work'])),
                           minimum_population=float(states[:, :-2].min()), maximum_population=float(states[:, :-2].max())))
        if steps == 320:
            trajectory = rows
    cell = ElectronicCell(catalog, .72, .2)
    thermal = cell.fermi_dirac(.25)
    direction = cell.energies*thermal*(1-thermal)
    direction *= power/(4*np.dot(cell.weights*cell.energies, direction))
    thermal_error = float(np.max(abs(cell.heating(thermal, power, bath)-direction)))
    return dict(duration=duration, power=power, tau=tau, bath_temperature=bath, fields='Delta=.72+.04t; Gamma=.1+.03sin(pi*t/2)',
                resolutions=result, source_energy_moment_max=largest_source_moment, thermal_direction_abs_error=thermal_error), trajectory


def run_self_consistent(catalog):
    mobility, tau, bath, duration = .2, .7, .08, 6.
    if not np.isfinite(mobility) or mobility <= 0:
        raise ValueError('synthetic mobility must be positive')
    initial = np.r_[.6, .08*np.exp(-catalog.count_nodes/.25), 0.]
    result = []; trajectory = []; balance = 0.; minimum_heat = np.inf

    def rhs(t, state):
        nonlocal balance, minimum_heat
        cell = ElectronicCell(catalog, state[0], 0.)
        p = state[1:-1]
        _, force, _ = cell.moments(p)
        adot, power = cell.condensate_relaxation(p, mobility)
        dp = cell.bgk(p, tau)[0]+cell.heating(p, power, bath)
        balance = max(balance, abs(force*adot+4*np.dot(cell.weights*cell.energies, dp)))
        minimum_heat = min(minimum_heat, power)
        return np.r_[adot, dp, power]

    u0 = ElectronicCell(catalog, initial[0], 0.).moments(initial[1:-1])[0]
    for steps in (80, 160, 320):
        times, states = rk4_trajectory(rhs, initial, duration, steps)
        rows = []
        for t, state in zip(times, states):
            cell = ElectronicCell(catalog, state[0], 0.)
            p = state[1:-1]; u = cell.moments(p)[0]
            rows.append(dict(time=float(t), amplitude=float(state[0]), temperature=cell.equivalent_temperature(p),
                             total_energy=u, vacuum_energy=cell.catalog.vacuum.evaluate(state[0], 0.)[0],
                             excitation_energy=cell.excitation_energy(p), converted_heat=float(state[-1]), energy_residual=u-u0))
        error = np.asarray([row['energy_residual'] for row in rows])
        result.append(dict(steps=steps, energy_ledger_scaled_max=float(np.max(abs(error))/max(1., abs(u0))),
                           energy_ledger_final_abs=float(abs(error[-1])), final_amplitude=float(states[-1, 0]),
                           final_temperature=rows[-1]['temperature'], converted_heat=float(states[-1, -1]),
                           minimum_population=float(states[:, 1:-1].min()), maximum_population=float(states[:, 1:-1].max())))
        if steps == 320:
            trajectory = rows
    return dict(initial_amplitude=.6, gamma=0., mobility=mobility, mobility_status='constant synthetic; not calibrated KWT',
                tau=tau, duration=duration, bath_temperature=bath, resolutions=result,
                instantaneous_balance_max=balance, minimum_dissipation=minimum_heat), trajectory


def run_transport(catalog):
    left, right = ElectronicCell(catalog, .43, .02), ElectronicCell(catalog, 1.2, .2)
    duration = 5.; rows = []; trajectory = []; time_rows = []
    for number in (129, 257, 513):
        common_upper = min(left.energies[-1], right.energies[-1])
        upper = max(left.energies[-1], right.energies[-1])
        energies = np.linspace(1e-6, common_upper, number)
        if upper > common_upper:
            energies = np.r_[energies, upper]  # private tail: zero interface conductance
        pair = EnergyFacePair(left, right, energies, .2)
        initial = np.asarray([expit(-energies/.4), expit(-energies/.12)])
        same = np.tile(expit(-energies/.25), (2, 1))
        checkpoints = []
        initial_count = pair.quasiparticle_count(initial)
        initial_u = float(np.sum(pair.reconstructed_energy(initial)))
        max_native_bias = 0.; min_population = 1.; max_population = 0.
        for t in np.linspace(0, duration, 101):
            state = pair.exact(initial, t)
            energy = pair.reconstructed_energy(state)
            native = pair.native_energy(state)
            max_native_bias = max(max_native_bias, float(np.max(abs(energy-native)/np.maximum(abs(native), 1e-100))))
            min_population = min(min_population, float(state.min())); max_population = max(max_population, float(state.max()))
            checkpoints.append(dict(nodes=number, time=float(t), left_energy=float(energy[0]), right_energy=float(energy[1]),
                                    left_native_R2_energy=float(native[0]), right_native_R2_energy=float(native[1]),
                                    energy_residual=float(np.sum(energy)-initial_u),
                                    quasiparticle_count=pair.quasiparticle_count(state)))
        final = pair.exact(initial, duration)
        native_fd = np.asarray([left.excitation_energy(left.fermi_dirac(.25)), right.excitation_energy(right.fermi_dirac(.25))])
        fd_reconstructed = pair.native_energy(same)
        remap = float(np.max(abs(fd_reconstructed-native_fd)/native_fd))
        count_error = abs(pair.quasiparticle_count(final)-initial_count)/initial_count
        rows.append(dict(nodes=number, reconstructed_energy_scaled_max=float(max(abs(row['energy_residual']) for row in checkpoints)/max(1., abs(initial_u))),
                         native_R2_reconstruction_relative_max=max_native_bias, FD_native_remap_relative_max=remap,
                         quasiparticle_count_relative_drift=count_error, quasiparticle_count_initial=initial_count,
                         quasiparticle_count_final=pair.quasiparticle_count(final), common_distribution_rhs_max=float(np.max(abs(pair.rhs(same)))),
                         energy_transferred_left_to_right=float(pair.reconstructed_energy(final)[1]-pair.reconstructed_energy(initial)[1]),
                         minimum_population=min_population, maximum_population=max_population,
                         active_energy_nodes=int(np.count_nonzero(pair.active)),
                         total_representation_nodes=len(energies), common_upper_energy=float(common_upper),
                         private_tail_native_energy_initial=[float(4*np.dot(cell.weights[cell.energies>common_upper]*cell.energies[cell.energies>common_upper],
                                                                                 native[cell.energies>common_upper]))
                                                             for cell, native in zip((left, right), pair.native_populations(initial))],
                         private_tail_population_change_max=float(np.max(abs(final[:, energies>common_upper]-initial[:, energies>common_upper]))),
                         common_upper_boundary_population_change_max=float(np.max(abs(final[:, energies==common_upper]-initial[:, energies==common_upper]))),
                         discarded_count_below_lower_energy=[float(1e-6*cell.catalog.count_nodes[0]/cell.energies[0]) for cell in (left, right)],
                         discarded_lower_energy_full_occupation_bound=[float(2e-12*cell.catalog.count_nodes[0]/cell.energies[0]) for cell in (left, right)],
                         represented_count_upper=[float(catalog.count_nodes[-1])]*2,
                         nominal_R2_count_support=float(np.sum(catalog.count_weights)),
                         unrepresented_count_interval_to_R2_cutoff=float(np.sum(catalog.count_weights)-catalog.count_nodes[-1])))
        if number == 513:
            trajectory = checkpoints
            for steps in (80, 160, 320):
                _, states = rk4_trajectory(lambda t, state: pair.rhs(state), initial, duration, steps)
                time_rows.append(dict(steps=steps, final_state_scaled_error=float(np.max(abs(states[-1]-final))),
                                      minimum_population=float(states.min()), maximum_population=float(states.max())))
    native_same_temperature = np.asarray([left.fermi_dirac(.25), right.fermi_dirac(.25)])
    wrong_difference = native_same_temperature[0]-native_same_temperature[1]
    negative = dict(equal_index_occupation_difference_max=float(np.max(abs(wrong_difference))),
                    incorrect_equal_index_left_power=float(4*.2*np.dot(left.weights*left.energies, wrong_difference)),
                    interpretation='Deliberately wrong equal-count-index exchange at equal physical temperature; not applied to any trajectory.')
    return dict(left_fields=[.43, .02], right_fields=[1.2, .2], initial_temperatures=[.4, .12], duration=duration,
                diffusion_over_length_squared=.2, lower_shared_energy=1e-6,
                method='linear E(x) reconstruction; exact common-domain hat moments; C=H/E; harmonic causal D_L; exponential two-state exchange; private tail constant and disconnected',
                resolutions=rows, time_resolutions=time_rows, negative_control=negative), trajectory


def figures(bgk, moving, coupled, transport, out):
    plt.rcParams.update({'font.size':12, 'axes.titlesize':13, 'axes.labelsize':12, 'legend.fontsize':11})
    fig, axes = plt.subplots(2, 1, figsize=(10, 6.8), constrained_layout=True)
    t = [row['time'] for row in coupled]
    axes[0].plot(t, [row['amplitude'] for row in coupled], label=r'$|\Delta|/\Delta_0$')
    axes[0].plot(t, [row['temperature'] for row in coupled], label=r'$k_BT_E/\Delta_0$')
    axes[0].set(title='Una celda: recuperación del condensado y poblaciones completas', ylabel='Valor normalizado')
    axes[0].legend(); axes[0].grid(alpha=.25)
    axes[1].plot(t, [row['excitation_energy']-coupled[0]['excitation_energy'] for row in coupled], label='Aumento de energía de excitaciones')
    axes[1].plot(t, [coupled[0]['vacuum_energy']-row['vacuum_energy'] for row in coupled], '--', label='Energía liberada por el fondo')
    axes[1].set(xlabel=r'Tiempo $t/\tau_{ref}$ (escala sintética)', ylabel=r'Energía / $(N_0\Delta_0^2)$')
    axes[1].legend(); axes[1].grid(alpha=.25)
    for suffix in ('png', 'pdf'): fig.savefig(out/f'cells_recovery.{suffix}', dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(2, 1, figsize=(10, 6.8), constrained_layout=True)
    tt = [row['time'] for row in transport]
    axes[0].plot(tt, [row['left_energy'] for row in transport], label='Celda izquierda, inicialmente caliente')
    axes[0].plot(tt, [row['right_energy'] for row in transport], label='Celda derecha, inicialmente fría')
    axes[0].plot(tt, [row['left_energy']+row['right_energy'] for row in transport], '--', label='Total reconstruido')
    axes[0].set(title='Dos gaps distintos: intercambio a energía compartida', ylabel=r'Energía / $(N_0\Delta_0^2)$')
    axes[0].legend(); axes[0].grid(alpha=.25)
    axes[1].plot([row['time'] for row in moving], [row['energy_residual'] for row in moving], label='Balance con trabajo de campos')
    axes[1].plot([row['time'] for row in moving], [row['incorrect_residual_without_field_work'] for row in moving], label='Control incorrecto: omitir ese trabajo')
    axes[1].set(xlabel=r'Tiempo $t/\tau_{ref}$ (escala sintética)', ylabel='Residuo energético normalizado', title='Espectro móvil: la contabilidad importa')
    axes[1].legend(); axes[1].grid(alpha=.25)
    for suffix in ('png', 'pdf'): fig.savefig(out/f'cells_transport_and_work.{suffix}', dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUT)
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    criteria = json.loads(CRITERIA.read_text())
    immutable = criteria['immutable_inputs']
    for path_key, hash_key in (('catalog_path','catalog_sha256'), ('query_source_path','query_source_sha256')):
        if sha(ROOT/immutable[path_key]) != immutable[hash_key]:
            raise ValueError('immutable R2 input differs from the frozen acceptance contract')
    catalog = OccupationEnergyCatalog.load(CATALOG)
    bgk, bgk_rows = run_bgk(catalog); print('BGK complete', flush=True)
    moving, moving_rows = run_moving(catalog); print('Moving fields complete', flush=True)
    coupled, coupled_rows = run_self_consistent(catalog); print('Self-consistent cell complete', flush=True)
    transport, transport_rows = run_transport(catalog); print('Transport complete', flush=True)
    result = dict(schema='pysnspd.stage1_closure.cells.v1', host=platform.node(), python=platform.python_version(),
                  numpy=np.__version__, scipy=scipy.__version__, runtime_seconds=time.perf_counter()-started,
                  criteria_sha256=sha(CRITERIA), catalog_sha256=sha(CATALOG),
                  query_code_sha256=sha(ROOT/immutable['query_source_path']),
                  cell_code_sha256=sha(ROOT/'pysnspd/experimental/cell_validation.py'), runner_sha256=sha(__file__),
                  units=dict(energy='N0 Delta0^2', temperature='kB T / Delta0', fields='Delta0',
                             time='synthetic tau_ref; no absolute material calibration'),
                  bgk=bgk, moving_spectrum=moving, self_consistent_cell=coupled, two_cell_transport=transport,
                  scope='Reduced electronic tests. No electron-phonon populations, spatial core, KWT calibration, circuit or full transient.',
                  status='PENDING_INDEPENDENT_ADJUDICATION', production_connected=False)
    (args.output/'cells_results.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    for name, rows in (('bgk', bgk_rows), ('moving', moving_rows), ('recovery', coupled_rows), ('transport', transport_rows)):
        write_csv(args.output/f'cells_{name}.csv', rows)
    figures(bgk_rows, moving_rows, coupled_rows, transport_rows, args.output)
    print(json.dumps(dict(runtime_seconds=result['runtime_seconds'], source_sha256=result['cell_code_sha256'],
                         bgk=bgk['resolutions'], moving=moving['resolutions'], coupled=coupled['resolutions'],
                         transport=transport['resolutions'], transport_time=transport['time_resolutions']), indent=2), flush=True)


if __name__ == '__main__':
    main()
