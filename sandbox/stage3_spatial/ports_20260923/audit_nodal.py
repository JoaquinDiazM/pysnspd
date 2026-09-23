"""Audit archived nodal 3A data; never import/evaluate a physical model.

Only hashes, saved arrays, finite-difference arithmetic, Fourier sums and
recorded 2x2 matrices are evaluated. Outputs are new, never overwritten.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
REG = ROOT/'docs/implementation/stage3/nodal_20260923/registration.json'
DEST = ROOT/'docs/implementation/stage3/ports_20260923'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def carray(value):
    return np.asarray(value['real'])+1j*np.asarray(value['imag'])


def audit(raw):
    reg, manifest, summary = read(REG), read(raw/'manifest.json'), read(raw/'summary.json')
    errors, records, rows, probes, fd_tasks, no_fd_tasks = [], [], [], [], [], []

    def check(value, label):
        if not value:
            errors.append(label)

    def same(a, b, label, atol=1e-13):
        check(np.allclose(a, b, rtol=1e-12, atol=atol), label)

    expected = {f"{c['id']}_n{n}": (c, n) for c in reg['cases'] for n in reg['cell_counts']}
    sources = manifest['source_hashes']
    source_checks = {p: dict(recorded_sha256=h, current_sha256=sha(ROOT/p)) for p, h in sources.items()}
    check(len(sources) == 8, 'eight source contract')
    for p, v in source_checks.items():
        check(v['recorded_sha256'] == v['current_sha256'], 'source '+p)
    check(sources == summary['source_hashes'], 'summary sources')
    check(set(manifest['tasks']) == set(expected) and len(manifest['tasks']) == 18, '18 registered tasks')
    check(manifest['fd_cells'] == reg['fd_cells'] == [16], 'FD registration')
    check(summary['completed_cases'] == 18, 'summary completed count')
    populations = {}
    for task in manifest['tasks']:
        path = raw/(task+'.json')
        r, receipt = read(path), read(path.with_suffix('.receipt.json'))
        case, cells = expected[task]
        check(receipt['json_sha256'] == sha(path), task+' result receipt')
        check(receipt['arrays_sha256'] == r['arrays_sha256'] == sha(path.with_suffix('.npz')), task+' arrays receipt')
        initial_path = raw/(task+'_initial.npz')
        check(receipt['initial_sha256'] == r['initial_sha256'] == sha(initial_path), task+' initial receipt')
        check(receipt['source_hashes'] == r['source_hashes'] == sources, task+' source contract')
        check(r['case'] == case and r['cells'] == cells, task+' registered parameters')
        with np.load(path.with_suffix('.npz'), allow_pickle=False) as arr, np.load(initial_path, allow_pickle=False) as init:
            check(all(np.isfinite(arr[k]).all() for k in arr.files), task+' finite arrays')
            check(all(np.isfinite(init[k]).all() for k in init.files), task+' finite initial arrays')
            p, z, x = arr['p_nodes'], arr['delta'], arr['count_nodes']
            check(p.shape == (cells, 630) and z.shape == (cells,) and x.shape == (630,), task+' shapes')
            check(np.all((p >= 0) & (p <= 1)) and np.all(np.diff(x) > 0), task+' population/count domain')
            check(np.all(p == p[0]), task+' globally fixed population')
            for key in ('delta', 'p_nodes', 'count_nodes'):
                check(np.array_equal(arr[key], init[key]), task+' initial/final '+key)
            if case['occupation'] == 'vacuum':
                check(not np.any(p), task+' vacuum population')
            else:
                if case['occupation'] in populations:
                    check(np.array_equal(p[0], populations[case['occupation']]), task+' same preparation across meshes/cases')
                populations[case['occupation']] = p[0].copy()
            coords = np.arange(cells)*r['h_bar']
            amp = case['amplitude_base']+case['amplitude_modulation']*np.cos(2*np.pi*coords/r['length_bar'])
            phase = case['q_bare_ell0']*coords+case['phase_modulation']*np.sin(2*np.pi*coords/r['length_bar'])
            same(z, amp*np.exp(1j*phase), task+' prescribed state')
            obs = r['observable_arrays']
            same(carray(obs['force_density']), (arr['gradient'][:, 0]+1j*arr['gradient'][:, 1])/r['h_bar'], task+' stored force density')
            same(obs['current_A'], arr['current_A'], task+' stored current', 1e-20)
            same(obs['amplitude'], abs(z), task+' stored amplitude')
        fd = r['derivative_checks']
        required = cells in reg['fd_cells']
        force_error = current_error = budget_fraction = None
        if fd['performed']:
            fd_tasks.append(task)
            check(required or task == 'weak_phase_thermal_n8', task+' FD scope')
            check(fd['status'] == 'PASS' and r['status'] == 'PASS_STATIC_IDENTITIES_AND_SAMPLED_SYMBOL', task+' FD status')
            for name in ('force', 'current'):
                v = fd[name]
                error = abs(v['predicted']-v['fd_estimates'][-1])
                tolerance = reg['thresholds']['derivative_absolute']+reg['thresholds']['derivative_relative']*abs(v['fd_estimates'][-1])
                uncertainty = abs(v['fd_estimates'][-1]-v['fd_estimates'][-2])
                same([error, tolerance, uncertainty], [v['absolute_error'], v['tolerance'], v['fd_uncertainty']], task+' FD arithmetic '+name, 1e-18)
                check(v['passed'] and v['status'] == 'PASS' and error <= tolerance and uncertainty <= reg['thresholds']['reference_budget_fraction']*tolerance, task+' FD decision '+name)
            force_error, current_error = fd['force']['absolute_error'], fd['current']['absolute_error']
            budget_fraction = max(fd[n]['absolute_error']/fd[n]['tolerance'] for n in ('force', 'current'))
        else:
            no_fd_tasks.append(task)
            check(not required and fd['status'] == 'NOT_REPEATED_UNIT_AND_N16_REFERENCES' and r['status'] == 'PASS_STATIC_STRUCTURE_WITHOUT_REPEATED_FD', task+' omitted FD scope')
        check(max(r['global_phase_checks'].values()) <= reg['thresholds']['global_phase_scaled_error'], task+' gauge')
        check(max(r['uniform_reference_errors'].values()) <= reg['thresholds'].get('uniform_oracle_scaled_error', 1e-4), task+' uniform oracle')
        symbols = r['principal_symbols']
        check(len(symbols) == cells, task+' symbol coverage')
        for i, symbol in enumerate(symbols):
            check(symbol['stable'] and min(symbol['eigenvalues']) > symbol['uncertainty'], task+f' symbol margin {i}')
            same(np.linalg.eigvalsh(symbol['matrix']), symbol['eigenvalues'], task+f' eigenvalues {i}')
        spec = r['spectral_sign_comparison']
        if cells == 32 and case['occupation'] != 'vacuum':
            selected = sorted(set((int(np.argmin([s['eigenvalues'][0] for s in symbols])), int(np.argmax(obs['gamma'])), int(np.argmin(obs['amplitude'])))))
            check([v['link'] for v in spec['probes']] == selected and spec['status'] == 'PASS_SAMPLED_SPECTRAL_SIGN', task+' spectral selection/status')
        elif case['occupation'] == 'vacuum':
            check(spec['status'] == 'VACUUM_NO_OCCUPIED_SPECTRAL_QUADRATURE' and not spec['probes'], task+' vacuum spectral scope')
        else:
            check(spec['status'] == 'NOT_RUN_PILOT_OR_COARSER_GRID' and not spec['probes'], task+' coarse spectral scope')
        for probe in spec['probes']:
            coarse, fine = np.asarray(probe['coarse_matrix']), np.asarray(probe['fine_matrix'])
            same(np.linalg.eigvalsh(coarse), probe['coarse_eigenvalues'], task+' spectral coarse eigenvalues')
            same(np.linalg.eigvalsh(fine), probe['fine_eigenvalues'], task+' spectral fine eigenvalues')
            matrix_shift = float(np.linalg.norm(fine-coarse, 2))
            same(matrix_shift, probe['matrix_operator_norm_difference'], task+' matrix shift', 1e-18)
            moments = probe['electronic_moments']
            delta = abs(np.asarray(moments['fine'])-np.asarray(moments['coarse']))
            same(delta, moments['absolute_differences'], task+' moment differences', 1e-18)
            same(delta/np.maximum(1, abs(np.asarray(moments['fine']))), moments['scaled_differences'], task+' moment scaling', 1e-18)
            margin = min(probe['coarse_min'], probe['fine_min'])-probe['coarse_uncertainty']-probe['fine_uncertainty']-abs(probe['fine_min']-probe['coarse_min'])
            same(margin, probe['remaining_positive_margin'], task+' spectral margin')
            check(margin > 0, task+' spectral sign')
            probes.append(dict(task=task, **probe))
        same((r['energy_bar']-r['uniform_reference']['energy_bar'])/r['length_bar'], r['excess_energy_density_bar'], task+' excess energy')
        # Reproduce the observable extraction, not the physical fields.
        # The archived arrays order the same five registered modes as -2..2;
        # the permutation relative to prose does not change their L2 norm.
        induced_force = carray(obs['force_density'])*np.exp(-1j*case['q_bare_ell0']*coords)-r['uniform_reference']['force_radial_density_bar']
        current_scale = 2*1.602176634e-19/1.054571817e-34*r['energy_scale_J']
        induced_current = np.asarray(obs['current_A'])/current_scale-r['uniform_reference']['link_current_bar']
        fm = [np.mean(induced_force*np.exp(-2j*np.pi*k*coords/r['length_bar'])) for k in (-2, -1, 0, 1, 2)]
        jm = [np.mean(induced_current*np.exp(-2j*np.pi*k*(np.arange(cells)+.5)/cells)) for k in (-2, -1, 0, 1, 2)]
        same(fm, carray(r['induced_force_modes']), task+' Fourier force extraction')
        same(jm, carray(r['induced_current_modes']), task+' Fourier current extraction')
        rows.append(dict(task=task, status=r['status'], performed_fd=fd['performed'], force_error=force_error, current_error=current_error,
            max_fd_budget_fraction=budget_fraction, gauge_error=max(r['global_phase_checks'].values()), uniform_oracle_error=max(r['uniform_reference_errors'].values()),
            min_eigenvalue=min(min(s['eigenvalues']) for s in symbols), symbol_uncertainty=max(s['uncertainty'] for s in symbols), spectral_status=spec['status']))
        records.append(r)
    check(fd_tasks == summary['cases_with_directional_fd'] and no_fd_tasks == summary['cases_without_repeated_fd'], 'summary FD scopes')
    check(len(fd_tasks) == 7 and len(no_fd_tasks) == 11, 'FD scope counts')
    metrics, gradient = [], None
    check(len(summary['grid_comparisons']) == 6, 'six grid comparisons')
    for row in summary['grid_comparisons']:
        rr = sorted((r for r in records if r['case']['id'] == row['case']), key=lambda r: r['cells'])
        if row['case'].startswith('uniform'):
            check(row['status'] == 'STATIC_IDENTITIES_ONLY', row['case']+' uniform scope')
            same([r['energy_bar']/r['length_bar'] for r in rr], row['energy_density'], row['case']+' uniform energies')
            continue
        for key, metric in row['metrics'].items():
            if key == 'isolated_gradient':
                gradient = []
                for r in rr:
                    v = r['analytic_gradient']
                    error = float(np.linalg.norm(np.asarray(v['measured'])-np.asarray(v['expected']))/np.linalg.norm(v['expected']))
                    same(error, v['relative_rms'], row['case']+' analytic RMS')
                    gradient.append(error)
                same(gradient, metric['errors'], 'analytic gradient comparison')
                check(gradient[2] <= reg['thresholds']['gradient_finest_relative_error'] and gradient[2] < gradient[1] < gradient[0] and metric['passed'], 'analytic gradient verdict')
                continue
            arrays = [np.atleast_1d(carray(r[key]) if isinstance(r[key], dict) else r[key]) for r in rr]
            d1, d2, norm = float(np.linalg.norm(arrays[1]-arrays[0])), float(np.linalg.norm(arrays[2]-arrays[1])), float(np.linalg.norm(arrays[2]))
            same([d1, d2, norm], [metric['coarse_medium_difference'], metric['medium_fine_difference'], metric['fine_norm']], row['case']+' '+key+' differences', 1e-17)
            denominator = reg['richardson_denominator_by_observable'][key]
            check(metric['richardson_denominator'] == denominator, row['case']+' '+key+' registered order')
            floor, target = reg['thresholds']['spatial_response_absolute_floor'], reg['thresholds']['nonlinear_mesh_response_target_relative']
            analytic_zero = key == 'induced_current_modes' and rr[0]['case']['q_bare_ell0'] == 0 and rr[0]['case']['phase_modulation'] == 0
            if analytic_zero:
                status = 'PASS_ANALYTIC_ZERO_ABSOLUTE' if max(float(np.linalg.norm(a)) for a in arrays) <= floor else 'FAIL_ANALYTIC_ZERO_ABSOLUTE'
                check(metric['relative_estimate'] is None and not metric['relative_certificate'], row['case']+' analytic zero scope')
            elif norm <= floor or min(d1, d2) <= floor:
                status = 'ROUND_OFF_LIMITED'
                check(metric['passed'] is None and metric['relative_estimate'] is None and not metric['relative_certificate'], row['case']+' unresolved scope')
            elif d1 <= d2:
                status = 'INCONCLUSIVE_SPATIAL_ESTIMATE'
            else:
                estimate = max(d2/denominator, d2/(d1/d2-1))
                same(estimate, metric['richardson_conservative_estimate'], row['case']+' '+key+' estimate', 1e-18)
                same(estimate/norm, metric['relative_estimate'], row['case']+' '+key+' relative', 1e-16)
                status = 'PASS_RELATIVE_ESTIMATE' if estimate <= target*norm else 'FAIL_PRECISION_TARGET'
                check(metric['relative_certificate'] == (status == 'PASS_RELATIVE_ESTIMATE'), row['case']+' relative scope')
            check(metric['status'] == status, row['case']+' '+key+' verdict')
            metrics.append(dict(case=row['case'], observable=key, **metric))
    negative = read(raw/'negative_control.json')
    check(negative == summary['negative_control'], 'negative provenance')
    same(max(abs(np.asarray(negative['measured_eigenvalues'])-np.asarray(negative['reference_eigenvalues']))), negative['error'], 'negative error', 1e-18)
    check(negative['status'] == 'PASS_EXPECTED_REJECTION' and not negative['admitted_for_evolution'] and min(negative['measured_eigenvalues']) < 0 and negative['error'] < 1e-8, 'negative rejection')
    inventory = {p.name: dict(sha256=sha(p), bytes=p.stat().st_size) for p in sorted(raw.iterdir()) if p.is_file()}
    check(len(inventory) == 76, '76 raw files')
    failed = [m for m in metrics if m['status'].startswith('FAIL') or m['status'] == 'INCONCLUSIVE_SPATIAL_ESTIMATE']
    limited = [m for m in metrics if m['status'] == 'ROUND_OFF_LIMITED']
    passed = [m for m in metrics if m['status'] == 'PASS_RELATIVE_ESTIMATE']
    check(summary['spatial_precision_fully_resolved'] == (not failed and not limited), 'summary resolution scope')
    check(not summary['stage3_closed'] and not summary['production_promotion'] and summary['new_trajectories'] == 0 and not summary['circuit_and_physical_boundaries_implemented'], 'summary physical scope')
    fds = [r for r in rows if r['performed_fd']]
    eligible = not errors and not failed and bool(gradient)
    return dict(schema='pysnspd.stage3.nodal_saved_audit.v1', status='PASS_SAVED_DATA_AUDIT' if not errors else 'FAIL_SAVED_DATA_AUDIT',
        scope='Only archived JSON/NPZ arithmetic and hashes; no model call, spectrum solve or trajectory.', script_sha256=sha(__file__), registration_sha256=sha(REG),
        raw_directory=str(raw.relative_to(ROOT)), raw_inventory=inventory, source_checks=source_checks, integrity_errors=errors,
        completed_results=len(rows), sampled_nodes=sum(r['cells'] for r in records), cases_with_fd=fd_tasks, cases_without_repeated_fd=no_fd_tasks,
        recorded_runtime_seconds=summary['runtime_seconds'], recorded_batch_status=summary['status'], case_checks=rows,
        maxima=dict(force_fd_error=max(r['force_error'] for r in fds), current_fd_error=max(r['current_error'] for r in fds), fd_budget_fraction=max(r['max_fd_budget_fraction'] for r in fds),
            gauge_error=max(r['gauge_error'] for r in rows), uniform_oracle_error=max(r['uniform_oracle_error'] for r in rows), symbol_fd_uncertainty=max(r['symbol_uncertainty'] for r in rows),
            spectral_matrix_shift=max(p['matrix_operator_norm_difference'] for p in probes), spectral_moment_scaled_shift=max(max(p['electronic_moments']['scaled_differences']) for p in probes),
            relative_spatial_estimate=max(m['relative_estimate'] for m in passed)),
        minimum_eigenvalue=min(r['min_eigenvalue'] for r in rows), spectral_probe_count=len(probes), spectral_probes=probes, negative_control=negative,
        spatial_metrics=metrics, relative_spatial_pass_count=len(passed), relative_spatial_unresolved_count=len(limited), spatial_fail_count=len(failed), analytic_gradient_errors=gradient,
        decision=dict(stage3a_development_closed_registered_scope=eligible, accepted_for_stage3b_boundary_development=eligible, spatial_precision_fully_resolved=False,
            universal_spatial_accuracy_certified=False, stage3_closed=False, production_promotion=False),
        limitations=['One vacuum amplitude energy response remains ROUND_OFF_LIMITED: below the declared absolute difference floor, not a demonstrated machine-roundoff limit or relative PASS.',
            'Largest conditional response estimate is the thermal amplitude current; its observed order is about 0.608. No fourth-order current claim follows.',
            'Directional force/current FD covers six N16 cases and one reused N8 pilot; no claim that FD was repeated in all 18 cases.',
            'D.36 is a local longitudinal constitutive condition, not the complete discrete Hessian or a proof of dynamical/2D stability.',
            '630/1260 comparisons cover preregistered representative nodes, not all fields or the exact spectral continuum.',
            'Periodic static fields with frozen populations do not certify physical terminals, reservoirs, coupled potential/circuit, time evolution, material NbN or detector observables.'],
        next_work=['Reuse this exact periodic nodal evidence without rerunning unchanged cases.',
            'Validate open-boundary work, charge/current continuity and reservoir exchange with their own registered controls.',
            'Validate circuit state and work/power balances before a spatial transient; refine space and time against the intended transient observables.'])


def notes(report):
    if report['integrity_errors'] or report['spatial_fail_count']:
        return '# Auditoría nodal sin cierre\n\nNo se acepta el cierre: véanse los errores y dictámenes originales en `nodal_audit.json`.\n'
    m = report['maxima']
    lines = ['# Cierre de desarrollo de 3A: ensayos nodales registrados', '',
        'La campaña nodal permite cerrar 3A en el alcance estático registrado y continuar el desarrollo de bordes. '
        'No cierra la etapa 3 completa ni certifica precisión espacial universal. Los resultados y las tolerancias originales se conservan.', '',
        f"La auditoría de datos guardados dio **{report['status']}**: ocho fuentes, 18 casos, {len(report['raw_inventory'])} archivos y {report['sampled_nodes']} nodos espaciales comprobados. "
        f"El lote registró {report['recorded_runtime_seconds']/60:.3f} minutos. No se evaluaron espectros, fuerzas ni trayectorias nuevas.", '',
        f"Las diferencias finitas de fuerza y corriente se ejecutaron en los seis casos de 16 celdas y en el piloto reutilizado de fase térmica de ocho celdas. "
        f"Sus errores absolutos máximos son {m['force_fd_error']:.3e} y {m['current_fd_error']:.3e}. "
        'Los otros 11 resultados conservan expresamente que no repitieron diferencias finitas. '
        f"La covariancia de fase y el oráculo uniforme discreto alcanzan diferencias máximas de {m['gauge_error']:.3e} y {m['uniform_oracle_error']:.3e}.", '',
        f"Las 11 respuestas no nulas tienen diez estimaciones relativas aceptadas y una sin certificado relativo. "
        f"La mayor estimación es {100*m['relative_spatial_estimate']:.6f}%, frente al objetivo registrado de 1%. "
        'Corresponde a la corriente inducida por la modulación térmica de amplitud; el orden observado es 0.608. '
        'El estimador conservador incorpora esa razón medida, pero sigue siendo una extrapolación condicional; no demuestra orden cuatro de la corriente.', '',
        '| Caso | Energía de exceso | Fuerza inducida | Corriente inducida |', '|---|---:|---:|---:|']
    for case in ('weak_amplitude_thermal', 'weak_phase_thermal', 'weak_phase_nonthermal', 'weak_amplitude_vacuum'):
        values=[]
        for key in ('excess_energy_density_bar', 'induced_force_modes', 'induced_current_modes'):
            v=next(x for x in report['spatial_metrics'] if x['case']==case and x['observable']==key)
            values.append(f"{v['relative_estimate']*100:.6g}%" if v['relative_estimate'] is not None else 'Bajo piso declarado; sin pase relativo' if v['status']=='ROUND_OFF_LIMITED' else 'Cero analítico; pase absoluto')
        lines.append('| '+case+' | '+' | '.join(values)+' |')
    lines += ['', 'En amplitud/vacío, la diferencia de energía entre 16 y 32 celdas es 1.835e-9, inferior al piso declarado de 1e-8. '
        'Se conserva `ROUND_OFF_LIMITED`; no significa que se haya demostrado el límite de redondeo de la máquina ni se convierte en un pase relativo. '
        'La corriente exactamente nula pasa sólo su control absoluto. Los uniformes se cotejan con su oráculo discreto, sin cocientes 0/0.', '',
        'El error del gradiente frente a la referencia analítica disminuye en 8, 16 y 32 celdas: '+', '.join(f'{v*100:.6g}%' for v in report['analytic_gradient_errors'])+'. '
        'El valor fino cumple el objetivo de 1%.', '',
        f"El menor valor propio local D.36 es {report['minimum_eigenvalue']:.8f}, frente a una incertidumbre de diferenciación máxima de {m['symbol_fd_uncertainty']:.3e}. "
        f"Los {report['spectral_probe_count']} puntos representativos comparados entre 630 y 1260 nodos ocupacionales conservan el signo positivo. "
        f"La variación máxima de matriz es {m['spectral_matrix_shift']:.3e}; la variación escalada máxima de momentos es {m['spectral_moment_scaled_shift']:.3e}. "
        'El control negativo conserva el valor propio -0.30661962 y su rechazo esperado. Estas comprobaciones no prueban estabilidad dinámica, todo el Hessiano nodal ni el espectro continuo.', '',
        'Se reutiliza la evidencia periódica íntegra. Lo siguiente requiere controles propios de bordes abiertos, intercambio con reservorios, continuidad de corriente y balance de trabajo/potencia del circuito. '
        'Después corresponde el ensayo espacial en el tiempo y su refinamiento sobre observables del dispositivo. No hace falta repetir esta campaña sin cambios para iniciar ese trabajo.', '',
        'Fuentes y trazabilidad: registro nodal `../nodal_20260923/registration.json`; resultados íntegros `raw/nodal_campaign/`; '
        'inventario SHA256, dictámenes y reproducción numérica en `nodal_audit.json`; script `sandbox/stage3_spatial/ports_20260923/audit_nodal.py`.', '']
    return '\n'.join(lines)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, default=DEST/'raw/nodal_campaign')
    parser.add_argument('--output-root', type=Path, default=DEST)
    parser.add_argument('--verify-only', action='store_true')
    args=parser.parse_args()
    report=audit(args.raw.resolve())
    if not args.verify_only:
        paths=[args.output_root/'nodal_audit.json', args.output_root/'nodal_closure.md']
        if any(p.exists() for p in paths):
            raise SystemExit('Preserve existing audit outputs; use --verify-only or a new output directory.')
        args.output_root.mkdir(parents=True, exist_ok=True)
        paths[0].write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')
        paths[1].write_text(notes(report), encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('status','integrity_errors','completed_results','spectral_probe_count','relative_spatial_pass_count','relative_spatial_unresolved_count','spatial_fail_count','maxima','decision')}, indent=2))
    raise SystemExit(bool(report['integrity_errors']))


if __name__=='__main__':
    main()
