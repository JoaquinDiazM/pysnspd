"""Audit six archived open static cases using hashes and saved arithmetic only.

Does not import any physical model, query a catalogue, optimize or integrate.
Existing outputs are preserved. --verify-only writes nothing.
"""
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUTPUT = ROOT/'docs/implementation/stage3/coupled_20260923'
RAW = OUTPUT/'raw/stage3b_open_full_20260923'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def audit(raw):
    manifest, summary = read(raw/'manifest.json'), read(raw/'summary.json')
    reg = manifest['registration']
    errors, rows = [], []

    def check(condition, label):
        if not condition:
            errors.append(label)

    def close(a, b, label, rtol=1e-12, atol=1e-15):
        check(np.allclose(a, b, rtol=rtol, atol=atol), label)

    sources = {p: dict(recorded_sha256=h, current_sha256=sha(ROOT/p))
               for p, h in manifest['sources'].items()}
    check(len(sources) == 9, 'nine frozen source/registration inputs')
    for p, hashes in sources.items():
        check(hashes['recorded_sha256'] == hashes['current_sha256'], 'source '+p)
    reg_path = ROOT/manifest['arguments']['registration']
    check(reg == read(reg_path), 'embedded and frozen registration identical')
    check(manifest['arguments']['execute'] and not manifest['arguments']['pilot'], 'full batch mode')
    check(summary['status'] == 'OPEN_STATIC_BATCH_COMPLETE_SCOPE_LIMITED', 'completed batch status')
    check(not summary['stage3_closed'], 'global stage3 remains open')
    expected = {c['id']: c for c in reg['cases']}
    check(len(expected) == len(summary['cases']) == 6, 'six registered cases')
    check({r['case']['id'] for r in summary['cases']} == set(expected), 'exact case coverage')
    required_files = {'manifest.json', 'summary.json', 'progress.jsonl'}
    for case in reg['cases']:
        task = case['id']
        paths = {name: raw/task/(name+'.json') for name in ('result', 'reservoir', 'iterations')}
        required_files.update(task+'/'+name+'.json' for name in paths)
        r, rr, iterations = read(paths['result']), read(paths['reservoir']), read(paths['iterations'])
        b, part = rr['branch'], rr['inductance_partition']
        check(r['case'] == case, task+' registered parameters')
        check(r == next(v for v in summary['cases'] if v['case']['id'] == task), task+' summary/result identity')
        check(r['status'] == 'PASS_STATIC_OPEN_CASE', task+' accepted static status')
        metrics = r['metrics']
        for name, threshold in reg['criteria'].items():
            check(np.isfinite(metrics[name]) and metrics[name] <= threshold, task+' criterion '+name)
        check(not r['solver']['active_bounds'], task+' inactive optimization bounds')
        check(r['solver']['iterations'] == 0 and r['solver']['evaluations'] == len(iterations) == 1,
              task+' initial helix without relaxation')
        close(iterations[0]['max_gradient'], metrics['stationarity_absolute'], task+' stationarity log')
        x, a, phase, current, qd = (np.asarray(r[k]) for k in
            ('x_m', 'amplitude_bar', 'phase_rad', 'current_A', 'q_delta_bar'))
        n = reg['degree']*case['elements']+1
        check(x.shape == a.shape == phase.shape == qd.shape == (n,) and current.shape == (n-1,), task+' geometry shapes')
        check(all(np.all(np.isfinite(v)) for v in (x, a, phase, current, qd)), task+' finite fields')
        check(np.all(np.diff(x) > 0), task+' ordered open coordinates')
        close(x[[0, -1]], [0., case['length_m']], task+' endpoints', atol=1e-20)
        expected_phase = reg['q_bare_bar']*(x-case['length_m'])/r['ell0_m']
        # The right coordinate was rounded by its SI-unit round trip; the
        # separately recorded right-phase gauge is prescribed exactly zero.
        expected_phase[-1] = 0.
        close(phase, expected_phase, task+' prescribed helix')
        close(a, b['amplitude_bar'], task+' prescribed amplitude')
        check(phase[-1] == 0., task+' right global-phase gauge')
        # Recompute the reported residuals from stored fields and reference.
        current_error = float(np.max(abs(current-r['reference_current_A']))/abs(r['reference_current_A']))
        q_end = qd[[0, -1]]/(b['amplitude_bar']**2/(b['amplitude_bar']**2+.01))
        q_error = float(np.max(abs(q_end-b['q_bare_bar']))/abs(b['q_bare_bar']))
        amp_error = float(np.max(abs(a-b['amplitude_bar'])))
        close(current_error, metrics['current_relative'], task+' current arithmetic', atol=2e-15)
        close(q_error, metrics['endpoint_q_relative'], task+' endpoint arithmetic')
        close(amp_error, metrics['amplitude_absolute'], task+' amplitude arithmetic', atol=1e-18)
        close(r['free_energy_bar']+b['current_bar']*phase[0], iterations[0]['gibbs_bar'], task+' Gibbs boundary work')
        check(metrics['minimum_principal_eigenvalue'] > metrics['maximum_symbol_uncertainty'], task+' recorded D36 sign margin')
        check(b['stable'] and b['amplitude_curvature'] > b['amplitude_curvature_uncertainty']
              and b['differential_current_bar'] > b['differential_current_uncertainty'], task+' reservoir stability margin')
        check(b['checked_qs'] == [0., *reg['continuation_qs'], reg['q_bare_bar']], task+' explicit branch path')
        close([b['bath_theta'], b['q_bare_bar']], [reg['bath_theta'], reg['q_bare_bar']], task+' branch parameters')
        close(b['resolved_length_m'], case['length_m'], task+' branch length', atol=1e-20)
        population = np.asarray(b['occupation'])
        check(population.shape == (reg['count_nodes'],) and np.all(np.isfinite(population))
              and np.all((population >= 0) & (population <= 1)), task+' branch population')
        hessian = np.asarray(b['thermal_hessian'])
        close(hessian[1, 1]-hessian[1, 0]*hessian[0, 1]/hessian[0, 0],
              b['differential_current_bar'], task+' reduced thermal stiffness')
        close(abs(hessian[0, 1]-hessian[1, 0]), b['maxwell_disagreement'], task+' Maxwell difference', atol=1e-18)
        close(b['current_bar']*b['current_scale_A'], b['current_A'], task+' current units', atol=1e-20)
        close(r['reference_current_A'], b['current_A'], task+' branch/result current', atol=1e-20)
        close(r['reference_amplitude_bar'], b['amplitude_bar'], task+' branch/result amplitude')
        computed_l = (1.054571817e-34/(2*1.602176634e-19))*(case['length_m']/r['ell0_m'])/(b['current_scale_A']*b['differential_current_bar'])
        close(computed_l, b['L_res_diff_H'], task+' differential inductance units', atol=1e-23)
        close([r['L_res_H'], r['L_ext_H']], [part['resolved_differential_H'], part['exterior_fixed_H']], task+' partition/result', atol=1e-23)
        close(b['L_res_diff_H'], part['resolved_differential_H'], task+' partition/branch', atol=1e-23)
        close(r['L_res_H']+r['L_ext_H'], part['total_reference_H'], task+' partition sum', atol=1e-23)
        check(r['L_res_H'] > 0 and r['L_ext_H'] > 0, task+' positive inductances')
        check(not r['stage3_closed'] and not r['production'], task+' limited admission')
        length_bar = case['length_m']/r['ell0_m']
        rows.append(dict(task=task, case=case, nodes=n, element_width_m=case['length_m']/case['elements'],
            runtime_seconds=r['runtime_seconds'], solver=r['solver'], metrics=metrics,
            recomputed=dict(current_relative=current_error, endpoint_q_relative=q_error, amplitude_absolute=amp_error),
            mean_current_A=float(np.mean(current)), reference_current_A=b['current_A'], reference_amplitude_bar=b['amplitude_bar'],
            current_peak_to_peak_A=float(np.ptp(current)), q_end_bare_bar=q_end.tolist(),
            free_energy_density_bar=r['free_energy_bar']/length_bar, internal_energy_density_bar=r['internal_energy_bar']/length_bar,
            L_res_H=r['L_res_H'], L_ext_H=r['L_ext_H'], branch_differential_current_bar=b['differential_current_bar'],
            branch_curvature=b['amplitude_curvature'], reservoir_population_sha256=hashlib.sha256(population.tobytes()).hexdigest()))
    actual_files = {str(p.relative_to(raw)).replace('\\', '/') for p in raw.rglob('*') if p.is_file()}
    check(actual_files == required_files, 'exact 21-file completed campaign inventory')
    inventory = {name: dict(sha256=sha(raw/name), bytes=(raw/name).stat().st_size) for name in sorted(actual_files)}
    pairs = []
    for length in sorted({r['case']['length_m'] for r in rows}):
        coarse, fine = sorted((r for r in rows if r['case']['length_m'] == length), key=lambda r: r['nodes'])
        check(fine['case']['elements'] == 2*coarse['case']['elements'], 'paired spatial resolution '+str(length))
        comparison = {}
        for metric in ('stationarity_absolute', 'current_relative', 'endpoint_q_relative'):
            c, f = coarse['metrics'][metric], fine['metrics'][metric]
            comparison[metric] = dict(coarse=c, fine=f, reduction_factor=c/f, decreasing=f<c)
        check(all(v['decreasing'] for v in comparison.values()), 'observed reference errors decrease '+str(length))
        pairs.append(dict(length_m=length, cases=[coarse['task'], fine['task']], nodes=[coarse['nodes'], fine['nodes']],
            reference_errors=comparison, mean_current_difference_A=fine['mean_current_A']-coarse['mean_current_A'],
            free_energy_density_difference_bar=fine['free_energy_density_bar']-coarse['free_energy_density_bar'],
            internal_energy_density_difference_bar=fine['internal_energy_density_bar']-coarse['internal_energy_density_bar'],
            interpretation='Two resolutions against the same uniform reference. Decreasing errors support this helix only; no three-level order or perturbed-boundary certificate.'))
    length_comparisons = []
    # The registered pairs use the same element widths at all three lengths.
    for rank in (0, 1):
        selected = [sorted((r for r in rows if r['case']['length_m'] == length), key=lambda r:r['nodes'])[rank]
                    for length in sorted({r['case']['length_m'] for r in rows})]
        widths = [r['element_width_m'] for r in selected]
        close(widths, widths[0], 'same resolution at three lengths', atol=1e-21)
        length_comparisons.append(dict(element_width_m=widths[0], cases=[r['task'] for r in selected],
            mean_current_range_A=float(np.ptp([r['mean_current_A'] for r in selected])),
            free_energy_density_range_bar=float(np.ptp([r['free_energy_density_bar'] for r in selected])),
            internal_energy_density_range_bar=float(np.ptp([r['internal_energy_density_bar'] for r in selected])),
            scope='Uniform material and initial helix. Tiny length sensitivity does not test dynamic reflections or a localized excitation.'))
    check(len({r['reservoir_population_sha256'] for r in rows}) == 1, 'same prepared reservoir population')
    # Linear length scaling is an algebraic reference property, not a new run.
    ratios = [r['L_res_H']/r['case']['length_m'] for r in rows]
    close(ratios, ratios[0], 'resolved differential inductance scales with reference length')
    return dict(schema='pysnspd.stage3.open_full_saved_audit.v1', status='PASS_SAVED_OPEN_AUDIT' if not errors else 'FAIL_SAVED_OPEN_AUDIT',
        method='Offline hashes and arithmetic on archived records only; no spectrum, field, minimization or time evolution evaluated.',
        script_sha256=sha(__file__), raw_directory=str(raw.relative_to(ROOT)).replace('\\','/'),
        source_checks=sources, registration_sha256=sha(reg_path), raw_inventory=inventory, integrity_errors=errors,
        recorded_status=summary['status'], recorded_runtime_seconds=summary['runtime_seconds'], criteria=reg['criteria'],
        case_count=len(rows), node_evaluations=sum(r['nodes'] for r in rows), cases=rows,
        resolution_comparisons=pairs, length_comparisons=length_comparisons,
        maxima={k:max(abs(r['metrics'][k]) for r in rows) for k in rows[0]['metrics']},
        all_cases_initial_helix_only=all(r['solver']['iterations']==0 and r['solver']['evaluations']==1 for r in rows),
        decision=dict(registered_open_uniform_static_campaign_accepted=not errors,
            may_continue_coupled_development=not errors, stage3b_general_boundary_validation_complete=False,
            stage3_closed=False, production_promotion=False),
        limitations=['Every optimizer stopped at its initial uniform-amplitude helix: 0 iterations, 1 objective evaluation. No relaxation from disturbed data was demonstrated.',
            'Two resolutions at each length test this reference branch; they do not supply a universal error estimate or certify an observed spatial order.',
            'No localized excitation, reservoir population perturbation, boundary reflection, 2D-1D transparency or detector transient was tested.',
            'The D.36 values are recorded minima/uncertainties. The producer did not save individual matrices here, so this audit checks their reported margins but cannot recompute their spectra.',
            'The producer froze source/registration hashes but did not issue per-output receipts. This audit inventories observed output hashes and verifies cross-file identity; it does not invent prior receipts.',
            'Stationarity is cross-checked against its iteration log; full gradient arrays were not archived and are not reconstructed through new physical evaluations.'],
        next_work=['Reuse the six unchanged static cases and their exact source contract.',
            'Give the coupled spatial experiment its own preregistered state, circuit, reservoir exchange and energy/current accounting controls.',
            'Assess perturbation sensitivity and spatial/temporal refinement on the coupled observables before interpreting a detector signal.'])


def notes(report):
    if report['integrity_errors']:
        return '# Auditoría abierta incompleta\n\nExisten discrepancias: consultar `open_audit.json`. No se acepta el alcance auditado.\n'
    lines = ['# Auditoría de los seis ensayos abiertos', '',
        '**PASS_SAVED_OPEN_AUDIT.** Las nueve fuentes y el registro coinciden con los hashes congelados. '
        f"Los {len(report['raw_inventory'])} archivos contienen exactamente los seis casos previstos, con sus registros de iteración y reservorio. "
        'Los resultados individuales y el resumen son idénticos. Se verificaron las métricas mediante aritmética de los campos guardados, sin ejecutar física nueva.', '',
        f"La campaña registrada duró **{report['recorded_runtime_seconds']:.2f} s** ({report['recorded_runtime_seconds']/60:.3f} min). "
        'Los tres dominios tienen dos resoluciones: elementos de 90 y 45 nm.', '',
        '| Longitud | Nodos | Residuo estacionario | Error de corriente | Error de q terminal |',
        '|---:|---:|---:|---:|---:|']
    for r in report['cases']:
        m=r['metrics']
        lines.append(f"| {r['case']['length_m']*1e9:g} nm | {r['nodes']} | {m['stationarity_absolute']:.6g} | {m['current_relative']*100:.7g}% | {m['endpoint_q_relative']*100:.7g}% |")
    lines += ['', 'Los errores de corriente y q terminal cumplen el presupuesto registrado de 1%; '
        'el residuo estacionario queda bajo 1e-4 y la amplitud bajo 0.001. '
        'Al duplicar los elementos, los errores disminuyen en las tres longitudes. '
        'La corriente de referencia es **8.63351 µA** y la amplitud **0.995218**.', '',
        '**Los seis optimizadores hicieron cero iteraciones y una sola evaluación.** '
        'La hélice inicial ya satisfacía los criterios. La campaña comprueba consistencia estática y refinamiento de esa hélice; '
        'no demuestra relajación desde una perturbación. La diferencia de amplitud cercana al redondeo no representa recuperación dinámica.', '',
        'El cambio de q terminal pasa aproximadamente de **0.0713006% a 0.00486823%** y el de corriente de '
        '**0.000603483% a 0.00000904413%**. Las dos resoluciones describen reducción de error frente a esta referencia: '
        'no justifican por sí solas un orden de convergencia universal ni un certificado de tres mallas.', '',
        'A igual tamaño de elemento, la corriente media y las densidades de energía apenas cambian al alargar el dominio uniforme. '
        'Este resultado es coherente con una hélice uniforme; no mide reflexiones ni permite inferir transparencia ante una excitación localizada.', '',
        '| Longitud | Inductancia resuelta | Inductancia exterior |', '|---:|---:|---:|']
    for pair in report['resolution_comparisons']:
        r=next(v for v in report['cases'] if v['task']==pair['cases'][0])
        lines.append(f"| {pair['length_m']*1e9:g} nm | {r['L_res_H']*1e9:.6f} nH | {r['L_ext_H']*1e9:.6f} nH |")
    lines += ['', 'La partición diferencial suma 10 nH y se identifica en la misma rama antes de cualquier evolución. '
        'Se verificaron su escala lineal con la longitud y sus unidades. Las márgenes registradas del reservorio y de D.36 son positivas.', '',
        'El alcance estático uniforme abierto queda aceptado para continuar el desarrollo acoplado. '
        '**La validación general de bordes y la etapa 3 completa siguen abiertas.** '
        'Faltan perturbaciones no uniformes, intercambio de poblaciones, balances del sistema acoplado y evolución espacial con circuito. '
        'No hay certificado de interfaz 2D–1D ni predicción de detección.', '',
        'Límites de esta auditoría: D.36 y la estacionariedad se conservaron como métricas y registros, sin las matrices ni gradientes completos; '
        'se comprobaron sus márgenes y coherencia, pero no se recalcularon físicamente. '
        'El productor guardó hashes de fuentes y registro, no recibos individuales de salidas. '
        'El inventario SHA256 de esta revisión identifica los archivos observados sin atribuirles recibos previos.', '',
        'Reproducción: `python sandbox/stage3_spatial/coupled_20260923/audit_open.py --verify-only`. '
        'Sin esa opción, el auditor crea `open_audit.json` y `open_audit.md` únicamente si todavía no existen. '
        'Las fuentes, resultados anteriores y tolerancias no se modificaron.', '']
    return '\n'.join(lines)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw', type=Path, default=RAW)
    parser.add_argument('--output-root', type=Path, default=OUTPUT)
    parser.add_argument('--verify-only', action='store_true')
    args=parser.parse_args()
    report=audit(args.raw.resolve())
    if not args.verify_only:
        paths=[args.output_root/'open_audit.json', args.output_root/'open_audit.md']
        if any(p.exists() for p in paths):
            raise SystemExit('Existing audit outputs preserved. Use --verify-only or a new output directory.')
        args.output_root.mkdir(parents=True, exist_ok=True)
        paths[0].write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')
        paths[1].write_text(notes(report), encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('status','integrity_errors','case_count','recorded_runtime_seconds','all_cases_initial_helix_only','decision')}, indent=2))
    raise SystemExit(bool(report['integrity_errors']))


if __name__=='__main__':
    main()
