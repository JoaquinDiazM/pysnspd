"""Audit saved stage3A outputs; no catalogue query or physical evaluation.

Only JSON/NPZ reading, hashes, algebra on recorded arrays, and an analytic
small-mesh expansion of the already registered vacuum formula are used.
Existing audit outputs are preserved. --verify-only writes nothing.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math

import numpy as np


ROOT = Path(__file__).resolve().parents[3]
REGISTRATION = ROOT/'docs/implementation/stage3/iteration_20260923/registration.json'
OUTPUT = ROOT/'docs/implementation/stage3/review_20260923'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def complex_array(value):
    return np.asarray(value['real'])+1j*np.asarray(value['imag'])


def vacuum_midpoint_expansion(records):
    """Algebra of D.8 plus the documented gapped vacuum; no spectrum solve.

    This is a postprocessing explanation, not an extra numerical experiment.
    All measured values below come from the saved uniform-vacuum cases.
    """
    cases=sorted((r for r in records if r['case']['id']=='uniform_vacuum'),key=lambda r:r['cells'])
    case=cases[0]['case'];a=case['amplitude_base'];q=case['q_bare_ell0']
    ratio=math.pi*math.exp(-0.5772156649015328606)
    kappa=math.pi/4;rho=a*a;s=rho+.1**2;m=rho/s
    gamma=(m*q)**2/ratio
    assert 0<=gamma<a
    ua=2*a*math.log(a)+math.pi*gamma/2
    ug=math.pi*a/2-2*gamma/3
    uag=math.pi/2;ugg=-2/3
    c=m/4-1/6
    b=-1/12-m*m*(2*c-1/4)
    energy_parts=dict(amplitude_midpoint=-a*ua*q*q/8,
                      spectral_flow_map=2*ug*gamma*q*q*c,
                      explicit_gradient=kappa*a*a*q**4*b)
    current_parts=dict(amplitude_midpoint=-a*q*(ua+uag*gamma)/4,
                       spectral_flow_map=4*c*gamma*q*(ugg*gamma+2*ug),
                       explicit_gradient=4*kappa*a*a*q**3*b)
    rows=[]
    for record in cases:
        h=record['h_bar'];oracle=record['uniform_reference']
        e=record['energy_bar']/record['length_bar'];ec=oracle['continuum_energy_density']
        j=oracle['link_current_bar'];jc=oracle['continuum_current_bar']
        energy_shift=sum(energy_parts.values())*h*h
        current_shift=sum(current_parts.values())*h*h
        rows.append(dict(cells=record['cells'],h_bar=h,
            actual_energy_density=e,continuum_energy_density=ec,
            measured_energy_shift=e-ec,predicted_leading_energy_shift=energy_shift,
            actual_current_bar=j,continuum_current_bar=jc,
            measured_current_shift=j-jc,predicted_leading_current_shift=current_shift,
            current_relative_difference=abs(j-jc)/abs(jc),
            current_shift_fraction_explained=current_shift/(j-jc),
            remaining_current_shift_after_leading_term=(j-jc)-current_shift))
    return dict(scope='Analytic expansion/postprocessing of registered vacuum data; no new physical calculation.',
        field_units='Delta0',coordinate_units='ell0',a=a,q_bare_ell0=q,
        gamma_continuum=gamma,u_a_continuum=ua,
        definitions=dict(m='a^2/(a^2+delta^2)',c='m/4-1/6',b='-1/12-m^2*(2*c-1/4)'),
        midpoint='a_link=a*cos(q*h/2); |d_link|=2*a*sin(q*h/2)/h',
        flow_expansion='q_delta,h=m*q*[1+c*q^2*h^2]+O(h^4)',
        energy_expansion='e_h=e+C_e*h^2+O(h^4)',
        current_expansion='J_h=J+C_J*h^2+O(h^4), C_J=partial_q(C_e) at fixed a,p',
        energy_coefficient_parts=energy_parts,current_coefficient_parts=current_parts,
        energy_coefficient_formula='C_e=-a*u_a*q^2/8+2*u_Gamma*Gamma*q^2*c+kappa*a^2*q^4*b',
        current_coefficient_formula='C_J=-a*q*(u_a+u_aGamma*Gamma)/4+4*c*Gamma*q*(u_GammaGamma*Gamma+2*u_Gamma)+4*kappa*a^2*q^3*b',
        vacuum_derivatives='u_a=2*a*ln(a)+pi*Gamma/2; u_Gamma=pi*a/2-2*Gamma/3; u_aGamma=pi/2; u_GammaGamma=-2/3',
        source_paths=['docs/modelo_v0_4/D_sintesis_plan_y_verificaciones_v0_4.md',
                      'pysnspd/experimental/energy_catalog.py',
                      'pysnspd/experimental/spatial_functional.py',
                      'sandbox/stage3_spatial/run_static_batch.py'],rows=rows,
        interpretation='Cartesian midpoint quadrature changes the sampled amplitude of a helical nodal field. At this off-stationary amplitude u_a is nonzero, so a variationally correct current can retain a sizeable spatial bias. This is not a sign or energy-derivative failure.')


def audit(raw):
    reg=read(REGISTRATION);manifest=read(raw/'manifest.json');summary=read(raw/'summary.json')
    errors=[];records=[];rows=[];probes=[]
    def check(condition,label):
        if not condition:errors.append(label)
    expected={f"{c['id']}_n{n}":(c,n) for c in reg['cases'] for n in reg['cell_counts']}
    check(set(manifest['tasks'])==set(expected) and len(manifest['tasks'])==len(expected),'declared coverage')
    sources=manifest['source_hashes']
    source_checks={path:dict(recorded_sha256=value,current_sha256=sha(ROOT/path)) for path,value in sources.items()}
    for path,value in source_checks.items():check(value['recorded_sha256']==value['current_sha256'],'source '+path)
    check(sources==summary['source_hashes'],'summary source contract')
    for task in manifest['tasks']:
        path=raw/(task+'.json');r=read(path);receipt=read(path.with_suffix('.receipt.json'))
        case,cells=expected[task]
        check(receipt['json_sha256']==sha(path),task+' result receipt')
        check(receipt['arrays_sha256']==r['arrays_sha256']==sha(path.with_suffix('.npz')),task+' arrays receipt')
        check(receipt['initial_sha256']==r['initial_sha256']==sha(raw/(task+'_initial.npz')),task+' initial receipt')
        check(receipt['source_hashes']==r['source_hashes']==sources,task+' source contract')
        check(r['case']==case and r['cells']==cells,task+' case parameters')
        with np.load(path.with_suffix('.npz'),allow_pickle=False) as arrays, np.load(raw/(task+'_initial.npz'),allow_pickle=False) as initial:
            check(all(np.isfinite(arrays[key]).all() for key in arrays.files),task+' finite arrays')
            p=arrays['p_links'];z=arrays['delta'];x=arrays['count_nodes']
            check(p.shape==(cells,630) and z.shape==(cells,) and x.shape==(630,),task+' array shapes')
            check(np.all((p>=0)&(p<=1)) and np.all(np.diff(x)>0),task+' populations/count')
            for key in ('delta','p_links','count_nodes'):
                check(np.array_equal(arrays[key],initial[key]),task+' static state '+key)
        derivatives=r['derivative_checks'];symbols=r['principal_symbols']
        check(r['status']=='PASS_STATIC_IDENTITIES_AND_SAMPLED_SYMBOL',task+' structural status')
        for name,value in derivatives.items():
            err=abs(value['predicted']-value['fd_estimates'][-1])
            tolerance=reg['thresholds']['derivative_absolute']+reg['thresholds']['derivative_relative']*abs(value['fd_estimates'][-1])
            uncertainty=abs(value['fd_estimates'][-1]-value['fd_estimates'][-2])
            check(value['passed'] and value['status']=='PASS' and err<=tolerance and
                  uncertainty<=reg['thresholds']['reference_budget_fraction']*tolerance,task+' derivative '+name)
        check(max(r['global_phase_checks'].values())<=reg['thresholds']['global_phase_scaled_error'],task+' gauge')
        check(max(r['uniform_reference_errors'].values())<=1e-4,task+' discrete uniform oracle')
        for symbol in symbols:
            check(symbol['stable'] and min(symbol['eigenvalues'])>symbol['uncertainty'],task+' symbol margin')
            check(np.allclose(np.linalg.eigvalsh(symbol['matrix']),symbol['eigenvalues'],rtol=1e-12,atol=1e-13),task+' recorded eigenvalues')
        for probe in r['spectral_sign_comparison']['probes']:
            margin=min(probe['coarse_min'],probe['fine_min'])-probe['coarse_uncertainty']-probe['fine_uncertainty']-abs(probe['fine_min']-probe['coarse_min'])
            check(margin>0 and abs(margin-probe['remaining_positive_margin'])<1e-12,task+' spectral sign margin')
            cd=np.asarray(probe['coarse_matrix']);fd=np.asarray(probe['fine_matrix'])
            check(np.allclose(np.linalg.eigvalsh(cd),probe['coarse_eigenvalues'],rtol=1e-12,atol=1e-13),task+' coarse spectral eigenvalues')
            check(np.allclose(np.linalg.eigvalsh(fd),probe['fine_eigenvalues'],rtol=1e-12,atol=1e-13),task+' fine spectral eigenvalues')
            probes.append(dict(task=task,**probe))
        rows.append(dict(task=task,status=r['status'],runtime_seconds=r['runtime_seconds'],
            force_error=derivatives['force']['absolute_error'],current_error=derivatives['current']['absolute_error'],
            max_derivative_budget_fraction=max(v['absolute_error']/v['tolerance'] for v in derivatives.values()),
            max_gauge_error=max(r['global_phase_checks'].values()),max_uniform_error=max(r['uniform_reference_errors'].values()),
            minimum_eigenvalue=min(min(s['eigenvalues']) for s in symbols),
            maximum_symbol_fd_uncertainty=max(s['uncertainty'] for s in symbols),
            spectral_status=r['spectral_sign_comparison']['status']))
        records.append(r)
    comparisons=[]
    for row in summary['grid_comparisons']:
        case_records=sorted((r for r in records if r['case']['id']==row['case']),key=lambda r:r['cells'])
        for key,metric in row.get('metrics',{}).items():
            if key=='isolated_gradient':
                values=[r['analytic_gradient']['relative_rms'] for r in case_records]
                check(np.allclose(values,metric['errors'],rtol=1e-13,atol=0),row['case']+' gradient reproduction')
                continue
            arrays=[np.atleast_1d(complex_array(r[key]) if isinstance(r[key],dict) else r[key]) for r in case_records]
            d1=float(np.linalg.norm(arrays[1]-arrays[0]));d2=float(np.linalg.norm(arrays[2]-arrays[1]));norm=float(np.linalg.norm(arrays[-1]))
            check(np.allclose([d1,d2,norm],[metric['coarse_medium_difference'],metric['medium_fine_difference'],metric['fine_norm']],rtol=1e-13,atol=1e-16),row['case']+' '+key+' grid differences')
            if metric['relative_estimate'] is not None:
                ratio=d1/d2;estimate=max(d2/3,d2/(ratio-1))
                check(abs(estimate/norm-metric['relative_estimate'])<=1e-13,row['case']+' '+key+' extrapolation')
            comparisons.append(dict(case=row['case'],observable=key,status=metric['status'],
                relative_estimate=metric['relative_estimate'],observed_order=metric['observed_order']))
    negative=read(raw/'negative_control.json')
    check(negative==summary['negative_control'],'negative control provenance')
    check(negative['status']=='PASS_EXPECTED_REJECTION' and not negative['admitted_for_evolution']
          and min(negative['measured_eigenvalues'])<0 and negative['error']<1e-8,'negative control')
    inventory={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(raw.iterdir()) if p.is_file()}
    maxima=dict(force_error=max(r['force_error'] for r in rows),current_error=max(r['current_error'] for r in rows),
        gauge_error=max(r['max_gauge_error'] for r in rows),uniform_oracle_error=max(r['max_uniform_error'] for r in rows),
        symbol_fd_uncertainty=max(r['maximum_symbol_fd_uncertainty'] for r in rows),
        spectral_matrix_shift=max(p['matrix_operator_norm_difference'] for p in probes),
        spectral_moment_scaled_shift=max(max(p['electronic_moments']['scaled_differences']) for p in probes))
    failed_spatial=sum(v['status']=='FAIL_PRECISION_TARGET' for v in comparisons)
    return dict(schema='pysnspd.stage3.saved_results_audit.v1',
        status=('AUDIT_INTEGRITY_OR_CONTRACT_FAILURE' if errors else
                'SAVED_RESULTS_VERIFIED_SPATIAL_PRECISION_NOT_MET' if failed_spatial else 'SAVED_RESULTS_VERIFIED'),
        scope='Audit of saved data only; no physical reevaluation, no new spectra or trajectories.',
        script_sha256=sha(__file__),registration_sha256=sha(REGISTRATION),raw_directory=str(raw.relative_to(ROOT)),
        raw_inventory=inventory,source_checks=source_checks,integrity_errors=errors,
        completed_results=len(records),sampled_links=sum(r['cells'] for r in records),spectral_probe_count=len(probes),
        recorded_runtime_seconds=summary['runtime_seconds'],recorded_batch_status=summary['status'],
        case_checks=rows,maxima=maxima,minimum_eigenvalue=min(r['minimum_eigenvalue'] for r in rows),
        spectral_probes=probes,negative_control=negative,spatial_metrics=comparisons,
        failed_relative_spatial_observables=failed_spatial,
        analytic_gradient_errors=next(r['metrics']['isolated_gradient']['errors'] for r in summary['grid_comparisons'] if 'isolated_gradient' in r.get('metrics',{})),
        midpoint_diagnosis=vacuum_midpoint_expansion(records),
        admitted_scope=['Discrete energy derivatives and their signs in the 18 saved configurations.',
            'Global-phase covariance and exact discrete uniform-helical oracle.',
            'Positive longitudinal D.36 on the candidate mesh and recorded representative spectral comparisons.',
            'Analytic amplitude-gradient control at N32 within 1%, and exact zero-current absolute control.'],
        not_admitted=['General 1% spatial precision of the registered nonuniform responses.',
            'Full stage3 completion, physical boundaries/reservoirs, potential/circuit or time evolution.',
            'Production, material NbN rates, photon threshold or detector latency.'],
        repetition_policy=['Conservar los 18 resultados y sus fallos espaciales; no repetir identidades sin cambios sólo para recrear evidencia.',
            'Esta evidencia no exige repetir trayectorias de etapa2 ni aumentar globalmente la resolución espectral.',
            'Una reconstrucción espacial nueva necesita comprobaciones variacionales y de respuesta espacial propias con los mismos perfiles físicos, además de D.36 en sus estados muestreados.',
            'Si sólo se añaden mallas espaciales más finas con el mismo operador, reutilizar esta evidencia bajo sus fuentes y entradas exactas y registrar previamente la ampliación.',
            'No atribuir precisión al siguiente ensayo espacial usando solamente conservación de energía o error relativo sobre una gran energía de fondo.'])


def notes(report):
    maxima=report['maxima'];fine=report['midpoint_diagnosis']['rows'][-1]
    lines=['# Auditoría de resultados guardados de 3A','',
        'Se verificaron los 18 resultados y sus archivos de estado, recibos y siete fuentes. '
        'La campaña completó los casos, pero no cumplió la precisión espacial registrada. '
        'Esta auditoría no evaluó espectros ni fuerzas nuevos, y no integró trayectorias.','',
        f"Integridad: {'PASS' if not report['integrity_errors'] else 'FAIL'}. "
        f"Archivos inventariados: {len(report['raw_inventory'])}. Enlaces comprobados: {report['sampled_links']}. "
        f"Duración registrada del lote: {report['recorded_runtime_seconds']/60:.3f} min.",'',
        '## Lo comprobado','',
        f"- Los 18 casos aprobaron las identidades discretas. Diferencias máximas de fuerza y corriente: "
        f"{maxima['force_error']:.3e} y {maxima['current_error']:.3e}; invariancia/covariancia de fase: {maxima['gauge_error']:.3e}.",
        f"- El menor valor propio es {report['minimum_eigenvalue']:.8f}; la incertidumbre de diferenciación máxima es {maxima['symbol_fd_uncertainty']:.3e}.",
        f"- Los {report['spectral_probe_count']} puntos representativos 630/1260 conservan signo positivo. "
        f"Desplazamiento máximo de la matriz: {maxima['spectral_matrix_shift']:.3e}; diferencia máxima escalada de momentos: {maxima['spectral_moment_scaled_shift']:.3e}. "
        'Esto describe los estados y resoluciones muestreados, no una cota global del espectro continuo.',
        '- El control negativo D.36 permanece rechazado: valor propio -0.30661962. No se interpreta su rechazo esperado como fallo del código.',
        f"- El error del gradiente analítico baja de {report['analytic_gradient_errors'][0]*100:.4f}% a "
        f"{report['analytic_gradient_errors'][1]*100:.4f}% y {report['analytic_gradient_errors'][2]*100:.4f}%. "
        'Pasa el control de 1% en 32 celdas; la corriente nula de este caso pasa sólo su control absoluto.','',
        '## Precisión espacial que falta','',
        'Los siguientes porcentajes son estimaciones condicionales de extrapolación, no errores frente a una solución exacta. '
        'El presupuesto original sigue siendo 1%.','',
        '| Caso | Energía de exceso | Fuerza inducida | Corriente inducida |','| --- | ---: | ---: | ---: |']
    cases=[]
    for item in report['spatial_metrics']:
        if item['case'] not in cases:cases.append(item['case'])
    for case in cases:
        values={v['observable']:v for v in report['spatial_metrics'] if v['case']==case}
        formatted=[]
        for name in ('excess_energy_density_bar','induced_force_modes','induced_current_modes'):
            value=values[name]['relative_estimate']
            formatted.append(f'{100*value:.4f}%' if value is not None else 'Cero analítico; control absoluto')
        lines.append('| '+case+' | '+' | '.join(formatted)+' |')
    lines += ['',f"Son {report['failed_relative_spatial_observables']} observables fuera del objetivo. "
        'Las diferencias decrecen al refinar, pero eso no basta para afirmar precisión del 1%. '
        'El caso de amplitud en vacío queda apenas por encima; su fallo se conserva.','',
        '## Por qué una identidad exacta puede coexistir con ese error','',
        'Para una hélice nodal de amplitud a y flujo q, el promedio cartesiano consulta '
        '`a_link = a cos(qh/2)` y `|d_link| = 2a sin(qh/2)/h`. '
        'La amplitud del punto medio depende así del flujo y del tamaño de celda. '
        'La fuerza local no es cero en el estado prescrito a=0.9: no es un equilibrio autoconsistente. '
        'Por ello, derivar correctamente esa energía discreta incluye un sesgo de corriente de orden h².','',
        'Con `m=a²/(a²+delta²)`, `c=m/4-1/6`, `b=-1/12-m²(2c-1/4)` y '
        '`Gamma=(mq)²/r`, la expansión es:', '',
        '`e_h=e+h²[-a*u_a*q²/8+2*u_Gamma*Gamma*q²*c+kappa*a²*q⁴*b]+O(h⁴)`.', '',
        'Su derivada a ocupaciones y amplitud fijas da `J_h=J+h²*C_J+O(h⁴)`, con:', '',
        '`C_J=-a*q*(u_a+u_aGamma*Gamma)/4+4*c*Gamma*q*(u_GammaGamma*Gamma+2*u_Gamma)+4*kappa*a²*q³*b`.', '',
        'Se usa aquí sólo la expresión analítica ya documentada del vacío con gap. '
        'Los tres aportes y las medidas originales quedan separados en el JSON de auditoría.','',
        f"En el uniforme de vacío a 32 celdas, la corriente discreta difiere del límite continuo en "
        f"{100*fine['current_relative_difference']:.3f}%. Su desplazamiento medido es "
        f"{fine['measured_current_shift']:.6g}; el término h² predice {fine['predicted_leading_current_shift']:.6g}, "
        f"es decir {100*fine['current_shift_fraction_explained']:.2f}% del desplazamiento. "
        'Este acuerdo identifica un sesgo de reconstrucción espacial, aunque no demuestra que todos los errores no uniformes provengan exclusivamente de ese término.','',
        'El error en respuestas inducidas puede ser mucho mayor que en la corriente total, '
        'porque se resta una referencia uniforme y la señal restante es pequeña. '
        'No corresponde ocultarlo normalizando con la energía total del fondo ni reinterpretarlo como fallo de la corriente física.','',
        '## Qué conservar y qué volver a comprobar','']
    lines += ['- '+value for value in report['repetition_policy']]
    lines += ['',
        'La evidencia sostiene el desarrollo del funcional, sus signos y la estabilidad local ensayada. '
        'Falta resolver o delimitar la precisión de la reconstrucción espacial antes de atribuir exactitud a la respuesta del siguiente ensayo. '
        'Bordes, reservorios, potencial y circuito conservan sus verificaciones propias; la etapa 3 completa sigue abierta.','',
        '## Reproducción y fuentes','',
        '`python sandbox/stage3_spatial/review_20260923/audit_saved.py --raw <directorio_del_lote> --verify-only`', '',
        'El programa lee exclusivamente JSON/NPZ, verifica hashes y reproduce comparaciones algebraicas. '
        'Para generar una copia de la auditoría, indicar un directorio nuevo mediante `--output-root`. '
        'El inventario por archivo y los hashes de las fuentes están en `saved_results_audit.json`.','',
        '- Modelo continuo: `docs/modelo_v0_4/D_sintesis_plan_y_verificaciones_v0_4.md`, D.8–D.10 y D.36.',
        '- Vacío analítico: `pysnspd/experimental/energy_catalog.py`, `vacuum_state`.',
        '- Reconstrucción discreta: `pysnspd/experimental/spatial_functional.py`.',
        '- Registro previo: `docs/implementation/stage3/iteration_20260923/registration.json`.',
        '- Resultados originales: manifiesto, 18 casos y `summary.json` del directorio del lote.','']
    return '\n'.join(lines)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw',type=Path,required=True)
    parser.add_argument('--output-root',type=Path,default=OUTPUT)
    parser.add_argument('--verify-only',action='store_true')
    args=parser.parse_args();raw=args.raw.resolve();result=audit(raw)
    if not args.verify_only:
        target=args.output_root.resolve();target.mkdir(parents=True,exist_ok=True)
        paths=[target/'saved_results_audit.json',target/'audit_notes.md']
        if any(path.exists() for path in paths):raise SystemExit('Audit output exists; use --verify-only or a fresh directory.')
        paths[0].write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
        paths[1].write_text(notes(result),encoding='utf-8')
    print(json.dumps({key:result[key] for key in ('status','integrity_errors','completed_results','spectral_probe_count','failed_relative_spatial_observables')},indent=2))
    if result['integrity_errors']:raise SystemExit(1)


if __name__=='__main__':
    main()
