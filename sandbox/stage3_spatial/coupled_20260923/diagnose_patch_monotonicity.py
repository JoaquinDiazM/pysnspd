"""Inspect archived Chebyshev coefficients only; no source spectrum calls.

Replays the registered probe/stencil locations, reconstructs training nodes,
and samples an 81x81 field grid. It does not admit or fit a replacement model.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time

import numpy as np
from numpy.polynomial.chebyshev import chebval2d, chebvander2d

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT/'docs/implementation/stage3/coupled_20260923'
RAW = OUT/'raw/patch_check'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def diagnose(raw):
    start = time.perf_counter()
    manifest = read(raw/'manifest.json')
    reg = manifest['registration']
    source_checks = {p: dict(recorded_sha256=h, current_sha256=sha(ROOT/p))
                     for p,h in manifest['sources'].items()}
    source_matches = all(v['recorded_sha256'] == v['current_sha256'] for v in source_checks.values())
    replay, scans, training = {}, {}, {}
    for degree in reg['degrees']:
        with np.load(raw/f'patch_degree{degree}.npz', allow_pickle=False) as data:
            coeff, count, weights = data['coefficients'], data['count_nodes'], data['count_weights']
            ab, gb = data['amplitude_bounds'], data['gamma_bounds']

        def evaluate(a, g):
            xx = 2*(a-ab[0])/(ab[1]-ab[0])-1
            yy = 2*(g-gb[0])/(gb[1]-gb[0])-1
            return chebval2d(xx, yy, coeff)

        def record(a, g, *, probe=None, offset=None):
            energy = evaluate(a, g)
            increments = np.diff(energy)
            bad = np.flatnonzero(increments <= 0)
            pairs=[]
            for i in bad:
                pairs.append(dict(index_left=int(i), index_right=int(i+1),
                    count_left=float(count[i]), count_right=float(count[i+1]),
                    weight_left=float(weights[i]), weight_right=float(weights[i+1]),
                    energy_left=float(energy[i]), energy_right=float(energy[i+1]),
                    signed_increment=float(increments[i]),
                    jump_over_energy_scale=float(increments[i]/max(abs(energy[i]), abs(energy[i+1]))),
                    jump_over_float64_ulp=float(abs(increments[i])/np.spacing(max(abs(energy[i]), abs(energy[i+1]))))))
            return dict(amplitude=float(a), gamma=float(g), probe_one_based=probe,
                gamma_stencil_offset=offset, finite=bool(np.all(np.isfinite(energy))),
                positive=bool(np.all(energy>0)), monotone=not len(bad),
                minimum_energy=float(np.min(energy)), minimum_increment=float(np.min(increments)), bad_pairs=pairs)

        locations=[]
        for index, (fa, fg) in enumerate(reg['fractional_probes']):
            a, g = ab[0]+fa*(ab[1]-ab[0]), gb[0]+fg*(gb[1]-gb[0])
            for offset in (0., -.0002, .0002, -.0001, .0001, -.00005, .00005):
                if not gb[0] <= g+offset <= gb[1]:
                    raise ValueError('replay requested outside archived patch domain')
                locations.append(record(a, g+offset, probe=index+1, offset=offset))
        failures=[row for row in locations if not row['monotone'] or not row['positive'] or not row['finite']]
        replay[str(degree)] = dict(query_count=len(locations), failure_count=len(failures),
            failing_probe_centres=[row['probe_one_based'] for row in failures if row['gamma_stencil_offset']==0],
            bad_left_indices=sorted({p['index_left'] for row in failures for p in row['bad_pairs']}),
            first_failure=failures[0] if failures else None, failures=failures)

        # Reconstruct the training abscissae from the saved interpolant. These
        # are not fresh source values; residual inverse-transform roundoff remains.
        nodes=np.cos(np.pi*np.arange(degree+1)/degree)
        xx, yy=np.meshgrid(nodes, nodes, indexing='ij')
        vander=chebvander2d(xx.ravel(), yy.ravel(), [degree,degree])
        values=np.asarray([chebval2d(x, y, coeff) for x,y in zip(xx.ravel(),yy.ravel())])
        increments=np.diff(values,axis=1)
        training[str(degree)]=dict(sample_count=len(values), finite=bool(np.all(np.isfinite(values))),
            all_positive=bool(np.all(values>0)), all_strictly_monotone=bool(np.all(increments>0)),
            minimum_energy=float(np.min(values)), minimum_increment=float(np.min(increments)),
            maximum_increment=float(np.max(increments)),
            tensor_vandermonde_condition_number=float(np.linalg.cond(vander)),
            sample_19_to_20_increments=increments[:,19].tolist(),
            scope='Training-node values reconstructed from archived coefficients, not newly queried source spectra.')

        failed_count=0; bad_indices=set(); worst=None; nonpositive=0; total=0
        # Modest deterministic coverage: sampled evidence, never a supremum proof.
        for a in np.linspace(ab[0], ab[1], 81):
            for g in np.linspace(gb[0], gb[1], 81):
                row=record(a,g); total+=1
                nonpositive += not row['positive']
                if not row['monotone']:
                    failed_count+=1
                    bad_indices.update(p['index_left'] for p in row['bad_pairs'])
                if worst is None or row['minimum_increment'] < worst['minimum_increment']:
                    worst=row
        scans[str(degree)]=dict(shape=[81,81], query_count=total, nonpositive_query_count=nonpositive,
            nonmonotone_query_count=failed_count, nonmonotone_sample_fraction=failed_count/total,
            bad_left_indices=sorted(bad_indices), worst_sample=worst,
            scope='Uniform grid evaluation of archived coefficients only; sampled failure coverage, not continuous-domain admission.')
    failure=read(raw/'failure.json')
    inventory={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(raw.iterdir()) if p.is_file()}
    return dict(schema='pysnspd.local_patch.saved_monotonicity_diagnosis.v1',
        status='ARCHIVED_INTERPOLANTS_NOT_ADMITTED', source_hashes_match=source_matches, source_checks=source_checks,
        input_inventory=inventory, script_sha256=sha(__file__), original_failure=failure,
        selected_degree=reg['selected_degree'], completed_original_probe_count=len(read(raw/'probes.json')),
        method='Only evaluation of saved polynomial coefficients; no physical source query, model RHS, trajectory, clipping, coefficient mutation or refit.',
        runtime_seconds=time.perf_counter()-start, registered_probe_and_stencil_replay=replay,
        reconstructed_training_nodes=training, sampled_field_box=scans,
        interpretation=[
            'The actual loop reaches degree4 at probe2 and rejects its Gamma+2e-4 stencil. Offline degree6 also fails at the centre of probe2; skipping the diagnostic degree would not repair candidate admission.',
            'Independent unconstrained polynomial interpolation of energies at each count node does not preserve their ordering between field samples, even when every training spectrum is ordered.',
            'Negative increments exceed floating-point spacing by many orders of magnitude. This is an interpolation-order failure, not a justification for roundoff clipping.',
            'The sampled failures concentrate near the very small-count spectral-edge transition. Their state-count weights are recorded, but small weight does not turn unordered energies into a physical spectrum.',
            'A higher polynomial degree can reduce some inversion magnitudes without removing them or monotonically improving every probe. Neither degree is admitted by these coefficients.'
        ],
        proposed_representation=dict(status='PROPOSAL_ONLY_NOT_IMPLEMENTED_OR_ADMITTED',
            training='At each existing field sample define positive increments d0=E0, dj=Ej-E(j-1); require all positive before taking logs.',
            interpolation='Fit a field interpolant Lj(a,Gamma) to log(dj) rather than fitting each Ej independently.',
            reconstruction='Ej=sum(k<=j) exp(Lk). For t=a or Gamma, dEj/dt=sum(k<=j) exp(Lk)*dLk/dt.',
            second_derivative='d2Ej/(ds dt)=sum(k<=j) exp(Lk)*(dLk/ds*dLk/dt+d2Lk/(ds dt)).',
            preservation='One differentiable reconstructed energy supplies both forces. Positivity and count ordering follow algebraically, without projection, sorting, clipping or population changes.',
            numerical_guards='Reject overflow, underflow to unresolved increments, or lost represented ordering; use stable summation if necessary. Keep all630 state-count nodes and original weights.',
            remaining_checks='Register the replacement before physical tests; compare E, Ea, EGamma, weighted moments and D36 against independent source probes and degree sensitivity under unchanged error budgets.',
            limitation='Log increments can vary sharply near the edge; positivity alone does not guarantee accurate derivatives. Smaller boxes or piecewise interpolation may still be required.'),
        admission=dict(original_failure_preserved=True, degree4_accepted=False, degree6_accepted=False,
                       stage3_closed=False, new_representation_accepted=False))


def notes(report):
    first=report['registered_probe_and_stencil_replay']['4']['first_failure']
    pair=first['bad_pairs'][0]
    selected=report['registered_probe_and_stencil_replay']['6']
    lines=['# Diagnóstico de la monotonía del catálogo local', '',
        'Los dos interpolantes archivados quedan **sin admisión**. Se evaluaron únicamente sus coeficientes guardados: '
        'no se consultó la fuente física, no se cambiaron coeficientes y no se ejecutó dinámica.', '',
        f"El primer fallo del recorrido es el grado4, probe{first['probe_one_based']}, a={first['amplitude']:.8g}, "
        f"Gamma={first['gamma']:.8g} (desplazamiento +2e-4 del stencil). "
        f"La pareja de índices {pair['index_left']}→{pair['index_right']} (base cero) tiene E={pair['energy_left']:.16g} "
        f"y {pair['energy_right']:.16g}: incremento **{pair['signed_increment']:.8g}**. "
        f"Los conteos correspondientes son {pair['count_left']:.8g} y {pair['count_right']:.8g}.", '',
        f"El grado6 también falla: sus centros de probe rechazados son {selected['failing_probe_centres']}. "
        'En particular ya pierde monotonía en el centro del probe2, antes de su stencil. Omitir el grado4 no arreglaría el candidato.', '',
        '| Grado | Probes + stencils no monótonos | Muestras de caja no monótonas | Peor incremento muestreado |',
        '|---:|---:|---:|---:|']
    for degree in ('4','6'):
        replay=report['registered_probe_and_stencil_replay'][degree]; scan=report['sampled_field_box'][degree]
        lines.append(f"| {degree} | {replay['failure_count']}/{replay['query_count']} | {scan['nonmonotone_query_count']}/{scan['query_count']} | {scan['worst_sample']['minimum_increment']:.8g} |")
    lines+=['', 'La malla de diagnóstico81×81 sólo mide fallos muestreados; no certifica toda la caja. '
        'Las energías permanecen positivas en este muestreo. La pérdida de orden se concentra cerca de la transición de borde a conteo muy pequeño. '
        'Los saltos negativos exceden por muchos órdenes de magnitud el espaciado de float64: no son redondeo que deba recortarse.', '',
        'La reconstrucción de los nodos de entrenamiento conserva energías ordenadas. El problema aparece entre nodos: '
        'ajustar cada energía con un polinomio independiente no preserva que las curvas vecinas estén ordenadas. '
        'Aumentar el grado reduce algunos errores, pero no aporta esa garantía.', '',
        'Una representación adecuada para la siguiente iteración es interpolar los logaritmos de los incrementos positivos: '
        '`d0=E0`, `dj=Ej−E(j−1)`, `Lj=log(dj)`, y reconstruir `Ej=sum(k≤j) exp(Lk)`. '
        'Las derivadas se obtienen de esa misma energía: `∂tEj=sum exp(Lk)∂tLk`, con t=a o Gamma. '
        'Esto conserva positividad y orden por construcción sin clipping, ordenamiento posterior ni cambio de poblaciones.', '',
        'La propuesta aún debe registrarse y contrastarse con la fuente física usando el mismo presupuesto de precisión. '
        'Sus derivadas pueden seguir necesitando una caja menor o una representación por tramos cerca del borde. '
        'Mantener la monotonía no basta para admitir momentos, curvaturas ni dinámica. '
        'También deben rechazarse desbordamientos o incrementos que pierdan resolución numérica al reconstruir.', '',
        'Trazabilidad completa, índices, pesos, todos los puntos fallidos, hashes y fórmulas: `patch_monotonicity_diagnosis.json`. '
        'Reproducción sin sobrescribir: `python sandbox/stage3_spatial/coupled_20260923/diagnose_patch_monotonicity.py --verify-only`.', '']
    return '\n'.join(lines)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw',type=Path,default=RAW)
    parser.add_argument('--output-root',type=Path,default=OUT)
    parser.add_argument('--verify-only',action='store_true')
    args=parser.parse_args(); result=diagnose(args.raw)
    if not args.verify_only:
        paths=[args.output_root/'patch_monotonicity_diagnosis.json', args.output_root/'patch_monotonicity_diagnosis.md']
        if any(p.exists() for p in paths):raise FileExistsError('Existing diagnostic outputs preserved')
        args.output_root.mkdir(parents=True,exist_ok=True)
        paths[0].write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
        paths[1].write_text(notes(result),encoding='utf-8')
    print(json.dumps(dict(status=result['status'],source_hashes_match=result['source_hashes_match'],runtime_seconds=result['runtime_seconds'],
        degrees={k:{'probe_failures':result['registered_probe_and_stencil_replay'][k]['failure_count'],
                    'grid_failures':v['nonmonotone_query_count'],'worst_increment':v['worst_sample']['minimum_increment']}
                 for k,v in result['sampled_field_box'].items()}),indent=2))


if __name__=='__main__':main()
