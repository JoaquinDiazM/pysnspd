"""Independent bounded postprocessing of self-consistent thermal core fields.

No solver is called. Gap fields and reconstructed face-current densities are
compared on common coordinates. Iteration histories are optimization histories,
never detector transients. An off-node vortex is located by bilinear zero finding.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from analyze_spatial_energy import nodal_current

DATA = ROOT/'docs/implementation/stage4/self_consistent_review_20260924'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf8'))


def norm(values, mass, mask):
    return float(np.sqrt(np.sum(mass[mask]*abs(values[mask])**2)))


def vortex_position(fields, gap):
    """Locate the complex-gap zero; preserve winding when zero is off nodes."""
    xy, d = fields['coordinates_bar'], fields['d']
    minimum = int(np.argmin(abs(d)))
    if abs(d[minimum]) < 1e-12*gap:
        return dict(coordinates_ell0=xy[minimum].tolist(), method='Resolved node zero within 1e-12 of reference gap',
                    interpolation_residual_relative=float(abs(d[minimum])/gap))
    xs, ys = np.unique(xy[:,0]), np.unique(xy[:,1])
    grid = d.reshape(len(xs),len(ys))
    roots = []
    for i in range(len(xs)-1):
        for j in range(len(ys)-1):
            corners = np.array([grid[i,j],grid[i+1,j],grid[i+1,j+1],grid[i,j+1]])
            winding = float(np.sum(np.angle(np.conj(corners)*np.roll(corners,-1)))/(2*np.pi))
            if abs(winding) < .5:
                continue
            a = corners[0]; b = corners[1]-a; c = corners[3]-a
            e = corners[2]-a-b-c
            uv = np.array([.5,.5])
            for _ in range(12):
                u,v = uv; z = a+b*u+c*v+e*u*v
                derivatives = [b+e*v,c+e*u]
                matrix = np.array([[z.real for z in derivatives],[z.imag for z in derivatives]])
                uv -= np.linalg.solve(matrix,np.array([z.real,z.imag]))
                if abs(z) < 1e-13*gap:
                    break
            u,v = uv; residual = abs(a+b*u+c*v+e*u*v)/gap
            if np.any(uv < -1e-8) or np.any(uv > 1+1e-8) or residual > 1e-10:
                raise ValueError('Bilinear vortex zero not resolved in winding cell')
            roots.append(dict(coordinates_ell0=[float(xs[i]+u*(xs[i+1]-xs[i])),float(ys[j]+v*(ys[j+1]-ys[j]))],
                method='Complex bilinear interpolant zero in a plaquette with unit winding; position is diagnostic, not sub-grid validation',
                interpolation_residual_relative=float(residual), cell_index=[i,j], cell_winding=winding))
    if len(roots) != 1:
        raise ValueError(f'Expected exactly one resolved vortex; got {len(roots)}')
    return roots[0]


def radial_half_gap(fields, gap):
    xy, d = fields['coordinates_bar'], fields['d']
    # Positive x-axis intersects an exactly centered radial vortex. Linear
    # interpolation removes binary area-threshold staircasing, not mesh error.
    line = (xy[:,1] == 0) & (xy[:,0] >= 0)
    x, amplitude = xy[line,0], abs(d[line])/gap
    index = np.flatnonzero(amplitude >= .5)
    if not len(index) or index[0] == 0:
        raise ValueError('No radial half-gap crossing')
    k = int(index[0])
    return float(x[k-1]+(.5-amplitude[k-1])*(x[k]-x[k-1])/(amplitude[k]-amplitude[k-1]))


def comparison(low, high):
    side_low = len(np.unique(low['coordinates_bar'][:,0]))
    side_high = len(np.unique(high['coordinates_bar'][:,0]))
    factor = (side_high-1)//(side_low-1)
    if factor < 1 or (side_high-1) % (side_low-1):
        raise ValueError('Only nested identical-domain Cartesian comparisons supported')
    def sample(value):
        return value.reshape(side_high,side_high)[::factor,::factor].ravel()
    sampled_xy = high['coordinates_bar'].reshape(side_high,side_high,2)[::factor,::factor].reshape(-1,2)
    if not np.array_equal(low['coordinates_bar'],sampled_xy):
        raise ValueError('Physical coordinates differ')
    mask = np.linalg.norm(low['coordinates_bar'],axis=1) <= 4
    mass = low['area_weights']
    hd, hj = sample(high['d']), sample(nodal_current(high))
    lj = nodal_current(low)
    return dict(gap_complex_relative_L2=norm(low['d']-hd,mass,mask)/norm(hd,mass,mask),
        gap_magnitude_relative_L2=norm(abs(low['d'])-abs(hd),mass,mask)/norm(hd,mass,mask),
        gap_maximum_absolute_change=float(np.max(abs(low['d'][mask]-hd[mask]))),
        current_density_relative_L2=norm(lj-hj,mass,mask)/norm(hj,mass,mask),
        current_density_maximum_absolute_change=float(np.max(abs(lj[mask]-hj[mask]))),
        domain='Area-weighted r<=4ell0, including origin, common coarse coordinates',
        current='Cartesian nodal diagnostic from adjacent oriented link currents divided by dual-face width; raw unequal-mesh link currents are not compared',
        uncertainty='Observed paired change of separately relaxed finite-cutoff states; not an asymptotic error bound, and stopping residual remains finite')


def process(raw):
    identity, summary, receipt, plan = (read(raw/name) for name in ('identity.json','summary.json','extraction_receipt.json','executed_plan.json'))
    if sha(raw/'identity.json') != receipt['identity_sha256'] or sha(raw/'summary.json') != receipt['summary_sha256']:
        raise ValueError('Downloaded run identity or summary changed')
    if sha(raw/'executed_plan.json') != identity['plan_sha256']:
        raise ValueError('Downloaded executed plan changed')
    for name, expected in identity['sources'].items():
        if sha(ROOT/name) != expected:
            raise ValueError('Source mismatch: '+name)
    if sha(HERE/'extract_self_consistent_results.py') != receipt['extraction_script_sha256']:
        raise ValueError('Extraction source differs')
    reference = ROOT/'docs/implementation/stage4/followup_20260923/raw/stage4_radial_reference_20260923/identity.json'
    if sha(reference) != identity['radial_reference_identity_sha256']:
        raise ValueError('Radial reference identity differs')
    gap = identity['gap_reference']
    histories = {case['id']:[] for case in plan['cases']}
    history_hashes = {}
    for path in sorted((raw/'histories').glob('*.json')):
        record = read(path); histories[record['case_id']].append(record)
        history_hashes[path.name] = sha(path)
    arrays, cases, history_table = {}, {}, {}
    actual_jobs = 0
    for case in plan['cases']:
        cid = case['id']; final = summary['cases'][cid]; extraction = receipt['final_cases'][cid]
        path = raw/extraction['compact_file']
        if sha(path) != extraction['compact_sha256']:
            raise ValueError('Compact file differs')
        with np.load(path) as values:
            fields = {name:values[name].copy() for name in values.files}
        if any(np.any(~np.isfinite(value)) for value in fields.values()):
            raise ValueError('Nonfinite compact fields')
        arrays[cid] = fields
        items = histories[cid]
        if [item['completed_sweeps'] for item in items] != list(range(1,final['completed_sweeps']+1)):
            raise ValueError('Incomplete sweep history')
        if items[-1] != final:
            raise ValueError('Summary/checkpoint differs')
        previous, history_table[cid] = None, []
        gap_drops, spectral_drops, energy_drops = [],[],[]
        max_mode_residual = 0.
        for record in items:
            metric = record['metrics']; sweep = record['completed_sweeps']
            modes = metric['modes']
            if len(modes) != case['matsubara_count'] or sorted(m['n'] for m in modes) != list(range(case['matsubara_count'])):
                raise ValueError('Missing or duplicate spectral mode')
            max_mode_residual = max(max_mode_residual,max(m['residual'] for m in modes))
            if sweep > 2:
                actual_jobs += len(modes)
            drop = metric['gap_block']['observed_energy_change']
            gap_drops.append(drop)
            if drop > 2e-11 or abs(drop-metric['gap_block']['predicted_energy_change']) > 2e-11:
                raise ValueError('Gap block failed its exact energy identity')
            if previous is not None:
                change = metric['energy']-previous['energy_after_gap_block']
                spectral_drops.append(change)
                energy_drops.append(metric['energy']-previous['energy'])
                if change > 2e-10 or abs(change-metric['spectral_block_energy_change']) > 2e-10:
                    raise ValueError('Spectral block failed energy descent')
            previous = metric
            history_table[cid].append(dict(sweep=sweep, energy=metric['energy'],
                energy_after_gap_block=metric['energy_after_gap_block'], core_mass_rms_relative=metric['core_mass_rms_relative'],
                global_free_maximum_relative=metric['global_free_maximum_relative'], gap_energy_change=drop,
                spectral_block_energy_change=metric['spectral_block_energy_change']))
        metric = final['metrics']
        if max_mode_residual > plan['spectral_tolerance'] or not metric['core_criterion_met']:
            raise ValueError('Declared residual admission failed')
        if metric['core_mass_rms_relative'] > plan['core_mass_rms_relative_tolerance']:
            raise ValueError('Declared core tolerance failed')
        winding = metric['observables']['winding_on_square_near_two_ell0']
        if winding is None or abs(winding-1) > 1e-12:
            raise ValueError('Vortex winding not preserved')
        current = nodal_current(fields)
        core = np.linalg.norm(fields['coordinates_bar'],axis=1) <= plan['core_radius_ell0']
        phase_mask = core & (abs(fields['d'])>1e-8*gap) & (abs(fields['f_lowest'])>1e-8)
        phase = np.angle(fields['f_lowest'][phase_mask]*np.conj(fields['d'][phase_mask]))
        cases[cid] = dict(completed_sweeps=final['completed_sweeps'], core_mass_rms_relative=metric['core_mass_rms_relative'],
            global_free_maximum_relative=metric['global_free_maximum_relative'], core_maximum_relative=metric['core_maximum_relative'],
            energy_initial=items[0]['metrics']['energy'], energy_final=metric['energy'],
            energy_change=metric['energy']-items[0]['metrics']['energy'],
            maximum_gap_energy_change=max(gap_drops), maximum_spectral_energy_change=max(spectral_drops),
            maximum_accepted_state_energy_change=max(energy_drops), monotonic_energy=True,
            maximum_spectral_residual_all_sweeps=max_mode_residual,
            maximum_recomputed_final_spectral_residual=extraction['checks']['maximum_recomputed_spectral_residual'],
            winding=winding, vortex_position=vortex_position(fields,gap),
            observables=metric['observables'],
            current_density_core_L2=norm(current,fields['area_weights'],core),
            gap_core_L2=norm(fields['d'],fields['area_weights'],core),
            lowest_spectral_phase_difference_max_rad=float(np.max(abs(phase))),
            coherence_checks=extraction['checks'])
        if case['profile']=='radial_vortex':
            cases[cid]['half_gap_radius_positive_x_ell0']=radial_half_gap(fields,gap)
    if actual_jobs != summary['completed_jobs']:
        raise ValueError('Completed job count does not match histories')
    comparisons = {}
    for label,low,high in (('cutoff_N128_to_N256','radial_65_N128','radial_65_N256'),
                           ('mesh_65_to_129','radial_65_N256','radial_129_N256')):
        comparisons[label]=comparison(arrays[low],arrays[high])
        comparisons[label].update(low=low, high=high,
            half_gap_radius_relative_change=abs(cases[high]['half_gap_radius_positive_x_ell0']-cases[low]['half_gap_radius_positive_x_ell0'])/cases[high]['half_gap_radius_positive_x_ell0'],
            energy_relative_change=abs(cases[high]['energy_final']-cases[low]['energy_final'])/abs(cases[high]['energy_final']))
    return dict(schema='pysnspd.stage4.self_consistent_postprocessing.v1',
        provenance=dict(identity_sha256=sha(raw/'identity.json'),summary_sha256=sha(raw/'summary.json'),
            extraction_receipt_sha256=sha(raw/'extraction_receipt.json'),analysis_script_sha256=sha(__file__),
            sources_matching_identity=identity['sources'], plan_sha256=identity['plan_sha256'],
            full_remote_checkpoints_verified=len(receipt['checkpoint_files_verified']),
            full_remote_checkpoint_bytes=sum(item['bytes'] for item in receipt['checkpoint_files_verified']),
            downloaded_compact_map_sha256={cid:item['compact_sha256'] for cid,item in receipt['final_cases'].items()},
            history_sha256=history_hashes,history_records=len(history_hashes),
            full_spectral_arrays='Remain on Geminga; hashes verified there, final coherent fields recalculated there without solves'),
        execution=dict(runtime_seconds=summary['runtime_seconds'],completed_jobs=summary['completed_jobs'],
            all_core_criteria_met=summary['all_core_criteria_met'],process_budget=identity['budget'],
            extraction_seconds=receipt['runtime_seconds'],physical_time_steps=0,new_postprocessing_solves=0),
        cases=cases,histories=history_table,comparisons=comparisons,
        decision=dict(static_control_accepted=True,additional_static_rerun_required=False,stage4_complete=False,production_changed=False,
            accepted='Finite-cutoff thermal stationary core on the declared fixed-boundary Cartesian graph, to its core RMS criterion',
            next='Proceed to explicit retarded spectral/kinetic and dissipation matching; another tighter static refinement is not justified by these paired changes',
            limitations=['Core mass RMS <=0.1% is not a rigorous error bound or full-domain maximum criterion',
                'Finite-N paired differences remain; no infinite-frequency remainder bound was inferred',
                'Pinned outer winding does not establish vortex stability, barrier or nucleation',
                'Bilinear off-node vortex location is a visualization diagnostic, not a separately validated sub-grid solver',
                'No nonthermal kinetic closure, retarded branch, photon transfer, physical time or detector circuit transient is admitted']))


def markdown(result):
    lines=['# Núcleo térmico autoconsistente: revisión de resultados','',
        'Los cuatro casos alcanzaron el criterio anunciado. Se acepta este control estático y no se solicita otra corrida estática más estricta. La etapa 4 continúa abierta para el enlace espectral retardado, cinético y disipativo.','',
        f"La continuación terminó en {result['execution']['runtime_seconds']:.2f} s y resolvió {result['execution']['completed_jobs']} consultas espectrales. Se comprobaron los SHA-256 de {result['provenance']['full_remote_checkpoints_verified']} checkpoints completos y la coherencia algebraica de los cuatro estados finales. El postproceso no volvió a resolver espectros.",'',
        '| Caso | Barridos totales | Residuo RMS del núcleo | Máximo libre | Energía final |','|---|---:|---:|---:|---:|']
    for cid,case in result['cases'].items():
        lines.append(f"| {cid} | {case['completed_sweeps']} | {100*case['core_mass_rms_relative']:.5f}% | {100*case['global_free_maximum_relative']:.5f}% | {case['energy_final']:.8f} |")
    lines += ['', 'Todos los bloques redujeron la misma energía y conservaron el arrollamiento de fase +1. Se verificó el estado final con su propio espectro: el campo next_d es una propuesta distinta y no se mezcló con el espectro anterior para fingir convergencia. Los máximos libres se muestran aparte porque el criterio era RMS en r≤4ℓ₀.','',
        '## Sensibilidad útil para continuar','',
        '| Cambio | Diferencia del gap (L2) | Diferencia de densidad de corriente (L2) | Cambio del radio de medio gap |','|---|---:|---:|---:|']
    for label,case in result['comparisons'].items():
        lines.append(f"| {label} | {100*case['gap_complex_relative_L2']:.4f}% | {100*case['current_density_relative_L2']:.4f}% | {100*case['half_gap_radius_relative_change']:.4f}% |")
    asym=result['cases']['asymmetric_65_N256']; pos=asym['vortex_position']['coordinates_ell0']
    lines += ['', 'Las comparaciones usan estados relajados separadamente, áreas duales y coordenadas comunes. La corriente de cada enlace se divide por el ancho de su cara dual antes de reconstruir una densidad nodal: comparar directamente corrientes integradas de mallas distintas produciría un cambio geométrico espurio. Las diferencias son sensibilidades observadas, no cotas rigurosas de error.','',
        f"El vórtice asimétrico se desplaza a aproximadamente ({pos[0]:.6f}, {pos[1]:.6f})ℓ₀. Su arrollamiento permanece +1. El mínimo nodal de |Δ| no debe interpretarse como desaparición del vórtice: el cero queda entre nodos; su ubicación se obtuvo del interpolante bilineal complejo.",'',
        '## Decisión','',
        'El problema estático que motivó esta corrida queda resuelto al nivel declarado. No se vuelve a ajustar el cierre local ni se exige convergencia arbitrariamente más dura. Se continúa con la correspondencia entre este oráculo térmico y la descripción espectral retardada/cinética que necesita el transiente. La suma de Matsubara representa equilibrio; no autoriza sustituir las poblaciones no térmicas por una temperatura ni extender esta parametrización compleja por sustitución ingenua de frecuencia.','',
        'No se ha validado aún el transiente de un fotón, la transferencia de energía inicial, la latencia de Korzh ni el circuito dinámico. Los índices de barrido son iteraciones de minimización, nunca tiempos físicos.','']
    return '\n'.join(lines)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw',type=Path,default=DATA/'raw')
    parser.add_argument('--output-directory',type=Path,default=DATA)
    args=parser.parse_args()
    result=process(args.raw)
    args.output_directory.mkdir(parents=True,exist_ok=True)
    for name,body in (('analysis.json',json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+'\n'),('analysis.md',markdown(result))):
        with (args.output_directory/name).open('x',encoding='utf8') as stream:
            stream.write(body)
    print(json.dumps(dict(status='STATIC_CONTROL_ACCEPTED_NO_EXTRA_RERUN',cases=len(result['cases']),
        comparisons=result['comparisons']),ensure_ascii=False))


if __name__=='__main__':
    main()
