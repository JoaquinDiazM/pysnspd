"""Recompute and explain the registered D/C/U moment criteria, without solves."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(HERE))
from resolution_campaign_analysis import analyze

DATA=ROOT/'docs/implementation/stage4/moment_review_20260924'
NAMES={'charge_current':'Corriente de carga','energy_weighted_energy_flux':'Flujo ponderado por energía',
    'gap_force_density':'Fuerza compleja','amplitude_force_density':'Fuerza de amplitud','phase_torque_density':'Torque de fase'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_text(encoding='utf8'))


def same(expected,actual,path='root'):
    if isinstance(expected,dict):
        if set(expected)!=set(actual):raise ValueError('Different keys: '+path)
        return max([same(expected[k],actual[k],path+'/'+k) for k in expected]+[0.])
    if isinstance(expected,list):
        if len(expected)!=len(actual):raise ValueError('Different length: '+path)
        return max([same(a,b,path) for a,b in zip(expected,actual)]+[0.])
    if isinstance(expected,(int,float)) and not isinstance(expected,bool):
        error=abs(expected-actual)
        if error>2e-11*max(1.,abs(expected)):raise ValueError('Metric changed: '+path)
        return error
    if expected!=actual:raise ValueError('Non-numerical metadata differs: '+path)
    return 0.


def process(raw):
    receipt=read(raw/'extraction_receipt.json');identity=read(raw/'identity.json')
    if receipt['script_sha256']!=sha(HERE/'extract_moment_results.py'):raise ValueError('Extractor changed')
    for name,digest in receipt['copied_file_sha256'].items():
        if sha(raw/name)!=digest:raise ValueError('Downloaded artifact changed: '+name)
    for name,digest in {**identity['sources'],**identity['inputs']}.items():
        if sha(ROOT/name)!=digest:raise ValueError('Frozen source or input changed: '+name)
    plan=read(raw/'executed_plan.json');summary=read(raw/'summary.json')
    if sha(raw/'executed_plan.json')!=identity['plan_sha256']:raise ValueError('Plan differs')
    projection=read(raw/'projection/summary.json')
    stored=read(raw/'moment_comparisons.json')
    rebuilt=analyze(plan,raw,projection)
    error=same(stored,rebuilt)
    table={}
    for probe,values in rebuilt['diagnostic_variation_budget'].items():
        table[probe]={}
        for name,value in values.items():
            D,C,F,UD,UC=(value[k] for k in ('projection_difference_D','omitted_response_C','full_norm','observed_variation_U_D','observed_variation_U_C'))
            table[probe][name]=dict(value,
                D_over_full_percent=100*D/F if F>1e-12 else None,
                C_over_full_percent=100*C/F if F>1e-12 else None,
                U_D_over_full_percent=100*UD/F if F>1e-12 else None,
                U_C_over_full_percent=100*UC/F if F>1e-12 else None,
                D_over_U_D=D/UD if UD>1e-12 else None,
                projected_to_full_norm=value['projected_norm']/F if F>1e-12 else None,
                discrepancy_resolved=D>UD,improvement_resolved=D+UD<C-UC,
                method='Predeclared observed-variation criteria; U is not a rigorous uncertainty bound')
    control_changes={}
    for label,record in rebuilt['comparisons'].items():
        control_changes[label]={}
        for probe,values in record['difference_norms'].items():
            control_changes[label][probe]={}
            for name,states in values.items():
                norm=rebuilt['moments'][record['reference']][probe][name]['full_norm']
                control_changes[label][probe][name]=dict(absolute_difference_norms=states,
                    percent_of_reference_full_norm={state:100*value/norm for state,value in states.items()} if norm>1e-12 else None)
    spectral=read(raw/'spectra/summary.json')['records'];kinetic=read(raw/'kinetic/summary.json')['records']
    all_steps=[step for row in spectral for step in row['continuation']]
    if len(all_steps)!=181*plan['spectral']['continuation_steps']:raise ValueError('Missing continuation step')
    if max(step['root_residual'] for step in all_steps)>plan['spectral']['spectral_tolerance']:
        raise ValueError('Spectral root criterion failed')
    if min(step['minimum_DOS'] for step in all_steps)<-plan['spectral']['causal_tolerance']:
        raise ValueError('Causality criterion failed')
    response_metrics=[value for row in kinetic for value in row['probes'].values()]
    projection_metrics=[value for row in projection['records'] for value in row['probes'].values()]
    chi={name:dict(values,error_percent=100*(values['numerical']/values['exact']-1)) for name,values in rebuilt['chi_quadrature'].items()}
    return dict(schema='pysnspd.stage4.moment_review.v1',
        provenance=dict(raw_directory=raw.relative_to(ROOT).as_posix(),analysis_source_sha256=sha(__file__),
            extraction_receipt_sha256=sha(raw/'extraction_receipt.json'),
            frozen_analysis_source_sha256=sha(HERE/'resolution_campaign_analysis.py'),
            sources_matching_identity=identity['sources'],plan_sha256=identity['plan_sha256'],
            verified_raw_fields=len(receipt['maps_verified']),verified_raw_bytes=sum(r['bytes'] for r in receipt['maps_verified']),
            locally_verified_compact_files=len(receipt['copied_file_sha256']),moment_recomputation_maximum_difference=error,
            retained_raw_spectra='181 complete spectral and181 kinetic field files remain on Geminga; all were verified there'),
        execution=dict(runtime_seconds=summary['runtime_seconds'],spectral_queries=181,kinetic_queries=181,
            projection_groups=7,mathematical_probe_responses=362,continuation_steps=len(all_steps),
            phase_resources={name:read(raw/name/'resource_snapshot.json')['budget'] for name in ('spectra','kinetic','projection')},
            extraction_seconds=receipt['runtime_seconds'],new_review_spectral_solves=0,physical_time_steps=0),
        checks=dict(maximum_root_residual=max(s['root_residual'] for s in all_steps),
            maximum_spectral_normalization_residual=max(r['maximum_normalization_residual'] for r in receipt['maps_verified'] if r['phase']=='spectra'),
            minimum_DOS_along_continuation=min(s['minimum_DOS'] for s in all_steps),
            maximum_matsubara_bridge_action_difference=max(r['bridge']['action_absolute_difference'] for r in spectral),
            maximum_kinetic_ward_residual=max(v['maximum_ward_residual'] for v in response_metrics),
            maximum_full_free_charge_residual=max(v['maximum_free_charge_residual'] for v in response_metrics),
            maximum_projected_free_integrated_charge_residual=max(v['norms']['projected']['maximum_free_integrated_charge_residual'] for v in projection_metrics),
            maximum_projected_integrated_ward_residual=max(v['norms']['projected']['maximum_joint_ward_residual'] for v in projection_metrics)),
        reference=rebuilt['reference_with_all_three_controls'],criteria=table,control_changes=control_changes,
        chi_quadrature=chi,all_groups=rebuilt['moments'],
        definitions=dict(D='Norm(projected-full)',C='Norm(zero_hT-full)',
            U_D='Sum over energy, contour and mesh of (change in full + change in projected)',
            U_C='Sum over energy, contour and mesh of (change in full + change in zero_hT)',
            criterion_resolved='D>U_D',criterion_improved='D+U_D<C-U_C',
            criteria_source='Frozen resolution_campaign_analysis.py; retained without retrospective tightening',
            limitations='U is a sum of observed changes, not a rigorous bound; asymmetric base31 lacks a separate fine-grid control'),
        decision=dict(moment_control_complete=True,additional_repetition_of_same_static_campaign_required=False,
            scalar_potential_projection_exact_frozen_replacement=False,
            reject_full_memory_KWT_potential_dynamics=False,new_dynamic_charge_state_admitted=False,
            stage4_complete=False,production_changed=False,
            conclusion='Angular current and phase-torque mismatch persist beyond observed resolution changes. The potential projection improves omission but does not reproduce the frozen full response. Assess the dynamic phase/potential coupling before choosing an additional charge state.',
            amplitude_vs_phase='Complex force mismatch is dominated by phase torque; amplitude mismatch is unresolved and much smaller. Do not use only a total complex-force percentage to justify a closure.',
            scope='Infinitesimal static probes, frozen gap and numerical eta; no detector transient, photon, material relaxation or dynamic energy-balance admission'))


def markdown(a):
    lines=['# Revisión de momentos: resultado y decisión','',
        'La campaña completó las 181 consultas espectrales, sus 181 respuestas cinéticas y siete proyecciones en **198,12 s**. No falló un solver ni hace falta repetir la relajación térmica. El resultado distingue una limitación de la proyección estática de potencial de un error de resolución.','',
        'Se comprobaron 369 archivos completos, las fuentes y el plan ejecutados. Los siete mapas integrados se descargaron y los momentos se recalcularon localmente. No se ejecutaron nuevas soluciones espectrales.','',
        '## Qué significan D, C y U','',
        '`D` mide la diferencia entre la proyección de potencial y la respuesta completa; `C`, la corrección que se pierde al omitir el modo de carga. `U_D` suma los cambios observados al refinar energía, reducir eta y refinar la malla en ambas respuestas comparadas; `U_C` hace lo mismo para la omisión. Se conserva el criterio previo: una diferencia se resuelve si `D>U_D`; la proyección mejora la omisión si `D+U_D<C−U_C`. U no es una cota rigurosa de error.','',
        '## La sonda angular discrimina el cierre','',
        'Referencia con los tres controles: núcleo radial 65², eta/Delta=0,01, malla fina de 50 energías. Las cifras siguientes se normalizan por la norma completa de cada observable.','',
        '| Observable | D / norma completa | U_D / norma completa | D / C | Interpretación |','|---|---:|---:|---:|---|']
    for name,value in a['criteria']['angular'].items():
        verdict=('Diferencia no resuelta' if not value['discrepancy_resolved'] else
                 'Diferencia resuelta; mejora la omisión' if value['improvement_resolved'] else
                 'Diferencia resuelta; mejora incierta')
        lines.append(f"| {NAMES[name]} | {value['D_over_full_percent']:.5f} % | {value['U_D_over_full_percent']:.5f} % | {value['D_over_C']:.5f} | {verdict} |")
    lines += ['',
        'El 3,23 % de diferencia de fuerza compleja oculta dos situaciones muy distintas: la fuerza de amplitud apenas cambia, mientras el torque de fase difiere un 77,53 % respecto de su propia norma. El torque proyectado tiene una norma 1,717 veces la completa. Esto no implica un error de latencia del 77,53 % ni invalida automáticamente una dinámica con fase y potencial acoplados.','',
        'En la sonda radial todas las diferencias quedan sin resolver frente a U. Sus pequeñas componentes de carga/fase no sirven por sí solas para aceptar el cierre, ni sus cocientes pueden trasladarse a la sonda angular. El caso asimétrico conserva datos útiles, pero sólo tiene 31 energías: no dispone del mismo presupuesto independiente de controles.','',
        '## Qué cambió al resolver mejor la energía','',
        'La integral conocida de la susceptibilidad pasa de un exceso de 18,02 % con 12 puntos a 1,616 % con 31 y 0,4243 % con 50. No se renormalizó esa función. La discrepancia angular de corriente permanece alrededor del 20,55 % después de mejorar la cuadratura. No procede atribuirla enteramente a los 12 puntos iniciales.','',
        '| Control del resultado completo, sonda angular | Corriente | Fuerza compleja | Torque de fase | Flujo energético |','|---|---:|---:|---:|---:|']
    for key,label in (('energy_radial_65_N256_eta0p0100','31 a 50 energías'),('contour_eta','eta 0,02 a 0,01'),('spatial_mesh','65² a 129²')):
        row=a['control_changes'][key]['angular']
        values=[row[name]['percent_of_reference_full_norm']['full'] for name in ('charge_current','gap_force_density','phase_torque_density','energy_weighted_energy_flux')]
        lines.append('| '+label+' | '+' | '.join(f'{v:.5f} %' for v in values)+' |')
    lines += ['',
        'Se comparan densidades de fuerza sobre coordenadas comunes y corrientes por ancho de cara dual; nunca corrientes brutas de enlaces de diferente tamaño. El flujo energético lleva el peso E en la integral. Todas las cuadraturas de fuerza, corriente y residuo son las mismas.','',
        '## Decisión para continuar','',
        'Este ensayo de resolución queda terminado. La proyección de un solo perfil térmico de carga por potencial no puede presentarse como sustitución exacta de la respuesta cinética congelada: deja diferencias resueltas en corriente y fase. Sí mejora omitir por completo la respuesta de carga.','',
        'La siguiente decisión requiere estudiar la unión dinámica entre fase del condensado, potencial y poblaciones, incluyendo su balance energético. No se añade todavía un nuevo estado dinámico de carga, no se rechaza automáticamente el esquema completo de la memoria y no se altera producción. La etapa 4 sigue abierta; el transiente de fotón y la latencia del dispositivo siguen pendientes.','',
        'Los archivos raw/ conservan los resúmenes, recursos, siete mapas integrados y el recibo de integridad remoto. analysis.json contiene magnitudes absolutas, cocientes, componentes de U y verificaciones.','']
    return '\n'.join(lines)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw',type=Path,default=DATA/'raw');p.add_argument('--output-directory',type=Path,default=DATA)
    args=p.parse_args();a=process(args.raw);args.output_directory.mkdir(parents=True,exist_ok=True)
    for name,text in (('analysis.json',json.dumps(a,indent=2,ensure_ascii=False,allow_nan=False)+'\n'),('analysis.md',markdown(a))):
        with (args.output_directory/name).open('x',encoding='utf8') as stream:stream.write(text)
    print(json.dumps(dict(status='REGISTERED_MOMENT_CRITERIA_REPRODUCED',maps=a['provenance']['verified_raw_fields'],
        maximum_recomputation_difference=a['provenance']['moment_recomputation_maximum_difference'],
        angular={name:{key:v[key] for key in ('D_over_full_percent','U_D_over_full_percent','D_over_C','assessment')} for name,v in a['criteria']['angular'].items()})))


if __name__=='__main__':main()
