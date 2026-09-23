"""Complete all127 range records with their coupled restrictions and use limits.

This is a documentary scenario/constraint ledger, never solver configuration.
"""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'docs/implementation/stage3_5'
OUT=BASE/'assessment_r2_20260923'

FAMILY_USE={
 'material':('Fuerza, corriente, difusión y energía','No mezclar escalas de películas o ajustes distintos.','material/material_assessment.md'),
 'bath':('Equilibrio inicial, poblaciones y movilidad','Cambiar el baño requiere recalcular la rama y los estados; no escalar una curva previa.','../research_20260923/material_kinetics.md'),
 'geometry':('Hotbelt, barrera, difusión y señal','Influencia observable de terminales o variación transversal incompatible con una continuación1D.','../research_20260923/geometry_numerics.md'),
 'photon':('Preparación, formación de hotbelt y latencia por color','Reparto, pérdidas, reloj o perfil no respaldados por una transferencia caracterizada.','cascade/cascade_assessment.md'),
 'electronic_state':('Fuerza y transporte no térmicos','Salida del soporte físico, del catálogo o de una rama espectral causal.','stage4_readiness.md'),
 'condensate':('Energía, núcleo, fuerzas y primeros eventos de fase','Símbolo principal indefinido, amplitud no cubierta o sensibilidad física no resuelta a delta.','stage4_readiness.md'),
 'kinetic':('Reparto de energía y velocidad de respuesta','Tratar tiempos efectivos como tasas microscópicas o usar soportes insuficientes.','material/material_assessment.md'),
 'boundaries':('Corriente, intercambio de energía y contaminación desde extremos','Borde sin trabajo coherente o interfaz sin transporte cinético a igual energía.','stage4_readiness.md'),
 'circuit':('Corrientes y voltaje de carga','Cambiar dominio sin actualizar partición inductiva o confundir CM con toda la lectura experimental.','../research_20260923/circuit_observables.md'),
 'catalogue':('Energía y derivadas electrónicas','Extrapolar estados o energía fuera de soporte, o usar aceleración lejos de su dominio probado.','stage4_readiness.md'),
 'spatial_numerics':('Núcleo, perfil inicial y observables espaciales','Estructura más pequeña que la resolución contrastada o cambio de evento al refinar.','../research_20260923/geometry_numerics.md'),
 'time_numerics':('Cruces y balances acumulados','Confundir muestreo con integración o adoptar paso sin control del operador acoplado.','stage4_readiness.md'),
 'observation':('Latencia relativa, hotbelt y recuperación','Equiparar cruce determinista, máximo de IRF y jitter o tratar censura como no detección.','../research_20260923/circuit_observables.md'),
 'provenance':('Reproducibilidad e interpretación de todos los resultados','Entrada sin unidades/procedencia o corrida que modifica resultados anteriores.','README.md'),
 'experiment_reference':('Comparación con una muestra y lectura identificables','Atribuir a Korzh una preparación/material/instrumento no identificado con esa muestra.','../research_20260923/circuit_observables.md'),
}

def main():
    inventory=json.loads((BASE/'parameter_inventory.json').read_text(encoding='utf8'))
    prior=json.loads((BASE/'research_20260923/range_ledger.json').read_text(encoding='utf8'))
    original={x['id']:x for x in inventory['parameters']}
    rows=[]
    for item in prior['parameters']:
        i=item['id'];p=original[i];obs,trigger,review=FAMILY_USE[item['family']]
        row=dict(item)
        row.update(sample_or_scenario='K20 fitted material reference with memory circuit; source-dependent alternatives stay separate',
          mathematical_domain=p['mathematical_domain'],
          coupled_constraints=p['couplings_and_restrictions'],
          affected_observables=obs,extrapolation_trigger=trigger,
          next_evidence=p['research_required'],
          uncertainty='Not assigned a statistical interval; declared reference/proposal/unknown and coupled restrictions apply.',
          margin_rationale=item['rationale'],review='../research_20260923/'+item['review'],assessment_review=review,
          implementation_activation=False)
        if i in ['retention','Ein','photon_sigma','photon_position','spectral_injection','handoff_time']:
            row['next_evidence']='Transfer state for each color: electron/phonon profiles and retained energies, prior losses and clock. A20 historical extraction is not a measured Korzh input.'
        if i=='photon_sigma':
            row['status']='OPEN_AFTER_DOCUMENTARY_CHARACTERIZATION'
            row['value_or_proposal']={'value':None,'units':'nm','range':None}
            row['margin_rationale']='User chose characterization first. No5–20nm or1.4–1.9nm scenario adopted; one total-energy percentile cannot define the phonon Gaussian.'
        if i=='kwt_times':
            row['implementation_constraint']='Experimental KWTMobility hardcodes0.50/2.47ps; a config-only edit does not select another pair. Stage4 must expose pairs consistently before comparison.'
        if i=='taukin':
            row['units']='normalized solver time; multiply by t_ref to obtain physical time'
            row['value_or_proposal']={'synthetic_normalized':0.7,'physical_ps':'0.7 * t_ref_ps','physical_range':None}
            row['margin_rationale']='Correction to r1 units: CoupledCellSystem evolves normalized time. The inherited 0.7 equals 0.7ps only for t_ref=1ps; no physical BGK time is identified.'
            row['assessment_review']='stage4_readiness.md'
        if i=='spatial_mesh':
            row['margin_rationale']+=' Compact-source mesh figures remain cost illustrations while physical width is open.'
        rows.append(row)
    assert len(rows)==127 and {r['id'] for r in rows}==set(original)
    result=dict(schema='pysnspd.stage3_5.use_domain.r2',status='REVIEWED_CONSTRAINTS_AND_EXPLICIT_UNKNOWN_INPUTS',
      base_inventory_sha256=hashlib.sha256((BASE/'parameter_inventory.json').read_bytes()).hexdigest(),
      r1_ledger_sha256=hashlib.sha256((BASE/'research_20260923/range_ledger.json').read_bytes()).hexdigest(),
      parameters=rows,count=len(rows),families=len(FAMILY_USE),all_solver_parameters_unchanged=True,
      code_parameters_changed=False,photon_width_nm=None,
      physical_confidence_box_admitted=False,new_transients=0,
      interpretation='The domain is a set of coupled restrictions per observable, not127 independent scalar intervals.',
      source_tables=[{'path':'../research_20260923/'+name} for name in ['material_kinetics_sources.json','photon_sources.json','gaussian_sources.json','circuit_sources.json']])
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'use_domain.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(json.dumps({'entries':len(rows),'families':len(FAMILY_USE),'new_physics':False}))

if __name__=='__main__':main()
