"""Publish the audited scoped result without issuing global stage admission."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'docs/implementation/stage2/time_pass_20260922'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
load=lambda p:json.loads(p.read_text(encoding='utf-8'))
write=lambda p,v:p.write_text(json.dumps(v,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
admission=ROOT/'docs/implementation/stage2/stage2_admission.json'
archive=DATA/'admission_before_time_pass.json'
if not archive.exists():archive.write_bytes(admission.read_bytes())
old=load(archive);current=dict(old)
audit=load(DATA/'one_cell_audit.json');assert audit['checks_passed']==audit['checks_total']==85
fields=load(DATA/'one_ssp_admitted_fields_20260922.json')
support=load(DATA/'one_ssp_admitted_support_20260922.json')
assert fields['status']==support['status']=='PASS'
inventory=[]
for line in (DATA/'remote_sha256.txt').read_text().splitlines():
    digest,remote=line.split(maxsplit=1)
    local=DATA/('raw/'+Path(remote).name if '/stage2_ssp_time_' in remote else Path(remote).name)
    assert sha(local)==digest,(remote,'remote/local mismatch')
    inventory.append(dict(remote=remote,path=local.relative_to(ROOT).as_posix(),bytes=local.stat().st_size,sha256=digest))
write(DATA/'import_inventory.json',dict(status='PASS',count=len(inventory),files=inventory))
current.update(status='ONE_CELL_TIME_ADMITTED_GLOBAL_STAGE2_PENDING',stage2_status='NOT_CLOSED',
    numerical_admission=False,production_promotion=False,
    recommendation='Retain the admitted one-cell SSP trajectories. Run the registered two-cell temporal block manually before selecting full energy-mesh calculations. Stage3 is prepared, not started.',
    latest_evidence='docs/implementation/stage2/time_pass_20260922/one_cell_audit.json',
    previous_checkpoint=dict(path=archive.relative_to(ROOT).as_posix(),sha256=sha(archive)))
for gate in current['gates']:
    if gate['id']=='SSP_temporal_admission':
        gate.clear();gate.update(id='SSP_one_cell_temporal_admission',status='PASS_ONE_CELL_SCOPE',
            admitted_measured_steps=[160,320,640],reference_steps=1280,
            maximum_error_by_steps={'160':5.415702718006263e-5,'320':1.5699838897606807e-5,'640':3.7418095755529397e-6},
            scope='630 electrons,1025 phonons,duration2,isolated driven one-cell scenario; no fine-grid or two-cell inference.',
            all_event_limiters_inactive=True,evidence=current['latest_evidence'])
    if gate['id']=='actual_final_trajectory_fields_and_support':
        gate.update(status='PASS_ONE_CELL_PENDING_OTHER_TRAJECTORIES',reason='SSP160 one cell checked at11 actual times; two-cell/fine-grid/equilibrium field and support gates remain pending.',
            one_cell_fields=fields['maxima'],one_cell_support_bound=support['worst_relative_upper_bound'])
current['gates'].insert(7,dict(id='SSP_two_cell_temporal_admission',status='PENDING_USER_LONG_RUN',
    plan='docs/implementation/stage2/time_pass_20260922/two_cell_time_plan.json'))
plan=ROOT/'docs/implementation/stage2/time_pass_20260922/two_cell_time_plan.json'
current['computation_handoff']=dict(mode='USER_FOREGROUND_COMMAND',plan=plan.relative_to(ROOT).as_posix(),plan_sha256=sha(plan),
    output_root='/home/jdiaz/pysnspd/tmp/stage2_ssp_two_time_20260922',expected_minutes=[45,65],threads=1,memory_reservation_GiB=2,
    memory_peak_measured=False,agent_executed=False,new_trajectories_this_review=0,
    scope='Two-cell temporal block only. Full remaining gate design is not yet an executable complete closure batch.')
current['latest_audit']=dict(path=current['latest_evidence'],sha256=sha(DATA/'one_cell_audit.json'),checks_passed=85,checks_total=85,
    new_postprocessing_checks=[dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p)) for p in (DATA/'one_ssp_admitted_fields_20260922.json',DATA/'one_ssp_admitted_support_20260922.json')])
current['historical_review_records']=[old['latest_review']]
current['latest_review']=dict(directory=DATA.relative_to(ROOT).as_posix(),status=current['status'],
    report='output/pdf/implementation/Informe_resultados_temporales_etapa_2_20260922.pdf',
    audit=current['latest_audit'],global_stage2_closed=False)
current['generator_sha256']=sha(Path(__file__))
current['stage3_preparation']=dict(status='PREPARADA_NO_INICIADA',entry_contract='docs/implementation/stage3/entry_contract.json',blocked_by='Stage2 global numerical admission')
current['publication']=dict(status='CLOSURE_PUSH_WITHHELD',reason='User conditioned the closure push on all stage2 gates being in order; mandatory gates are still missing. Working evidence synchronized to Geminga.')
write(admission,current)
notebook=ROOT/'docs/GEMINGA_COMMANDS.md';oldnotebook=DATA/'GEMINGA_COMMANDS_before_time_pass.md'
if not oldnotebook.exists():oldnotebook.write_bytes(notebook.read_bytes())
assert sha(oldnotebook)=='8647d223b12364fd8474c45dd54ff3d103ee7814a672367616ff041c68a01303'
notebook.write_text('''# Geminga: comandos vigentes

Actualizado: 22 de septiembre de 2026. Cuenta `jdiaz`; repositorio
`/home/jdiaz/pysnspd`; base Git `031991e` más la revisión de trabajo temporal.

## Estado comprobado

El lote `tmp/stage2_ssp_time_20260922/` terminó: **85/85 controles PASS**.
SSP160, SSP320 y SSP640 de una celda cumplen precisión y energía frente a 1280.
Campos y soporte de SSP160 también pasan. No repetir ninguna de estas corridas.
La etapa 2 global sigue abierta; la etapa 3 está preparada pero no iniciada.

La continuación ejecutable de abajo cubre **tiempo de DOS celdas**, no el cierre
global. Después quedan mallas electrónicas/fonónicas, equilibrio, fronteras y
campos/soporte de las trayectorias aceptadas. El diseño completo está en
`docs/implementation/stage2/time_pass_20260922/remaining_gates_plan.md`.
Sus estimaciones de 5 a 7,5 h para los bloques restantes, si bastan 160 pasos, NO son
una orden de ejecución: todavía faltan evaluadores y la selección del paso.

Los comandos anteriores quedan conservados en el
[histórico de la revisión anterior](/home/jdiaz/pysnspd/docs/implementation/stage2/time_pass_20260922/GEMINGA_COMMANDS_before_time_pass.md),
que enlaza a su vez con el histórico original completo. No relanzar esos lotes.

## Preparar la terminal

```bash
cd /home/jdiaz/pysnspd
NOTEBOOK_PY=/home/jdiaz/.conda/envs/snspd/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
```

## Siguiente cálculo largo: tiempo de DOS celdas

Propósito: comprobar SSP160/320/640 contra una referencia separada SSP1280
en el mismo problema TWO, incluyendo transporte, escape 15 y calentamiento 0.01.
Después se evalúa la pareja 320/640 con el presupuesto suplementario 2.5e-5 para
orientar las futuras comparaciones de malla. Ese diagnóstico no sustituye
la admisión temporal de tres niveles ni decide por sí solo el cierre global.

Recursos estimados: **45 a 65 minutos**, un hilo; reservar **2 GiB**. El pico de memoria
no está medido; la estimación corresponde sólo a 630 estados y 1025 fonones.
No repetir ONE, DOP853 ni la malla 2520 de los lotes anteriores.

Vista previa opcional; sólo comprueba el plan y los hashes:

```bash
"$NOTEBOOK_PY" sandbox/stage2_cells/time_pass_20260922/run_two_cell_time_batch.py --plan docs/implementation/stage2/time_pass_20260922/two_cell_time_plan.json --output-root tmp/stage2_ssp_two_time_20260922 --dry-run
```

Ejecución manual en esta terminal:

```bash
"$NOTEBOOK_PY" -u sandbox/stage2_cells/time_pass_20260922/run_two_cell_time_batch.py --plan docs/implementation/stage2/time_pass_20260922/two_cell_time_plan.json --output-root tmp/stage2_ssp_two_time_20260922 --execute
```

Salidas: cuatro JSON de trayectoria `two_ssp_{160,320,640,1280}.json`, sus NPZ y
estados iniciales, evaluación temporal, evaluación pareada, recibos y registros
en `tmp/stage2_ssp_two_time_20260922/`. Parada al primer fallo, sin reintentos,
recortes ni sobrescritura. Conservar todo y escribir «reinicia» al terminar o
detenerse; basta indicar la carpeta. No hace falta copiar todo el output al chat.

## Diagnósticos ligeros reproducibles

Auditoría de las corridas ya importadas, sin nueva dinámica (unos 5 segundos),
guardando la repetición fuera de la evidencia congelada:

```bash
timeout --signal=TERM --kill-after=5s 240s "$NOTEBOOK_PY" sandbox/stage2_cells/time_pass_20260922/audit_one_cell.py --output tmp/one_cell_time_audit_repeat.json
```

Pruebas del orquestador (segundos, sin RHS) y comprobación de la entrega:

```bash
timeout --signal=TERM --kill-after=5s 240s "$NOTEBOOK_PY" sandbox/stage2_cells/time_pass_20260922/run_two_cell_time_batch.py --self-test
"$NOTEBOOK_PY" sandbox/stage2_cells/verify_delivery.py
```

Campos y soporte de la trayectoria aprobada (menos de 3 s observados por control).
Usar una ruta nueva de salida para conservar la evidencia original:

```bash
timeout --signal=TERM --kill-after=5s 240s "$NOTEBOOK_PY" sandbox/stage2_cells/check_trajectory_fields.py tmp/stage2_ssp_time_20260922/one_ssp_160.json --output tmp/one_ssp_fields_repeat.json
timeout --signal=TERM --kill-after=5s 240s "$NOTEBOOK_PY" sandbox/stage2_cells/check_trajectory_support.py tmp/stage2_ssp_time_20260922/one_ssp_160.json --output tmp/one_ssp_support_repeat.json
```

Un trabajo previsto por encima de cinco minutos lo ejecuta el usuario. Si un
diagnóstico acotado agota su tiempo, queda incompleto; no se divide ni reintenta
para eludir ese límite. El push de cierre y el inicio espacial permanecen
condicionados a la admisión global de la etapa 2.
''',encoding='utf-8')
print(json.dumps(dict(status=current['status'],imported_files=len(inventory),notebook_sha256=sha(notebook)),indent=2))
