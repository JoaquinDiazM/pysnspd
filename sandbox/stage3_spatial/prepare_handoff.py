"""Publish the measured pilot and the unexecuted user-run command notebook."""
from pathlib import Path
import hashlib
import json
import re

ROOT=Path(__file__).resolve().parents[2]
DATA=ROOT/'docs/implementation/stage3/iteration_20260923'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
notebook=ROOT/'docs/GEMINGA_COMMANDS.md'
archive=DATA/'GEMINGA_COMMANDS_before_stage3.md'
expected='fc9940a313f85af5a42a5d679387990f80142e8c149b2aeb620f5724ccef008c'
if not archive.exists():
    assert sha(notebook)==expected,'Notebook changed: preserve new user entries before updating'
    archive.write_bytes(notebook.read_bytes())
assert sha(archive)==expected
pilot=json.loads((DATA/'pilot/weak_phase_thermal_n8.json').read_text(encoding='utf-8'))
tests_passed=int(re.search(r'(\d+) passed', (DATA/'checks/pytest.log').read_text()).group(1))
estimated=pilot['runtime_seconds']*((4*sum((8,16,32))-8)/8)+25
command=r'''cd /home/jdiaz/pysnspd
NOTEBOOK_PY=/home/jdiaz/.conda/envs/snspd/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$NOTEBOOK_PY" -u sandbox/stage3_spatial/run_static_batch.py \
  --registration docs/implementation/stage3/iteration_20260923/registration.json \
  --reuse tmp/stage3_start_20260923/pilot_reviewed \
  --output-root tmp/stage3a_static_20260923 \
  --execute'''
contents=f'''# Geminga: comandos vigentes

Actualizado el 23 de septiembre de 2026. Cuenta `jdiaz`, sin administrador.

## Pendiente de ejecución: etapa 3A, ensayo espacial estático

**Estado: preparado, NO lanzado por el agente.** El piloto de ocho celdas pasó
en {pilot['runtime_seconds']:.1f} s, y las {tests_passed} pruebas breves pasaron. El lote completo
supera cinco minutos y queda para ejecución del usuario.

Propósito: comparar energía, fuerza y corriente del funcional espacial en seis
casos y mallas de 8/16/32 celdas; controlar el signo de D.36 y contrastar
630/1260 estados electrónicos en puntos representativos. Es estático: todavía
no integra el detector con bordes ni circuito.

Duración estimada por el piloto: **unos {estimated/60:.0f} minutos**; margen orientativo
12-25 minutos. Recursos estimados: un proceso CPU, menos de 1 GB RAM y
decenas de MB de salida. El costo no incluye una ejecución del transitorio completo.

Ejecutar este bloque directamente en la terminal:

```bash
{command}
```

El piloto se reutiliza únicamente si coinciden fuentes, caso, malla, condiciones
iniciales y hashes. La terminal muestra barras del caso y del lote, unidades
terminadas, tiempo transcurrido y **ETA aproximada**. Al comenzar puede decir
`unknown` hasta tener mediciones; `100%` se reserva para finalización exitosa.
Las tareas tienen costos distintos, por lo que la ETA puede corregirse.

Salidas esperadas en `tmp/stage3a_static_20260923/`:

- `progress.jsonl`: avance y tiempos guardados automáticamente.
- `manifest.json`: fuentes, entradas, entorno y lista de casos.
- `*_initial.npz`, `*.json`, `*.npz`, `*.receipt.json`: condiciones, resultados e integridad.
- `negative_control.json`: estado que debe rechazarse por inestabilidad.
- `summary.json`: comparaciones y estado del lote; no cierra automáticamente toda la etapa 3.
- `failure.json`: motivo y diagnóstico si ocurre un fallo de ejecución.

El directorio debe ser nuevo. Si falla, conservarlo y comunicarlo; no repetir
el bloque sobre la misma carpeta ni relanzar automáticamente. Basta indicar
«terminó» o «reinicia» para que el agente lea los archivos de Geminga.

## Consultas ligeras opcionales

Ver la lista de casos sin calcular:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_spatial/run_static_batch.py
```

Ver las últimas líneas guardadas mientras tu ejecución está activa:

```bash
tail -n 8 /home/jdiaz/pysnspd/tmp/stage3a_static_20260923/progress.jsonl
```

Verificar archivos de esta entrega y del piloto, sin física:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_spatial/verify_delivery.py
```

Informe de inicio:
`output/pdf/implementation/Informe_inicio_etapa_3_20260923.pdf`.

## Historial preservado y política de cálculo

La libreta del cierre de etapa 2 está íntegra en
[GEMINGA_COMMANDS_before_stage3.md](/home/jdiaz/pysnspd/docs/implementation/stage3/iteration_20260923/GEMINGA_COMMANDS_before_stage3.md).
Es un antecedente; no una cola de tareas pendientes. Su etapa permanece cerrada
en el alcance declarado y sus datos no se modifican.

Todo cálculo previsto de más de cinco minutos se prepara aquí para ejecución
del usuario. No se ejecutan trabajos de fondo ni se programa sondeo automático
para eludir esa entrega. Las comprobaciones breves inciertas usan un límite
de 240 segundos y un timeout se conserva como resultado incompleto.
'''
notebook.write_text(contents,encoding='utf-8',newline='\n')
handoff=dict(status='AWAITING_USER_LONG_STATIC_CAMPAIGN',date='2026-09-23',
    implementation_started=True,stage3_closed=False,production_promotion=False,
    implemented_scope='3A periodic 1D static functional, Cartesian forces, conjugate current and local D.36',
    current_primary_document='docs/implementation/stage3/iteration_20260923/README.md',
    preserved_proposal='docs/implementation/stage3/entry_contract.json',
    tests_passed=tests_passed,pilot_status=pilot['status'],pilot_runtime_seconds=pilot['runtime_seconds'],
    expected_long_runtime_seconds=estimated,estimated_runtime_range_minutes=[12,25],
    command=command,long_run_launched_by_agent=False,
    user_output_root='/home/jdiaz/pysnspd/tmp/stage3a_static_20260923',
    notebook='/home/jdiaz/GEMINGA_COMMANDS.md',
    next_dependency='Review the user-run spatial refinement and spectral sign comparisons before continuing 3B.',
    source_hashes=pilot['source_hashes'],
    full_stage3_sequence=['3A static spatial','3B boundaries/reservoirs','3C charge/potential',
                          '3D thesis three-state circuit','3E admitted weak dynamics'])
(DATA/'handoff.json').write_text(json.dumps(handoff,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
print(json.dumps(dict(status=handoff['status'],estimated_minutes=estimated/60,notebook_sha256=sha(notebook))))
