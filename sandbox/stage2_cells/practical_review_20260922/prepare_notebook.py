"""Retire the costly automatic continuation and retain all prior instructions."""
from pathlib import Path
import hashlib
ROOT=Path(__file__).resolve().parents[3];DATA=ROOT/'docs/implementation/stage2/practical_review_20260922'
notebook=ROOT/'docs/GEMINGA_COMMANDS.md';archive=DATA/'GEMINGA_COMMANDS_before_scope_review.md'
if not archive.exists():archive.write_bytes(notebook.read_bytes())
assert hashlib.sha256(archive.read_bytes()).hexdigest()=='cd3ef2bb32b579b77fe845a736aa0cc4b0ee98b71fd8c77fc060d0edcc6edc6d'
notebook.write_text('''# Geminga: comandos vigentes

Actualizado el 22 de septiembre de 2026. Cuenta `jdiaz`, sin administrador.

## No hace falta ejecutar otro lote largo ahora

El lote `tmp/stage2_guarded_acceptance_20260922` completó 21 tareas y 13
trayectorias en unos 94,5 minutos. Las trayectorias tienen poblaciones físicas
y balances válidos. El fallo fue una diferencia temporal fonónica de
**0,004565 %**, frente a un presupuesto auxiliar de **0,0025 %**.
El resultado original queda como FAIL; no se relajan sus tolerancias.

Ese presupuesto no es un requisito físico fundamental ni una dependencia de
un ensayo espacial estático con ocupaciones congeladas. La ruta vigente permite
preparar 3A estática con sus controles propios y conserva abierta la certificación
dinámica de mallas. No se propone repetir ahora las 35 tareas ni completar
automáticamente las mallas más caras. Los datos ya obtenidos se reutilizarán
cuando sean compatibles con la pregunta siguiente.

El modelo previsto adopta el circuito de tres variables de la memoria.
Documento vigente: `docs/implementation/MODELO_VIGENTE.md`.
Ecuaciones: `docs/modelo_v0_4/actualizaciones/circuito_memoria_20260922.md`.
Todavía no se implementa el nuevo acoplamiento espacial/circuital.

Histórico íntegro de la libreta anterior:
[GEMINGA_COMMANDS_before_scope_review.md](/home/jdiaz/pysnspd/docs/implementation/stage2/practical_review_20260922/GEMINGA_COMMANDS_before_scope_review.md).
Sus órdenes son históricas; no constituyen una cola de trabajos pendientes.

## Comprobaciones ligeras opcionales

Preparación:

```bash
cd /home/jdiaz/pysnspd
NOTEBOOK_PY=/home/jdiaz/.conda/envs/snspd/bin/python
```

Ver resumen y comprobar hashes de datos recibidos (segundos; sin dinámica):

```bash
timeout --signal=TERM --kill-after=5s 240s "$NOTEBOOK_PY" sandbox/stage2_cells/practical_review_20260922/show_checkpoint.py --verify
```

Verificar la entrega documental y numérica completa por hashes (segundos):

```bash
"$NOTEBOOK_PY" sandbox/stage2_cells/verify_delivery.py
```

Informe con resultados y plots:
`output/pdf/implementation/Informe_revision_criterios_y_circuito_20260922.pdf`.
La comparación circuital usa una resistencia prescrita, no un transitorio
calculado del detector; su CSV y parámetros están en
`docs/implementation/stage2/practical_review_20260922/`.

## Política para cálculos posteriores

Conservar unidades, signos, balances y dominios físicos. Antes de exigir una
precisión dinámica nueva, fijar el observable del dispositivo y la diferencia
que queremos resolver. Registrar el presupuesto antes de nuevas mediciones.
Repetir trayectorias sólo ante cambios pertinentes o una incertidumbre que
impida decidir. No transformar un fallo histórico en PASS cambiando su umbral.

Todo cálculo estimado por encima de cinco minutos quedará aquí con propósito,
comando, salidas y recursos para ejecución del usuario. No hay ningún comando
largo nuevo recomendado en esta entrega. No se han lanzado trabajos de fondo.
''',encoding='utf-8')
print(hashlib.sha256(notebook.read_bytes()).hexdigest())
