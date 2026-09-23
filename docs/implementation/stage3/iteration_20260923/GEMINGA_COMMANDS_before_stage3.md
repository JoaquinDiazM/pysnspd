# Geminga: comandos vigentes

Actualizado el 22 de septiembre de 2026. Cuenta `jdiaz`, sin administrador.

## Cierre de desarrollo de etapa 2

La etapa 2 queda cerrada como desarrollo por decisión del usuario. El informe
final está en `docs/implementation/stage2/closure_20260922/Informe_cierre_etapa_2.md`
y `output/pdf/implementation/Informe_cierre_etapa_2.pdf`; la decisión verificable
está en `docs/implementation/stage2/closure_20260922/closure_decision.json`.
El certificado numérico estricto de mallas permanece incompleto y los fallos
históricos conservan su dictamen. No hay promoción a producción.

La etapa 3 queda preparada, no iniciada: energía espacial y estabilidad,
bordes/reservorios, conservación de carga, circuito de tres variables y después
dinámica débil y deposición localizada sintética. Su contrato está en
`docs/implementation/stage3/entry_contract.json`.

## No hace falta ejecutar otro lote largo ahora

El lote `tmp/stage2_guarded_acceptance_20260922` completó 21 tareas y 13
trayectorias en unos 94,5 minutos. Las trayectorias tienen poblaciones físicas
y balances válidos. El fallo fue una diferencia temporal fonónica de
**0,004565 %**, frente a un presupuesto auxiliar de **0,0025 %**.
El resultado original queda como FAIL; no se relajan sus tolerancias.

Ese presupuesto no es un requisito físico fundamental ni una dependencia de
un ensayo espacial estático con ocupaciones congeladas. La ruta vigente prepara
la secuencia de etapa 3, empezando por ese ensayo y sus controles propios,
y conserva abierta la certificación dinámica de mallas. No se propone repetir
ahora las 35 tareas ni completar automáticamente las mallas más caras.
Los datos ya obtenidos se reutilizarán
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

Ver el cierre vigente y comprobar sus referencias (segundos; sin dinámica):

```bash
timeout --signal=TERM --kill-after=5s 240s "$NOTEBOOK_PY" sandbox/stage2_cells/closure_20260922/show_closure.py --verify
```

Salida esperada: resumen del cierre de desarrollo, certificado numérico pendiente
y etapa 3 preparada, con verificación de sus referencias. No crea trayectorias;
usa un proceso ligero de lectura y memoria de orden de decenas de MB.

Consultar la revisión anterior y sus hashes (segundos; antecedente histórico):

```bash
timeout --signal=TERM --kill-after=5s 240s "$NOTEBOOK_PY" sandbox/stage2_cells/practical_review_20260922/show_checkpoint.py --verify
```

Verificar la entrega documental y numérica completa por hashes (segundos):

```bash
"$NOTEBOOK_PY" sandbox/stage2_cells/verify_delivery.py
```

Informe anterior de revisión, conservado como antecedente:
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
