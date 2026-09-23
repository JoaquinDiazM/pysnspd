# Geminga: comandos vigentes

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
