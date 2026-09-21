
## 2026-09-21 - Etapa 2: referencia temporal pendiente y continuación

Estado: la etapa 2 NO está cerrada. La implementación experimental y sus pruebas
estáticas avanzaron; falta admitir la dinámica con refinamientos temporales y
de malla. Una referencia DOP853 de una sola celda, duración 0,1 en unidades de
ensayo, permaneció INCOMPLETE bajo una ejecución limitada a 240 s. No se volvió
a ejecutar ni se dividió ese cálculo. El registro está en
docs/implementation/stage2/reference_handoff.json. No se dispone de evidencia
que atribuya la interrupción a un fallo físico.

La configuración candidata se afinó hasta 630 estados electrónicos y 1025
fonónicos. El catálogo R2 y v1.0.0 siguen intactos. Los 513 nodos fonónicos
anteriores no superaron todos los controles de distribución, aunque sus
balances y trayectorias preliminares fueran estables.

### Ejecución manual necesaria antes de continuar

Este comando NO fue lanzado por el agente. Ejecuta el diagnóstico de referencia
en la configuración candidata, con mensajes cada 100 evaluaciones del RHS. La
ecuación física y el integrador no se modifican por el observador. `t_rhs` es el
tiempo de una evaluación de prueba, no necesariamente el último paso aceptado:
puede retroceder al rechazar un paso. Si tarda mucho, esas líneas permiten
distinguir avance físico de un control adaptativo excesivamente restrictivo.

Recursos: un hilo de CPU; reservar 2 GiB de RAM. Tiempo: más de cuatro minutos
es posible; no hay una estimación fiable de finalización para el integrador
adaptativo porque el sondeo acotado no terminó. No es un transiente completo de
un detector. La ejecución queda a cargo del usuario, sin límite artificial de
240 s ni trabajo en segundo plano.

```bash
cd /home/jdiaz/pysnspd
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
STAGE2_PY=/home/jdiaz/.conda/envs/snspd/bin/python
STAGE2_RUN=/home/jdiaz/pysnspd/tmp/stage2_user_$(date +%Y%m%d_%H%M%S)
mkdir -p "$STAGE2_RUN"
"$STAGE2_PY" -u sandbox/stage2_cells/run_reference_progress.py --case one --duration 0.1 --steps 2 --escape inf --heating 0 --method dop853 --rtol 1e-9 --phonon-nodes 1025 --electron-refinement 1 --reaction-layout resolved_panels --reaction-order 2 --output "$STAGE2_RUN/one_reference_probe.json"
```

Salidas: JSON de parámetros/resultados, NPZ inicial, NPZ de trayectoria si
completa y `one_reference_probe.progress.jsonl`. Conservar todos los archivos.
Enviar la salida del terminal y la ruta STAGE2_RUN junto con `reinicia`.
Si se interrumpe manualmente, enviar también las últimas líneas de progreso;
una trayectoria incompleta no se contabiliza como aprobada.

Para reproducir exactamente el sondeo anterior, cambiar sólo
`--phonon-nodes 1025` por `--phonon-nodes 513` y usar una ruta nueva. Esa variante
es histórica: no es la candidata de aceptación.

### Diagnósticos ligeros que resultaron útiles

Pruebas del repositorio: aproximadamente un minuto, un hilo, menos de 2 GiB.
No sustituyen la referencia temporal pendiente.

```bash
cd /home/jdiaz/pysnspd
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout --signal=TERM --kill-after=5 240 /home/jdiaz/.conda/envs/snspd/bin/python -m pytest -q
```

Orden temporal de BGK/escape frente a soluciones exponenciales: menos de un
segundo de cálculo. Escribe `docs/implementation/stage2/time_suboperator_order.json`;
la reproducción modifica ese archivo de evidencia y puede cambiar su hash.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout --signal=TERM --kill-after=5 240 /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage2_cells/check_suboperator_order.py
```

La matriz estática de colisiones en 630/1260/2520 estados terminó en 229,70 s,
con picos de unos 4,5 GiB. Está conservada en `complementary_event_final.json`.
Su proximidad al límite desaconseja volver a lanzarla automáticamente. Para
reproducirla manualmente, se conserva el runner exacto en
`sandbox/stage2_cells/history/projected_event_diagnostics_uniform_v1.py`; copiarlo
a un nombre nuevo directamente dentro de `sandbox/stage2_cells/` mantiene su
resolución relativa de rutas. Sus parámetros fueron:
`--method projected --refinements 1 2 4 --orders 2 --integration-layout resolved_panels`.
Usar una salida nueva en tmp y revisar las convenciones fonónicas registradas.

Una trayectoria larga con 2520 estados y cientos de pasos se estima por encima
de veinte minutos. No fue ejecutada ni es un reemplazo de la referencia temporal
pendiente. Su necesidad debe decidirse después de diagnosticar el integrador.

Al continuar: revisar primero el progreso y el coste de DOP853; después fijar el
control de error temporal, completar los refinamientos sobre la configuración
admitida, los controles de soporte en las trayectorias y el informe final. No
avanzar al sistema espacial mientras stage2_admission.json siga pendiente.

Interpolación fonónica aislada con Lobatto real (15-25 s, un hilo, menos de 2 GiB;
salida en tmp, sin alterar el resultado publicado):

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout --signal=TERM --kill-after=5 240 /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage2_cells/check_projected_phonons.py --nodes 129 257 513 1025 --infrared .005 --output tmp/phonon_check.json
```
