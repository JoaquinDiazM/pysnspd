# Geminga: comandos vigentes

Actualizado el 23 de septiembre de 2026. Cuenta `jdiaz`, sin administrador.

## Pendiente de ejecución: etapa 3A, ensayo espacial estático

**Estado: preparado, NO lanzado por el agente.** El piloto de ocho celdas pasó
en 33.9 s, y las 74 pruebas breves pasaron. El lote completo
supera cinco minutos y queda para ejecución del usuario.

Propósito: comparar energía, fuerza y corriente del funcional espacial en seis
casos y mallas de 8/16/32 celdas; controlar el signo de D.36 y contrastar
630/1260 estados electrónicos en puntos representativos. Es estático: todavía
no integra el detector con bordes ni circuito.

Duración estimada por el piloto: **unos 16 minutos**; margen orientativo
12-25 minutos. Recursos estimados: un proceso CPU, menos de 1 GB RAM y
decenas de MB de salida. El costo no incluye una ejecución del transitorio completo.

Ejecutar este bloque directamente en la terminal:

```bash
cd /home/jdiaz/pysnspd
NOTEBOOK_PY=/home/jdiaz/.conda/envs/snspd/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$NOTEBOOK_PY" -u sandbox/stage3_spatial/run_static_batch.py \
  --registration docs/implementation/stage3/iteration_20260923/registration.json \
  --reuse tmp/stage3_start_20260923/pilot_reviewed \
  --output-root tmp/stage3a_static_20260923 \
  --execute
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
