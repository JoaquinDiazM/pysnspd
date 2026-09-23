# Geminga: comandos vigentes

Actualizado el 23 de septiembre de 2026. Cuenta `jdiaz`, sin administrador.

## Pendiente: etapa 3A, esquema nodal de cuarto orden

**Listo para ejecución del usuario; no lanzado por el agente.**

La campaña anterior terminó sus 18 casos en 15,15 minutos. Conserva sus
11 fallos de precisión espacial: el punto medio cartesiano introduce un sesgo
de amplitud y corriente. Los controles de energía, signos y espectro pasaron.
No repetir el lote anterior ni las trayectorias de etapa 2.

El nuevo esquema conserva amplitudes nodales, usa operadores covariantes
compatibles y obtiene fuerza y corriente de una sola energía. Mantiene el
modelo continuo, los seis perfiles, las mallas 8/16/32 y el objetivo del 1 %.
Ya pasaron 102 pruebas breves y un piloto de 33.2 s.

Tiempo estimado: **unos 7 minutos**; margen orientativo 6-12 minutos.
Recursos estimados: un proceso CPU, menos de 1 GB RAM y decenas de MB de salida.
Este es un cálculo estático; no un transitorio completo del detector.

Ejecutar directamente en la terminal:

```bash
cd /home/jdiaz/pysnspd
NOTEBOOK_PY=/home/jdiaz/.conda/envs/snspd/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$NOTEBOOK_PY" -u sandbox/stage3_spatial/run_nodal_batch.py \
  --registration docs/implementation/stage3/nodal_20260923/registration.json \
  --reuse tmp/stage3_nodal_20260923/pilot \
  --output-root tmp/stage3a_nodal_full_20260923 \
  --execute
```

La terminal muestra barras por caso y lote, tiempo transcurrido y ETA aproximada.
`progress.jsonl` guarda el avance. El piloto se reutiliza solo si coinciden
fuentes, estado, malla y hashes. El directorio de salida debe ser nuevo.

El lote mantiene los controles de fase, referencias uniformes, estabilidad y
respuesta espacial en las tres mallas. Las 24 evaluaciones extra de diferencias
de energía se hacen por caso en N=16 y en el piloto N=8. Las otras mallas declaran
explícitamente que ese control no se repitió. La energía/fuerza tiene orden formal
cuatro; para la corriente conjugada de caras se usa una estimación de orden dos,
además del orden observado. El presupuesto final sigue siendo 1 %.

Salidas en `tmp/stage3a_nodal_full_20260923/`:

- `manifest.json`, `progress.jsonl`: fuentes, entorno, avance y tiempos.
- `*_initial.npz`, `*.npz`, `*.json`, `*.receipt.json`: estados, resultados e integridad.
- `negative_control.json`: rechazo esperado del estado físicamente inestable.
- `summary.json`: precisión, controles ejecutados y alcance; no cierra automáticamente la etapa 3.
- `failure.json`, si corresponde: diagnóstico de una interrupción de cálculo.

Un resumen con precisión insuficiente puede aparecer tras completar los 18 casos.
Conservarlo; no relanzar ni sobrescribir automáticamente. Indicar «terminó» o
«reinicia» es suficiente: el agente leerá los archivos directamente de Geminga.

## Consultas ligeras opcionales

Listar casos sin calcular física:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_spatial/run_nodal_batch.py
```

Verificar entrega y compatibilidad del piloto antes de lanzar el lote:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_spatial/verify_nodal_delivery.py
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_spatial/check_nodal_handoff.py
```

Últimas líneas del progreso guardado:

```bash
tail -n 8 /home/jdiaz/pysnspd/tmp/stage3a_nodal_full_20260923/progress.jsonl
```

## Historial y secuencia pendiente

Libreta anterior íntegra:
`docs/implementation/stage3/review_20260923/GEMINGA_COMMANDS_before_review.md`.
Resultados anteriores íntegros:
`docs/implementation/stage3/review_20260923/raw/static_campaign/`.
Informe:
`output/pdf/implementation/Informe_revision_espacial_etapa_3_20260923.pdf`.

Tras revisar la convergencia nodal siguen bordes y reservorios, potencial,
el circuito de tres estados de la memoria y dinámica débil con depósito
sintético. La etapa 3 no se cierra con esta prueba periódica estática.

Todo cálculo previsto de más de cinco minutos se deja al usuario. No se usan
trabajos de fondo ni sondeos automáticos para eludir esa entrega. Las pruebas
breves inciertas tienen un límite de 240 s; un timeout se conserva incompleto.
