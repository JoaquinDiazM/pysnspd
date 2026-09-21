# Cierre de etapa 1: celdas y datos fonónicos

La entrada de lectura es el [informe final](../../docs/implementation/stage1_closure/Informe_cierre_etapa_1.md).
El [dictamen](../../docs/implementation/stage1_closure/closure_admission.json)
identifica exactamente qué queda aceptado y para qué siguiente trabajo.
Los [criterios](../../docs/implementation/stage1_closure/acceptance_criteria.json)
se fijaron antes de los resultados. El contrato de la siguiente implementación
está en [NEXT_STAGE.md](../../docs/implementation/stage1_closure/NEXT_STAGE.md).

El módulo `pysnspd.experimental.cell_validation` utiliza el catálogo R2 guardado,
sin volver a resolverlo ni modificar su energía. Distingue la cuadratura nativa
del catálogo de la reconstrucción continua usada por el prototipo de transporte.
Los balances y errores de ambas representaciones se publican por separado.
El tiempo de las celdas es sintético; no representa una calibración cinética NbN.

El preprocesado fonónico genera una tabla derivada con procedencia y operaciones
explícitas. Recortar valores en una entrada material autorizada no equivale a
recortar poblaciones o energía durante una trayectoria. Los operadores de celdas
rechazan poblaciones no físicas y salidas de soporte.

## Entorno de Geminga

```bash
cd /home/jdiaz/pysnspd
export CLOSURE_PY=/home/jdiaz/.conda/envs/snspd/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p tmp/stage1_closure_reproduction
```

Las recetas y sus tiempos medidos quedan en `/home/jdiaz/GEMINGA_COMMANDS.md`,
con copia en `docs/GEMINGA_COMMANDS.md`. Los cálculos inciertos se limitan a
240 segundos. Si el límite se agota, la ejecución queda incompleta y debe
entregarse al usuario; no se reinicia en partes.

## Comprobación del catálogo retenido y de la entrega

```bash
"$CLOSURE_PY" sandbox/stage1_catalog_r2/verify_delivery.py
timeout --signal=TERM --kill-after=5 240 "$CLOSURE_PY" \
  docs/implementation/stage1_r2/review/assess_candidate.py \
  docs/implementation/stage1_r2/catalogs/occupation_catalog.npz \
  --output tmp/stage1_closure_reproduction/catalog_recheck.json
"$CLOSURE_PY" sandbox/stage1_closure/verify_delivery.py
```

El primer manifiesto sigue verificando la evidencia R2 original. El último
verifica los archivos de este cierre y los hashes de los insumos retenidos.
Reproducir un cálculo no autoriza a sobrescribir su dictamen de aceptación.

## Pruebas y documento

Repetir celdas y revisión de la tabla material en una carpeta nueva (~2 s cada
diagnóstico, un hilo; la consulta pública depende de la conexión):

```bash
timeout --signal=TERM --kill-after=5 240 "$CLOSURE_PY" \
  sandbox/stage1_closure/run_cells.py \
  --output tmp/stage1_closure_reproduction/cells
timeout --signal=TERM --kill-after=5 240 "$CLOSURE_PY" \
  sandbox/stage1_closure/review_material.py \
  --nbn-path /home/jdiaz/scratch/big_data/catalogs/simon_2025/nbn-a2f-ph.dat \
  --output-root tmp/stage1_closure_reproduction/material \
  --source-dir tmp/stage1_closure_reproduction/material_sources --online
```

La segunda receta verifica el cuerpo numérico frente a la revisión pública
fijada. Sin `--online` necesita las fuentes exactas ya descargadas en
`--source-dir`. La tabla CSV derivada de la entrega permite estudiar la forma
sin la tabla original, pero no sustituye ésta en la auditoría del preprocesado.
Los resultados finales fueron adjudicados con
`docs/implementation/stage1_closure/review/adjudicate.py`; el certificado se
refiere a esos bytes y a su fuente, no a una reproducción nueva.

```bash
timeout --signal=TERM --kill-after=5 240 "$CLOSURE_PY" -m pytest -q
PYTHONPATH=/home/jdiaz/pysnspd/tmp/stage1_report_deps \
  "$CLOSURE_PY" sandbox/stage1_closure/build_report.py
```

La regresión de cierre pasó 317 pruebas en 34,69 s. El registro con los hashes
antes/después de módulos y pruebas es `regression_result.json`; su generador es
`run_regression.py`. Invocarlo escribe ese registro en la ruta de entrega, por
eso la receta ordinaria anterior ejecuta pytest sin reemplazarlo.

El constructor del informe sólo combina las figuras y el contenido ya revisados;
no ejecuta una trayectoria. Cambia el PDF publicado y su Markdown. Deben revisarse
todas las páginas renderizadas antes de generar un manifiesto nuevo con
`verify_delivery.py --write`. Las fuentes disponibles pueden cambiar la
maquetación entre Windows y Linux: el PDF entregado queda fijado por hash.

La etiqueta `v1.0.0`, los resultados históricos y los módulos de producción
permanecen fuera de los cambios de este cierre.
