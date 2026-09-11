# Etapa 1, iteración R2: reproducción y evidencia

Este directorio prepara y verifica el catálogo electrónico experimental. El estado de la entrega se consulta en el [certificado externo](../../docs/implementation/stage1_r2/catalog_admission.json), la [evaluación final](../../docs/implementation/stage1_r2/review/final_assessment.json) y el [informe](../../docs/implementation/stage1_r2/Informe_etapa_1_r2.md). Si alguno todavía no existe, esa parte del cierre está pendiente. Un `PASS` de un piloto o del constructor no sustituye al certificado.

La versión de producción basada en la memoria quedó en la etiqueta `v1.0.0`, commit `5ea0cd6a08d2cae414c82944f8230469ca5aa7d6`. La evidencia de la primera iteración permanece en `docs/implementation/stage1/` y en el commit `d17d7c3fc3943699631126eaa636732bba56aabe`. R2 conserva sus propios criterios, referencias, pilotos rechazados y puntos de regresión. Estos scripts no conectan el catálogo a PRE, SS, gTDGL ni a un transiente de fotón.

## Qué compara cada figura

| Script | Dato y comparación | Evidencia |
|---|---|---|
| `material_comparison.py` | DOS **fonónica** y función de interacción α²F del archivo público de Simón (2025), frente al tratamiento real del cargador de producción. Distingue ordenación, abscisas duplicadas y recorte de valores negativos. El tratamiento experimental conserva el rechazo; no inventa una tercera DOS. | [Diagnóstico material](../../docs/implementation/stage1_r2/material_comparison.json), [procedencia](../../docs/implementation/stage1_r2/material_source_notice.md), [curvas](../../docs/implementation/stage1_r2/material_curves.csv). |
| `electronic_comparison.py` | DOS **electrónica** normalizada de Usadel, a campos fijos. Separa la comparación de implementaciones con el mismo η del efecto de reducir ese regulador numérico. | [Diagnóstico electrónico](../../docs/implementation/stage1_r2/electronic_comparison.json). |

La DOS fonónica cuenta modos vibracionales; la DOS electrónica describe estados electrónicos y usa la referencia normal por espín `N0`. No son dos versiones de una misma curva. La cabecera local del archivo de NbN no certifica por sí sola la unidad o la normalización de la ordenada pública. Las fuentes primarias y sus hashes están registrados en el diagnóstico material; una aprobación numérica electrónica no admite ese material.

## Representación y soporte

El esquema de archivo electrónico `v5` almacena una energía BCS de referencia más una corrección energética separada. El eje interno es `r = Gamma / absDelta`; la interpolación usa `log1p(r / gamma_coordinate_scale)`. El archivo conserva los valores y las derivadas nodales de energía, el eje `r`, la escala de transformación, η y la cuadratura común de conteo de estados `x`.

Las consultas reciben los campos físicos normalizados `absDelta/Delta0` y `Gamma/Delta0`. Las fuerzas se devuelven respecto a esos campos después de aplicar la regla de la cadena. Al variar los campos se mantiene fija la población `p(x)` en la cuadratura común: el catálogo admite distribuciones no térmicas.

El contrato R2 exige `absDelta/Delta0` en `[0.08, 1.5]` y `Gamma/Delta0` en `[0, 1.2]`, además del límite normal exacto `absDelta=0`. La extensión del eje interno `r` permite cubrir ese rectángulo; **no amplía el soporte físico de consulta**. Las consultas fuera de soporte, energías no positivas o pérdida del orden de los estados se rechazan, sin reparación de los datos.

## Entorno y política de ejecución

En Geminga se usa la cuenta no administrativa `jdiaz`. Preparación de una terminal Bash:

```bash
cd /home/jdiaz/pysnspd
export PY_R2=/home/jdiaz/.conda/envs/snspd/bin/python
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R2_RUN=tmp/stage1_r2_reproduction
```

Los comandos siguientes son recetas; sus tiempos finales se leen en los JSON de cada ejecución. Un trabajo que se sabe que supera cinco minutos se anota en `/home/jdiaz/GEMINGA_COMMANDS.md` para que lo ejecute el usuario. Si el coste es incierto, el límite de 240 segundos evita dejar un cálculo largo en marcha. No se encadenan fragmentos para eludir ese límite. Las nuevas reproducciones deben tener una ruta de salida propia.

### Comparaciones y pruebas breves

```bash
timeout 240 "$PY_R2" sandbox/stage1_catalog_r2/material_comparison.py \
  --nbn-path /home/jdiaz/scratch/big_data/catalogs/simon_2025/nbn-a2f-ph.dat \
  --output-root "$R2_RUN/material"
timeout 240 "$PY_R2" sandbox/stage1_catalog_r2/electronic_comparison.py --eta 1e-8
timeout 240 "$PY_R2" -m pytest \
  tests/test_experimental_material_admission.py \
  tests/test_experimental_energy_catalog.py \
  tests/test_experimental_energy_catalog_r2.py -q
```

`material_comparison.py --offline` omite la consulta nueva al repositorio público y usa la evidencia fijada por hash; el JSON deja constancia de esa elección. `electronic_comparison.py` escribe sus salidas en `docs/implementation/stage1_r2/`; al repetirlo se actualizan esos archivos, por lo que se debe conservar la versión que corresponda al manifiesto que se desea verificar.

### Construcción y refinamiento

Ejemplo explícito de construcción de un **candidato**, no de una aprobación:

```bash
timeout 240 "$PY_R2" sandbox/stage1_catalog_r2/build_catalog.py \
  --delta-nodes 17 --gamma-low-nodes 49 --gamma-high-nodes 17 \
  --count-order 10 --count-cutoff 12 --eta 1e-8 --ratio-coordinate \
  --output "$R2_RUN/base"
```

El constructor escribe `occupation_catalog.npz`, `vacuum_catalog.npz` y `candidate_diagnostics.json`. La opción `--catalog RUTA.npz` permite elegir el destino del catálogo de ocupaciones. Para reproducir un candidato concreto, se toman sus dimensiones, ejes, cuadratura y η de su `candidate_diagnostics.json`; el ejemplo anterior no declara ser la malla final. Cambiar la representación o la malla exige repetir la evaluación del nuevo archivo.

`refine_catalog.py` inserta únicamente los valores de `r` indicados mediante `--ratios`. La API `refine_occupation_catalog` exige conservar todos los nodos anteriores y ambos extremos. Los valores y derivadas de esos nodos se copian exactamente; sólo los nodos nuevos resuelven el espectro y, donde hace falta, su pequeña corrección energética. El informe de refinamiento registra el hash del catálogo de partida.

```bash
# Ejemplo de sintaxis: sustituir 0.45 y 0.6 por los nodos diagnosticados.
timeout 240 "$PY_R2" sandbox/stage1_catalog_r2/refine_catalog.py \
  "$R2_RUN/base/occupation_catalog.npz" --ratios 0.45 0.6 \
  --output "$R2_RUN/refined/occupation_catalog.npz" \
  --report "$R2_RUN/refined/refinement.json"
```

No se eligen nodos nuevos para desplazar los casos difíciles del contrato. Se conservan los casos fijos y los [puntos de regresión de soporte](../../docs/implementation/stage1_r2/review/support_regression_points.json), y se añaden controles en las nuevas celdas.

### Evaluación del archivo exacto

```bash
export R2_CATALOG="$R2_RUN/refined/occupation_catalog.npz"
timeout 240 "$PY_R2" docs/implementation/stage1_r2/review/assess_candidate.py \
  "$R2_CATALOG" --output "$R2_RUN/assessment.json"
timeout 240 "$PY_R2" sandbox/stage1_catalog_r2/check_shape_cells.py \
  "$R2_CATALOG" --output "$R2_RUN/shape_cells.json"
timeout 240 "$PY_R2" sandbox/stage1_catalog_r2/check_catalog_quadrature.py \
  "$R2_CATALOG" --output "$R2_RUN/count_rule_check.json"
```

La [evaluación independiente](../../docs/implementation/stage1_r2/review/assess_candidate.py) aplica los [criterios congelados](../../docs/implementation/stage1_r2/acceptance_criteria.json) al NPZ guardado: energía y fuerzas totales, referencias causales y térmicas, soporte, poblaciones válidas, derivadas, factores SI, paridad y persistencia. `check_shape_cells.py` calcula el mínimo del polinomio cúbico de las diferencias entre energías en cada intervalo interno, a amplitudes muestreadas. Esa comprobación no constituye una prueba uniforme sobre toda amplitud ni sustituye la comprobación de energías positivas. Sus controles de Bernstein y sus mínimos cúbicos tienen significados distintos; un control negativo de Bernstein aislado no equivale a una consulta fallida.

`check_catalog_quadrature.py` contrasta la cuadratura realmente guardada con la [referencia resuelta al mismo η](../../docs/implementation/stage1_r2/quadrature_diagnostics.json), sin mezclar ese error con interpolación de campos. Acepta `--reference OTRA_REFERENCIA.json`. El [cuaderno de revisión](../../docs/implementation/stage1_r2/review/REVIEW.md) explica las referencias independientes y los fallos históricos.

| Generador de referencia o diagnóstico | Salida R2 |
|---|---|
| `docs/implementation/stage1_r2/review/causal_reference.py` | [Referencia causal](../../docs/implementation/stage1_r2/review/causal_reference.json). |
| `docs/implementation/stage1_r2/review/gamma_zero_reference.py` | [Referencia de Γ=0 y regulador](../../docs/implementation/stage1_r2/review/gamma_zero_reference.json). |
| `docs/implementation/stage1_r2/review/thermal_reference.py` | [Referencia térmica y Matsubara](../../docs/implementation/stage1_r2/review/thermal_reference.json). |
| `docs/implementation/stage1_r2/review/spectral_crosscheck.py` | [Causalidad e inversión de conteo](../../docs/implementation/stage1_r2/review/spectral_crosscheck.json). |
| `sandbox/stage1_catalog_r2/diagnose_quadrature.py` | [Cuadratura, corte y regulador](../../docs/implementation/stage1_r2/quadrature_diagnostics.json); permite `--output`. |

Los generadores sin opciones escriben junto a su código. Sus resultados incluyen hashes y deben corresponder a la implementación y los criterios utilizados en la evaluación. No basta con reutilizar un `PASS` anterior después de cambiar el módulo.

## Certificado y entrega

El [certificado de admisión](../../docs/implementation/stage1_r2/catalog_admission.json) es externo y referencia los hashes del NPZ, los criterios y la evidencia de cierre. El catálogo conserva sus metadatos de construcción; no se reescribe después de medirlo sólo para cambiar una etiqueta de estado. La admisión se limita al dominio, poblaciones y tolerancias comprobados, y es independiente de la admisión material de NbN.

Los destinos de entrega son [catálogo de ocupaciones](../../docs/implementation/stage1_r2/catalogs/occupation_catalog.npz), [vacío](../../docs/implementation/stage1_r2/catalogs/vacuum_catalog.npz), [evaluación final](../../docs/implementation/stage1_r2/review/final_assessment.json) y [PDF](../../output/pdf/implementation/Informe_etapa_1_r2.pdf). Para generar los gráficos del dictamen y el documento:

```bash
"$PY_R2" sandbox/stage1_catalog_r2/plot_acceptance.py \
  docs/implementation/stage1_r2/review/final_assessment.json
"$PY_R2" sandbox/stage1_catalog_r2/build_report.py
```

`build_report.py` usa `report_content.json`; requiere revisar visualmente el PDF antes de cerrar. Tras completar la revisión numérica y visual se puede crear el [manifiesto de entrega](../../docs/implementation/stage1_r2/delivery_manifest.json). Después, la ejecución sin `--write` sólo verifica los archivos contra ese manifiesto:

```bash
# Crear únicamente al cerrar una entrega nueva.
"$PY_R2" sandbox/stage1_catalog_r2/verify_delivery.py --write
# Verificar una entrega existente sin reemplazar su manifiesto.
"$PY_R2" sandbox/stage1_catalog_r2/verify_delivery.py
```

No se deduce de estos controles la duración de un transiente, la latencia del detector ni la aceptación del modelo completo en producción.
