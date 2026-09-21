# Comandos recurrentes en Geminga

Los ejemplos asumen una terminal interactiva en Geminga. Antes de trabajar,
ejecutar una vez:

```bash
source /opt/miniforge3/etc/profile.d/conda.sh
conda activate snspd
cd /home/jdiaz/pysnspd
CFG=configs/geminga_local_v3.yaml
PRE=pre_oe6_v4_L360nm_mesh4p0nm_k101T121D141Q_phase200T31D41Q2400W_power200Tph_01
NEW_PRE=pre_qwide_01
SS=ss_phasecg_I30uA_200ps_circuitthermal_stiffness_01
PHOTON=photon_phasecg_I30uA_0p8eV_sigma10nm_t50ps_1500ps_stiffness_01
PHOTON_EDGE=photon_phasecg_I30uA_0p8eV_sigma10nm_t50ps_1500ps_stiffness_01_edge
SWEEP_FAILED=ss_phasecg_sweep_I20to100uA_200ps_isothermal_nocircuit_01
SWEEP_PREFIX=ss_phasecg_sweep_I20to100uA_150ps_isothermal_nocircuit_01
```


`CFG` contiene los parámetros de producción y `PRE` identifica el catálogo
compartido; `SS` y `PHOTON` mantienen enlazada la cadena temporal. Los comandos
omiten flags que ya coinciden con sus defaults. Usar siempre nombres nuevos
para no mezclar resultados.

ToDo: Verificar que paso con "photon_phasecg_I20uA_0p8eV_sigma10nm_t50ps_1500ps_01"

## Validación ligera

```bash
python -m compileall -q pysnspd pipelines plot_pipelines tests
python -m pytest -q
```

La ayuda de cualquier entrada se consulta con, por ejemplo:

```bash
python pipelines/03_photon_run_template.py --help
```

## Simulaciones

### PRE

Esta PRE reemplaza obligatoriamente los catálogos que sólo guardaban
`js_A_m2`: SS y photon requieren ahora la rigidez Matsubara en
`(Te, |Delta|^2, |q|)`. Respecto de la última producción, los ejes de rigidez
pasan de `81x101x121` a `101x121x141`; fase/potencia conservan
`200x31x41x2400` y `200` temperaturas de fonones. Los 16 procesos coinciden
con la asignación paralela normal de Geminga.

```bash
python pipelines/01_prerun_template.py \
  --config "$CFG" --run-name "$PRE" --workers 16 \
  --stiffness-n-Te 101 --stiffness-n-delta 121 --stiffness-n-q 141 \
  --phase-n-Te 200 --phase-n-delta 31 --phase-n-q 41 --phase-n-omega 2400 \
  --power-n-Tph 200
```

### PRE temporal con rango q ampliado (pendiente de ejecución por el usuario)

D3 encontró que el catálogo actual termina en `q*xi=1.221`, con clipping de q
en hasta `28.8%` de los nodos centrales durante `41.2 ps`. El 99% de los
máximos temporales queda por debajo de `q*xi=2.792`; el máximo aislado de
`9.832` no se usa para dimensionar esta primera ampliación porque haría crecer
el catálogo cinético casi ocho veces. `NEW_PRE` lleva el rango nominal a
`q*xi~=3.05` (`2.5x`) y conserva la densidad nominal con 64 nodos q en
DOS/fase y 351 en la tabla uniforme de rigidez. Se espera que el catálogo de
fase crezca aproximadamente de 3.3 GB a 8.1 GB.

Este comando queda preparado para que lo ejecute jdiaz; Codex no lanza la
corrida de producción. No usar `NEW_PRE` en SS/photon hasta comprobar
`Status: OK`, los ejes realizados y el clipping residual con D3.

```bash
python pipelines/01_prerun_template.py \
  --config "$CFG" --run-name "$NEW_PRE" --workers 16 \
  --gamma-max-fraction 5.0 --dos-n-q 64 \
  --stiffness-n-Te 101 --stiffness-n-delta 121 --stiffness-n-q 351 \
  --phase-n-Te 200 --phase-n-delta 31 --phase-n-q 64 --phase-n-omega 2400 \
  --power-n-Tph 200
```

No lanzar SS hasta que PRE termine con `Status: OK`. Para una sensibilidad
mesoscópica distinta, `--allmaras-diffusion-factor` sigue siendo la opción
física que merece variarse explícitamente.

### SS de 30 uA y hasta 200 ps

El circuito, la térmica, el paso adaptativo, la condición de estacionariedad
total, la tolerancia Allmaras y unos 10 snapshots/ps ya están habilitados por
defecto. Cambiar normalmente el nombre, la corriente y el horizonte:

```bash
python pipelines/02_ss_run_template.py \
  --config "$CFG" --run-name "$SS" \
  --pre-run-name "$PRE" --ss-target-current-uA 30 --ss-time-ps 200
```

Para forzar el horizonte completo aunque la estacionariedad total se cumpla,
añadir `--ss-no-early-stop`. Para aislar física, existen
`--circuit-disable` y `--thermal-disable`.

### Barrido SS isotermal sin circuito: 15 corrientes y 150 ps

Este barrido conserva el PRE, la malla, el cierre Usadel-Poisson, Allmaras y el
paso adaptativo del caso de referencia. La parte termica y el circuito se
desactivan; se fuerza el horizonte completo de 150 ps para comparar todas las
corrientes aunque `photon_ready` se cumpla antes. A las 13 corrientes originales
se agregan 40 y 44 uA para resolver mejor la transicion inmediatamente por
encima del `Ic` Usadel del PRE (aproximadamente 38.8 uA).

Diagnostico de la run de `code3`: la desconexion no fue la causa del fallo. La
maquina no reinicio y, entre 06:20 y 06:32 UTC, el uso de RAM subio de 76.4% a
98.1%, la swap llego a 100%, la carga a 51.5 y el commit de memoria a 103.3%.
Esto coincide con la escritura simultanea de los NPZ: el caso base produjo unos
26 GB y los historiales/snapshots parciales de 25 y 30 uA siguieron creciendo.
El contador del kernel registra OOM kills, aunque el usuario no tiene permiso
para fechar el PID exacto en el journal. La alerta `DatasourceNoData` de HOME y
SCRATCH es consistente con la saturacion de RAM/I/O; un SSH desconectado no
rompe por si solo un `ProcessPoolExecutor` dentro de `screen`. Al morir un hijo,
el pool marco los doce futures como `BrokenProcessPool`.

Geminga tiene 16 nucleos fisicos (32 CPUs logicas). Hay 14 workers extra y el
caso base corre en el proceso principal: son 15 casos/hilos simultaneos y queda
un nucleo fisico nominalmente libre. Los 301 snapshots conservan 2 muestras/ps
y 11 muestras en la cola fisica de 5 ps, pero reducen la retencion al 15% de la
run fallida de 2000 snapshots. El volumen esperado es del orden de 60 GB para
los 15 casos, no cientos de GB residentes/escritos a la vez. Cada proceso limita
BLAS/OpenMP a un hilo y `MALLOC_ARENA_MAX=2` acota arenas de memoria. Usar este
nombre nuevo: la run `_01` es parcial y no se debe reanudar ni sobrescribir.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
MALLOC_ARENA_MAX=2 PYTHONUNBUFFERED=1 \
python pipelines/02_ss_run_template.py \
  --config "$CFG" \
  --run-name "$SWEEP_PREFIX" \
  --pre-run-name "$PRE" \
  --ss-target-current-uA 20 \
  --extra-currents-uA 25 30 33 36 38 39 40 42 44 48 58 74 92 100 \
  --ss-time-ps 150 --ss-snapshots 301 \
  --thermal-disable --circuit-disable --ss-no-early-stop \
  --ss-sweep-workers 14 --ss-threads-per-case 1
```

### Photon-run de hasta 1.5 ns

Los defaults representan el caso central de `0.8 eV`, `sigma=10 nm`,
`t_gamma=50 ps`, evolución térmica/circuital, umbral de `100 uV` y parada tras
recuperación eléctrica sostenida durante `10 ps`. El estado circuital se
hereda de la SS.

```bash
python pipelines/03_photon_run_template.py \
  --config "$CFG" --run-name "$PHOTON" \
  --pre-run-name "$PRE" --ss-run-name "$SS"
```

Las opciones que normalmente se barren son `--photon-energy-eV`,
`--photon-x-nm`, `--photon-y-nm`, `--photon-sigma-nm`,
`--detection-threshold-uV` y `--total-time-ps`. Usar
`--early-stop-mode none` sólo cuando sea necesario conservar todo el horizonte.

## Diagnósticos desde resultados existentes

Estos comandos sólo leen datos ya calculados.

### PRE

```bash
python plot_pipelines/E1_plot_prerun.py --config "$CFG" --pre-run-name "$PRE"
```

Añadir `--with-gap-plot` a E1 cuando se necesite reconstruir la curva térmica
del gap.

### SS individual

E2 usa el PRE declarado en el resumen de la SS y selecciona por defecto los
snapshots más cercanos a `0, 5, 10, 15, 20, 30, 50 ps`.

```bash
python plot_pipelines/E2_phasecg_ss_diagnostics.py \
  --config "$CFG" --run-name "$SS"
```

Para otra ventana, añadir por ejemplo `--times-ps 0 10 25 50 100 200`.

### Photon-run individual

La E3 individual genera evolución escalar y mapas con
`|Delta|/Delta0`, potencial, `|q|xi`, `Te` y `Tph`.

```bash
python plot_pipelines/E3_photon_run_diagnostics.py \
  --config "$CFG" --pre-run-name "$PRE" \
  --run-name "$PHOTON"
```

Los tiempos por defecto de E3 individual son `50, 55, 60, 100 ps`; se cambian con
`--times-ps`.

### Comparación central/lateral

```bash
python plot_pipelines/E3_photon_position_comparison.py \
  --config "$CFG" --pre-run-name "$PRE" \
  --center-run-name "$PHOTON" \
  --edge-run-name "$PHOTON_EDGE"
```

La comparacion usa por defecto `50, 51, 52, 53 ps` y genera mapas emparejados
de campos (incluido `|q|*xi`), potencias, energias/capacidades calorificas,
respuesta circuital y, cuando corresponde, recuperacion censurada. Usar
`--times-ps` para cambiar la ventana y `--output-run-name` cuando se quieran
conservar varias comparaciones.

### Barrido de corriente

```bash
python plot_pipelines/Z2_current_sweep_analysis.py \
  --config "$CFG" --run-prefix "$SWEEP_PREFIX" --stage ss
```

`--delta-inset-currents-uA` y
`--terminal-delta-inset-currents-uA` controlan únicamente los casos mostrados
en los insets.

Z2 usa por defecto un inventario superficial y solo carga el resumen y el
estado final de cada run completa. No abre los historiales ni los snapshots
multigigabyte. Genera las dos curvas I-V y
`Z2_sweep_regime_summary.pdf`/`.yaml`, donde separa runs incompletas de fallos
físicos y resume los rangos muestreados de SS estricta, SS dinámica,
`photon_ready` y comportamiento aproximadamente óhmico (10% por defecto).
`--deep-inventory` queda reservado para una auditoría explícita de claves NPZ;
no usarlo en el barrido normal.


## Memoria

```bash
cd /home/jdiaz/memoria/main
tectonic main.tex --keep-logs --keep-intermediates
grep -nE "Overfull|Underfull|undefined references|LaTeX Error|^!" main.log
pdfinfo main.pdf
cd /home/jdiaz/pysnspd
```

## Screen

```bash
screen -ls
screen -r code1
screen -r code2
```

Dentro de una sesión, `Ctrl-a d` la desconecta sin detener el proceso.


## pySNSPD 1.0.0 y etapa 1 - 2026-09-11

Política vigente: Codex usa la cuenta no administradora `jdiaz`. No lanza
cálculos que previsiblemente superen cinco minutos. Los deja en este archivo
con propósito y salidas, y espera a que el usuario ejecute el comando, copie
su salida y diga `reinicia` o equivalente. No usar sondeos repetidos ni
trabajos en segundo plano para eludir esa entrega. Los comandos históricos
anteriores se conservan; no son una autorización para ejecutarlos ahora.

La etiqueta `v1.0.0` congela el solver de la memoria y la documentación 0.4
con el cuaderno E-r02. La admisión de datos y el catálogo experimental irán
en commits posteriores. No se activa el nuevo condensado ni la cinética.

### Comprobar la entrega documental (ligero)

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/model_v0_4/verify_distribution.py
```

Salida esperada: `status: verified` y hashes del manifiesto y ZIP. No regenera
figuras ni sobrescribe respuestas del cuaderno.

### Pruebas del código publicado

La validación automática de la release usa una copia aislada del código y
un límite de 240 s. Si ese límite se alcanza, la prueba queda incompleta;
el comando sin límite siguiente se reserva al usuario. No es una orden de
ejecución automática ni genera un PRE o transitorio de producción.

```bash
cd /home/jdiaz/pysnspd
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/home/jdiaz/.conda/envs/snspd/bin/python -m pytest -q
```

Para comprobar el paquete sin instalar o descargar sus dependencias:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -m pip wheel . --no-deps --no-build-isolation --wheel-dir /home/jdiaz/pysnspd/tmp/release_wheels
```

Los comandos del catálogo se incorporarán a continuación, con los tamaños
realmente ensayados y sus tiempos medidos. No usar los comandos PRE de arriba
para validar automáticamente el nuevo catálogo experimental.

## Etapa 1 ejecutada y dictamen - 2026-09-11

La etiqueta v1.0.0 quedó publicada en GitHub antes de implementar la etapa:
`5ea0cd6a08d2cae414c82944f8230469ca5aa7d6`. Los siguientes módulos son
experimentales y no se conectan al solver de producción.

Resultado: 232 pruebas aprobadas (57 nuevas), catálogo construido y auditado.
**No promover todavía a etapa 2:** NbN está rechazado y el regulador del piloto
mantiene un sesgo de 8.61% en la pendiente de corriente del caso de baja energía.
El informe está en `output/pdf/implementation/Informe_etapa_1_catalogo.pdf`.
El dictamen verificable está en `docs/implementation/stage1/catalog_admission.json`.

Todos estos comandos se midieron con un hilo en Geminga, usando la cuenta jdiaz.
No hay un cálculo largo pendiente o ejecutándose para esta entrega. Un futuro
catálogo más fino debe estimarse antes de lanzarlo: si se prevén más de cinco
minutos, dejar aquí el comando concreto y esperar la salida manual con `reinicia`.
No relanzar ni fragmentar una corrida para eludir ese límite.

### Preparación común

```bash
cd /home/jdiaz/pysnspd
SNSPD_PY=/home/jdiaz/.conda/envs/snspd/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
```

### Regresión completa - 22.41 s

```bash
timeout --signal=TERM --kill-after=5s 240s "$SNSPD_PY" -m pytest -q
```

La configuración limita la colección a `tests`; las copias de fuentes en `tmp`
no se recogen como otra suite. El resultado esperado de esta entrega es 232 passed.

### Auditar los datos NbN - aproximadamente 0.4 a 2 s

```bash
timeout --signal=TERM --kill-after=5s 240s "$SNSPD_PY" sandbox/stage1_catalog/audit_materials.py --nbn-path /home/jdiaz/scratch/big_data/catalogs/simon_2025/nbn-a2f-ph.dat --verify-upstream
```

Comprueba seis fuentes públicas por SHA-256, conserva los datos y devuelve NbN
REJECTED con las causas; admite únicamente los fixtures sintéticos declarados.
Las referencias de red están fijadas a una revisión. Un fallo de red no es una
verificación de procedencia aprobada.

### Reconstruir pilotos y diagnósticos - 27.37 s

```bash
timeout --signal=TERM --kill-after=5s 240s "$SNSPD_PY" sandbox/stage1_catalog/build_catalog.py
timeout --signal=TERM --kill-after=5s 240s "$SNSPD_PY" sandbox/stage1_catalog/check_boundary_comparison.py
timeout --signal=TERM --kill-after=5s 240s "$SNSPD_PY" sandbox/stage1_catalog/benchmark_queries.py
```

El primer programa genera NPZ en `docs/implementation/stage1/catalogs` y resultados
en `docs/implementation/stage1/electronic`. Los otros tardan menos de un segundo.
Medianas observadas: 0.048 ms por consulta con población y 0.009 ms en el vacío.
Estos comandos sobrescriben los resultados de diagnóstico, no los datos externos
ni un transitorio. Se puede conservar otra corrida con las opciones `--output`
y `--catalogs` del constructor, eligiendo directorios nuevos.

### Referencias independientes - menos de 2 s en conjunto

```bash
"$SNSPD_PY" docs/implementation/stage1/review/independent_vacuum_reference.py
timeout 60 "$SNSPD_PY" docs/implementation/stage1/review/verify_query_contract.py
"$SNSPD_PY" docs/implementation/stage1/review/low_energy_regulator_reference.py
```

El último comando es especialmente útil para comprobar el bloqueo actual:
separa el sesgo del regulador del error de interpolación y de cuadratura. En
la población de baja energía probada, eta/Delta0=1e-3, 1e-4 y 1e-5 dejan 8.61%,
2.99% y 0.986% de sesgo en la pendiente j/q, respectivamente. Bajar eta sobre
la misma malla global de 64 nodos no certifica ese límite.

### Componer y verificar el informe - pocos segundos

Las dependencias de PDF están aisladas en `tmp/stage1_report_deps` (ReportLab y
pypdf); no modifican el entorno de simulación. Para restaurarlas si se elimina
esa carpeta: `"$SNSPD_PY" -m pip install --target tmp/stage1_report_deps reportlab pypdf`.

```bash
"$SNSPD_PY" sandbox/stage1_catalog/prepare_report.py
PYTHONPATH=/home/jdiaz/pysnspd/tmp/stage1_report_deps "$SNSPD_PY" sandbox/stage1_catalog/build_report.py
"$SNSPD_PY" sandbox/stage1_catalog/verify_delivery.py
```

El verificador detecta cualquier archivo distinto del manifiesto entregado.
Regenerar resultados exige revisar el dictamen y las tres páginas del PDF antes
de crear un manifiesto nuevo con `verify_delivery.py --write`. No sobrescribir
un dictamen de admisión únicamente para hacer pasar su verificación.

Siguiente trabajo recomendado: cuadratura concentrada cerca del borde, malla de
Gamma que resuelva su escala estrecha y convergencia conjunta del regulador;
además, obtener datos fonónicos positivos, unívocos y documentados. Esa ampliación
no se ha implementado ni ejecutado como un transitorio en esta etapa.


## Etapa 1 R2 — catálogo electrónico aceptado en el contrato muestreado (2026-09-11)

R2 pasa 278 pruebas y admite la siguiente validación en una o dos celdas
sintéticas; la entrada fonónica NbN continúa REJECTED. El informe y el certificado
están en `docs/implementation/stage1_r2/`. No se ejecutó un transiente completo.
La etiqueta v1.0.0 permanece intacta. Los comandos siguientes son breves y
reproducibles; cada ejecución incierta conserva el límite de 240 segundos. Si
vence, conservar el resultado incompleto y ejecutar sólo por decisión del usuario
la versión sin límite: no repetir ni subdividir para eludir la entrega de cómputo.

Preparar una terminal propia (los resultados nuevos van a `tmp`, sin reemplazar
el certificado entregado):

```bash
cd /home/jdiaz/pysnspd
export PY_R2=/home/jdiaz/.conda/envs/snspd/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export R2_RUN=tmp/stage1_r2_reproduction
mkdir -p "$R2_RUN"
export R2_CATALOG=docs/implementation/stage1_r2/catalogs/occupation_catalog.npz
```

Verificar archivos entregados, regresión completa (~34 s), física (~6 s), orden
cúbico (~2 s), interpolación (~8 s) y cuadratura realmente guardada (~1 s):

```bash
"$PY_R2" sandbox/stage1_catalog_r2/verify_delivery.py
timeout --signal=TERM --kill-after=5 240 "$PY_R2" -m pytest -q
timeout --signal=TERM --kill-after=5 240 "$PY_R2" docs/implementation/stage1_r2/review/assess_candidate.py "$R2_CATALOG" --output "$R2_RUN/assessment.json"
timeout --signal=TERM --kill-after=5 240 "$PY_R2" sandbox/stage1_catalog_r2/check_shape_cells.py "$R2_CATALOG" --output "$R2_RUN/shape.json"
timeout --signal=TERM --kill-after=5 240 "$PY_R2" sandbox/stage1_catalog_r2/check_catalog_interpolation.py "$R2_CATALOG" --output "$R2_RUN/interpolation.json"
timeout --signal=TERM --kill-after=5 240 "$PY_R2" sandbox/stage1_catalog_r2/check_catalog_quadrature.py "$R2_CATALOG" --output "$R2_RUN/quadrature.json"
timeout --signal=TERM --kill-after=5 240 "$PY_R2" sandbox/stage1_catalog/benchmark_queries.py --vacuum docs/implementation/stage1_r2/catalogs/vacuum_catalog.npz --occupation "$R2_CATALOG" --output "$R2_RUN/query_benchmark.json"
```

Reconstrucción numérica del catálogo en tres ejecuciones con resultados útiles
propios: base ~144 s; insertar seis nodos ~5 s; insertar el último ~1 s. Los nodos
anteriores se reutilizan exactamente. No es una partición de un cálculo que haya
agotado un límite; surgió del diagnóstico de interpolación de cada candidato.
No ejecutar el constructor con una malla mayor por defecto sin estimar su coste.

```bash
timeout --signal=TERM --kill-after=5 240 "$PY_R2" sandbox/stage1_catalog_r2/build_catalog.py --delta-nodes 17 --gamma-low-nodes 49 --gamma-high-nodes 17 --count-order 10 --count-cutoff 12 --eta 1e-8 --ratio-coordinate --output "$R2_RUN/base"
timeout --signal=TERM --kill-after=5 240 "$PY_R2" sandbox/stage1_catalog_r2/refine_catalog.py "$R2_RUN/base/occupation_catalog.npz" --ratios 1.463897380883855e-6 5.500371809126641e-6 1.3293713801755863e-5 4.9949108470627004e-5 2.9176708501341483e-4 0.2335109071328366 --output "$R2_RUN/refined.npz" --report "$R2_RUN/refined.json"
timeout --signal=TERM --kill-after=5 240 "$PY_R2" sandbox/stage1_catalog_r2/refine_catalog.py "$R2_RUN/refined.npz" --ratios 0.0002612940447849904 --output "$R2_RUN/final.npz" --report "$R2_RUN/final.json"
timeout --signal=TERM --kill-after=5 240 "$PY_R2" docs/implementation/stage1_r2/review/assess_candidate.py "$R2_RUN/final.npz" --output "$R2_RUN/rebuilt_assessment.json"
```

Los metadatos de procedencia pueden dar un hash NPZ diferente al reconstruir;
hay que comparar contenido numérico y volver a evaluar el archivo nuevo. El
certificado de R2 admite exclusivamente el archivo entregado y su fuente por hash.

Comparar el tratamiento legado del archivo fonónico de Simón (~1 s):

```bash
timeout --signal=TERM --kill-after=5 240 "$PY_R2" sandbox/stage1_catalog_r2/material_comparison.py --nbn-path /home/jdiaz/scratch/big_data/catalogs/simon_2025/nbn-a2f-ph.dat --output-root "$R2_RUN/material"
```

Las siguientes recetas regeneran las figuras y el PDF EN sus rutas publicadas.
Guardar primero cualquier versión que se quiera conservar. Cambiar el archivo
publicado invalida su manifiesto anterior; no recrear éste sin revisar los datos
y renderizar todas las páginas. La reproducción visual puede variar por las
fuentes instaladas en Windows/Linux. No cambia el dictamen físico.

```bash
timeout --signal=TERM --kill-after=5 240 "$PY_R2" sandbox/stage1_catalog_r2/electronic_comparison.py --eta 1e-8
"$PY_R2" sandbox/stage1_catalog_r2/plot_acceptance.py docs/implementation/stage1_r2/review/final_assessment.json
PYTHONPATH=/home/jdiaz/pysnspd/tmp/stage1_report_deps "$PY_R2" sandbox/stage1_catalog_r2/build_report.py
```

Próximo paso recomendado: conservación de energía, relajación al equilibrio y
convergencia temporal en una o dos celdas sintéticas (D.4.2), manteniendo el
rechazo material y el límite de soporte del catálogo. No hay un cálculo pesado
pendiente de ejecución por el usuario al cerrar R2.


## Documento F-r01: cambios respecto de la memoria (2026-09-13)

Entrega PDF de seis paginas, con 22 bloques de ecuaciones viejas/nuevas,
motivos y una seccion de cambios menores, numericos y de codigo. Resume el
candidato 0.4 frente a la memoria; diferencia propuesta, catalogo R2 y produccion.
Solo se agrega F al directorio de A-E; no se crea una version Markdown del documento.
Las seis paginas se inspeccionaron visualmente; ecuaciones F.1-F.22 completas.
A-E conservan sus hashes anteriores. No se ejecutaron simulaciones.

Archivo: /home/jdiaz/pysnspd/output/pdf/modelo_v0_4/F_resumen_de_cambios_respecto_a_la_memoria_v0_4.pdf
SHA-256: cfc6315a60d9913f35715f967d15944a7cbd47107cbb0f23afb45277a9cc3b34
Tamano: 105305 bytes. Memoria cotejada: /home/jdiaz/memoria/main/memoria_02.pdf,
SHA-256 50a75f1bd84f06f32820dbf9477c0fdcaf809fb502485e84249ebb46b286c1e8.

Diagnostico reproducible (lectura, 1 CPU, menos de 30 s, memoria despreciable):
```bash
timeout 30s sha256sum /home/jdiaz/pysnspd/output/pdf/modelo_v0_4/F_resumen_de_cambios_respecto_a_la_memoria_v0_4.pdf
```
Salida esperada: el SHA-256 anterior y la ruta del PDF. No genera nuevos archivos.


## Cierre de etapa 1: celdas y NbN derivado (2026-09-21)

El cierre admite el uso electrónico reducido del catálogo R2: 53 puertas pasan,
317 pruebas de regresión pasan en 34,69 s. El derivado NbN tiene admisión de forma
experimental condicionada; no acredita unidades, normalización ni tasas absolutas.
La cinética electrón-fonón simultánea sigue siendo la siguiente implementación.
No se lanzó un transiente completo ni hay un cálculo pesado pendiente del usuario.

Preparación (cuenta jdiaz, sin sudo):

```bash
cd /home/jdiaz/pysnspd
export CLOSURE_PY=/home/jdiaz/.conda/envs/snspd/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export CLOSURE_RUN=tmp/stage1_closure_reproduction
mkdir -p "$CLOSURE_RUN"
```

Integridad de las dos entregas y regresión (~35 s, un hilo, memoria estimada
menor de 1 GiB; salida: hashes verificados y resumen pytest):

```bash
"$CLOSURE_PY" sandbox/stage1_catalog_r2/verify_delivery.py
"$CLOSURE_PY" sandbox/stage1_closure/verify_delivery.py
timeout --signal=TERM --kill-after=5 240 "$CLOSURE_PY" -m pytest -q
```

Reevaluar el catálogo retenido (~1 s, un hilo, menos de 1 GiB estimado):

```bash
timeout --signal=TERM --kill-after=5 240 "$CLOSURE_PY" docs/implementation/stage1_r2/review/assess_candidate.py docs/implementation/stage1_r2/catalogs/occupation_catalog.npz --output "$CLOSURE_RUN/catalog_recheck.json"
```

Repetir celdas (~2,15 s, un hilo, menos de 1 GiB estimado): BGK, calentamiento,
trabajo espectral, recuperación del condensado y transporte a energía compartida.
Salidas: cells_results.json, cuatro CSV y dos figuras PNG/PDF. Los resultados
nuevos quedan en tmp; no reemplazan los artefactos aceptados.

```bash
timeout --signal=TERM --kill-after=5 240 "$CLOSURE_PY" sandbox/stage1_closure/run_cells.py --output "$CLOSURE_RUN/cells"
```

Revisar fuente y recorte fonónico (~1,79 s medidos, más latencia de Internet;
un hilo, menos de 1 GiB estimado). Descarga la revisión pública fijada, compara
el cuerpo numérico original, deriva soporte común y calcula sensibilidades.
Salidas: material_results.json, material_nbn_shape_v1.csv, manifiesto de filas,
material_thermal.csv y figures/material_{support,impact}_review en PNG/PDF.
No renormaliza la DOS; la hipótesis THz queda explícita en los resultados.

```bash
timeout --signal=TERM --kill-after=5 240 "$CLOSURE_PY" sandbox/stage1_closure/review_material.py --nbn-path /home/jdiaz/scratch/big_data/catalogs/simon_2025/nbn-a2f-ph.dat --output-root "$CLOSURE_RUN/material" --source-dir "$CLOSURE_RUN/material_sources" --online
```

Ver sólo las pruebas nuevas de operadores/preprocesado (~0,3 s de pruebas,
más arranque de pytest; un hilo, menos de 1 GiB estimado):

```bash
timeout --signal=TERM --kill-after=5 240 "$CLOSURE_PY" -m pytest tests/test_experimental_cell_validation.py tests/test_experimental_material_preprocessing.py -q
```

Para regenerar el informe (segundos, un hilo, menos de 1 GiB estimado), la receta
siguiente escribe EN las rutas publicadas el PDF y Markdown. Guardar la versión
anterior si se desea conservar. Deben revisarse todas las páginas antes de
actualizar un manifiesto: fuentes de Windows/Linux pueden cambiar la maquetación.

```bash
PYTHONPATH=/home/jdiaz/pysnspd/tmp/stage1_report_deps "$CLOSURE_PY" sandbox/stage1_closure/build_report.py
```

Resultado: output/pdf/implementation/Informe_cierre_etapa_1.pdf y
docs/implementation/stage1_closure/Informe_cierre_etapa_1.md. El paquete entregado
se comprueba con verify_delivery.py; --write sólo corresponde a una entrega
nueva ya revisada, no a la reproducción ordinaria.

La siguiente secuencia está en docs/implementation/stage1_closure/NEXT_STAGE.md.
Toda prueba nueva prevista de más de cinco minutos se dejará aquí con sus
recursos y salidas para ejecución del usuario. Si un diagnóstico limitado a
240 segundos agota el tiempo, queda incompleto; no se reinicia en fragmentos.

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


## 2026-09-21 - Etapa 2: referencia recuperada y lote de convergencia RK4

La referencia manual de code_000 terminó: 1164,55 s y 6843 RHS. Sus archivos
originales siguen en tmp/stage2_user_20260921_193252 y tienen copia verificada
en docs/implementation/stage2/resume_20260921/manual_reference. No repetirla.
RK4 con 10/20/40 pasos pasó la comparación corta: error máximo fino 1,57e-6,
reducciones 3,70/4,15 y 24,94 s frente a 1164,55 s. Los kernels y umbrales
siguen iguales. La etapa completa aún está pendiente.

### Ejecutar manualmente: lote único de validaciones pendientes

El agente NO ejecutó este lote. Está preregistrado en
docs/implementation/stage2/resume_20260921/manual_validation_plan.json.
Primero compara tres resoluciones temporales con una cuarta trayectoria de
paso reducido a la mitad, para una y dos celdas hasta t=2. Si pasan, evalúa
las mallas electrónicas 630/1260/2520 y fonónicas 1025/2049/4097 en el caso de
dos celdas, más equilibrio, fronteras de ocupación, campos y soporte.

Recursos estimados: un hilo de CPU, reservar 8 GiB de RAM, unas 2-3 horas
(aproximadamente 150 minutos; puede excederse). La red de 2520 estados domina
el coste. Los datos de cada trayectoria se guardan antes de continuar. No es
un transiente espacial de SNSPD ni una calibración de NbN.

Ejecutar en la terminal de screen, como jdiaz, desde el repositorio:

```bash
cd /home/jdiaz/pysnspd
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
STAGE2_BATCH=/home/jdiaz/pysnspd/tmp/stage2_validation_$(date +%Y%m%d_%H%M%S)
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage2_cells/resume_20260921/run_acceptance_batch.py --plan docs/implementation/stage2/resume_20260921/manual_validation_plan.json --output-root "$STAGE2_BATCH" --execute
```

El proceso permanece en primer plano; screen permite desconectar la terminal.
Hay progreso cada 100 evaluaciones y mensajes por tarea. Reutiliza el piloto
TWO40 existente sólo tras verificar parámetros, fuentes, versiones y hashes.
Se detiene ante el primer fallo; no reintenta ni sobrescribe resultados. Una
interrupción deja evidencia incompleta y requiere diagnóstico antes de seguir.

Salidas: batch_manifest.json, batch_events.jsonl, recibos de tareas, JSON,
NPZ y registros .console.log/.progress.jsonl en STAGE2_BATCH. Al terminar o
detenerse, escribir «reinicia» e indicar la screen si cambió. El agente puede
leer los archivos y la screen directamente, sin pegar el output al chat.
BATCH_COMPLETE sólo significa que las tareas listadas pasaron: todavía se
revisan actividad, estacionariedad y el dictamen final antes de cerrar etapa 2.

### Vista previa sin cálculos

No crea carpetas ni lanza subprocesos/RHS; valida hashes y muestra las tareas:

```bash
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage2_cells/resume_20260921/run_acceptance_batch.py --plan docs/implementation/stage2/resume_20260921/manual_validation_plan.json --output-root tmp/stage2_preview --dry-run
```

### Diagnósticos ligeros reutilizables

Suavidad del RHS: 17 evaluaciones estáticas, 2,65 s medidos, menos de 1 GiB.
No integra una trayectoria. Elegir una salida nueva:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout --signal=TERM --kill-after=5 240 /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage2_cells/resume_20260921/diagnose_rhs_smoothness.py --initial tmp/stage2_user_20260921_193252/one_reference_probe_initial.npz --output tmp/rhs_smoothness_recheck.json
```

Comparación continua seleccionada: 3,65 s medidos; sólo integra las referencias,
sin reconstruir redes ni trayectorias. La reproducción reemplaza su JSON de
evidencia y cambia su hash por los tiempos medidos; conservar antes una copia.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout --signal=TERM --kill-after=5 240 /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage2_cells/resume_20260921/close_selected_continuum.py
```
