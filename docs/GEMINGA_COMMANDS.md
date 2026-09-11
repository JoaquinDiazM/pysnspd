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
