## Actualización RHS y balance — 5 de octubre de 2026

La continuación necesaria está ahora en `/home/jdiaz/GEMINGA_COMMANDS_PRESENT.md`. Hay cuatro comandos manuales de auditoría; el estado anterior sin corrida pendiente queda como historia. No repetir DC ni ejecutar recuperación eléctrica total.

# Geminga: continuación vigente del modelo — 5 de octubre de 2026

## Estado y siguiente paso

**DC_PREPARATION_ADMITTED revalidado.** La campaña del 24–25 de septiembre
ya terminó: cuatro referencias, cinco holds de 10 ps y comparaciones de
malla, longitud y paso admitidas. No repetir el comando DC que se conserva
más abajo como historia. El conjunto tomó unas 4,50 h según recibos.

Se implementaron un observador de energía interna de la acción térmica finita,
el controlador gatillo más margen y la guarda de dirección DC global.
**92 pruebas y 11 subpruebas pasaron en 0,34 s**. Esto no completa el integrador
no térmico ni admite aún calor finito, tasas NbN o preparación fotónica.

**No hay una corrida pesada pendiente definida en esta entrega.** El próximo
paso es completar el RHS y su balance independiente, no repetir DC ni iniciar
el solver térmico heredado como si fuera el modelo mejorado.

Informes vigentes del repo `/home/jdiaz/pysnspd`:

- `docs/implementation/stage5/dc_admission_20261005/README.md`
- `docs/implementation/stage5/continuation_20261005/README.md`
- `docs/implementation/stage5/continuation_20261005/thermal_internal_energy.md`
- `docs/implementation/stage5/continuation_20261005/experimental_trigger_scale.md`

## Diagnósticos ligeros opcionales, sin simulación

Revalidar datos existentes sin cambiar informes originales (4,58 s medidos):

```bash
cd /home/jdiaz/pysnspd
timeout 240s /home/jdiaz/.conda/envs/snspd/bin/python \
  sandbox/stage5_prephoton/verify_completed_dc.py \
  --output-root /home/jdiaz/scratch/stage5_prephoton_dc_20260924
```

Verificar piezas nuevas y compatibilidad de las configuraciones de tiempo:

```bash
cd /home/jdiaz/pysnspd
timeout 240s /home/jdiaz/.conda/envs/snspd/bin/python \
  sandbox/stage5_prephoton/check_continuation.py
```

Ambos comandos consultan CPU/afinidad/cuotas/RAM a cada ejecución y admiten
un CPU y un hilo numérico dentro del 90 % disponible. La prueba reserva 1 GiB;
coordinador y subproceso de pruebas comparten ese CPU. Los checks terminan
sin trabajos en fondo y los tests internos tienen timeout de 120 s. No generan
nuevos estados físicos ni modifican los recibos versionados sin pedir una
ruta nueva explícita con `--receipt`.

## Contrato futuro del transiente corto

- Máximo **200 ps sin gatillo confirmado**, desde la transferencia al modelo.
  Un cruce confirmado permite completar su margen aunque termine después del
  techo. No se espera el pico ni la recuperación eléctrica completa.
- Vout en SI con presentación mV. Rango provisional de estudio **0,1–0,5 mV**,
  inspirado por escalas de dispositivo/lectura, no calibración experimental.
  Umbral definitivo, confirmación, histéresis y margen siguen abiertos.
- Guardar Ib e Is globales firmados y detener la propuesta antes de aceptarla
  si invierte la dirección DC. Preservar último estado aceptado y diagnóstico;
  no recortar ni rectificar corrientes ni confundir inversión local con global.
- Conservar física, circuito CM completo, material, depósito, malla y bordes.
  Fijar Lext una sola vez; no usar Lk(t) ni acortar constantes de recuperación.
- El futuro runner debe guardar estado físico completo, reloj y depósito,
  poblaciones y espectros, circuito, procedencia y controlador, más diagnóstico
  de cada paso aceptado. El checkpoint heredado no satisface ese contrato.

Cuando haya una campaña válida de más de cinco minutos, su comando exacto,
recursos, propósito, salidas y tiempo estimado se agregarán aquí y al chat
para ejecución manual del usuario. Esta entrega no prepara comandos de
transiente completo con reset eléctrico total.

## Historial íntegro de entregas previas

La primera entrada histórica dice «pendiente» porque fue escrita antes de su
corrida. Su comando DC ya fue ejecutado y admitido; no es la cola vigente.
Se conserva el contenido previo completo a continuación.

# Geminga: preparación DC de una sección interior antes del fotón

## 2026-09-24 — Comando vigente: equilibrios intrínsecos y conservación DC

**Pendiente de ejecución manual.** Prepara cuatro referencias con bordes del
mismo superconductor portador de corriente, y cinco trayectorias de 10 ps con
fuente DC y circuito CM completo. No añade sonda AC ni fotón. Lext se calcula
una sola vez por referencia: 96 nH del inductor añadido más el tramo activo
omitido, sin inductancia circuital dependiente del tiempo. Compara corriente,
longitud, malla y paso temporal. No hace falta configurar ninguna screen.

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_prephoton/run_prephoton_dc.py \
  --plan docs/implementation/stage5/prephoton_dc_20260924/campaign.json \
  --output-root /home/jdiaz/scratch/stage5_prephoton_dc_20260924 --execute
```

**Tiempo orientativo: 10–24 horas para el conjunto**, sujeto al coste de las
mallas mayores y a cuántas veces deba renovarse el Jacobiano. Es una proyección,
no un tiempo medido de esta campaña. El piloto de 1712 nodos/256 frecuencias
completó 50 pasos en 7,26 s con 12 trabajadores. La primera versión reconstruía
Newton en cada paso (~4,7 s/paso) y se detuvo al límite de 180 s; la versión
actual reutiliza las factorizaciones y verifica el mismo residuo no lineal.
La ETA se actualiza con el progreso real; la de preparación estática indica
el presupuesto de iteraciones restante, que puede terminar mucho antes.

**Recursos:** hasta 28 CPU lógicas en el inventario actual de 16 núcleos/32
hilos; dos núcleos completos libres. Dos casos simultáneos, con 13 y 12
trabajadores espectrales, dos coordinadores de caso y uno global. Un hilo
BLAS/OpenMP por proceso, afinidades disjuntas y reducción de concurrencia si
RAM/cuotas lo exigen. Reserva estimada por el planificador ≈34 GiB para dos
casos; los factores LU residentes están dentro de ese presupuesto. Los datos
se guardan en scratch; las referencias no se duplican en Git.

**Salidas:** `references/cases/*/reference.npz`, `summary.json`; después
`holds/cases/*/{summary.json,history.json,fields.npz}`. Al finalizar:
`analysis.md`, `analysis.json`, `figures/*.png` y `workflow_result.json`.
El informe incluye perfiles del gap, mapas, deriva temporal, Vdev y Vout en
paneles separados, además de la partición inductiva. El destino debe ser nuevo.
Los casos independientes continúan si otro falla; sus dependientes se bloquean.
No se reintenta, sobrescribe ni ajustan tolerancias silenciosamente.

**Alcance:** verifica el fondo DC térmico y el circuito. No certifica por sí
solo el balance de calor fuera del equilibrio ni la transferencia del fotón.
Los detalles y márgenes están en
`docs/implementation/stage5/prephoton_dc_20260924/README.md`.

Comprobación útil y ligera de entradas, versiones y recursos, sin simulación:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage5_prephoton/run_prephoton_dc.py \
  --plan docs/implementation/stage5/prephoton_dc_20260924/campaign.json \
  --output-root /home/jdiaz/scratch/stage5_prephoton_dc_20260924
```

Rehacer sólo las figuras y el análisis de datos ya terminados, segundos:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage5_prephoton/analyze_prephoton_dc.py \
  --output-root /home/jdiaz/scratch/stage5_prephoton_dc_20260924
```

Avisar al terminar; no es necesario pegar los registros extensos en el chat.

## Historial: cierre de etapa 4, previo a las observaciones sobre el bulk DC

## Entrega actual: cierre físico y preparación de etapa 5

Dos referencias y cinco respuestas completadas en 48,57 min. El informe final
admite desarrollo térmico y respuesta reactiva condicionada; disipación absoluta
y continuidad interior siguen pendientes antes de un fotón. No hay otro comando
largo requerido ahora. Informe: docs/implementation/stage4/closure_20260924/README.md.

Posprocesamiento opcional, segundos, sin simulaciones nuevas:

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage4_core/analyze_final_closure.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage4_core/plot_final_closure.py
```

Los datos completos siguen en /home/jdiaz/scratch/stage4_final_coupling_20260924.
El informe PDF vigente es output/pdf/implementation/Informe_cierre_etapa_4.pdf.
Los PDF intermedios y 15 copias NPZ exactas se retiraron del repo; véase
docs/implementation/ARCHIVOS_ARCHIVADOS.md. Recuperación opcional segura:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/maintenance/restore_archived_artifacts.py --list
```


## Historial: unión dinámica polarizada completada el 24 de septiembre

Comando histórico ya completado; no repetir para este cierre. No repetir Euler dual ni las campañas anteriores.
Construye dos referencias con corriente y luego cinco respuestas acopladas de
amplitud, fase, poblaciones y potencial sobre la malla dual, con el circuito
completo de la memoria. La comparación incluye un refinamiento independiente.

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/run_final_coupling.py \
  --plan docs/implementation/stage4/coupled_closure_20260924/final_coupling_plan.json \
  --output-root /home/jdiaz/scratch/stage4_final_coupling_20260924 --execute
```

**Tiempo orientativo:** una a varias horas; no es un plazo garantizado. Una
consulta de espectro/cinética con 13 perturbaciones tardó 1,92 s en la malla
completa sobre fondo uniforme. La corriente, las proximidades del gap y el
refinamiento aumentan ese coste. La ETA usa el avance de la fase activa.

**Recursos:** hasta 28 hilos de los 32 detectados, con dos núcleos físicos
reservados. Dos casos simultáneos comparten ese límite; nunca se asigna 90 % a
cada caso por separado. Un hilo BLAS/OpenMP por proceso. Reserva de RAM
aproximada de 34 GiB para dos casos, revisada contra la disponibilidad real;
disco esperado del orden de unos GB, dependiente de los registros espectrales.

**Salidas:** `references/cases/*/reference.npz`,
`responses/cases/*/results.json`, `coupled_fields.npz`,
`integrated_moments.npz`, registros por energía y `workflow_result.json`.
Todo queda bajo el destino del comando. Barras/ETA aparecen en la terminal,
y la salida detallada permanece en los archivos `logs/*.console.log` de cada
campaña. No hace falta copiar el registro al chat: basta avisar que terminó.

El destino debe ser nuevo. Si falla un caso, se guardan sus resultados ya
obtenidos, se cancelan sus propios trabajadores y los casos independientes
siguen. No hay reintento ni reinicio automático a partir de sumas parciales.
Una referencia incompleta sólo bloquea sus respuestas dependientes.

La integral de energía incluye las colas semiinfinitas mediante `x=C/E`;
el valor C separa las dos cuadraturas y no elimina estados de alta energía.
Se mantiene el margen práctico del 2 % en observables relevantes. El estado
de ejecución `COMPLETED` no declara por sí solo aprobado el cierre físico.

**Preflight opcional, segundos y sin cálculo físico:**

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage4_core/run_final_coupling.py \
  --plan docs/implementation/stage4/coupled_closure_20260924/final_coupling_plan.json \
  --output-root /home/jdiaz/scratch/stage4_final_coupling_20260924
```

**Diagnóstico ligero reproducible opcional:** una consulta sobre los 1712
nodos y una tangente térmica, sin recalcular el lote. Usar un destino nuevo.

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 120 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/pilot_biased_coupling.py \
  --output /home/jdiaz/scratch/stage4_coupling_pilot_reproduction
```

El piloto original ya terminó en 2,21 s y no necesita repetirse. El modelo
de producción y v1.0.0 no cambian; aún no se simula un fotón. El informe
de alcance y criterios está en
`docs/implementation/stage4/coupled_closure_20260924/README.md`.

También se completó el flujo real sobre una malla pequeña en 31,92 s, con
referencia aceptada y 1112 energías. Su comando y los planes exactos están en
`docs/implementation/stage4/coupled_closure_20260924/pipeline_smoke/integration_receipt.json`.
Es un diagnóstico de integración ya terminado; no es necesario repetirlo
antes del lote principal. Para una reproducción opcional debe elegirse otro
`--output-root` y conservar la afinidad indicada en ese recibo.

---

## Registro anterior conservado

# Geminga: Euler dual completado; no hay otra ejecución pendiente

Actualizado el 24 de septiembre de 2026. Cuenta jdiaz, sin administrador.

La campaña optimizada completó 1 ps en 153,81 s dentro del intento acotado de 240 s.
Sondas de amplitud/fase, dos pasos y balance integrado cumplen el margen del 2 %.
No repetir: los resultados están en
`/home/jdiaz/scratch/stage4_dual_kwt_time_20260924`.
La campaña usó 27 trabajadores + coordinador, 28 CPU lógicas y dejó dos núcleos
físicos libres; BLAS/OpenMP a un hilo. No hubo correcciones Newton adicionales.

El control longitudinal también terminó, en 0,078 s, en
`/home/jdiaz/scratch/stage4_dual_longitudinal_20260924`.
La investigación cuasiclásica está en
`docs/implementation/stage4/final_kwt_20260924/quasiclassical_assessment.md`.
Recomienda conservar las distribuciones dinámicas y revisar por separado la
aproximación espectral adiabática. Sigue pendiente su unión compatible con
trabajo espectral, fase y potencial. No se solicita otro barrido numérico.

Diagnóstico algebraico opcional (segundos, un CPU; no simula transientes):

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage4_core/quasiclassical_scales.py
```

Salida: `docs/implementation/stage4/final_kwt_20260924/scales/scales.json` y
`scale_hierarchy.png`. Recalcula escalas declaradas; no ajusta constantes.
Los artefactos ya están guardados; no es una ejecución requerida.

## Consultas ligeras y reproducción opcional

Lectura de resultados completados, sin nueva simulación:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -m json.tool docs/implementation/stage4/final_kwt_20260924/thermal/summary.json
```

Comando completo conservado sólo para reproducir el ensayo voluntariamente.
Coste observado 2,56 min, variable con carga del equipo; máximo 90 % CPU/RAM
detectados, progreso/ETA y fuentes registradas. No modificar constantes
físicas para reducir su duración. Salidas: observaciones, historia temporal,
balance, refinamiento y checkpoints, dentro de una carpeta nueva.

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/dual_kwt_predicted.py \
  --plan docs/implementation/stage4/final_kwt_20260924/dual_kwt_predicted_plan.json \
  --mesh docs/implementation/stage4/practical_time_review_20260924/dual_mesh/resampled/mesh.npz \
  --output-root /home/jdiaz/scratch/stage4_dual_kwt_manual_repeat \
  --execute
```

No usar nuevamente la carpeta de una ejecución existente. Los pilotos sin
predictor y con predictor permanecen conservados; no son tareas pendientes.

## Historial anterior preservado

Las peticiones de ejecución contenidas más abajo son históricas.

# Geminga: no hay otra corrida larga solicitada

Actualizado el 24 de septiembre de 2026. Cuenta jdiaz, sin administrador.

## Resultado de la ejecución ETD2

Ambas trayectorias llegaron a 1 ps en 1089,49 s. La interrupción fue posterior:
una comparación del torque de fase secundario de la sonda de amplitud no pasó.
Pasan 95/96 comparaciones; en tiempos con pasos realmente diferentes, 23/24.
La diferencia que falla es 8,30465×10⁻⁸ en la norma adimensional declarada:
1,500 % de su pequeño valor inicial y 0,004566 % del torque angular principal.
Ese contexto no cambia el certificado original ni acredita precisión de la cola.

No repetir la campaña ni sus pilotos. Se conserva todo en
`/home/jdiaz/scratch/stage4_nonlinear_time_20260924`. El ensayo térmico se cierra
como desarrollo con ese límite y se avanza a conectar la acción con la malla
dual y el paso KWT de la memoria. La etapa4 completa sigue abierta.

## Lecturas ligeras útiles

```bash
cd /home/jdiaz/pysnspd
tail -n 2 /home/jdiaz/scratch/stage4_nonlinear_time_20260924/progress.jsonl
/home/jdiaz/.conda/envs/snspd/bin/python -m json.tool docs/implementation/stage4/practical_time_review_20260924/decision.json
```

El [informe vigente](https://github.com/JoaquinDiazM/pysnspd/blob/c80c0f8612d8fd5a2b390938ca5ed4b9d58d3caa/output/pdf/implementation/Informe_etapa_4_revision_practica.pdf)
declara qué campo, norma, unidad, sonda y referencia se representan. Las
[fuentes consultadas](implementation/stage4/practical_time_review_20260924/published_methods.md)
explican qué código publicado y heredado se reutiliza. No hay otro cálculo
largo necesario para revisar este resultado. El próximo comando se registrará
sólo cuando corresponda al nuevo ensayo acoplado, no a repetir esta cola.

## Reproducción opcional de las comprobaciones ligeras

Ya completadas; no son tareas pendientes. Usan un solo proceso, bibliotecas a
un hilo y normalmente tardan segundos. La primera imprime el resultado de las
tres pruebas del adaptador dual. La segunda compara el paso KWT heredado con la
ley continua a tres pasos decrecientes y guarda el recibo indicado; no evoluciona
el dispositivo ni genera un depósito fotónico.

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -m unittest discover -s tests -p 'test_thermal_mesh_adapter.py'
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -m sandbox.stage4_core.kwt_bridge_check \
  --output tmp/kwt_bridge_repeat_20260924.json
```

La malla preparada para el siguiente ensayo está en
`docs/implementation/stage4/practical_time_review_20260924/dual_mesh/resampled/mesh.npz`.
Tiene 1712 nodos y contactos sólo en los extremos; no hace falta regenerarla.

## Historial preservado — comandos anteriores completados o sustituidos

Los comandos siguientes son históricos. Sus títulos de «ejecución vigente»
no constituyen una solicitud nueva. No hay que relanzarlos.

# Geminga: ejecución vigente — trayectoria térmica no lineal

Actualizado el 24 de septiembre de 2026. Cuenta jdiaz, sin administrador.

## Resultados ya recuperados: no repetir

La trayectoria térmica afín terminó en 71,10 s y pasa su refinamiento. El contraste
no lineal de campos guardados terminó en 22,74 s con 2048 raíces y pasa el margen
registrado. Se corrigió cancelación de redondeo en la comparación energética
de Newton, conservando la misma ecuación y tolerancia. El fallo previo permanece
en `/home/jdiaz/scratch/stage4_nonlinear_snapshots_20260924`; el resultado corregido
en `/home/jdiaz/scratch/stage4_nonlinear_stable_20260924`.

## Siguiente ejecución: ETD2 térmico no lineal hasta 1 ps

El piloto de 0,001 ps ya terminó en 48,95 s con 4096 raíces espectrales;
no repetirlo. Comprueba el arranque del mismo sistema y sus contactos, sin
acreditar todavía la trayectoria completa ni su refinamiento temporal.

Ejecutar una vez desde la terminal habitual:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/thermal_nonlinear_time.py \
  --plan docs/implementation/stage4/time_review_20260924/nonlinear_time/plan.json \
  --operator-root /home/jdiaz/scratch/stage4_thermal_weak_20260924 \
  --output-root /home/jdiaz/scratch/stage4_nonlinear_time_20260924 \
  --execute
```

Propósito: verificar si la evolución de referencia y las sondas débiles se
mantienen al evaluar el espectro y la fuerza no lineales en cada paso. El
integrador ETD2 reutiliza J0 para la parte rígida y añade el resto no lineal.
Conserva malla, contactos, 256 frecuencias, temperatura, tiempos KWT y potencial.
No cambia el sistema físico para acelerarlo.

Tiempo provisional: 5–30 min; depende del control adaptativo. Máximo 27 trabajadores
más coordinador, un hilo BLAS/OMP por proceso, presupuesto 90 % CPU/RAM, dejando
dos núcleos físicos completos libres. Los casos usan el mismo pool; cada
evaluación reparte las frecuencias independientes. Consola: barra, avance físico,
pasos aceptados/rechazados y ETA. Reserva 4 GiB en scratch; no llena home.

Salidas: identity.json, executed_plan.json, progress.jsonl, checkpoints aceptados,
observaciones, refinement.json y summary.json. Si falla conserva la evidencia
y se detiene; no reintenta ni cambia tolerancias automáticamente. El contraste
final usa señales iniciales además del remanente, evitando exigir precisión
relativa de una cola extinguida. No acredita todavía detector ni poblaciones.

Diagnósticos ligeros útiles después de la corrida:

```bash
tail -n 8 /home/jdiaz/scratch/stage4_nonlinear_time_20260924/progress.jsonl
/home/jdiaz/.conda/envs/snspd/bin/python -m json.tool /home/jdiaz/scratch/stage4_nonlinear_time_20260924/summary.json
```

## Historial de instrucciones anteriores — completadas o sustituidas

Se preserva a continuación el cuaderno anterior íntegro. Sus títulos de
«siguiente ejecución» son históricos; el único comando pendiente es ETD2 arriba.

# Geminga: siguiente ejecución vigente

Actualizado el 24 de septiembre de 2026 UTC. Cuenta jdiaz, sin administrador.

## Estado de las campañas anteriores

La campaña de momentos solicitada terminó correctamente en 198,12 s: 181
espectros, 362 respuestas cinéticas y siete proyecciones. No debe repetirse.
Su cuadratura fina reduce el defecto del peso térmico a 0,424 %, pero deja
una discrepancia real del potencial único: 20,55 % en corriente y 77,53 % en
torque de fase. Son respuestas de prueba, no errores de latencia.

El diagnóstico armónico también terminó: 24 trabajos en 21,19 s, sin nuevas
consultas espectrales. La eliminación algebraica que conserva cada energía
reproduce bien las anclas lentas: 0,0258 % de diferencia de corriente en
ν=0,001. Las frecuencias altas son exploratorias, no una validación de la
respuesta de pocos picosegundos. Véanse la
[revisión actual](implementation/stage4/moment_review_20260924/README.md)
y el [informe](https://github.com/JoaquinDiazM/pysnspd/blob/c80c0f8612d8fd5a2b390938ca5ed4b9d58d3caa/output/pdf/implementation/Informe_etapa_4_respuesta_de_carga.pdf).

## Siguiente ejecución: trayectoria térmica débil con condensado móvil

El operador completo ya se comprobó en 32,94 s: 256 frecuencias, 2304 raíces
espectrales y cinco tests nuevos. Se conservó la corrección de movilidad y
potencial debida al residuo de la referencia, que aporta aproximadamente 1,1 %
a la respuesta. **No repetir la preparación del operador.** Sus matrices están
en `/home/jdiaz/scratch/stage4_thermal_weak_20260924` y se reutilizan directamente.

Ejecutar una vez en la terminal habitual:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/thermal_time_campaign.py \
  --plan docs/implementation/stage4/moment_review_20260924/thermal_time/plan.json \
  --operator-root /home/jdiaz/scratch/stage4_thermal_weak_20260924 \
  --output-root /home/jdiaz/scratch/stage4_thermal_time_20260924 \
  --execute
```

**Propósito:** seguir hasta 1 ps una referencia térmica y dos perturbaciones
débiles de amplitud/fase, conservando todos los nodos y la derivada del mismo
sistema KWT con potencial normal. La amplitud inicial de prueba es 0,001;
se detiene si alguna trayectoria abandona el vecindario débil declarado del
2 % del gap de referencia. Ese límite no modifica ni recorta campos. No hay
fotón, circuito sustituto ni cambios de constantes para acelerar el sistema.

Las frecuencias se distribuyen en un único grupo residente de trabajadores:
cada matriz se factoriza una vez y se reutiliza en las aplicaciones del operador
temporal. Las trayectorias comparten ese presupuesto. Se programa desde el
inicio una comparación con mayor dimensión de Krylov e intervalos más cortos;
se observa su efecto sobre desplazamiento, fuerza, corriente y torque.

**Tiempo previsto:** reservar **5–45 min**, estimación conservadora todavía no
medida para la trayectoria completa. Reutilizar matrices evita raíces no lineales
por paso, pero cada aplicación temporal resuelve las 256 respuestas tangentes
para tres campos. La rigidez y los reinicios numéricos determinan el coste.
Por poder superar cinco minutos, la campaña física completa queda para ejecución
manual. No se ha iniciado ni se mantiene una espera automática del agente.

**Recursos:** hasta 27 trabajadores y un coordinador, 28 de 32 hilos lógicos
(87,5 %), dejando dos núcleos físicos completos libres. BLAS/OMP usa un hilo
por proceso. Se comprueban afinidad, memoria y almacenamiento disponibles;
la planificación reserva hasta 27 GiB de trabajadores más 3 GiB del coordinador,
y reduce el número de trabajadores si supera el 90 % de la memoria disponible.
Son reservas de planificación, no límites de memoria impuestos por el sistema.

**Progreso y salidas:** imprime preparación de factores, pasos intentados y
aceptados, tiempo físico alcanzado y ETA. Se reservan 4 GiB en scratch.
Conserva `identity.json`, plan ejecutado, `progress.jsonl`, cada estado aceptado,
observaciones, resúmenes de los pases y `refinement.json`. El resultado final es
`summary.json`; un fallo deja `failure.json` y todos los checkpoints previos.
La ETA temprana puede fluctuar con la rigidez.

No sobrescribe un destino existente ni relanza campañas fallidas. El ajuste
interno del paso temporal es parte declarada del integrador y queda registrado;
no cambia ecuaciones ni reduce tolerancias en silencio. El defecto de Arnoldi
no se presenta como cota rigurosa de error: se contrasta con el segundo pase.

**Interpretación:** se registra la energía libre cuadrática y el defecto de
truncación de la ley térmica linealizada. No se afirma conservación no lineal
exacta ni validación del dispositivo. El trabajo espectral no térmico y el
acoplamiento al dispositivo/circuito siguen siendo requisitos de etapa 4.
Al terminar, basta avisar en el chat; no hace falta copiar el registro completo.

Para comprobar únicamente fuentes, entradas y recursos, omitir `--execute`.
Esa variante no inicia trabajadores ni escribe salidas. Requiere un destino
todavía libre para confirmar que la ejecución posterior no sobrescribirá datos.


## Comandos ligeros útiles

Verificar esta entrega, sus fuentes y la cadena histórica sin resolver problemas físicos:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage4_core/verify_moment_review.py
```

Los comandos de reproducción del diagnóstico armónico están en
[su registro](implementation/stage4/moment_review_20260924/charge_frequency/README.md).
Ya se ejecutaron; no forman parte de la cola pendiente. Toda repetición usa
un destino nuevo y conserva los resultados anteriores.

## Horizonte posterior del dispositivo

Se mantiene el circuito completo de la memoria. El futuro fotón se observará
hasta el gatillo confirmado en Vout más margen, con techo si no aparece.
Sólo cambia la duración observada: ecuaciones, material, geometría, bordes,
precisión y constantes físicas son las del mismo dispositivo.

## Historial conservado: momentos ya ejecutados

El contenido siguiente se conserva como registro de la entrega anterior.
Sus instrucciones no son la cola vigente: las campañas térmica, espectral
y de momentos allí solicitadas ya terminaron. No configurar nuevas sesiones.


# Geminga: siguiente ejecución vigente

Actualizado el 24 de septiembre de 2026 UTC. Cuenta jdiaz, sin administrador.

## Resolver los momentos y comparar un potencial por nodo

La continuación térmica terminó bien: cuatro núcleos aceptados en 8,08 min.
También terminaron 72 consultas espectrales (88,48 s), 144 respuestas cinéticas
(5,13 s) y seis comparaciones con un potencial electroquímico (2,20 s).
No hace falta repetir la relajación térmica ni los comandos históricos de abajo.

La comparación preliminar deja diferencias de aproximadamente 20 % en corriente,
3,1 % en fuerza y 0,42 % en flujo de energía, pero la malla energética actual
sobrestima una integral térmica conocida en 18 %. El siguiente cálculo resuelve
esa limitación antes de decidir el cierre físico. No se descarta el potencial
de la memoria ni se introduce una variable dinámica nueva a partir de datos
insuficientemente resueltos. Véanse el
[informe vigente](implementation/stage4/self_consistent_review_20260924/README.md)
y la [campaña](implementation/stage4/self_consistent_review_20260924/resolution_campaign/README.md).

Ejecutar una vez en la terminal habitual, sin configurar nuevas sesiones:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/resolution_campaign.py \
  --plan docs/implementation/stage4/self_consistent_review_20260924/resolution_campaign/plan.json \
  --output-root /home/jdiaz/scratch/pysnspd_stage4_moment_resolution_20260924 \
  --execute
```

**Propósito:** comparar corriente, fuerza y flujo de energía integrados con
31 y 50 energías anidadas, dos desplazamientos del contorno y dos mallas
espaciales. Son 181 consultas únicas; la respuesta cinética y las siete
proyecciones reutilizan esos espectros. El núcleo térmico y sus bordes se
conservan. No se simula un fotón ni se ejecuta recuperación larga.

**Tiempo previsto:** reservar 5–15 min. La campaña anterior de 72 consultas
duró 88,48 s; ahora hay más energías en regiones difíciles y un contorno más
cercano al eje real. La extrapolación es orientativa. Al poder superar cinco
minutos, esta campaña completa queda para ejecución manual según la política
acordada. El agente no la inicia ni espera sondeando el proceso.

**Recursos:** hasta 27 trabajadores y un coordinador, 28 de 32 CPU lógicas
(87,5 %), reservando dos núcleos físicos completos en la topología actual.
Un hilo BLAS/OMP por proceso. El presupuesto se vuelve a comprobar por fase
y no supera el 90 % de la memoria disponible; reserva 1 GiB por trabajador
espectral/cinético y 2 GiB por trabajador de proyección, más el coordinador.
Son reservas de planificación, no límites impuestos por el sistema operativo.
Las proyecciones son siete trabajos paralelos, no 27 copias innecesarias.

**Almacenamiento y salidas:** se utiliza scratch propio, que tenía 982 GiB
libres. La reserva estimada es 1,16 GiB; el preflight la comprueba antes de crear el
destino. Se guardan `spectra/`, `kinetic/`, `projection/`, identidad, plan,
barra/ETA por fase, `progress.jsonl`, `summary.json` y `moment_comparisons.json`.
El programa conserva cada consulta terminada, no sobrescribe una ruta existente
y se detiene si hay un fallo. No hace reintentos ni cambia el método en silencio.
Una ETA temprana puede fluctuar porque el coste depende de la energía y la malla.

**Lectura de resultados:** se contrastará el cambio físico con la variación
numérica observada; no se exige una precisión universal en cada punto de DOS.
Un diagnóstico estático favorable prepara el siguiente ensayo dinámico débil,
pero no cierra automáticamente la etapa 4 ni acredita una latencia del detector.
Cuando termine, basta avisar en el chat: no hace falta copiar todo el registro.

## Comandos ligeros útiles

Verificar esta entrega y su cadena histórica, sin resolver problemas físicos:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage4_core/verify_self_consistent_delivery.py
```

Para revisar sólo fuentes, rutas y recursos del nuevo cálculo, usar el comando
principal sin `--execute`, con el destino todavía libre. No crea salidas ni
inicia la campaña. Si el destino ya existe se rechaza para evitar sobrescritura.

## Horizonte futuro con fotón

Se conserva el dispositivo completo y el circuito de la memoria. Sólo se
limitará el tiempo observado hasta un gatillo confirmado en Vout más margen,
con un techo finito si no aparece. No se modifican ecuaciones, geometría,
material, precisión ni constantes físicas para acelerar la recuperación.

## Historial conservado: continuación térmica ya ejecutada

Lo siguiente se conserva íntegro como registro de la entrega anterior.
Sus instrucciones de ejecución no son la cola vigente: la continuación
térmica y los pilotos mencionados ya terminaron.

# Geminga: comandos vigentes

Actualizado el 24 de septiembre de 2026 UTC. Cuenta jdiaz, sin administrador.

## Siguiente ejecución: núcleo térmico autoconsistente

La referencia espacial ya recupera la fuerza radial de Usadel con diferencia
de malla de 0,18 % al mismo corte espectral. El piloto siguiente completó dos
barridos en 79,60 s, con descenso de energía en los cuatro casos, pero su
residuo de autoconsistencia todavía es 1,85–2,44 %. Corresponde continuar esa
minimización de la misma energía, sin introducir otro cierre ni tiempo físico.
La etapa 4 permanece abierta. Véase el
[informe vigente](implementation/stage4/spatial_energy_20260924/README.md).

Ejecutar **una vez** desde la terminal; no requiere modificar screen ni activar
otro entorno. El destino debe estar libre: el programa no sobrescribe resultados.

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/run_self_consistent_core.py \
  --plan docs/implementation/stage4/spatial_energy_20260924/self_consistent_plan.json \
  --resume-from tmp/stage4_self_consistent_pilot_20260924 \
  --output-root /home/jdiaz/scratch/pysnspd_stage4_self_consistent_20260924 \
  --execute
```

**Propósito:** obtener núcleo, corriente y energía autoconsistentes para dos
mallas y dos cortes espectrales, con un caso inicialmente asimétrico. El criterio
inicial es un residuo RMS ponderado de 0,1 % en el núcleo; también se registra
el máximo global. No es una afirmación de estacionariedad en todo el dominio
ni un certificado del detector. Llegar al máximo de 150 barridos sin cumplirlo
produce un resultado incompleto, nunca una aprobación automática.

**Duración estimada:** reservar aproximadamente 15–100 min, extrapolando el
piloto; el coste por barrido y la cantidad necesaria pueden cambiar. La barra y
ETA estiman el trabajo presupuestado restante, no garantizan la convergencia.
El agente no lanza esta campaña ni espera sondeando su ejecución.

**Recursos:** un solo grupo compartido por frecuencias y casos, hasta 27
trabajadores y un coordinador (28 de 32 CPU lógicas, 87,5 %). Se dejan dos de
los 16 núcleos físicos completos libres; cada proceso usa un hilo BLAS/OMP.
Reserva conservadora de RAM: 30 GiB, revisada contra el 90 % de la memoria
disponible. No es un límite de memoria impuesto por el sistema operativo.
Los casos que cumplen el criterio salen del grupo y liberan capacidad.

**Almacenamiento:** se usa el scratch propio de jdiaz, con unos 982 GiB libres
en la comprobación. La proyección conservadora sin compresión es 8,51 GiB
(reservar 10 GiB). `/home` sólo tenía unos 4,3 GiB libres; no se borran ni mueven
resultados previos. El preflight verifica espacio, plan, fuentes y hashes antes
de iniciar. No se debe volver a utilizar esta ruta si la corrida ya existe.

**Salidas:** `identity.json`, `progress.jsonl`, `summary.json` y `checkpoints/`
dentro del destino indicado. Hay métricas JSON en cada barrido y campos NPZ
en el primero, cada cinco, al aceptar un caso o agotar su presupuesto. Cada
checkpoint distingue el condensado y sus espectros del siguiente condensado
propuesto. Los `resume_index_*.json` permiten una continuación explícita en una
carpeta nueva; como máximo se repetirían cuatro barridos sin checkpoint completo.
El resumen puede ser incompleto aunque el proceso finalice correctamente: se
revisarán su estado y sus residuos antes de interpretar los perfiles.
Una interrupción durante la escritura de un índice puede requerir recuperar
manualmente el último índice completo; no hay sustitución silenciosa.

Cuando termine, basta avisar en el chat; no hace falta copiar todo el registro.

## Comprobaciones ligeras opcionales

Verificar entrega y cadena histórica (segundos, sin cálculos físicos):

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage4_core/verify_spatial_energy.py
```

Para repetir sólo el preflight del comando principal, omitir `--execute`.
Esa modalidad consulta recursos, comprueba fuentes y lee checkpoints; no crea
el destino ni inicia la minimización. Puede usarse también después de la corrida
como revisión de la preparación, pero no acredita sus resultados nuevos.

Las 36 pruebas analíticas enfocadas pasaron en 0,15 s. La referencia prescrita
de 2.048 modos terminó en 60,35 s y no necesita repetirse antes de esta campaña.
Su comando reproducible opcional (menos de cinco minutos según la medición)
usa una carpeta nueva:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/run_spatial_energy.py \
  --plan docs/implementation/stage4/spatial_energy_20260924/campaign_plan.json \
  --output-root /home/jdiaz/scratch/pysnspd_spatial_reference_reproduction --execute
```

## Horizonte futuro con fotón

Se simulará hasta un cruce confirmado de Vout más un margen, con techo finito
si no se observa disparo. Sólo cambia la ventana: mismas ecuaciones, material,
depósito, geometría, malla, bordes, precisión y circuito completo de la memoria.
No se adelanta la recuperación ni se acortan constantes de tiempo. El estado
final permitirá continuar. El criterio se definirá para el ensayo; la falta
de cruce no recibirá una latencia ficticia. Esta política todavía no reemplaza
el modo heredado que espera un máximo del pulso.

## Historial íntegro

La [libreta anterior](implementation/stage4/spatial_energy_20260924/previous_delivery_exact/docs/GEMINGA_COMMANDS.md)
conserva íntegros los comandos y enlaces de entregas previas. No hace falta
repetir las campañas terminadas.

## Política de cómputo

Cálculos razonablemente previstos de más de cinco minutos se entregan al usuario
con comando copiable en el chat, propósito, salidas y recursos. No se lanzan
desde el agente ni mediante trabajos de fondo o sondeos. Los chequeos inciertos
llevan timeout de 240 s; si vence, quedan incompletos y se entrega el comando
sin reinicios repetidos ni subdivisiones para eludir el límite. No se requieren
nuevas sesiones ni comandos screen.


## Registro de comandos preparados: RHS y balance, 5 de octubre de 2026

Este cuaderno contiene la campaña pendiente de RHS y balance. El historial
completo permanece en `/home/jdiaz/GEMINGA_COMMANDS.md`. La preparación DC del
24–25 de septiembre ya está admitida: no se vuelve a ejecutar.

## Evidencia ya obtenida

- Cuatro referencias DC y cinco holds de10ps revalidados, sin repetirlos.
- Observador U_N real:1712nodos,256modos,0.9K,referencia21.5µA,13.40s.
- Derivada U independiente: errorrelativo8.6425e-7, defecto1.5490e-18W,
  diferencias centrales±1e-4 sobre amplitud;22.38s con2casos×4workers.
- Cuatro RHS reales de transporte: conservacióninstantánea~1e-17rel.
- Un RHS matricial a1.5meV: ambosmodos evolucionan; residuoproyectado4.24e-16.
  El residuo Nambu completo sin fuentes anómalas queda0.274por masa: no se
  admite por ello un cierre completo.

Las sondas de amplitud/población son controles estáticos del acoplamiento,
no una preparación fotónica ni una selección de ancho de absorción.
Las unidades de entrada internas son SI; plots usanmeV,ps,mV,µA,J,Wsegúnmagnitud.

## Reglas de ejecución

Ejecutar las campañas una por vez. Cada ejecutor paraleliza sus casos/consultas
internas, reevalúa afinidad/cuota/RAM y reserva almenos10% de los recursos.
No lanzar estos cuatro comandos simultáneamente: no comparten un presupuesto.
Cada ruta de salida debe ser nueva. No hay reintentos,resume,clipping ni
reparación de energía mediante una fuente residual. Si falla, conservar
failure.json/logs/últimoestadoaceptado y volver al chat con los resultados.
No cambiar tolerancias,eta,fracción de excitación ni ecuaciones para conseguir
un pase. Las estimaciones son de planificación; no un speedup supuesto.

## 1. Balance direccional de ambos modos con espectro espacial móvil

Necesita el observador ya calculado en
`/home/jdiaz/scratch/stage5_rhs_energy_reference_20261005` y la referenciaDC215.
40energías entre0y6meV×2reguladores=80consultas. Se contrae la RHS temporal
con el momento de energía, la tangenteDOS y U_N por separado; se publica el
trabajo candidato y su defecto sin cerrar la cuenta por definición.

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_rhs/run_matrix_rhs_audit.py --plan docs/implementation/stage5/rhs_balance_20261005/matrix_rhs_snapshot_audit.json --output /home/jdiaz/scratch/stage5_matrix_rhs_audit_20261005 --workers 4 --execute
```

Salida:identity.json,queries/*.npz y*.json,summary.json. Revisar
independent_candidate_work_defect_W, flujodeborde, residuoNambu, relaciónR1/R0,
comparación de gauge y dependenciaeta. 4workers+coordinador,1hilo/worker,
reserva7GiB como techo de planificación. Única consulta medida:1.30s en1CPU.
Los estados subgap/bordegap y3000iteracionesLSMR no están medidos: el conjunto
puede superar5min; ejecución del usuario. Un fallo físico no autoriza un fotón.

## 2. Corte espectral, malla, longitud y depresión finita del gap

Ocho casos estáticos reutilizan las referencias DC y comparan256/512/1024modos,
las mallas5/3.5nm,los cortes160/320nm y depresión deamplitud10%.
No relajan otraDC ni evolucionanunhotspot. Doscasosenparalelo×4workers,
11procesosincluidoscoordinadores, reservaaprox17GiB. Baseline256medido13.40s;
los mayores estados no lineales no están medidos. Plan conservador total
aprox10–20min, incierto; usuarioejecuta.

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_rhs/run_energy_campaign.py --suite convergence --output-root /home/jdiaz/scratch/stage5_rhs_energy_convergence_20261005 --execute
```

Salida:campaign.json,logs/,cases/*/{summary.json,fields.npz,spectral_receipts/},
analysis.json,provenance.json,integrity.json. Revisar cambios enU_J y su
potencia direccional W; ambos cortes modifican un observador finito del mismo
estado guardado. No interpretar el fondo normal ni la colaUV como admitidos.
No se certifica convergencia térmica global por un único Nc256.

## 3. Evolución conservativa y refinamientos del transporte neutral

Cuatro casos de2ps:dt=.001/.0005ps, malla5/3.5nm,longitud160/320nm.
Usa el catálogo heredado,volúmenes y caras de la malla dual real, con espectro
bulk congelado. Los cortes son de flujo cero para esta identidad aislada;
no sustituyen las continuaciones superconductivas del detector.
La mezclaFD de control aportaaprox.003047eV: no corresponde a un fotón óptico.
No activa colisiones,escape,BGK,fuente,condensado,potencial ni circuito.

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_rhs/run_nodal_transport_audit.py --workers 3 --output /home/jdiaz/scratch/stage5_nodal_transport_20261005 --execute
```

Salida:campaign_result.json,results.json,cases/*/{result.json,trajectory.npz,
progress.json},y last_accepted_state.npz/failure.json si el caso se detiene.
3casosenparalelo×1hilo+coordinador, reserva6GiB. ProxyRHSsinIO:75–186s/caso;
planconservador600–855s/caso por evolución,observación,auditoría ycontención.
Revisar energía y número globales, p∈[0,1], momentoinstantáneo y diferencias
entre dt,dt/2,malla ylongitud. No asignar latencia ni reset a este control.

## 4. Refinamiento de energía real y regulador del balance móvil

Ejecutar tras leer el resultado1; si falla su contrato de inputs/código,
corregir primero y no usar un lote refinado como reparación implícita.
96energías en0–12meV×2eta=192consultas: separa error de cuadratura/cola,
regulador y falta de cierre físico. Mismo gap,dirección y referencia.

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_rhs/run_matrix_rhs_audit.py --plan docs/implementation/stage5/rhs_balance_20261005/matrix_rhs_snapshot_refined.json --output /home/jdiaz/scratch/stage5_matrix_rhs_refined_20261005 --workers 4 --execute
```

Salidas/recurso como1. La consulta aislada no estima los estados cercanos al
gap ni el lotecompleto; su tiempo puede exceder5min. Registrar ambosdefectos
originalyrefinado. Si el defecto persiste, revisar fuerzaa poblacionesfijas,
self-energyR/A/K,neutralidad y doblecontabilizaciónKWT: no inventar una tasaNbN.

## Puertas posteriores

Estos resultados verifican observadores y sectores ejecutables. Antes de
admitir el nuevo detector falta cerrar fuentes R/A/K a energía finita,
neutralidad/fuerza en gauge común y trabajo de cada canal sin doble conteo;
además,documentar normalización materialNbN y transferencia fotónica.
No hay todavía comando que declare resuelto lo que está abierto.

Después viene un fotón individual y transiente hasta cruceVout confirmado
más margen. Techo sin gatillo:200ps. Umbral se estudia enmV, rangoexploratorio
0.1–0.5mV del estudio experimental previo; confirmación,histéresis,margen y
transferencia/ancho requieren contrato explícito, no valores inventados.
Se conserva el sentidoDC global y el circuitoCM completo con Lext congelada.
No se prepara aquí ningún transiente de reset eléctrico total.
Latencia absoluta y jitter sólo se analizan tras admitir ese detector y su
transferencia; jitter necesita además variabilidad física documentada.

## Diagnósticos ligeros útiles

```bash
cd /home/jdiaz/pysnspd
timeout 240s /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage5_rhs/check_rhs.py --receipt /home/jdiaz/scratch/stage5_rhs_light_validation_repeat_20261005.json
```

El nombre de recibo debe ser nuevo. No vuelve a correrDC ni un transiente.
Los ejecutores1,2y4 sin--execute hacen sólo dry-run; transporte con--plan-only
escribe un plan en una ruta nueva, sin evolución. No reusar las salidas finales
para inspeccionar el plan. Traslas corridas,volver al chat con `reinicia` y la
salida terminal; los artefactos completos quedanen geminga para análisis.


### Comandos usados para evidencia ligera (no repetir por defecto)

```bash
cd /home/jdiaz/pysnspd
timeout 240s /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage5_rhs/audit_internal_energy.py --reference /home/jdiaz/scratch/stage5_prephoton_dc_20260924/references/cases/dc215_l160_h5/reference.npz --count 256 --output /home/jdiaz/scratch/stage5_rhs_energy_reference_20261005 --workers 4 --probe-span-nm 80 --direction-rate-per-s 1e10 --execute
timeout 180s /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_rhs/run_energy_campaign.py --suite direction --output-root /home/jdiaz/scratch/stage5_rhs_direction_check_20261005 --execute
```

El segundo terminó ambas físicas con éxito pero exit1 por cambio de un driver ajeno incluido en el conjunto amplio de hashes. `integrity.json` original se preserva; `integrity_assessment.json` identifica que no cambió un productor numérico usado. La derivada U pasó; no se repitió cálculo físico. El ejecutor ahora fija sólo dependencias realmente cargadas.


## Actualización final del cuaderno vigente: cinco auditorías independientes

Se separaron cuadratura a 6 meV y corte UV a 12 meV. Las entradas previas se preservan como historia y esta versión del cuaderno vigente prevalece.

# Geminga: comandos necesarios ahora — 5 de octubre de 2026

La preparación DC del 24–25 de septiembre está admitida: **no repetirla**.
Este cuaderno contiene sólo las auditorías necesarias para continuar el RHS
espacial y su balance. `/home/jdiaz/GEMINGA_COMMANDS.md` conserva el historial
completo de comandos posibles, preparados y usados.

## Evidencia disponible

- Cuatro referencias DC y cinco holds de 10 ps revalidados sin repetirlos.
- Observador $U_N$ en 1712 nodos y 256 modos, a 0,9 K y referencia de
  21,5 µA: 13,40 s. La energía es relativa al metal normal del mismo corte.
- Diferencia finita independiente: defecto de 1,5490e-18 W y error relativo
  8,6425e-7; dos sondas de amplitud de ±1e-4, en 22,38 s.
- Cuatro RHS completos de transporte: conservación instantánea de energía
  y número con residuos relativos del orden de 1e-17.
- RHS matricial a 1,5 meV: evolucionan ambos modos; residuo proyectado
  4,24e-16. El residuo Nambu completo sin fuentes anómalas queda en 0,274
  por masa nodal. Eso no acredita el cierre completo del detector.
- Validación integrada: **140 pruebas y 13 subpruebas**, en 0,89 s.

Las sondas de amplitud y población son controles del acoplamiento, no una
preparación fotónica ni una selección del ancho de absorción. Las entradas
físicas son SI; los plots usan meV, ps, mV, µA, J y W según la magnitud.

## Ejecución y conservación de datos

Ejecutar las campañas **una por vez**. Cada ejecutor paraleliza sus consultas
o casos internos y reevalúa afinidad, cuota y RAM; reserva al menos el 10 %.
Estos ejecutores no comparten un presupuesto entre campañas simultáneas.
Cada ruta de salida debe ser nueva. Se conservan fallos, logs y estados
aceptados; no hay reintentos, resume, clipping ni calor residual para cerrar
una cuenta. Si falla, volver al chat con los artefactos y la salida terminal.
Las estimaciones son de planificación, sin suponer un speedup no medido.

## 1. Balance direccional de ambos modos y espectro espacial móvil

Reutiliza el observador de energía ya calculado en
`/home/jdiaz/scratch/stage5_rhs_energy_reference_20261005` y la referencia
DC215. Son 40 energías entre 0 y 6 meV por dos reguladores: 80 consultas.
Contrae por separado el RHS temporal, la tangente DOS y $U_N$; publica el
trabajo candidato y su defecto sin imponer que el balance sea cero.

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_rhs/run_matrix_rhs_audit.py --plan docs/implementation/stage5/rhs_balance_20261005/matrix_rhs_snapshot_audit.json --output /home/jdiaz/scratch/stage5_matrix_rhs_audit_20261005 --workers 4 --execute
```

Salidas: `identity.json`, `queries/*.npz`, `queries/*.json`, `summary.json`.
Revisar `independent_candidate_work_defect_W`, flujo de borde, residuo Nambu,
relación $R_1/R_0$, comparación de gauge y dependencia del regulador eta.
Usa 4 workers más coordinador, un hilo por worker y reserva de 7 GiB.
Una consulta a 1,5 meV tomó 1,30 s en una CPU. El subgap, el borde del gap y
LSMR de hasta 3000 iteraciones no están medidos: el conjunto puede exceder
cinco minutos y corresponde al usuario. Un fallo físico no autoriza un fotón.

## 2. Corte espectral, malla, longitud y depresión finita del gap

Ocho casos estáticos comparan 256/512/1024 modos, mallas de 5/3,5 nm,
cortes de 160/320 nm y depresión de amplitud del 10 % sobre las referencias
existentes. No relajan otra DC ni evolucionan un hotspot. Hay dos casos en
paralelo con cuatro workers cada uno: 11 procesos contando coordinadores,
reserva de unos 17 GiB. El baseline de 256 modos tomó 13,40 s; los mayores
estados no lineales no están medidos. Plan conservador aproximado de
10–20 min para el conjunto, sujeto a revisión con tiempos medidos.

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_rhs/run_energy_campaign.py --suite convergence --output-root /home/jdiaz/scratch/stage5_rhs_energy_convergence_20261005 --execute
```

Salidas: `campaign.json`, `logs/`, `cases/*/summary.json`, `fields.npz`,
`spectral_receipts/`, `analysis.json`, `provenance.json`, `integrity.json`.
Revisar cambios de energía interna en J y de su potencia direccional en W.
Los distintos cortes evalúan el mismo estado guardado; no certifican por sí
solos el fondo normal, la cola UV o el balance térmico dinámico completo.

## 3. Evolución conservativa y refinamientos del transporte neutral

Cuatro casos de 2 ps: dt de 0,001/0,0005 ps, malla de 5/3,5 nm y longitud de
160/320 nm. Usa el catálogo heredado, volúmenes y caras de la malla dual
real con espectro bulk congelado. Los cortes tienen flujo cero para esta
identidad aislada; no sustituyen las continuaciones superconductivas del
SNSPD. La mezcla FD del control aporta unos 0,003047 eV: no es un fotón
óptico. Colisiones, escape, BGK, fuente, condensado, potencial y circuito
están desactivados y no se consideran validados por este control.

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_rhs/run_nodal_transport_audit.py --workers 3 --output /home/jdiaz/scratch/stage5_nodal_transport_20261005 --execute
```

Salidas: `campaign_result.json`, `results.json`, `cases/*/result.json`,
`trajectory.npz`, `progress.json`; también `last_accepted_state.npz` y
`failure.json` si hay una parada. Tres casos en paralelo con un hilo cada
uno, más coordinador; reserva de 6 GiB. Proxy sin IO: 75–186 s por caso;
plan conservador de 600–855 s por caso para evolución, auditoría,
observación y contención. Revisar energía y número globales, soporte
$p\in[0,1]$, momentos instantáneos y refinamiento temporal y espacial.
No se asigna latencia ni reset a este control de transporte.

## 4. Refinamiento de cuadratura a igual corte físico

Leer primero el resultado 1. Si falla el contrato de inputs o código,
corregir antes de ampliar el lote. Se sube de orden Gauss 4 a 8 manteniendo
el corte de 6 meV, mismos paneles, reguladores, gap, dirección y referencia:
80 energías por dos eta, 160 consultas. Esta comparación aísla cuadratura.

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_rhs/run_matrix_rhs_audit.py --plan docs/implementation/stage5/rhs_balance_20261005/matrix_rhs_snapshot_refined.json --output /home/jdiaz/scratch/stage5_matrix_rhs_refined_20261005 --workers 4 --execute
```

Salidas y recursos como 1. Puede superar cinco minutos. Registrar ambos
defectos y todos los residuos; el lote no repara un cierre físico faltante.

## 5. Corte UV con cuadratura ya refinada

Comparar con 4: mismo orden 8 y reguladores, corte ampliado de 6 a 12 meV.
Son 96 energías por dos eta, 192 consultas. Así se distingue la sensibilidad
al corte de la sensibilidad a cuadratura y regulador.

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_rhs/run_matrix_rhs_audit.py --plan docs/implementation/stage5/rhs_balance_20261005/matrix_rhs_snapshot_uv.json --output /home/jdiaz/scratch/stage5_matrix_rhs_uv_20261005 --workers 4 --execute
```

Salidas y recursos como 1; usuario ejecuta. Si el defecto persiste al refinar,
revisar fuerza a poblaciones fijas, self-energy R/A/K, neutralidad y posible
doble conteo KWT. No inventar tasas NbN ni inyectar calor igual al defecto.

## Puertas posteriores

Estos resultados verifican observadores y sectores ejecutables. Antes de
admitir el nuevo detector falta cerrar fuentes R/A/K a energía finita,
neutralidad y fuerza en el mismo gauge, y trabajo de cada canal sin doble
conteo. También falta documentar la normalización material NbN y la
transferencia fotónica. No hay un comando que declare resuelto ese cierre.

Después viene un fotón y transiente hasta cruce Vout confirmado más margen.
Techo sin gatillo: **200 ps**. Umbral en mV; rango exploratorio previo de
0,1–0,5 mV, sin atribuirlo a una calibración experimental del dispositivo.
Confirmación, histéresis, margen, transferencia y ancho necesitan contrato
explícito. Se conserva el sentido DC global y el circuito CM completo con
Lext congelada. No hay comandos de reset eléctrico total en este cuaderno.
Latencia absoluta y jitter se analizan tras admitir detector y transferencia;
jitter requiere además una variabilidad física documentada.

## Diagnóstico ligero reproducible

```bash
cd /home/jdiaz/pysnspd
timeout 240s /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage5_rhs/check_rhs.py --receipt /home/jdiaz/scratch/stage5_rhs_light_validation_repeat_20261005.json
```

El nombre del recibo debe ser nuevo. No repite DC ni un transiente. Los
comandos 1, 2, 4 y 5 sin `--execute` son dry-run. Transporte con `--plan-only`
escribe sólo un plan en una ruta nueva; no usar la salida final para ello.
Tras las corridas, volver al chat con `reinicia` o equivalente y la salida
terminal. Los datos completos quedan en geminga para el análisis.


### Validación final con guardas del observador

145 pruebas y 13 subpruebas pasaron en 1,03 s. Se añadieron guardas que rechazan datos de gap, dirección, gradiente o conteo incompatibles con el observador. Recibo: `/home/jdiaz/pysnspd/docs/implementation/stage5/rhs_balance_20261005/validation_final.json`. Las entradas anteriores de 140 pruebas se conservan como historia.


## Actualización 8 de octubre de 2026: cinco campañas anteriores USADAS

Las cinco campañas del5deoctubre ya fueron ejecutadas y revisadas. No repetir DC, refinamiento UV ni los dos transportes h5 que pasaron. El cuaderno vigente está en `/home/jdiaz/GEMINGA_COMMANDS_PRESENT.md`. Sólo quedan las dos campañas de transporte con avance reparado y los tres planes de borde siguientes. Ejecutar uno por vez; el usuario realiza el cómputo largo.

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_rhs/run_nodal_transport_audit.py --workers 1 --case l160_h35_dt1 --output /home/jdiaz/scratch/stage5_nodal_transport_h35_pauli_20261008 --execute
```

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_rhs/run_nodal_transport_audit.py --workers 1 --case l320_h5_dt1 --output /home/jdiaz/scratch/stage5_nodal_transport_l320_pauli_20261008 --execute
```

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_rhs/run_matrix_rhs_audit.py --plan docs/implementation/stage5/rhs_balance_20261008/matrix_rhs_snapshot_edge.json --output /home/jdiaz/scratch/stage5_matrix_rhs_edge_20261008 --workers 4 --execute
```

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_rhs/run_matrix_rhs_audit.py --plan docs/implementation/stage5/rhs_balance_20261008/matrix_rhs_snapshot_edge_refined.json --output /home/jdiaz/scratch/stage5_matrix_rhs_edge_refined_20261008 --workers 4 --execute
```

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_rhs/run_matrix_rhs_audit.py --plan docs/implementation/stage5/rhs_balance_20261008/matrix_rhs_snapshot_edge_panels_halved.json --output /home/jdiaz/scratch/stage5_matrix_rhs_edge_panels_halved_20261008 --workers 4 --execute
```

Las salidas son nuevas. El cierre fotónico sigue pendiente de balance acoplado; ventana sólo gatillo más margen y200ps sincruce, sinreset eléctrico total.


## Preparación física corta ejecutada el 8 de octubre de 2026 - USADA

La ampliación siguiente tomó 2,773 s y calculó sólo 32 raíces R0 nuevas de 12 a 24 meV. No repitió DC, energías anteriores, R1 ni RHS; no avanzó tiempo físico. Está USADA y no pertenece al cuaderno de campañas pendientes.

```bash
cd /home/jdiaz/pysnspd
timeout 240s /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage5_photon/extend_static_photon_spectrum.py --cache /home/jdiaz/scratch/stage5_matrix_rhs_uv_20261005 --output /home/jdiaz/scratch/stage5_photon_static_uv24_20261008 --workers 2 --execute
```

Las seis preparaciones físicas actuales y sus hashes están seleccionadas por `docs/implementation/stage5/photon_readiness_20261008/CURRENT_PREPARATIONS.json`: dos longitudes de onda por tres varianzas. Hay mapas y fuentes de energía, seis estados FD/Bose a gap DC congelado y dos contratos nominales. Las dos hipótesis nuevas fueron autorizadas por el usuario: material efectivo exploratorio y ancho gaussiano de referencia con sensibilidad. No se supone resuelta la cascada óptica ni calibradas las tasas de la película. El source ya contenido en las poblaciones no debe depositarse por segunda vez.

Quedan pendientes el backend común y su balance no térmico. Los cinco comandos largos necesarios están en `GEMINGA_COMMANDS_PRESENT.md`; no existe un comando de transiente fotónico admitido que se pueda sustituir por producción 0.4. Los resultados futuros de formación de hotbelt deben ser estados temporales aceptados y llevar reloj físico, campos completos, energía y circuito. Ventana de 200 ps sin gatillo o gatillo confirmado más margen; no hay reset eléctrico total.


## Revisión final del 8 de octubre de 2026: campañas anteriores USADAS; único RHS fotónico inicial pendiente

Las cinco entradas anteriores de transporte/borde están ejecutadas y revisadas. No repetirlas. El piloto inicial real detectó R1 de orden uno cerca del gap; integrar sus momentos es el único lote nuevo necesario antes de unir la trayectoria. No cambia F256, no corrige calor, no avanza tiempo y no hace reset total.

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_photon/run_initial_rhs_audit.py --plan docs/implementation/stage5/final_review_20261008/initial_photon_rhs_plan.json --output /home/jdiaz/scratch/stage5_photon_initial_rhs_20261008 --workers 4 --execute
```

Planificación: cuatro workers con BLAS/OpenMP=1, presupuesto compartido <=90 % de CPU/RAM disponible, reserva 7 GiB, orden 1-2 GB en scratch. Coste nominal de minutos (1-5 min) aún incierto por picos espectrales: ejecución manual. Salida nueva; preservar summary, per_energy, spectral_fields, recursos y fallos. Cuaderno vigente: /home/jdiaz/GEMINGA_COMMANDS_PRESENT.md. Informe final y tres documentos matemáticos: output/pdf/implementation/.


## Revisión del lote fotónico inicial completo del 8 de octubre - 224 consultas USADAS

El lote terminó en45,808s. No repetirlo: todavía requiere cierre constitutivo conjunto y no admite photon-runs ni push del modelo concluido. No hay otra campaña pesada pendiente; el trabajo siguiente es de implementación. Pipelines de revisión y visual ligeras, OPCIONALES:

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/04_initial_photon_rhs_review.py --input /home/jdiaz/scratch/stage5_photon_initial_rhs_20261008 --run-name stage5_initial_photon_rhs_review_user_20261008 --big-data-root /home/jdiaz/scratch/big_data
```

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python -u plot_pipelines/04_plot_initial_photon_rhs_review.py --input /home/jdiaz/scratch/big_data/raw/stage5_initial_photon_rhs_review_user_20261008/photon --run-name stage5_initial_photon_rhs_review_user_20261008 --big-data-root /home/jdiaz/scratch/big_data
```

Son replay inicial cero pasos, no nueva simulación. Medidos~0,76/~4,4s, CPU1/BLAS1 y512MiBreservados, todoen scratch/big_data, runnamesfresh. Productos usadosya disponiblesen stage5_initial_photon_rhs_review_20261008_v2. Decisión/hipótesis/límites en docs/implementation/stage5/completed_photon_rhs_20261008/; cuaderno actual /home/jdiaz/GEMINGA_COMMANDS_PRESENT.md.


## Corrección de comandos opcionales: cachés en scratch (8 de octubre de 2026)

Las dos pipelines de revisión inicial ya se ejecutaron; no son photon-runs y no requieren cálculo pesado. Sus comandos de reproducción usan cachés en scratch y un nombre nuevo. No hay comando de transiente admitido.

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache /home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/04_initial_photon_rhs_review.py --input /home/jdiaz/scratch/stage5_photon_initial_rhs_20261008 --run-name stage5_initial_photon_rhs_review_20261008_user01 --big-data-root /home/jdiaz/scratch/big_data
```

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib /home/jdiaz/.conda/envs/snspd/bin/python -u plot_pipelines/04_plot_initial_photon_rhs_review.py --input /home/jdiaz/scratch/big_data/raw/stage5_initial_photon_rhs_review_20261008_user01/photon --run-name stage5_initial_photon_rhs_review_20261008_user01 --big-data-root /home/jdiaz/scratch/big_data
```


<!-- stage51_unified_rhs_20261009 -->
# Geminga: etapa 5.1 — único sistema acoplado — 9 de octubre de 2026

La tarea activa es 5.1: resolver conjuntamente gap, hL/hT por energía/nodo, neutralidad electroquímica, contactos superconductores con amplitud/q/fase móviles y circuito CM. Se implementaron también la tripleta electrón-fonón Born R/A/K y el RHS fonónico independiente del mismo diagrama; su entrada común exige consultas energéticas desplazadas y cuadraturas emparejadas. No se reutiliza BGK neutro como fuente Nambu general.

**5.1 está en pausa por esta única corrida larga necesaria. Etapa 5 no concluida.** No hay informes nuevos ni comandos de transiente completo/reset. Primero se resuelve el sistema conjunto en el estado real inicial de 775 nm; el material en esta consulta usa únicamente el balance detallado FD/BE local inicial. No acredita después una evolución sin colisiones ni sustituye el proveedor general. No hay cálculo pesado de 1550 nm solicitado antes de revisar éste.

Pruebas: 79 aprobadas en 2,55 s; dos pruebas de pipeline repetidas tras las últimas correcciones pasan en 0,37 s. Benchmark completo: 112 energías, 1712 nodos, 4940 incógnitas; construcción 12,58 s, evaluación 9,36 s, memoria residente 5,46 GB por proceso. Se revalidan las raíces guardadas y se preservan identidades de datos, circuito y código. Cero nuevas raíces o depósitos en esta corrida.

## NECESARIO — sistema real conjunto de 775 nm

La ejecución ensambla las 4940 columnas por eliminaciones espectrales/cinéticas, comprueba el único gauge global antes de fijar referencia, realiza SVD y reevalúa independientemente todas las filas del candidato. Conserva residuos Nambu completos, corriente interior y en ambos puertos, normalización, soporte Pauli, proporciones R1/R0 y K1/K0 y el balance eléctrico CM. No elimina otras singularidades ni introduce fricción, corriente óhmica añadida o calor-residuo.

```bash
cd /home/jdiaz/pysnspd
set -o pipefail
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_photon/run_joint_snapshot.py \
  --initial /home/jdiaz/pysnspd/docs/implementation/stage5/photon_readiness_20261008/fd_be_775_variance_1p0_v2_uv24 \
  --output /home/jdiaz/scratch/stage51_joint_775_20261009 \
  --workers 8 --execute \
  2>&1 | tee -a /home/jdiaz/scratch/stage51_joint_775_20261009.console.log
```

Estimación orientativa: aproximadamente 1,6 h de ensamblado ideal con ocho trabajadores a partir del benchmark serial; reservar 2–4 h por concurrencia/memoria y SVD, sin garantía de duración. Ocho trabajadores y un coordinador, BLAS/OpenMP=1, afinidad distribuida por núcleos físicos. El programa vuelve a leer cuotas/afinidad/RAM, respeta <=90 % y rechaza un presupuesto insuficiente. La reserva conservadora actual ronda 52 GB y se calcula realmente al ejecutar.

Salidas: `/home/jdiaz/scratch/stage51_joint_775_20261009/`. Guarda identidad completa, cada columna con hash, matriz, valores singulares, candidato con campos SI, recibos espectrales y `summary.json`. El estado de 5.1/5 sigue explícitamente pendiente: resolver las filas proyectadas no admite por sí solo el resto de Nambu, el orden de gradiente, energía ni tiempo.

## REANUDAR — únicamente si se interrumpe antes de terminar

Usar el mismo comando anterior añadiendo `--resume` junto a `--execute`. No modificar código, entradas, metadatos de referencia ni material entre ejecuciones. Verifica hashes antes de reutilizar columnas; si difieren, rechaza la mezcla. Un resultado terminado no se sobrescribe. El log usa `tee -a`.

## Secuencia global conservada

5.1: acoplamiento único y su solución consistente; luego proveedor material general emparejado y avance del mismo estado. 5.2: balance energético independiente al mismo orden, incluyendo energía/cortes superconductores y exterior CM. 5.3: soporte físico, precisión temporal/espacial y admisión de las corridas hasta detección más margen. No se repiten por defecto DC, UV ni campañas antiguas. Los nominales Allmaras 5/24,7 ps son calibradores abiertos, no valores universales ni centro de un rango impuesto.

El historial general se preserva en `/home/jdiaz/GEMINGA_COMMANDS.md`. La copia previa de este archivo queda en `/home/jdiaz/scratch/stage5_joint_51_development_20261009/GEMINGA_COMMANDS_PRESENT_before_20261009.md`.


---

# Geminga: etapa 5.1 — primer paso físico acoplado — 9 de octubre de 2026

**NECESARIO: ejecutar una sola campaña, el comando de abajo. Etapa 5.1/5 pendiente.** No hay comandos de reset completo ni informes nuevos. La campaña anterior ya terminó y NO debe repetirse.

## Resultado recibido

`/home/jdiaz/scratch/stage51_joint_775_20261009`: las 4940 columnas terminaron en unas dos horas. El rango 1838 de la SVD truncada estaba dominado por escala global; equilibrado diagonal y QR sobre la MISMA matriz recuperan 4939/4939 columnas físicas tras fijar únicamente el gauge. Esto no salva el cierre instantáneo: residuo original ≈0,0554 frente a tolerancia 1e-7, tasas enormes y potenciales de escala voltios. Se conserva el fallo; no se añade fricción para ocultarlo.

## Qué hace el nuevo comando

Resuelve conjuntamente el primer paso Rosenbrock–Euler linealmente implícito: rigidez del gap, cinética completa hL/hT con almacenamiento más difusión `M+dt Lh`, neutralidad, ambos puertos, continuación superconductora móvil y CM implícito. Compara 1 fs y 0,5 fs desde la misma preparación física; no los encadena ni vuelve a depositar el fotón. Reconstruye R0 y la referencia Matsubara en cada estado finito, guardando gap, enlaces, poblaciones, corrientes y energía observada en SI.

La reconstrucción puede observar un estado formal con exceso Pauli <=1e-7, sin alterar sus campos: el soporte físico sigue evaluándose por separado con 32 eps y permanece rechazado si corresponde. Se cuantifica la energía asociada al exceso con el DOS reconstruido. Ningún paso se declara aceptado por esa observación. Se impide invertir Ib/Is globales.

El intercambio electrón–fonón explícito es cero SOLO en el estado local FD/BE inicial. El escape fonónico se resuelve implícitamente; `--tau-escape-ps 15` es la entrada efectiva heredada de D.19, libre y sin atribuirle universalidad ni origen Allmaras. Los nominales Allmaras τee(Tc)=5 ps y τep(Tc)=24,7 ps siguen siendo calibradores libres. La anchura Born física no se ha resumado en R0: su régimen y el proveedor material general siguen pendientes. Este ensayo no es backward Euler no lineal convergido ni una trayectoria hasta detección.

```bash
cd /home/jdiaz/pysnspd
set -o pipefail
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage5_photon/run_joint_first_step.py \
  --initial /home/jdiaz/pysnspd/docs/implementation/stage5/photon_readiness_20261008/fd_be_775_variance_1p0_v2_uv24 \
  --snapshot /home/jdiaz/scratch/stage51_joint_775_20261009 \
  --output /home/jdiaz/scratch/stage51_first_step_775_20261009 \
  --dt-ps 0.001 0.0005 --tau-escape-ps 15 --workers 8 --execute \
  2>&1 | tee -a /home/jdiaz/scratch/stage51_first_step_775_20261009.console.log
```

Recursos: 8 trabajadores + coordinador, un hilo BLAS/OpenMP cada uno, núcleos reservados y presupuesto <=90 % de cuotas/afinidad/RAM reales. Piloto completo: 112 energías, 1712 nodos, 4940 incógnitas y 256 modos térmicos; 99,03 s totales, 24,49 s para dos evaluaciones, error de gauge 1,54e-13. Reserva conservadora medida 85.610.993.280 bytes (~79,7 GiB); el programa verifica el presupuesto al ejecutar. Reservar orientativamente 5–8 h para las dos matrices y reevaluaciones, sin garantía de duración; no son 200 ps de transiente.

Salidas: `/home/jdiaz/scratch/stage51_first_step_775_20261009/`: identidad y hashes, columnas y matrices por paso, QR equilibrado, candidatos completos, reconstrucción espectral independiente, balances fonónicos y energía neutral adiabática parcial; `summary.json` de campaña. Esa energía no incluye todavía completado K1/carga/trabajo de cortes y no acredita 5.2.

## Reanudación y decisión

Sólo si se interrumpe: mismo comando, mismos parámetros, añadiendo `--resume`. No cambiar código ni entradas; se verifican hashes. Un resultado final no se sobrescribe. No repetir la campaña anterior.

La comparación de incrementos y reconstrucciones decidirá si hay límite temporal regular o una capa/proyección algebraica inicial. Un salto independiente de dt o primer gradiente descontrolado rechaza ese arranque; no se prolonga hasta detección. Después corresponde el proveedor material general acoplado, 5.2 balance independiente completo y 5.3 convergencia/admisión. No repetir por defecto DC ni UV. El modelo de producción conserva su estado y no se declara cierre por una solución algebraica aislada.


## Verificacion y publicacion del avance

Push confirmado a main: abbeb7a678f9cedc2ffc526c5edc430b790b5482 (GitHub). Pasaron 182 pruebas y 11 subpruebas; recibo del nucleo en /home/jdiaz/scratch/stage51_first_step_development_20261009/quick_test_receipt.json, con una correccion posterior de espacio en blanco verificada por AST equivalente. El piloto real esta en /home/jdiaz/scratch/stage51_first_step_pilot_20261009/pilot.json. La etapa 5.1/5 sigue pendiente; el push publica la implementacion experimental, no el cierre fisico del detector. El unico comando necesario permanece arriba, sin cambios.

## Reanudación integral autorizada — 11 de octubre de 2026

El usuario autoriza corridas largas directas, cierre completo, informe final PDF y push final. Siguen vigentes scratch, jdiaz y presupuesto global de recursos <=90%.

La campaña stage51_first_step_775_20261009 terminó: los pasos de 1 fs y 0,5 fs fallan soporte Pauli y control del primer gradiente, especialmente en la región fría. No se repite ni se prolonga hasta detección. Auditoría: /home/jdiaz/scratch/stage51_failure_audit_20261011/audit.json.

Decisión de desarrollo: cierre neutral efectivo de estados BCS nativos, energía y corriente comunes, rigidez dirty radial, movilidad KWT calibrable, cinética electrón–fonón conservativa y circuito CM. No se presenta como evolución microscópica Usadel exacta. Requiere admisión DC y fotónica nueva. Las preparaciones DC nuevas sobre las mallas heredadas están en /home/jdiaz/scratch/stage5_effective_bcs_development_20261011.

Comando de diagnóstico usado (no repetir por defecto):

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_coupled_first_pilot.py

Los comandos de las pipelines finales se registrarán después de validar la trayectoria. Ningún comando autorizado de esta actualización apunta al reset eléctrico completo.

## 2026-10-11 - evolucion fotonica efectiva autorizada directamente

Dos casos concurrentes, cuatro nucleos fisicos cada uno, BLAS1, datos y logs en scratch. Criterios fijados en sandbox/stage5_effective_validation/admission_policy.json. Ninguna corrida completada se admite automaticamente.

Hasta gatillo+5ps o censura200ps, sin reset. Estimacion inicial: varias horas por caso, actualizada con progreso real.

cd /home/jdiaz/pysnspd && env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache /home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_dt002_20261011 --length-nm 320 --spacing-nm 5 --wavelength-nm 775 --dt-ps 0.02 --threads 4 --cpu-start 0 > /home/jdiaz/scratch/stage5_effective_L320_h5_775_dt002_20261011.console.log 2>&1

Hasta gatillo+5ps o censura200ps, sin reset. Estimacion inicial: varias horas por caso, actualizada con progreso real.

cd /home/jdiaz/pysnspd && env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache /home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_1550_dt002_20261011 --length-nm 320 --spacing-nm 5 --wavelength-nm 1550 --dt-ps 0.02 --threads 4 --cpu-start 4 > /home/jdiaz/scratch/stage5_effective_L320_h5_1550_dt002_20261011.console.log 2>&1

### Revision numerica de la primera evolucion efectiva

Las tres trayectorias dt0.02/dt0.01 fueron interrumpidas controladamente a ~0.14ps por colapso adaptativo; sus checkpoints y logs se conservan. No reanudar tras cambios de fuentes. Se reparo la discontinuidad del limitador min(1,0.95A/L), y el criterio poblacional ahora usa defecto L1 de energia fisica con proteccion independiente de fuerza sobre gap. No se cambia el RHS ni las tolerancias gap/CM. Evidencia /home/jdiaz/scratch/stage5_solver_diagnostic_20261011/.

### Campana v3 con norma de energia fisica y guarda gap

Eref=energia retenida, fija al comparar malla/dominio. Limiter continuo. Hasta gatillo+5ps o censura200ps. Dos procesos de4cores fisicas, BLAS1. Resultadosv1 interrumpidos preservados.

cd /home/jdiaz/pysnspd && env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache /home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_dt002_v3_20261011 --length-nm 320 --spacing-nm 5 --wavelength-nm 775 --dt-ps 0.02 --threads 4 --cpu-start 0 > /home/jdiaz/scratch/stage5_effective_L320_h5_775_dt002_v3_20261011.console.log 2>&1

cd /home/jdiaz/pysnspd && env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache /home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_1550_dt002_v3_20261011 --length-nm 320 --spacing-nm 5 --wavelength-nm 1550 --dt-ps 0.02 --threads 4 --cpu-start 4 > /home/jdiaz/scratch/stage5_effective_L320_h5_1550_dt002_v3_20261011.console.log 2>&1

### Contraste temporal completo v3

cd /home/jdiaz/pysnspd && env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache /home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_dt001_v3_20261011 --length-nm 320 --spacing-nm 5 --wavelength-nm 775 --dt-ps 0.01 --threads 4 --cpu-start 8 > /home/jdiaz/scratch/stage5_effective_L320_h5_775_dt001_v3_20261011.console.log 2>&1

### Control de malla fotonico v3

Malla3.5nm, film7nm, L320nm,775nm, dt.02ps. Dos coresfisicas12,13; totalcampana14/16cores87.5%, BLAS1. Misma ventana gatillo+5ps/censura200ps. Estimacion: varias horas, actualizada con avance.

cd /home/jdiaz/pysnspd && env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache /home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h35_775_dt002_v3_20261011 --length-nm 320 --spacing-nm 3.5 --wavelength-nm 775 --dt-ps 0.02 --threads 2 --cpu-start 12 > /home/jdiaz/scratch/stage5_effective_L320_h35_775_dt002_v3_20261011.console.log 2>&1


## 2026-10-11 — ensayo aislado de checkpoints (ejecutado, CPU13)

Diagnóstico numérico sin promoción: copias atómicas de los checkpoints 775 nm L320 h5 y h3.5, fuentes congeladas y verificadas por SHA. B: un paso de 0.02 ps con tolerancias 2e-9/2e-8/2e-7 y 120 iteraciones. A: 2 ps de continuación por hmax=0.02/0.05/0.1 ps, tolerancia2e-9. Una CPU física, Numba1 y BLAS1; no modifica ni reinserta resultados en las campañas vigentes. Corrida larga autorizada directamente por usuario2026-10-11.

```bash
cd /home/jdiaz/pysnspd
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/pycache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib TMPDIR=/home/jdiaz/scratch/tmp/jdiaz NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache PYTHONPATH=/home/jdiaz/pysnspd taskset -c 13 /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage5_effective_validation/probe_frozen_checkpoint.py --output /home/jdiaz/scratch/stage5_frozen_solver_probe_20261011 --phase all
```

Salida: experiment_receipt.json, fuentes/checkpoints congelados, B_single_steps/*_receipt.json y trazasJSONL, A_continuations/*/checkpoint_latest.npz/step_receipts.json/receipt.json. El directorio existente se preserva: para repetir utilizar otro nombre scratch.


## 2026-10-11 — continuación numérica práctica desde frozenA (autorizada, CPU13)

Finalidad: comparar durante los mismos2ps la tolerancia2e-7 frente a baseline2e-9 en hmax0.02ps, y después hmax0.05/0.1 frente a0.02 con tolerancia idéntica2e-7. Mismo estado inicial, circuito, fuente ya aplicada, ledger y controlador; mismo RHS congelado. No modifica ni reinserta datos en las campañas vigentes. UnaCPU13, Numba1, BLAS1; tiempo esperado10–20min. Datos nuevos exclusivamente scratch. Criterios predeclarados: curvas/gap≤2% (señal real, no normalización al gatillo), balance integrado≤1% del denominador físico original, Pauli/Bose y corrientes positivas estrictos.

```bash
cd /home/jdiaz/pysnspd
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/pycache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib TMPDIR=/home/jdiaz/scratch/tmp/jdiaz NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache PYTHONPATH=/home/jdiaz/pysnspd taskset -c 13 /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage5_effective_validation/probe_frozen_checkpoint.py --output /home/jdiaz/scratch/stage5_frozen_solver_probe_20261011 --phase practical
```

Salidas: /home/jdiaz/scratch/stage5_frozen_solver_probe_20261011/A_practical_norm_continuations/{hmax_0p02,hmax_0p05,hmax_0p1}/, checkpoint completo, trajectory/snapshots, step_receipts.json y receipt.json; la procedencia del helper3 está en practical_helper_receipt.json. Baseline estricto preservado en A_continuations/hmax_0p02. Ninguna admisión implícita.

## 2026-10-11 — malla fina: candidato de tolerancia, continuación explícita

Captura coherente v3 a 0.73 ps; se preserva estado físico, fuente única, circuito y balance. Candidato solver 2e-7 frente a 2e-9; no es admisión. Mismos 2 núcleos físicos 12,13; BLAS 1, memoria <4 GB; estimación hasta 20 ps 1–3 h, techo 200 ps sólo si sigue activo. Fuente congelada y recibo: /home/jdiaz/scratch/stage5_effective_L320_h35_775_dt002_tol2e7_v4_candidate_20261011/provenance/transfer_receipt.json. Se detiene sólo PID 160427 y se preserva su campaña histórica.

```bash
bash /home/jdiaz/scratch/stage5_effective_L320_h35_775_dt002_tol2e7_v4_candidate_20261011/provenance/resume_command.txt > /home/jdiaz/scratch/stage5_effective_L320_h35_775_dt002_tol2e7_v4_candidate_20261011.console.log 2>&1
```

Comando expandido:
```bash
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python /home/jdiaz/scratch/stage5_effective_L320_h35_775_dt002_tol2e7_v4_candidate_20261011/runtime_sources/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h35_775_dt002_tol2e7_v4_candidate_20261011 --resume --cpu-start 12 --length-nm 320.0 --spacing-nm 3.5 --xi-cutoff 32.0 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 775.0 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.02 --ceiling-ps 200.0 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5.0 --tau-ee-Tc-ps 5.0 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15.0 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1.0 --checkpoint-wall-s 60.0 --step-tolerance 2e-07 --max-iterations 40 --max-rejections 50 --threads 2
```


Actualización de recursos del ensayo práctico: el proceso privado se trasladó deCPU13 aCPU7 por autorización del coordinador, para liberar el trabajador de malla fina. Conserva un solo núcleo y BLAS/Numba1, conjunto global0..13. Recibo affinity_segment_13_to_7.json; las comparaciones de coste no atribuyen toda diferencia a tolerancia. Próximas ejecuciones de este diagnóstico/perfil: taskset -c7. Las corridas físicas y sus pins permanecen intactos.

## 2026-10-11 — control temporal: continuación explícita a tolerancia práctica

Ancla v3 coherente a 12.42 ps; mismo dt máximo 0.01 ps, física/fuente/CM/ledger exactos; tol 2e-9→2e-7, threads 4→2 declarados. Candidato de validación, no admisión. Cores físicos8,9; BLAS1; memoria<4GB; estimación hasta20ps0.5–2h. Sólo PID156492 se detiene; originales preservados. Comando expandido y SHA: /home/jdiaz/scratch/stage5_effective_L320_h5_775_dt001_tol2e7_v4_candidate_20261011/provenance/transfer_receipt.json y resume_command.txt.

```bash
bash /home/jdiaz/scratch/stage5_effective_L320_h5_775_dt001_tol2e7_v4_candidate_20261011/provenance/resume_command.txt > /home/jdiaz/scratch/stage5_effective_L320_h5_775_dt001_tol2e7_v4_candidate_20261011.console.log 2>&1
```

## 2026-10-11 — control del dominio espacial efectivo, 0–20 ps

Referencia interior L320nm frente a L480nm, mismo W80nm y film7nm, photon775nm, nueva preparación DC a igual corriente. Sin reset completo. Se prueba sensibilidad a reservorios/cortes; el alcance de20ps se ampliará desde checkpoint sólo si el evento físico lo requiere. Runtime v4 congelado con mismos13pins; tol2e-7 candidata; dt máximo.02ps.2núcleos físicos10,11, BLAS1; memoria<4GB; estimación1–3h. Salidas: manifiesto, historial SI, mapas, checkpoint, DC/material; no admisión automática.

```bash
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_runtime_v4_20261011/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L480_h5_775_dt002_tol2e7_v4_control20_20261011 --length-nm 480 --spacing-nm 5 --wavelength-nm 775 --dt-ps 0.02 --step-tolerance 2e-7 --duration-ps 20 --threads 2 --cpu-start 10 > /home/jdiaz/scratch/stage5_effective_L480_h5_775_dt002_tol2e7_v4_control20_20261011.console.log 2>&1
```

## 2026-10-11:1550nm v3→v4 desde27.42ps, lanzamiento pendiente de coordinación
Tolerancia no lineal2e−7 admitida por comparación A/B; dt0.02ps unchanged. Recursos2cores físicos CPU4/5,BLAS1. No reset completo; gatillo+0.5psconfirmación+5psextra,techo200ps. Fuente ya depositada y todos los balances/historia preservados. 
Ancla: /home/jdiaz/scratch/stage5_effective_L320_h5_1550_dt002_v3_anchor_20261011_PRIMARY_PENDING
Recibo: /home/jdiaz/scratch/stage5_effective_L320_h5_1550_dt002_tol2e7_v4_primary_candidate_20261011_PENDING/provenance/transfer_receipt.json
Original detenido156224; preservación final: /home/jdiaz/scratch/stage5_effective_L320_h5_1550_dt002_tol2e7_v4_primary_candidate_20261011_PENDING/provenance/post_stop_original/preservation_receipt.json
Comando expandido preparado; se ejecutará sólo después del acuse root.
```bash
set -o pipefail
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_L320_h5_1550_dt002_tol2e7_v4_primary_candidate_20261011_PENDING/runtime_sources/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_1550_dt002_tol2e7_v4_primary_candidate_20261011_PENDING --resume --cpu-start 4 --length-nm 320.0 --spacing-nm 5.0 --xi-cutoff 32.0 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 1550.0 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.02 --ceiling-ps 200.0 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5.0 --tau-ee-Tc-ps 5.0 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15.0 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1.0 --checkpoint-wall-s 60.0 --step-tolerance 2e-07 --max-iterations 40 --max-rejections 50 --threads 2 2>&1 | tee /home/jdiaz/scratch/stage5_effective_L320_h5_1550_dt002_tol2e7_v4_primary_candidate_20261011_PENDING/execution_v4.console.log
```

Lanzado después del acuse root: UTC 2026-10-11T14:55:56.937426+00:00, PID 178314, sesión foreground19088, afinidad4/5,2threads. Recibo /home/jdiaz/scratch/stage5_effective_L320_h5_1550_dt002_tol2e7_v4_primary_candidate_20261011_PENDING/provenance/launch_receipt.json. Admisión física pendiente.


## 2026-10-11 — Control de paso tardío1550, ejecución autorizada directa

Mismo checkpoint completo, dos continuaciones de2ps con hmax0.02 y0.1ps, tolerancia no lineal2e-7. CPU1,1hilo, BLAS1; no modifica corridas activas ni reinicia depósito. Se estima5–15min; salidas y fuentes congeladas en /home/jdiaz/scratch/stage5_1550_late_step_control_20261011. Criterios2%observables/1%energía; no admite automáticamente etapa5.

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache taskset -c 1 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_1550_late_step_control_20261011/run_control.py > /home/jdiaz/scratch/stage5_1550_late_step_control_20261011/execution.console.log 2>&1
```

Control1550 tardío: primer arranque falló antes de evolución por nombres de captura *_original. Se añadieron copias compatibles byteidénticas; salida fallida preservada, receta física/checkpoint intactos. Recibo provider_alias_receipt.json. Se repite el mismo execute.sh.


### 2026-10-11 — Perfil caliente congelado de coste, después del control A

Finalidad: cinco RHS de punto medio calientes y un paso implícito diagnóstico de 0.02 ps; medir CPU y pared por separado, nodos activos y buffers de eventos. Sin cambiar física, fuente, corridas primarias ni resultados de A. El paso de perfil no se reintegra en ninguna campaña.

Recursos: 1 CPU física7, Numba1, BLAS1; el resto del presupuesto compartido permanece intacto. Se esperan menos de5min. Salidas: scratch/stage5_frozen_solver_probe_20261011/warm_profile/{summary.json,midpoint.prof,five_evaluate.prof,*_top.txt}. El microbenchmark de buffers sólo es un proxy de asignación/materialización; no es una descomposición exacta del dispatcher Numba.

```bash
cd /home/jdiaz/pysnspd
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/pycache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib TMPDIR=/home/jdiaz/scratch/tmp/jdiaz NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache PYTHONPATH=/home/jdiaz/pysnspd taskset -c 7 /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage5_effective_validation/profile_frozen_checkpoint.py /home/jdiaz/scratch/stage5_frozen_solver_probe_20261011
```


## 2026-10-11 — 775 nm v5, guarda energética local desde ancla segura12,4ps
Estado: PREPARADO; todavía sin lanzar. Continuación explícita del parent A estricto; no se reaplica el fotón. Se preservan estado p/n/Delta/CM, historial, fuente y ledger físico completos. Cambios numéricos declarados: hmax=0,02ps, tolerancia no lineal2e-7, guarda local1e-5. Guarda: |DeltaE-DeltaWmid| <= 1e-5*(Eret+sum|DeltaWmid aceptado|+|DeltaWmid prueba|). El tramo previo sin guarda queda etiquetado, no certificado retrospectivamente. Recursos previstos:4 núcleos físicos0..3, Numba4, BLAS/OpenMP1, datos/cachés scratch. Objetivo: gatillo+.1mV, confirmación.5ps, margen5ps, techo200ps; sin reset completo. ETA desde el avance de esta sesión, no desde t=0 absoluto.

Parent inalterado: /home/jdiaz/scratch/stage5_effective_L320_h5_775_A_strict_12p4_parent_20261011
Checkpoint parent SHA256: d5fc7358dfda9fd0f4e41361bd80d193adf2cc4a88870921340940d704f2b99a
Receipt: /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v5_20261011/provenance/transfer_receipt.json
Receipt SHA256: 05307fe65eeb69cff2a1b7fa01b9f5e683b9e34bf4f7ed5204b2631ba4450be8
Wrapper: /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v5_20261011/provenance/launch_command_expanded.sh
Salida: /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v5_20261011

```bash
set -o pipefail
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=4 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v5_20261011/runtime_sources/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v5_20261011 --resume --cpu-start 0 --length-nm 320.0 --spacing-nm 5.0 --xi-cutoff 32.0 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 775.0 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.02 --ceiling-ps 200.0 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5.0 --tau-ee-Tc-ps 5.0 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15.0 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1.0 --checkpoint-wall-s 60.0 --step-tolerance 2e-07 --max-iterations 40 --max-rejections 50 --threads 4 --local-energy-tolerance 1e-05 2>&1 | tee /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v5_20261011/execution_v5.console.log
```

775 v5 lanzado en foreground: PID184461, sesión72323; afinidad física[0,1,2,3], Numba4, BLAS/OpenMP1. Recibo: /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v5_20261011/provenance/launch_receipt.json.


## 2026-10-11: control espectral completo v5, preparado

Refina la trayectoria775nm de32/6/96/48 a64/8/144/96; misma física, malla, guarda energética y fuente. Ventana inicial20ps, ampliable por checkpoint si la región crítica lo exige. No implica cierre de etapa5. Recursos:2 núcleos físicos CPU6/7, BLAS1, variosGB de RAM; duración prevista de varias horas. Esperar liberación coordinada deCPU6/7. Todo el dato y caché en scratch.

```bash
bash /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_v5_control20_20261011/launch_spectral_control.sh
```

Comando expandido:

```bash
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=2 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v5_20261011/runtime_sources/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_v5_control20_20261011 --cpu-start 6 --threads 2 --length-nm 320 --spacing-nm 5 --xi-cutoff 64 --xi-order 8 --phonon-nodes 144 --transport-order 96 --wavelength-nm 775 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.02 --local-energy-tolerance 1e-5 --step-tolerance 2e-7 --duration-ps 20 --ceiling-ps 200 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5 --tau-ee-Tc-ps 5 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1 --checkpoint-wall-s 60 --max-iterations 40 --max-rejections 50 2>&1 | tee /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_v5_control20_20261011/execution.console.log
```


## 2026-10-11 — continuación v5 fine, PREPARADA sin lanzamiento
Se conserva el proceso v4 original todavía activo, sus fuentes, checkpoint coherente, historial/fuente/ledger y log. stop_intent.json identifica UID/PID/start_ticks antes de cualquier señal. El siguiente proceso sólo podrá iniciar tras salir el anterior y tras mostrar este comando en chat. No se reaplica el fotón; RHS idéntico y control de gatillo/margen/techo conservado. NL2e-7, hcap=0.02ps, guarda local=1e-05, Numba2 y BLAS/OpenMP1; núcleos físicos[12, 13].
Parent elegido: /home/jdiaz/scratch/stage5_effective_L320_h35_775_v4_to_v5_anchor_20261011
Salida: /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v5_20261011
Receipt SHA256: b439f664275cb5984dd4d38aa092bab408d7b182b8cc98f7f68f609551b2f579
Wrapper SHA256: 2c4aa693c34860f127b276f4efaf28f73ea8a8febbc0760807f0528499e96cf9

```bash
set -o pipefail
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=2 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v5_20261011/runtime_sources/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v5_20261011 --resume --cpu-start 12 --length-nm 320.0 --spacing-nm 3.5 --xi-cutoff 32.0 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 775.0 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.02 --ceiling-ps 200.0 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5.0 --tau-ee-Tc-ps 5.0 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15.0 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1.0 --checkpoint-wall-s 60.0 --step-tolerance 2e-07 --max-iterations 40 --max-rejections 50 --threads 2 --local-energy-tolerance 1e-05 2>&1 | tee /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v5_20261011/execution_v5.console.log
```


## 2026-10-11 — continuación v5 domain, PREPARADA sin lanzamiento
Se conserva el proceso v4 original todavía activo, sus fuentes, checkpoint coherente, historial/fuente/ledger y log. stop_intent.json identifica UID/PID/start_ticks antes de cualquier señal. El siguiente proceso sólo podrá iniciar tras salir el anterior y tras mostrar este comando en chat. No se reaplica el fotón; RHS idéntico y control de gatillo/margen/techo conservado. NL2e-7, hcap=0.02ps, guarda local=1e-05, Numba2 y BLAS/OpenMP1; núcleos físicos[10, 11].
Parent elegido: /home/jdiaz/scratch/stage5_effective_L480_h5_775_v4_to_v5_anchor_20261011
Salida: /home/jdiaz/scratch/stage5_effective_L480_h5_775_adaptive_v5_control20_20261011
Receipt SHA256: 0bb03c18e56738cb00658ab9214987984156093108bc251cf01e6175d81a80a7
Wrapper SHA256: 1b092f5124d8475f05bdb96b970a9e4dd1ac6cd4d224be77313c89321d7f6669

```bash
set -o pipefail
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=2 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_L480_h5_775_adaptive_v5_control20_20261011/runtime_sources/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L480_h5_775_adaptive_v5_control20_20261011 --resume --cpu-start 10 --length-nm 480.0 --spacing-nm 5.0 --xi-cutoff 32.0 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 775.0 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.02 --ceiling-ps 200.0 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5.0 --tau-ee-Tc-ps 5.0 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15.0 --current-uA 21.29720614489798 --bath-K 0.9 --duration-ps 20.0 --snapshot-every-ps 0.1 --checkpoint-every-ps 1.0 --checkpoint-wall-s 60.0 --step-tolerance 2e-07 --max-iterations 40 --max-rejections 50 --threads 2 --local-energy-tolerance 1e-05 2>&1 | tee /home/jdiaz/scratch/stage5_effective_L480_h5_775_adaptive_v5_control20_20261011/execution_v5.console.log
```


## 2026-10-11 — continuación v5 temporal, PREPARADA sin lanzamiento
Se conserva el proceso v4 original todavía activo, sus fuentes, checkpoint coherente, historial/fuente/ledger y log. stop_intent.json identifica UID/PID/start_ticks antes de cualquier señal. El siguiente proceso sólo podrá iniciar tras salir el anterior y tras mostrar este comando en chat. No se reaplica el fotón; RHS idéntico y control de gatillo/margen/techo conservado. NL2e-7, hcap=0.01ps, guarda local=1.25e-06, Numba2 y BLAS/OpenMP1; núcleos físicos[8, 9].
Refinamiento temporal explícito: desde ancla segura v3_B12,42ps; hcap0,01ps y guarda1,25e-6 (también refina el h efectivo en slips). El estado v4 posterior se conserva sólo como diagnóstico; no se usa de parent.
Parent elegido: /home/jdiaz/scratch/stage5_effective_L320_h5_775_dt001_v3_anchor_20261011_B
Salida: /home/jdiaz/scratch/stage5_effective_L320_h5_775_temporal_adaptive_v5_20261011
Receipt SHA256: 77bf633134ed18fd30c0775ce2f8e10fafe5adc33db06b836e7fea528e2c2a97
Wrapper SHA256: b99d979a3ce3274c246d5251c589c89c2a19e2334f418afe78eb3c792361ede3

```bash
set -o pipefail
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=2 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_L320_h5_775_temporal_adaptive_v5_20261011/runtime_sources/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_temporal_adaptive_v5_20261011 --resume --cpu-start 8 --length-nm 320.0 --spacing-nm 5.0 --xi-cutoff 32.0 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 775.0 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.01 --ceiling-ps 200.0 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5.0 --tau-ee-Tc-ps 5.0 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15.0 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1.0 --checkpoint-wall-s 60.0 --step-tolerance 2e-07 --max-iterations 40 --max-rejections 50 --threads 2 --local-energy-tolerance 1.25e-06 2>&1 | tee /home/jdiaz/scratch/stage5_effective_L320_h5_775_temporal_adaptive_v5_20261011/execution_v5.console.log
```


## 2026-10-11 — continuación v5 1550, PREPARADA sin lanzamiento
Se conserva el proceso v4 original todavía activo, sus fuentes, checkpoint coherente, historial/fuente/ledger y log. stop_intent.json identifica UID/PID/start_ticks antes de cualquier señal. El siguiente proceso sólo podrá iniciar tras salir el anterior y tras mostrar este comando en chat. No se reaplica el fotón; RHS idéntico y control de gatillo/margen/techo conservado. NL2e-7, hcap=0.1ps, guarda local=1e-05, Numba2 y BLAS/OpenMP1; núcleos físicos[4, 5].
hcap0,1ps tardío admitido por control C desde29,62ps, decision JSON en stage5_1550_late_step_control_20261011. Parent vigente capturado a38.540000000000084ps, sin repetir desde29,62ps.
Parent elegido: /home/jdiaz/scratch/stage5_effective_L320_h5_1550_v4_to_v5_anchor_20261011
Salida: /home/jdiaz/scratch/stage5_effective_L320_h5_1550_adaptive_v5_20261011
Receipt SHA256: d5a978fcc30c7ffc8af51fb18b9c8e47f8395822801e3e553a35c4fb06910547
Wrapper SHA256: baf88f68ed98c49b07755001e7e2a5f49d2ee628c5a88cccd0d6e0d27dd09d37

```bash
set -o pipefail
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=2 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_L320_h5_1550_adaptive_v5_20261011/runtime_sources/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_1550_adaptive_v5_20261011 --resume --cpu-start 4 --length-nm 320.0 --spacing-nm 5.0 --xi-cutoff 32.0 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 1550.0 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.1 --ceiling-ps 200.0 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5.0 --tau-ee-Tc-ps 5.0 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15.0 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1.0 --checkpoint-wall-s 60.0 --step-tolerance 2e-07 --max-iterations 40 --max-rejections 50 --threads 2 --local-energy-tolerance 1e-05 2>&1 | tee /home/jdiaz/scratch/stage5_effective_L320_h5_1550_adaptive_v5_20261011/execution_v5.console.log
```


### 2026-10-11: decisión del control tardío1550

Control desde el mismo estado29.62ps hasta31.62ps: hcap0.1ps frente0.02ps pasa, Vdev cambia0.170136%, Vout0.000570771%, mínimo del gap1.7429e-8 del DC y defecto energético global máximo2.8319e-8. Se admite sólo para1550 desde29.62ps en el caso físico actual, con guarda local v5 activa. Evidencia: /home/jdiaz/scratch/stage5_1550_late_step_control_20261011/late_step_admission_decision.json. No autoriza cap grande para775ni arranque inicial.


### 2026-10-11 — Candidato aislado de eventos por nodo: equivalencia exacta y coste

Dos pruebas secuenciales acotadas a300s cadauna, CPU7 exclusiva, Numba1/BLAS1. Primero oráculo congelado, luego candidato scratch. DC, fotón inicial, estado vorticial12.4ps, malla3.5nm a.405ps y soporte exactamente0/1/0fonones/sin caras/calor agotado. Salidas completas de RHS en results/{oracle,candidate}, recibos de fuentes/checkpoints, tiemposCPUcalientes y RSSmáximo real. No main, no jobs actuales, no cambio de ecuaciones ni tolerancias físicas.

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/pycache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib TMPDIR=/home/jdiaz/scratch/tmp/jdiaz NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache taskset -c 7 timeout --signal=TERM 300s /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_node_fusion_candidate_20261011/verify_candidate.py oracle
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/pycache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib TMPDIR=/home/jdiaz/scratch/tmp/jdiaz NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache taskset -c 7 timeout --signal=TERM 300s /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_node_fusion_candidate_20261011/verify_candidate.py candidate
```

v5 fine lanzado foreground: PID186958, sesión73040, afinidad[12, 13], Numba2/BLAS1; original v4 salió tras SIGTERM intencional y sus finales quedan en provenance/post_stop_original. Recibo: /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v5_20261011/provenance/launch_receipt.json.

v5 domain lanzado foreground: PID187172, sesión17590, afinidad[10, 11], Numba2/BLAS1; original v4 salió tras SIGTERM intencional y sus finales quedan en provenance/post_stop_original. Recibo: /home/jdiaz/scratch/stage5_effective_L480_h5_775_adaptive_v5_control20_20261011/provenance/launch_receipt.json.

v5 temporal lanzado foreground: PID187012, sesión32692, afinidad[8, 9], Numba2/BLAS1; original v4 salió tras SIGTERM intencional y sus finales quedan en provenance/post_stop_original. Recibo: /home/jdiaz/scratch/stage5_effective_L320_h5_775_temporal_adaptive_v5_20261011/provenance/launch_receipt.json.

v5 1550 lanzado foreground: PID187118, sesión65866, afinidad[4, 5], Numba2/BLAS1; original v4 salió tras SIGTERM intencional y sus finales quedan en provenance/post_stop_original. Recibo: /home/jdiaz/scratch/stage5_effective_L320_h5_1550_adaptive_v5_20261011/provenance/launch_receipt.json.


## 2026-10-11 — equivalencia de almacenamiento por nodo, NumBa2/4

Finalidad: comparar RHS completo, un midpoint775 desde12.4ps, soporte exacto y fallback genérico; candidato aislado, sin cambiar main ni corridas. Cada proceso limitado300s. Num2 usaCPU6/7 yBLAS1. Num4 usaCPU0..3 sólo bajo controlador que pausa/reanuda nominal775 verificado; sin solapar Num2. Resultados: /home/jdiaz/scratch/stage5_node_fusion_candidate_20261011/results.

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=2 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/pycache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib TMPDIR=/home/jdiaz/scratch/tmp/jdiaz NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache taskset -c 6,7 timeout --signal=TERM 300s /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_node_fusion_candidate_20261011/verify_candidate.py oracle --threads 2
```

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=2 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/pycache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib TMPDIR=/home/jdiaz/scratch/tmp/jdiaz NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache taskset -c 6,7 timeout --signal=TERM 300s /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_node_fusion_candidate_20261011/verify_candidate.py candidate --threads 2
```

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=4 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/pycache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib TMPDIR=/home/jdiaz/scratch/tmp/jdiaz NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache taskset -c 0,1,2,3 timeout --signal=TERM 300s /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_node_fusion_candidate_20261011/verify_candidate.py oracle --threads 4
```

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=4 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/pycache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib TMPDIR=/home/jdiaz/scratch/tmp/jdiaz NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache taskset -c 0,1,2,3 timeout --signal=TERM 300s /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_node_fusion_candidate_20261011/verify_candidate.py candidate --threads 4
```


Wrapper secuencial Num4 (sólo bajo controlador foreground con pausa/reanudación verificada, timeout total300s): /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_node_fusion_candidate_20261011/run_num4_pair.py

### Microtest Num4 con pausa acotada del nominal775 (2026-10-11)

Root autorizó SIGSTOP/SIGCONT del PID184461/start_ticks34576706, tras terminar Num2; no SIGTERM ni cambio de fuentes/estado. Timeout total300s, CPU0..3/Numba4/BLAS1, otros8workers activos. Siempre SIGCONT en finally. Wrapper SHA 7f6986570b6d8a27f530710849ffefb20e30a5d3ed079404052c363248cf8ad6.

```bash
env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 0-3 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_controls_20261011/bounded_num4_pause_controller.py --command-file /home/jdiaz/scratch/stage5_node_fusion_candidate_20261011/run_num4_pair_reviewed.sh --command-sha256 7f6986570b6d8a27f530710849ffefb20e30a5d3ed079404052c363248cf8ad6 --output /home/jdiaz/scratch/stage5_node_fusion_candidate_20261011/num4_pause_receipt --timeout-s 300
```

Num4 terminó PASS; PID184461 reanudado con SIGCONT, UID/start_ticks/cmd/afinidad intactos y ticks CPU avanzando. Checkpoint no cambió durante la pausa. 15 NPZ completos más1excepción iguales bit a bit. Recibo: /home/jdiaz/scratch/stage5_node_fusion_candidate_20261011/num4_pause_receipt/resume_and_outputs_verification.json SHA bbf2517d8e1f4c2e55a5925e330cc7800f139b7960246aec95f61854fa1675e5.


## 2026-10-11 — regresión efectiva tras optimización exacta

Suite test_effective* completa+API/imports/circuito/gatillo/soporte. Límite300s, CPU6/7, NumBa2, BLAS1; logs, JUnit, temporales y recibo en scratch/publication_checks; sin evolución de campañas.

```bash
cd /home/jdiaz/pysnspd && env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=2 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/pycache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib TMPDIR=/home/jdiaz/scratch/tmp/jdiaz NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache taskset -c 6,7 timeout --signal=TERM 300s /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_node_fusion_candidate_20261011/run_publication_checks.py
```

## 2026-10-11 — Comparación espectral v6, fuentes equivalentes optimizadas

Objetivo: evolución completa de control espectral hasta 20ps inicialmente; ampliar sólo si falta cubrir dinámica crítica. Mismo modelo y guarda1e-5, NL2e-7, hcap0.02ps; cuadratura64/8/144/96. Conserva la preparación v5 sin ejecutarla. Datos y cachés scratch, CPU6,7/NumBa2/BLAS1, presupuesto global14/16 núcleos físicos. Duración prevista: varias horas, checkpoint60s. Autorizada ejecución directa en este chat.

Comando reproducible:

```bash
bash /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_v6_control20_20261011/launch_spectral_control.sh
```

Comando expandido:

```bash
set -o pipefail
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=2 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_v6_control20_20261011 --cpu-start 6 --threads 2 --length-nm 320 --spacing-nm 5 --xi-cutoff 64 --xi-order 8 --phonon-nodes 144 --transport-order 96 --wavelength-nm 775 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.02 --local-energy-tolerance 1e-5 --step-tolerance 2e-7 --duration-ps 20 --ceiling-ps 200 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5 --tau-ee-Tc-ps 5 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1 --checkpoint-wall-s 60 --max-iterations 40 --max-rejections 50 2>&1 | tee /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_v6_control20_20261011/execution.console.log
```

Recibo: /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_v6_control20_20261011/preparation_receipt.json

### Continuación v6 preparada 775_primary (2026-10-11)

Sin señales ni lanzamiento todavía. Misma física, dt/tolerancias/hilos/guarda/controlador/historia; únicos2paths cinéticos de almacenamiento bajo certificadoSHA a7e0170f5c4b295a89ede835d7c234f4ebda8919f4271d2098b60fee638f937f. Ancla /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v6_20261011_parent_anchor; t=15.429999999999954ps; recibo /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v6_20261011/provenance/transfer_receipt.json SHA 68521bd13cb7a047025a870230d7c0aab96a6c24c3b06ce2e0f800cdee3b84e2. Cápsula común /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011.

```bash
set -o pipefail
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=4 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v6_20261011 --resume --cpu-start 0 --length-nm 320.0 --spacing-nm 5.0 --xi-cutoff 32.0 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 775.0 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.02 --ceiling-ps 200.0 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5.0 --tau-ee-Tc-ps 5.0 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15.0 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1.0 --checkpoint-wall-s 60.0 --step-tolerance 2e-07 --local-energy-tolerance 1e-05 --max-iterations 40 --max-rejections 50 --threads 4 2>&1 | tee /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v6_20261011/execution_v6.console.log
```

### Continuación v6 preparada 775_fine (2026-10-11)

Sin señales ni lanzamiento todavía. Misma física, dt/tolerancias/hilos/guarda/controlador/historia; únicos2paths cinéticos de almacenamiento bajo certificadoSHA a7e0170f5c4b295a89ede835d7c234f4ebda8919f4271d2098b60fee638f937f. Ancla /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v6_20261011_parent_anchor; t=7.859999999999918ps; recibo /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v6_20261011/provenance/transfer_receipt.json SHA 1babb9b458e7175688cde56b4209ff50a2d1b73d2265ff7c089dddfa249c019f. Cápsula común /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011.

```bash
set -o pipefail
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=2 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v6_20261011 --resume --cpu-start 12 --length-nm 320.0 --spacing-nm 3.5 --xi-cutoff 32.0 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 775.0 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.02 --ceiling-ps 200.0 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5.0 --tau-ee-Tc-ps 5.0 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15.0 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1.0 --checkpoint-wall-s 60.0 --step-tolerance 2e-07 --local-energy-tolerance 1e-05 --max-iterations 40 --max-rejections 50 --threads 2 2>&1 | tee /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v6_20261011/execution_v6.console.log
```

### Continuación v6 preparada 775_domain (2026-10-11)

Sin señales ni lanzamiento todavía. Misma física, dt/tolerancias/hilos/guarda/controlador/historia; únicos2paths cinéticos de almacenamiento bajo certificadoSHA a7e0170f5c4b295a89ede835d7c234f4ebda8919f4271d2098b60fee638f937f. Ancla /home/jdiaz/scratch/stage5_effective_L480_h5_775_adaptive_v6_control20_20261011_parent_anchor; t=12.359999999999822ps; recibo /home/jdiaz/scratch/stage5_effective_L480_h5_775_adaptive_v6_control20_20261011/provenance/transfer_receipt.json SHA d0d720c0864f1f332f92e753b9222bf720b6ee7f0b5ce5fe344b77a1f8da723e. Cápsula común /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011.

```bash
set -o pipefail
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=2 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L480_h5_775_adaptive_v6_control20_20261011 --resume --cpu-start 10 --length-nm 480.0 --spacing-nm 5.0 --xi-cutoff 32.0 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 775.0 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.02 --ceiling-ps 200.0 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5.0 --tau-ee-Tc-ps 5.0 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15.0 --current-uA 21.29720614489798 --bath-K 0.9 --duration-ps 20.0 --snapshot-every-ps 0.1 --checkpoint-every-ps 1.0 --checkpoint-wall-s 60.0 --step-tolerance 2e-07 --local-energy-tolerance 1e-05 --max-iterations 40 --max-rejections 50 --threads 2 2>&1 | tee /home/jdiaz/scratch/stage5_effective_L480_h5_775_adaptive_v6_control20_20261011/execution_v6.console.log
```

### Continuación v6 preparada 775_temporal (2026-10-11)

Sin señales ni lanzamiento todavía. Misma física, dt/tolerancias/hilos/guarda/controlador/historia; únicos2paths cinéticos de almacenamiento bajo certificadoSHA a7e0170f5c4b295a89ede835d7c234f4ebda8919f4271d2098b60fee638f937f. Ancla /home/jdiaz/scratch/stage5_effective_L320_h5_775_temporal_adaptive_v6_20261011_parent_anchor; t=13.919999999999984ps; recibo /home/jdiaz/scratch/stage5_effective_L320_h5_775_temporal_adaptive_v6_20261011/provenance/transfer_receipt.json SHA a512e11fe9f9da1d734a71250e7ddc03616afa9b7c1000ad186e5dde4941dd29. Cápsula común /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011.

```bash
set -o pipefail
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=2 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_temporal_adaptive_v6_20261011 --resume --cpu-start 8 --length-nm 320.0 --spacing-nm 5.0 --xi-cutoff 32.0 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 775.0 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.01 --ceiling-ps 200.0 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5.0 --tau-ee-Tc-ps 5.0 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15.0 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1.0 --checkpoint-wall-s 60.0 --step-tolerance 2e-07 --local-energy-tolerance 1.25e-06 --max-iterations 40 --max-rejections 50 --threads 2 2>&1 | tee /home/jdiaz/scratch/stage5_effective_L320_h5_775_temporal_adaptive_v6_20261011/execution_v6.console.log
```

### Continuación v6 preparada 1550_primary (2026-10-11)

Sin señales ni lanzamiento todavía. Misma física, dt/tolerancias/hilos/guarda/controlador/historia; únicos2paths cinéticos de almacenamiento bajo certificadoSHA a7e0170f5c4b295a89ede835d7c234f4ebda8919f4271d2098b60fee638f937f. Ancla /home/jdiaz/scratch/stage5_effective_L320_h5_1550_adaptive_v6_20261011_parent_anchor; t=47.299999999999656ps; recibo /home/jdiaz/scratch/stage5_effective_L320_h5_1550_adaptive_v6_20261011/provenance/transfer_receipt.json SHA 2f96e58233bbc237f9ebc2b10f395a3f08d1fc032b7d4a9aa8cb261a5651861c. Cápsula común /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011.

```bash
set -o pipefail
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=2 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_1550_adaptive_v6_20261011 --resume --cpu-start 4 --length-nm 320.0 --spacing-nm 5.0 --xi-cutoff 32.0 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 1550.0 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.1 --ceiling-ps 200.0 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5.0 --tau-ee-Tc-ps 5.0 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15.0 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1.0 --checkpoint-wall-s 60.0 --step-tolerance 2e-07 --local-energy-tolerance 1e-05 --max-iterations 40 --max-rejections 50 --threads 2 2>&1 | tee /home/jdiaz/scratch/stage5_effective_L320_h5_1550_adaptive_v6_20261011/execution_v6.console.log
```

### Corrección de carpeta de arranque espectral v6

El primer intento se detuvo antes de preparar DC o evolucionar: la carpeta de preparación ya existía. Se preservan wrapper y log originales. Se separa preparación de los nuevos datos:

```bash
bash /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_v6_control20_20261011/launch_spectral_control_v2.sh
```

Expandido:

```bash
set -o pipefail
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=2 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_v6_control20_data_20261011 --cpu-start 6 --threads 2 --length-nm 320 --spacing-nm 5 --xi-cutoff 64 --xi-order 8 --phonon-nodes 144 --transport-order 96 --wavelength-nm 775 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.02 --local-energy-tolerance 1e-5 --step-tolerance 2e-7 --duration-ps 20 --ceiling-ps 200 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5 --tau-ee-Tc-ps 5 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1 --checkpoint-wall-s 60 --max-iterations 40 --max-rejections 50 2>&1 | tee /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_v6_control20_20261011/execution_v2.console.log
```

Datos: /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_v6_control20_data_20261011

Las cinco continuaciones v6 fueron autorizadas por root, tras mostrar wrappers. Padres v5 SIGTERM intencional con checkpoints pre/finales/logs preservados (INTERRUPTED, sin fallo físico). Recibo de paradas: /home/jdiaz/scratch/stage5_effective_controls_20261011/five_v5_intentional_stops_preserved_20261011.json. Lanzamientos foreground verificados: /home/jdiaz/scratch/stage5_effective_controls_20261011/five_v6_launches_verified_20261011.json SHA 18dda084033389dd8e0c8cd5d2e18d2a06b8c73515a6acf9bec2cefbf88af27c.

- 775_primary: PID199034, start_ticks34764557, session63637, CPU[0, 1, 2, 3], NumBa4.

- 775_fine: PID199354, start_ticks34766220, session92229, CPU[12, 13], NumBa2.

- 775_domain: PID199521, start_ticks34767920, session80985, CPU[10, 11], NumBa2.

- 775_temporal: PID199584, start_ticks34769673, session33111, CPU[8, 9], NumBa2.

- 1550_primary: PID199693, start_ticks34771345, session52517, CPU[4, 5], NumBa2.

Primera auditoría de los cinco nuevos checkpoints v6 PASS: continuidad exacta de historia/snapshots/fuente/guarda; recibo all_five_v6_first_post_launch_checkpoints_20261011.json SHA8bf9cd55f0c01760d271712226d187e3da75b7f3a2860dd01787a6aa8fc8deb8. Las cinco sesiones v5 antiguas cerraron exit1 por KeyboardInterrupt SIGTERM intencional, con manifiestos INTERRUPTED coherentes; no fallo físico.

GEMINGA_COMMANDS_PRESENT.md actualizado a seis casos reales v6, incluyendo espectral PID198919/sesión50866/CPU6,7 y su directorio de datos separado. Backup /home/jdiaz/scratch/stage5_effective_controls_20261011/GEMINGA_COMMANDS_PRESENT_before_v6_active.md.


## 2026-10-11 — Sensibilidad exploratoria de anchura775nm, cápsulav6 (PREPARADO; NO EJECUTADO)

Finalidad: comparar sigma_ref/sqrt(2) y sigma_ref*sqrt(2), con sigma_ref=7.549834435nm; misma física y calibradores nominales. Es sensibilidad de hipótesis postcascada, no calibración. Duracióninicial20ps; no declarar censura200ps o margen de un cruce tardío sin completar su horizonte. Recursos:2hilosNumba+BLAS1 porcaso, CPU_START pendiente de coordinación; noreserva ni lanzamiento. Datos/cache/logsscratch.13fuentesverificadas al inicio, salida nueva no precreada. Fuente y recibo: /home/jdiaz/scratch/stage5_effective_width_controls_v6_prepare_20261011/preparation_receipt.json.

### narrow — sigma=5.33853912582nm

```bash
bash /home/jdiaz/scratch/stage5_effective_width_controls_v6_prepare_20261011/launch_narrow.sh CPU_START_A
```

Comando interno exacto (CPU_START lo asigna el coordinador):

```bash
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib TMPDIR=/home/jdiaz/scratch/tmp/jdiaz OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=2 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_sigma_narrow_v6_control20_data_20261011 --threads 2 --sigma-nm 5.3385391258242061 --length-nm 320 --spacing-nm 5 --xi-cutoff 32 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 775 --retained-fraction 0.667 --dt-ps 0.02 --local-energy-tolerance 1e-5 --step-tolerance 2e-7 --duration-ps 20 --ceiling-ps 200 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5 --tau-ee-Tc-ps 5 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1 --checkpoint-wall-s 60 --max-iterations 40 --max-rejections 50 --cpu-start "$CPU_START"
```

### wide — sigma=10.6770782516nm

```bash
bash /home/jdiaz/scratch/stage5_effective_width_controls_v6_prepare_20261011/launch_wide.sh CPU_START_B
```

Comando interno exacto (CPU_START lo asigna el coordinador):

```bash
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib TMPDIR=/home/jdiaz/scratch/tmp/jdiaz OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=2 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_sigma_wide_v6_control20_data_20261011 --threads 2 --sigma-nm 10.677078251648414 --length-nm 320 --spacing-nm 5 --xi-cutoff 32 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 775 --retained-fraction 0.667 --dt-ps 0.02 --local-energy-tolerance 1e-5 --step-tolerance 2e-7 --duration-ps 20 --ceiling-ps 200 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5 --tau-ee-Tc-ps 5 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1 --checkpoint-wall-s 60 --max-iterations 40 --max-rejections 50 --cpu-start "$CPU_START"
```



##2026-10-11 — Verificaciónfinal wheeloptimizado v6, offline y sin instalación

Finalidad: empaquetarworkingtree1.1.0 consólo fuentes/copiasbuild/extracciónscratch; verificarorígenes y SHAde2kernels aprobados, APIacelerada y fallbacksinNumba. Recursos:CPU4compartidotemporalmentecon1550,1hiloNumba/BLAS, timeout300s; nootrocore, sinred/pipinstall/evolución. Wheelsanterioresintactos. Recibo: /home/jdiaz/scratch/stage5_wheel_optimized_v6_20261011/wheel_verification_receipt.json.

```bash
timeout --signal=TERM 300s taskset -c 4 env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib TMPDIR=/home/jdiaz/scratch/tmp/jdiaz OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_wheel_optimized_v6_20261011/build_and_verify.py
```

Resultado de la verificación wheel optimizado v6: PASS, 12.019s, CPU4 liberado. Versión1.1.0; ambas API desde wheel extraído, con NumBa y fallback sin NumBa; orígenes y dos SHA certificados correctos. Sin instalación/red/evolución ni modificación del repo. Recibo /home/jdiaz/scratch/stage5_wheel_optimized_v6_20261011/wheel_verification_receipt.json SHA256 673993531fa95545ccb4049776130342f67a3a094b2ef713ecc1fbb12a103e1b. Wheel /home/jdiaz/scratch/stage5_wheel_optimized_v6_20261011/dist/pysnspd-1.1.0-py3-none-any.whl SHA256 a645ce8f53e0b8ebaa728ce53736002483283b570e15465a8f0dd61ccff70cbf.


## Plan preparado: ampliar controles dominio/espectro después de20ps —2026-10-11

Estado: NO capturado, NO transferido, NO lanzado. Se conserva la ventana declarada inicialmente 0–20 ps. Tras DEVELOPMENT_DURATION_REACHED/ACTIVE sin cruce, duration_ps pasa únicamente de 20 ps a None: continúan gatillo +0,1 mV, confirmación de 0,5 ps, margen de 5 ps y techo de 200 ps. No hay reset eléctrico total. Las 13 fuentes fijadas por hash, física, paso máximo de 0,02 ps, tolerancia NL 2e−7, guarda 1e−5 y 2 hilos permanecen idénticos.

### domain: preparación y lanzamiento futuro separado

Padre: /home/jdiaz/scratch/stage5_effective_L480_h5_775_adaptive_v6_control20_20261011

Ancla nueva: /home/jdiaz/scratch/stage5_effective_L480_h5_775_v6_control20_completed_anchor_20261011

Destino nuevo: /home/jdiaz/scratch/stage5_effective_L480_h5_775_control_trigger_v6_20261011

Recursos: 2 hilos Numba, BLAS/OpenMP 1; CPU 10,11, sustituyendo al control completado. Techo global 14/16 núcleos físicos; 14/15 quedan reservados. Reconstrucción/preparación con un hilo y timeout 180 s; continuación hasta evento+margen o censura a 200 ps, tiempo de pared dependiente de dinámica adaptativa.

Preparar SÓLO cuando el padre complete 20 ps y root revise la transición:

```bash
bash /home/jdiaz/scratch/stage5_control_window_extension_plan_20261011/prepare_domain_extension.sh
```

Después de revisar el recibo y mostrar el comando, lanzamiento foreground validado con pipefail y log preservado:

```bash
bash /home/jdiaz/scratch/stage5_control_window_extension_plan_20261011/launch_domain_extension.sh
```

### spectral: preparación y lanzamiento futuro separado

Padre: /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_v6_control20_data_20261011

Ancla nueva: /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_v6_control20_completed_anchor_20261011

Destino nuevo: /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_trigger_v6_20261011

Recursos: 2 hilos Numba, BLAS/OpenMP 1; CPU 6,7, sustituyendo al control completado. Techo global 14/16 núcleos físicos; 14/15 quedan reservados. Reconstrucción/preparación con un hilo y timeout 180 s; continuación hasta evento+margen o censura a 200 ps, tiempo de pared dependiente de dinámica adaptativa.

Preparar SÓLO cuando el padre complete 20 ps y root revise la transición:

```bash
bash /home/jdiaz/scratch/stage5_control_window_extension_plan_20261011/prepare_spectral_extension.sh
```

Después de revisar el recibo y mostrar el comando, lanzamiento foreground validado con pipefail y log preservado:

```bash
bash /home/jdiaz/scratch/stage5_control_window_extension_plan_20261011/launch_spectral_extension.sh
```

Plan completo/hashes: /home/jdiaz/scratch/stage5_control_window_extension_plan_20261011/preparation_receipt.json. El launcher lee resume_command.txt, valida igualdad exacta con el recibo/argv, cadena fresh/heredada, arrays/ledger/guard/recursos y 13 fuentes, y rechaza logs de ejecución ya existentes. No se han creado las anclas ni los destinos numéricos.


##2026-10-11 — Current-v6 DC baseline hold (bounded)

Fresh equilibrium check of current v6 after a historical hold used older source pins. No photon; no inherited state.10ps,100 steps at0.1ps;1Numba thread on physicalCPU4 shared with1550, BLAS1,14occupied physical cores globally, scratch only. Expected <=300s; timeout preserves checkpoint. Outputs manifest, traces/snapshots, complete state and preparation under /home/jdiaz/scratch/stage5_effective_DC_v6_10ps_20261011. Historical hold preserved and not relabeled as current.

~~~bash
bash /home/jdiaz/scratch/stage5_effective_DC_v6_10ps_prepare_20261011/launch.sh
~~~

Expanded command:

~~~bash
#!/usr/bin/env bash
set -euo pipefail
cd /home/jdiaz/pysnspd
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache
export NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache
export MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
timeout --signal=TERM --kill-after=30s 300s /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011/pipelines/05_effective_photon_run.py \
 --output-dir /home/jdiaz/scratch/stage5_effective_DC_v6_10ps_20261011 \
 --dc-only --duration-ps 10 --wavelength-nm 775 --length-nm 320 --spacing-nm 5 \
 --xi-cutoff 32 --xi-order 6 --phonon-nodes 96 --transport-order 48 \
 --sigma-nm 7.549834435 --retained-fraction 0.667 \
 --current-uA 21.29720614489798 --bath-K 0.9 \
 --tau-ee-Tc-ps 5 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15 \
 --dt-ps 0.1 --step-tolerance 2e-7 --local-energy-tolerance 1e-5 --max-iterations 40 \
 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5 --ceiling-ps 200 \
 --snapshot-every-ps 1 --checkpoint-every-ps 1 --checkpoint-wall-s 60 \
 --threads 1 --cpu-start 4 2>&1 | tee /home/jdiaz/scratch/stage5_effective_DC_v6_10ps_prepare_20261011/execution.console.log
~~~


### Result — current-v6 DC baseline,2026-10-11

Completed10ps/100steps in135.97s with zero rejections. Independent saved-data audit PASS; Vout0V, max|Vdev|1.4558046545474175e-17V; no saved gap/temperature/current drift; max energy defect2.3700227555395377e-41J. All13 runtime pins match currentv6. Evidence: /home/jdiaz/scratch/stage5_effective_DC_v6_10ps_20261011/current_v6_DC_evidence.json (SHA256 f4c26bb1b75abb61baa47227ec0c195fc57b9ebf2a7330e1114c0df2cdd6d721). Independent audit: /home/jdiaz/scratch/stage5_effective_DC_v6_10ps_audit_20261011/audit.json (SHA256 ba2cb6ee065ed86ad2752bed33d0a6b9fc6b08d1339611ae214a2597330b6a7d). This is basal validation only, not stage5global admission. Historical hold preserved separately.


##2026-10-11 — Variantes sigma775 compartidas1hilo (PREPARADAS; NO LANZADAS)

Conserva wrappers originales2hilos en /home/jdiaz/scratch/stage5_effective_width_controls_v6_prepare_20261011. Sólo cambia asignación runtime: narrowNumba1CPU5 compartida con1550(4/5), wideNumba1CPU7 compartida conespectral(6/7).13pins iguales, parámetros/cap0.02ps/guarda1e-5/NL2e-7/ventanainicial20ps intactos. No crea salidas de runner. La sensibilidad no es calibración. Presupuesto verificado:6corridas/14hilos vigentes; propuesto16hilos sobre unión14/16núcleos físicos (87.5%), reserva completa14/15 y hermanos30/31,16<28lógicosadmitidos, BLAS1. Memoria disponible 112.46GiB; reserva conservadora adicional4GiB. Hay contención explícitaCPU5/7; no se afirma mejora sin medir y sigue siendo cálculo largo. Root decide y revalida presupuesto antesde lanzar.

Recibo /home/jdiaz/scratch/stage5_effective_width_controls_v6_shared1_prepare_20261011/preparation_receipt.json SHA256 332a5a913ce2ea99f190b999a25c8e92e5a9e31aa477f2ecf6220cafb1711c26.

### narrow

```bash
bash /home/jdiaz/scratch/stage5_effective_width_controls_v6_shared1_prepare_20261011/launch_narrow_shared1.sh
```

Comando interno exacto:

```bash
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib TMPDIR=/home/jdiaz/scratch/tmp/jdiaz OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_sigma_narrow_v6_control20_data_20261011 --threads 1 --sigma-nm 5.3385391258242061 --length-nm 320 --spacing-nm 5 --xi-cutoff 32 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 775 --retained-fraction 0.667 --dt-ps 0.02 --local-energy-tolerance 1e-5 --step-tolerance 2e-7 --duration-ps 20 --ceiling-ps 200 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5 --tau-ee-Tc-ps 5 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1 --checkpoint-wall-s 60 --max-iterations 40 --max-rejections 50 --cpu-start 5
```

### wide

```bash
bash /home/jdiaz/scratch/stage5_effective_width_controls_v6_shared1_prepare_20261011/launch_wide_shared1.sh
```

Comando interno exacto:

```bash
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib TMPDIR=/home/jdiaz/scratch/tmp/jdiaz OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_effective_runtime_v6_20261011/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_sigma_wide_v6_control20_data_20261011 --threads 1 --sigma-nm 10.677078251648414 --length-nm 320 --spacing-nm 5 --xi-cutoff 32 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 775 --retained-fraction 0.667 --dt-ps 0.02 --local-energy-tolerance 1e-5 --step-tolerance 2e-7 --duration-ps 20 --ceiling-ps 200 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5 --tau-ee-Tc-ps 5 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1 --checkpoint-wall-s 60 --max-iterations 40 --max-rejections 50 --cpu-start 7
```



Resultado de activación de sensibilidades compartidas: STARTED0 confirmado para estrecha PID212853/sesión55315/CPU5 y ancha PID212907/sesión35086/CPU7, un hilo cada una. Misma cápsula v6 archivada con13pins+173transitivos; recibo /home/jdiaz/scratch/stage5_effective_controls_20261011/two_width_v6_runtime_archives_20261011.json SHA3e224020229fc3a20514bd3272200ba27e174cf392188aa2237d1cf9eecbd408. PRESENT actualizado a8trayectorias activas:16hilos Numba sobre unión14/16núcleos físicos, BLAS/OpenMP1; CPU14/15 libres. Verificación /home/jdiaz/scratch/stage5_effective_controls_20261011/eight_active_v6_status_20261011.json. Sin señales, cambios de RHS ni relanzamientos durante archivo/verificación.


## 2026-10-11 — precisión física puntual775 en captura v6, Num1 compartidaCPU4

Finalidad: desde una captura atómica del checkpoint actual, comparar un paso h vigente y dos h/2; NL2e-7, máximo40iter, sin reintentos, guarda, reinyección ni cambios de la corrida. Evidencia local, sin admisión global. Límite global300s; Num1/CPU4 compartida1550, BLAS1; reserva física global14/15 intacta; memoria adicional estimada2GiB. Todos los datos en scratch.

```bash
taskset -c 4 env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python /home/jdiaz/scratch/stage5_nominal775_single_step_precision_20261011/bounded_controller.py
```

Salidas: /home/jdiaz/scratch/stage5_nominal775_single_step_precision_20261011/{recipe,receipt,execution_receipt}.json, anchor, estados y diagnostic.log. El recibo distingue fallo/timeout incompleto; no se relanza.

Resultado del diagnóstico puntual anterior: COMPLETED57.27s, checkpoint22.220000ps, h=.005 frente2×.0025. Recibo físico SHAe0cc5fe2b09c211d50c44214aea65b6a081a32d36759f2bba0beabcb8c6cc0e1. ΔVout1.199e−10V (0.000405% señal real), ΔVdev1.873e−6V (0.01733%), L2 magnitud/DC central1.618e−6. Errorlocal9.4389e−6 grueso frente1.1707/1.1742e−6 por subpaso fino. Soporte/13pins íntegros. Evidencia local, sin admisión global; corrida y política intactas. Coste secuencial con compilación inicial/contención: no inferir speedup de los tiempos medidos.


## Contingencia guard2 nominal775 — preparada, NO ejecutada (2026-10-11)

La decisión prospectiva de root precedió la parada intencional. Original final coherente: t*=22.43500000000006ps, global0.6634625866%, ACTIVE; todos los archivos y logs originales quedan íntegros. Nueva guarda nominal2.5e-6; temporal3.125e-7 (razón8). Los cuatro controles775 deben alcanzar la banda[22.18500000000006,22.43500000000006]ps antes de sus transferencias; la diferencia acotada se declara y no acredita aislamiento exacto de eje. Sin cambio RHS/13pins, circuito, depósito, tiempo máximo200ps o gatillo+5ps.

Preparación futura, sólo tras revisión de root:
```bash
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/prepare_nominal_guard2.sh
```

Lanzamiento futuro foreground, sólo tras revisión de recibos reales y comando mostrado en chat:
```bash
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/launch_nominal_guard2.sh
```

Expansión prevista que el launcher exige idéntica a resume_command.txt generado por el helper:
```bash
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v6_guard2_20261011/runtime_sources/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v6_guard2_20261011 --resume --cpu-start 0 --length-nm 320.0 --spacing-nm 5.0 --xi-cutoff 32.0 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 775.0 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.02 --ceiling-ps 200.0 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5.0 --tau-ee-Tc-ps 5.0 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15.0 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1.0 --checkpoint-wall-s 60.0 --step-tolerance 2e-07 --local-energy-tolerance 2.5e-06 --max-iterations 40 --max-rejections 50 --threads 4
```

Recursos:4hilos Numba,CPU0–3;BLAS/OpenMP1. Salida futura:/home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v6_guard2_20261011;log execution.console.log. Tiempo físico pendiente hasta gatillo positivo0.1mV confirmado0.5ps más5ps, o censura200ps; no reset completo. Preparación/validación auxiliares ≤180/120s, un núcleo; simulación foreground puede exceder5min y está autorizada para esta tarea. Helper52tests PASS;launcher6tests sintéticos PASS;AST/bash -n PASS. Recibo de plan:/home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/preparation_review_receipt.json SHA8e8902cedbdf9f0ebbced1296ef465e9e3b9c0e14558b516f349e7d916ebf085. No ancla ni salida numérica futura creadas, no transferencia ni lanzamiento ejecutados.


## Guard2b — corrección de trazabilidad, pendiente de ejecución (2026-10-11)

Root ejecutó la preparación guard2 anterior sin evolución y sin lanzamiento. Todos27arrays,13pins y ledger eran exactos; el nuevo segmento heredaba incorrectamente el campo descriptivo guard_active_since_time_ps=12.4. Se conserva íntegra esa preparación como antecedente NO LANZADO y el helper94006d. Guard2b usa auxiliar nuevo SHA567d90a84c9c89083b37d759e8a07ee04377c69d57810cc19f402850ddf3de41 y asigna22.435ps únicamente al nuevo segmento; la guarda física/numérica2.5e-6 y el resto del modelo no cambian.53tests del auxiliar y6tests del wrapper PASS. Incluye copias iniciales inmutables de checkpoint y manifest para auditar una corrida avanzada.

Wrappers futuros, sólo después de revisión root:
```bash
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/prepare_nominal_guard2b.sh
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/launch_nominal_guard2b.sh
```

Expansión prevista, que debe coincidir con el resume_command generado:
```bash
env PYTHONDONTWRITEBYTECODE=1 PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v6_guard2b_20261011/runtime_sources/pipelines/05_effective_photon_run.py --output-dir /home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v6_guard2b_20261011 --resume --cpu-start 0 --length-nm 320.0 --spacing-nm 5.0 --xi-cutoff 32.0 --xi-order 6 --phonon-nodes 96 --transport-order 48 --wavelength-nm 775.0 --sigma-nm 7.549834435 --retained-fraction 0.667 --dt-ps 0.02 --ceiling-ps 200.0 --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5.0 --tau-ee-Tc-ps 5.0 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15.0 --current-uA 21.29720614489798 --bath-K 0.9 --snapshot-every-ps 0.1 --checkpoint-every-ps 1.0 --checkpoint-wall-s 60.0 --step-tolerance 2e-07 --local-energy-tolerance 2.5e-06 --max-iterations 40 --max-rejections 50 --threads 4
```

Destino nuevo:/home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v6_guard2b_20261011. Ancla nueva:/home/jdiaz/scratch/stage5_effective_L320_h5_775_adaptive_v6_guard2b_anchor_20261011. Recursos4hilosCPU0–3,BLAS1;foreground con execution.console.log y pipefail. Gatillo/extra/censura originales conservados. Recibo de revisión:/home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/guard2b_preparation_review_receipt.json SHA60db2ea220979ccfc420af4ce317e89d5133e7396407f5e21b2663de73c26a35. Ninguna captura, transferencia o corrida guard2b ejecutada por este plan todavía.


## Nominal guard2b ACTIVO — continuidad real verificada (2026-10-11)

Root lanzó la nueva continuación después de revisar53tests,27arrays y13pins. PID222187,sesiónROOT47902,CPU0–3/Numba4,STARTED22.43500000000006ps. Guarda2.5e-6 sólo desde esa ancla;E0,Σtrabajo firmado/absoluto,Σresiduos,fuente única,controlador e historia completos conservados. TransferSHA9b1c375520fb52718a98db781461087149a16c6b9c07b6ed0e41aeac1a6409c0. RootapprovalhelperSHAca185b92256360d92275681e008fb2e36e518bc796d929cfc3de748aa22c2db5;rootreviewtransferSHAfe9a8834498043dd4d43382a897a1f739e56cf3777ad355a140b818e3db3f62a. Ocho procesos activos y16threads sobre14núcleos físicos verificados;reservados14/15. Recibo:/home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/eight_active_guard2b_status_receipt.json SHAd8824f7508176a81d03ac2b00c551bbf55e9a4cbca94d716e2d5921d364aebee. Guard2 previo NO LANZADO preservado. Las cuatro transferencias de controles esperan su banda física22.185–22.435ps; no fueron ejecutadas. PRESENT actualizado con el directorio/log vigente; copia anterior:/home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/GEMINGA_COMMANDS_PRESENT_before_nominal_guard2b_active.md.


## Plan futuro mínimo: cuatro controles 775 nm hacia guard2b — NO EJECUTADO

Fecha UTC: 2026-10-11T17:09:14.543733+00:00. Nominal guard2b ya activo desde t*=22,43500000000006 ps. Los controles conservarán sus guardas originales hasta una ancla física en [22,18500000000006; 22,43500000000006] ps, nunca posterior. Temporal pasará 1,25e-6→3,125e-7; malla/dominio/espectro pasarán 1e-5→2,5e-6. La razón nominal/temporal es 8 fuera de la banda; los prefijos y la diferencia transitoria quedan explícitos, sin afirmar aislamiento exacto del eje.

Estos comandos son sólo PREPARACIÓN FUTURA, no autorizan detener ni lanzar. Root deberá revisar el final INTERRUPTED coherente dentro de la banda. Dominio/espectro deben completar primero 20 ps y ampliarse por el auxiliar de ventana aprobado. Cada preparador valida los mismos 13 pins, parámetros, fuente, arrays, controlador, historia, contadores, ledger global y recursos; preserva el checkpoint y manifest iniciales inmutables. Si la ancla queda fuera de la banda, rechaza la operación. CPUs respectivas 8/9, 12/13, 10/11 y 6/7; preparación de un hilo, timeout 180 s. No hay captura, transferencia, señal ni ejecución realizadas por este plan.

```bash
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/prepare_temporal_guard2b.sh
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/prepare_fine_guard2b.sh
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/prepare_domain_guard2b.sh
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/prepare_spectral_guard2b.sh
```

Plan: /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/plan.json; SHA b295aaa4918fe76f30f17783e328c3830eda5e90189c6142a86b8e20a398ed1a. Recibo: /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/preparation_review_receipt.json; SHA 5640454ccb62694db33f8536418766dfeac8769266ca07b53883a224b0959d02. Pruebas puras: 23 PASS en 0,06 s; AST y bash -n PASS. Los ocho futuros directorios de ancla/salida siguen ausentes. Los launchers sólo se prepararán y mostrarán después de revisar los recibos reales, conservando los originales.

## 2026-10-11 — Política espectral prospectiva, regresiones de análisis

Finalidad: conservar el fallo nodal original del 2 % y verificar la capa adoptada de 3 % nodal más ambos observables espaciales al 2 %, con todos los demás controles intactos. Sin evolución física ni admisión. CPU4 compartida, un proceso, BLAS1, límite60s. Quince pruebas con subcasos: PASS0,012s. Recibo: /home/jdiaz/scratch/stage5_spectral_policy_candidate_20261011/candidate_receipt.json; SHA256 f73de50fe427044a4fe6fcbfedada1054eedd78602a7adb99cc06aafed7005d5.

```bash
env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 60s taskset -c 4 /home/jdiaz/.conda/envs/snspd/bin/python /home/jdiaz/scratch/stage5_spectral_policy_candidate_20261011/test_spectral_observable_policy.py
```


## 2026-10-11 — lanzadores futuros de los cuatro controles guard2b (NO EJECUTADOS)

Preparación documental para revisión de raíz. No detiene, captura, transfiere ni lanza corridas. Sólo después de una parada coherente autorizada dentro de [22.18500000000006,22.43500000000006] ps, preparación real del caso y revisión de sus recibos puede ejecutarse su lanzador. Domain y spectral proceden de las extensiones explícitas a duration=None; no de los pilotos control20. Conserva 13 fuentes físicas, prefijos completos, reloj, energía inicial, trabajo firmado/absoluto, contadores y circuito. Nueva guarda nominal/espacial/espectral2.5e-6; temporal3.125e-7, NL2e-7, cap.02/.01ps,2hilos por control y afinidades existentes. No certifica etapa5.

Plan: /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/plan.json

Validador: /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/validate_control_guard2b_launch.py SHA256 9f47936c3370ae1054674a8baf2154cd68846e9f336e7b575dc829a1ef254577

Cada lanzador verifica SHA del validador; éste verifica plan/prepare/helper y todas sus dependencias, recibos reales, banda exacta, archivos iniciales inmutables, arrays/historia/ledger,13fuentes y argv. Rechaza destinos faltantes, avanzados o procesos existentes. Sólo imprime comando tras validación real. Foreground, pipefail, tee sin sobrescribir log y cachés en scratch. Recursos estimados:2núcleos/hilos por control, memoria heredada con presupuesto global≤90%; duración depende del evento+5ps/techo200ps y puede exceder5min. Root revisa/coordina antes de cada lanzamiento.

### temporal

Parent: /home/jdiaz/scratch/stage5_effective_L320_h5_775_temporal_adaptive_v6_20261011

Output futuro: /home/jdiaz/scratch/stage5_effective_L320_h5_775_temporal_adaptive_v6_guard2b_20261011

Afinidad: 8,9;2hilos.

~~~bash
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/launch_temporal_guard2b.sh
~~~

### fine

Parent: /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v6_20261011

Output futuro: /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v6_guard2b_20261011

Afinidad: 12,13;2hilos.

~~~bash
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/launch_fine_guard2b.sh
~~~

### domain

Parent: /home/jdiaz/scratch/stage5_effective_L480_h5_775_control_trigger_v6_20261011

Output futuro: /home/jdiaz/scratch/stage5_effective_L480_h5_775_control_trigger_v6_guard2b_20261011

Afinidad: 10,11;2hilos.

~~~bash
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/launch_domain_guard2b.sh
~~~

### spectral

Parent: /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_trigger_v6_20261011

Output futuro: /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_trigger_v6_guard2b_20261011

Afinidad: 6,7;2hilos.

~~~bash
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/launch_spectral_guard2b.sh
~~~

Pruebas sintéticas:41PASS; cuatro bash -n PASS. Ningún destino de corrida creado ni evolución ejecutada. Los recibos y fuentes se encuentran junto a plan.json.

## 2026-10-11 — Auditoría futura de sensibilidad física de anchura, sólo preparada

No ejecutado. Primer comando exige que narrow ywide terminen realmente20ps; audita energía1%, soporte y eventos originales sincomparación de convergencia. Segundo comando exige además nominal realmentecompleto yusa protocolo histórico en0--20ps conservando toda historia posterior. No asigna latencia ficticia nicensura200ps a20ps. CPU4compartida,1proceso/BLAS1,timeout300s; previsión2GiB. Salidas nuevas sóloal ejecutarse: /home/jdiaz/scratch/stage5_width_independent_audit_completed20_20261011. Recibo: /home/jdiaz/scratch/stage5_width_independent_audit_prepare_20261011/preparation_receipt.json SHA256 ec68bcd11aa5d7c95d4d9c3f8b4af754fde3fe7904f6cd037598a11c1ec6bdb4.

audit_widths_after20.sh

```bash
bash /home/jdiaz/scratch/stage5_width_independent_audit_prepare_20261011/audit_widths_after20.sh
```

analyze_sensitivity_after_nominal.sh

```bash
bash /home/jdiaz/scratch/stage5_width_independent_audit_prepare_20261011/analyze_sensitivity_after_nominal.sh
```



## Lectura terminal de dominio y revisión temporal guard2b — 18:05 UTC

Temporal original se detuvo intencionalmente dentro de banda y quedó INTERRUPTED a 22,252500000000108 ps. Revisión independiente de la preparación guard2b: 27 arrays exactos, metadata completa salvo campos prospectivos declarados, mismos 13 pins, ledger global intacto. Recibo temporal_independent_transfer_review.json en /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan; SHA 4bfac4523506967675c3d5b262b906d279c66e9bafbef9b2e3968e0580fff82c. El dominio original completó su ventana de 20 ps de forma natural, ACTIVE y sin cruce, sesión80985 exit0. Auditoría completa PASS sin modificación de sus productos originales: /home/jdiaz/scratch/stage5_domain_control20_terminal_audit_20261011/terminal_receipt.json; SHA c1eb2c581dadcf7ed935a625924fef15a806e3ec46986cab1f0228bdf3c87c8a. Este diagnóstico no capturó, extendió ni lanzó ninguna continuación. No se infiere etapa5 ni latencia por finalizar un corte de desarrollo.


## 2026-10-11 — corrección de arranque de los cuatro launchers guard2b

La primera llamada temporal con taskset externo8,9 abortó en resource_budget antes de construir sistema/evolucionar. Checkpoint y manifest permanecen byteexactos a sus iniciales; el log se preserva en execution_failed_outer_affinity.console.log. Recibo: /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/outer_affinity_failure_preservation_receipt.json SHA256 f42391b63302be1323749f778c2efafafb2ed98bf4df3a39919d071ae01caad7. La versión anterior de las fuentes se conserva en history_before_remove_outer_runner_affinity_20261011; no usar ese comando impreso con taskset externo.

Los launchers actuales conservan la afinidad anfitrión del runner: primero su constructor admite el presupuesto global90%, después --cpu-start selecciona exactamente2núcleos físicos. El validador de lectura sigue en1CPU. No cambiar la afinidad externa del runner a sus2CPUs antes del constructor. Constructor real verificado sin runner: temporal[8,9],fine[12,13],domain[10,11],spectral[6,7].42testsPASS+4bash-nPASS. Recibo fuentes corregidas: /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/corrected_host_affinity_launchers_receipt.json SHA256 63558e890c63c3fd208c4535be5c6bb574e254ac71f283cfc68db07705c38504. Nueva aprobación raíz de fuentes y transferencia real sigue requerida antes de lanzar cada caso; este registro no ejecuta nada.

~~~bash
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/launch_temporal_guard2b.sh
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/launch_fine_guard2b.sh
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/launch_domain_guard2b.sh
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/launch_spectral_guard2b.sh
~~~


## 2026-10-11 — CANDIDATO condicional resource-only del fino, Num2→Num4

Fuente revisable: /home/jdiaz/scratch/stage5_fine_resource_transition_candidate_20261011. Ninguna corrida preparada, detenida o lanzada. Las fuentes aún NO están promovidas.

# Candidato condicional: fino mesh2→4hilos

**NO EJECUTADO, NO APLICADO.** Requiere revisión de root, término real del nominal775 y liberación de los núcleos0–3. El fino continúa con2hilos hasta decisión explícita. No hay señal de parada en estos comandos.

La preparación sólo se permite después de que root haya detenido coherentemente el fino guard2b en estado INTERRUPTED/ACTIVE pretrigger y haya verificado el recibo de guarda real. Un evento pendiente o previo excluye esta transición. Se conservan los27arrays, toda metadata salvo threads, el siguiente paso, los contadores, la ventana de guarda y los ledgers desde t=0.

El protocolo indicado en APPROVED_FINE_RESOURCE_GUARD_PROTOCOL debe contener la entrada mesh real del padre guard2b. No se usa el protocolo nominal-only de las anchuras. El destino sigue sin existir hasta ejecutar expresamente la preparación aprobada.

```bash
export TMPDIR=/home/jdiaz/scratch/tmp/pysnspd
export PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache
export NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache
export MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/mpl
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH=/home/jdiaz/pysnspd
cd /home/jdiaz/pysnspd
: "${APPROVED_FINE_RESOURCE_GUARD_PROTOCOL:?root debe revisar el protocolo con mesh real}"
/home/jdiaz/.conda/envs/snspd/bin/python /home/jdiaz/scratch/stage5_fine_resource_transition_candidate_20261011/fine_resource_transition.py preflight
taskset -c4 /home/jdiaz/.conda/envs/snspd/bin/python /home/jdiaz/scratch/stage5_effective_controls_20261011/tighten_local_energy_guard_v6_guard2b.py capture --run /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v6_guard2b_20261011 --anchor /home/jdiaz/scratch/stage5_fine_resources4_future_anchor_20261011
/home/jdiaz/.conda/envs/snspd/bin/python /home/jdiaz/scratch/stage5_fine_resource_transition_candidate_20261011/fine_resource_transition.py prepare --anchor /home/jdiaz/scratch/stage5_fine_resources4_future_anchor_20261011 --output /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v6_guard2b_resources4_20261011 --guard-protocol "$APPROVED_FINE_RESOURCE_GUARD_PROTOCOL"
```

Parar aquí para revisar la transferencia real y sus SHA. Sólo después de esa revisión, el launcher exige APPROVED_RESOURCE_TRANSFER_SHA256; revalida el estado inicial y vuelve a verificar que nominal/parent ya no tienen proceso, CPU0–3 estén libres, la ocupación de pipelines quede≤14de16núcleos físicos y14/15permanezcan reservados.

```bash
: "${APPROVED_RESOURCE_TRANSFER_SHA256:?SHA del recibo real revisado por root}"
bash /home/jdiaz/scratch/stage5_fine_resource_transition_candidate_20261011/launch_resources4_REVIEW_ONLY.sh
```

El inventario automático identifica pipelines de jdiaz en ese instante; no reserva núcleos ni inventaría cualquier proceso del servidor. Root debe confirmar el presupuesto global al autorizar el launcher.

El runner se lanza foreground con tee/pipefail y checkpoint. No se recorta su afinidad antes de resource_budget; usa --threads4 --cpu-start0 y fija CPU0–3 internamente. No se promete identidad bit a bit de una evolución futura entre hilos. El certificado a7e0170 acredita la optimización y los casos comprobados, con pequeños redondeos entre algunas reducciones escalares.

Auditoría futura: añadir únicamente --resource-transition-receipt <hijo>/provenance/transfer_receipt.json a la comparación kind=mesh, después de promover expresamente el selector candidato. Padre guard2b completo e hijo completo se verifican antes de cualquier ventana finita; Vout/latencia, L2/D, soporte, energía global1%, guarda2.5e−6 y NL2e−7 conservan sus criterios.



## 2026-10-11 — FE DE ERRATAS: candidato resource-only, esquema real28keys

La propuesta candidata anterior0b345743/patch2f15 NO se aplica: usó incorrectamente27keys totales; el productor real escribe27arrays numéricos+metadata_json=28keys. Los tests anteriores eran sintéticos y no contrastaban ese esquema. Ninguna corrida se preparó, detuvo o lanzó. El nuevo candidato exige nombres exactos desde HISTORY_KEYS/SNAPSHOT_KEYS del runner pinned, contrastados con h.decode y un NPZ inicial inmutable. Se corrige además taskset -c4 a taskset -c 4 y el launcher verifica SHA/archivo/segmentos ancestrales antes de exec.

Fuente única revisable: /home/jdiaz/scratch/stage5_fine_resource_transition_candidate_20261011/candidate_receipt.json. Patchfinal SHA e423789c54de15cf0f98cc75c3aeecfa36039c3da973ef1663e93e98d8391c40. Sigue bloqueado hasta revisión root y transferencia real aprobada.

# Candidato condicional: fino mesh2→4hilos

**NO EJECUTADO, NO APLICADO.** Requiere revisión de root, término real del nominal775 y liberación de los núcleos0–3. El fino continúa con2hilos hasta decisión explícita. No hay señal de parada en estos comandos.

La preparación sólo se permite después de que root haya detenido coherentemente el fino guard2b en estado INTERRUPTED/ACTIVE pretrigger y haya verificado el recibo de guarda real. Un evento pendiente o previo excluye esta transición. Se conservan los27 arrays numéricos + metadata_json (28 claves), toda metadata salvo threads, el siguiente paso, los contadores, la ventana de guarda y los ledgers desde t=0.

El protocolo indicado en APPROVED_FINE_RESOURCE_GUARD_PROTOCOL debe contener la entrada mesh real del padre guard2b. No se usa el protocolo nominal-only de las anchuras. El destino sigue sin existir hasta ejecutar expresamente la preparación aprobada.

```bash
export TMPDIR=/home/jdiaz/scratch/tmp/pysnspd
export PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache
export NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache
export MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/mpl
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH=/home/jdiaz/pysnspd
cd /home/jdiaz/pysnspd
: "${APPROVED_FINE_RESOURCE_GUARD_PROTOCOL:?root debe revisar el protocolo con mesh real}"
/home/jdiaz/.conda/envs/snspd/bin/python /home/jdiaz/scratch/stage5_fine_resource_transition_candidate_20261011/fine_resource_transition.py preflight
taskset -c 4 /home/jdiaz/.conda/envs/snspd/bin/python /home/jdiaz/scratch/stage5_effective_controls_20261011/tighten_local_energy_guard_v6_guard2b.py capture --run /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v6_guard2b_20261011 --anchor /home/jdiaz/scratch/stage5_fine_resources4_future_anchor_20261011
/home/jdiaz/.conda/envs/snspd/bin/python /home/jdiaz/scratch/stage5_fine_resource_transition_candidate_20261011/fine_resource_transition.py prepare --anchor /home/jdiaz/scratch/stage5_fine_resources4_future_anchor_20261011 --output /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v6_guard2b_resources4_20261011 --guard-protocol "$APPROVED_FINE_RESOURCE_GUARD_PROTOCOL"
```

Parar aquí para revisar la transferencia real y sus SHA. Sólo después de esa revisión, el launcher exige APPROVED_RESOURCE_TRANSFER_SHA256; revalida el estado inicial y vuelve a verificar que nominal/parent ya no tienen proceso, CPU0–3 estén libres, la ocupación de pipelines quede≤14de16núcleos físicos y14/15permanezcan reservados.

```bash
: "${APPROVED_RESOURCE_TRANSFER_SHA256:?SHA del recibo real revisado por root}"
bash /home/jdiaz/scratch/stage5_fine_resource_transition_candidate_20261011/launch_resources4_REVIEW_ONLY.sh
```

El inventario automático identifica pipelines de jdiaz en ese instante; no reserva núcleos ni inventaría cualquier proceso del servidor. Root debe confirmar el presupuesto global al autorizar el launcher.

El runner se lanza foreground con tee/pipefail y checkpoint. No se recorta su afinidad antes de resource_budget; usa --threads4 --cpu-start0 y fija CPU0–3 internamente. No se promete identidad bit a bit de una evolución futura entre hilos. El certificado a7e0170 acredita la optimización y los casos comprobados, con pequeños redondeos entre algunas reducciones escalares.

Auditoría futura: añadir únicamente --resource-transition-receipt <hijo>/provenance/transfer_receipt.json a la comparación kind=mesh, después de promover expresamente el selector candidato. Padre guard2b completo e hijo completo se verifican antes de cualquier ventana finita; Vout/latencia, L2/D, soporte, energía global1%, guarda2.5e−6 y NL2e−7 conservan sus criterios.



## 2026-10-11 — Candidato de recursos cerrado: dominio/espectro (NO ejecutado)

Estado: REVIEW_ONLY. Fuentes y launchers en `/home/jdiaz/scratch/stage5_control_resource_transition_candidate_20261011`. El puente fino ya promovido no cambia. Este bloque NO autoriza por sí solo una transferencia real ni acredita el modelo. Root debe revisar las fuentes, detener coherentemente fuera de estos scripts, verificar el padre final todavía ACTIVE/pretrigger y aprobar cada recibo real antes de lanzar.

Casos únicos: dominio Num2→4 CPU10..13 tras traslado real del fino a0..3; espectro Num2→4 CPU4..7 tras1550 y wide realmente terminales; espectro Num4→8 CPU4..11 tras temporal+dominio realmente terminales. Reservar físicos14/15 y revalidar inventario GLOBAL≤90%; el inventario automático sólo identifica pipelines de jdiaz y no reserva CPUs. Si el padre ya tiene cualquier evento/pending, estos candidatos se bloquean y conserva sus recursos.

Los28keys del productor real (27arrays numéricos+metadata) conservan estado, historia/snapshots, E0, contador, next_h, origen global y toda la ventana de guarda; sólo cambia parameter_identity.threads. NL2e-7/guarda2.5e-6/cap.02ps, fuentes13, campos físicos, puertos, gatillo+5ps y techo200ps permanecen iguales. No se promete identidad futura entre hilos. Todos los temporales/cachés/logs siguen en scratch; no usar taskset exterior antes del runner.

Preparación FUTURA: el helper exige una captura inmutable ya revisada de un padre estable INTERRUPTED, sin dueño de proceso. Los destinos reales no existen por este plan. Reemplazar las variables exclusivamente por evidencias reales revisadas; el caso8 además exige certificadoNum8 independiente y su aprobaciónSHA. Comandos condicionales de preparación (ninguno ejecutado):

```bash
cd /home/jdiaz/pysnspd
export TMPDIR=/home/jdiaz/scratch/tmp/pysnspd PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python /home/jdiaz/scratch/stage5_control_resource_transition_candidate_20261011/proposed/sandbox/stage5_effective_validation/control_resource_transition.py --case domain_2to4 prepare --anchor "$ANCHOR_DOMAIN_REVIEWED" --output /home/jdiaz/scratch/stage5_effective_L480_h5_775_control_trigger_v6_guard2b_resources4_20261011 --guard-protocol "$GUARD_PROTOCOL_REVIEWED"
/home/jdiaz/.conda/envs/snspd/bin/python /home/jdiaz/scratch/stage5_control_resource_transition_candidate_20261011/proposed/sandbox/stage5_effective_validation/control_resource_transition.py --case spectral_2to4 prepare --anchor "$ANCHOR_SPECTRAL4_REVIEWED" --output /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_trigger_v6_guard2b_resources4_20261011 --guard-protocol "$GUARD_PROTOCOL_REVIEWED"
/home/jdiaz/.conda/envs/snspd/bin/python /home/jdiaz/scratch/stage5_control_resource_transition_candidate_20261011/proposed/sandbox/stage5_effective_validation/control_resource_transition.py --case spectral_4to8 prepare --anchor "$ANCHOR_SPECTRAL8_REVIEWED" --output /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_trigger_v6_guard2b_resources8_20261011 --guard-protocol "$GUARD_PROTOCOL_REVIEWED" --prior-resource-receipt /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_trigger_v6_guard2b_resources4_20261011/provenance/transfer_receipt.json --num8-evidence "$NUM8_CERTIFICATE_REVIEWED" --num8-approval "$ROOT_NUM8_CERTIFICATE_APPROVAL"
```

Lanzamientos FUTUROS foreground, pipefail y tee: cada launcher vuelve a verificar el recibo original/linaje/arrays/guarda, inicial aún inmutable, afinidad y aprobación GLOBAL de root. Ninguno envía señales, captura, prepara ni crea destinos.

```bash
bash /home/jdiaz/scratch/stage5_control_resource_transition_candidate_20261011/launch_domain_2to4_CONDITIONAL.sh "$GUARD_PROTOCOL_REVIEWED" "$ROOT_TRANSFER_APPROVAL_DOMAIN"
bash /home/jdiaz/scratch/stage5_control_resource_transition_candidate_20261011/launch_spectral_2to4_CONDITIONAL.sh "$GUARD_PROTOCOL_REVIEWED" "$ROOT_TRANSFER_APPROVAL_SPECTRAL4"
bash /home/jdiaz/scratch/stage5_control_resource_transition_candidate_20261011/launch_spectral_4to8_CONDITIONAL.sh "$GUARD_PROTOCOL_REVIEWED" "$ROOT_TRANSFER_APPROVAL_SPECTRAL8"
```

Num8 sigue PENDING/BLOCKED. El certificadoa7e0170 cubre implementación1/2/4, noNum8 ni una trayectoriaL480 completa. Su futuro checker usa cpu_start6→físicos6..13: deben estar todos libres tras parar espectral4 y completar temporal/dominio/1550/wide. El launcher adicional exige aprobación específica de recursos conSHAchecker, sourcepins actuales y CPU6..13; usa timeout300 por oráculo/candidato, salidasn8 nuevas y un certificado SEPARADO. Expone los11otros hashes históricos y el puente exacto de2frontends hacia la cápsula actual; no finge13fuentesdriver idénticas. No usar build_certificate.py ni sobrescribir a7e0170. Esta evidencia es de implementación por nodo/espectro/midpoint congelado, no de convergencia física futura.

```bash
bash /home/jdiaz/scratch/stage5_control_resource_transition_candidate_20261011/run_num8_oracle_CONDITIONAL.sh "$GUARD_PROTOCOL_REVIEWED" "$ROOT_NUM8_TEST_RESOURCES_APPROVAL" "$NEW_NUM8_CERTIFICATE_PATH"
```

Tras auditoría futura, usar únicamente el selector nuevo `--control-resource-transition-receipt` para domain/spectral; nunca alias mesh ni ignorar sources/metadata. Antes de cada ejecución root revisa el comando concreto, los8cores si corresponde y sus recibos. Estado actual de este bloque: sin captura, transferencia, señal, lanzamiento, PDF ni stage5admission.


### Corrección de revisión del candidato de recursos (2026-10-11, NO ejecutado)

El candidato inicial8573928a queda preservado como antecedente NO aplicable. Usar sólo la revisión `/home/jdiaz/scratch/stage5_control_resource_transition_candidate_20261011/candidate_receipt.json`, SHA `7a4067cd99794c70fa191bb77f73c3c3218b19feac6bbe4dc12e284f805c8113` y patch `4c73318f8da5243fc1d88f1935a1e25fc4de678fc1d9347f4a97581a14e26645`, aún pendiente de aprobación de fuentes/transferencias reales.

Las mismas rutas de comandos anteriores ahora contienen tres correctivos: resume_argv conserva su contrato de argumentos, mientras runner_argv añade Python actual y el runner de la cápsula congelada (SHA verificada, PYTHONPATH=runtime_sources antesexec); Num8 compara también las excepciones contra la clase/texto históricoNum4; el checker NO fija afinidad por cpu_start. Por eso su preflight GLOBAL sigue sin taskset y, únicamente después, cada invocación oráculo/candidato usa `taskset -c 6-13 timeout 300 ... --threads 8`. Los launchers de simulación continúan SIN taskset exterior: el runner admite global y fija su afinidad internamente. 20 pruebas sintéticas y vectores de shell PASS; sinRHS/Num8/captura/transferencia/señales/lanzamientos/PDF/admisión.


## Ejecución real 2026-10-11T21:05:36.850794+00:00 — extensión espectral desde20ps

Preparación real revisada:27arraysdtype/shape/Cbytes idénticos; metadatos sólo duration20.0→None; fuente13/controlador/E0/ledgers/counters/next_h conservados. Recibo raíz /home/jdiaz/scratch/stage5_control_window_extension_plan_20261011/root_spectral_extension_prelaunch_review.json SHA9c2ba3ac2418807d000a055a05eb3844a7b424b67f6dd078a8e906ec58cea878. Globalhostrechecked:16núcleos físicos,14máximo tras reanudar, reserva14/15,111GiBdisponibles,quotas unlimited;2hilosNumBaCPU6/7/BLAS1. Horizonte hasta cruce confirmado+5ps o200ps; nunca reset completo. Previsión:varias horas según adaptación, sin promesa de speedup.

~~~bash
bash /home/jdiaz/scratch/stage5_control_window_extension_plan_20261011/launch_spectral_extension.sh
~~~

Salida: /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_trigger_v6_20261011 (foreground/checkpoints).


## Ejecución real 2026-10-11T21:07:29.134813+00:00 — fino guard2b desde22.26ps

Parada coherente SIGTERM raíz enbanda[22.185,22.435]ps; checkpointfinal22.2599999999998psSHA963ba2a7dc223e69e577404c16dfeed7943221af9ad992e1d241bf124c10eba8. Cambio local1e-5→2.5e-6,NL2e-7/cap.02sin cambio.27arraysdtype/shape/Cbytes exactos y todosmetadatos salvo4clavesguarda;E0/ledgers/counters/controller/next_h/fuente13/historialmáximospreservados. Revisión raíz /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/root_fine_initial_prelaunch_review.json SHA dd4ece8a96f1e2171a8910640ae7f255ed81b662ac777de675385fdb53991954. Endpointglobal0.563286%,límiteglobal1%independiente. Globalhostrechecked:2hilosCPU12/13,14núcleosfísicosasignados/reserva14/15,BLAS1,scratch. Duración esperadavariashoras,8ospeedupsNOadmitidos aún.

~~~bash
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/launch_fine_guard2b.sh
~~~

Salida /home/jdiaz/scratch/stage5_effective_L320_h35_775_adaptive_v6_guard2b_20261011. Horizonte cruce confirmado+5ps o200ps,sinresetcompleto.


###2026-10-11 continuación real de malla fina Num4, misma física

Nominal775 terminó detección+margen; fino guard2b preservado en22.4625ps. Los27 arrays, todo el ledger y guardas son exactos; sólo threads2→4. Auditoría derivada autentica hashes que el runner no declara al interrumpirse; originales intactos. CPU0–3, BLAS1,12núcleos físicos usados globalmente de16 (reserva14/15). Corrida larga autorizada en este chat; salida scratch/checkpoints. RootreviewSHA ed059053ffa4896e973a014d7e95a4bfb5b093d6f987b0df4c14e602a017d33e.

```bash
APPROVED_RESOURCE_TRANSFER_SHA256=fd8289953ef1d6a8bdeacbd60c3ee8150b541d6677d43cfd4768752cadf6dc54 bash /home/jdiaz/scratch/stage5_fine_resource_transition_candidate_20261011/launch_resources4_REVIEW_ONLY.sh
```


### 2026-10-11 continuación real control de dominio480nm Num4

Fino ya Num4 en0–3. Control480nm guard2b conserva27arrays y todo ledger; sólo threads2→4. Ancla recurso 26.34000000000002ps; ancla guarda22.235ps intacta. CPU10–13, BLAS1,14/16núcleos físicos globales y reserva14/15; larga autorizada. Rootaprobación a85246bd6b5e112cc7aa8727a98fd91e790b703c3779aa8add686f12bd6ccdef.

```bash
bash /home/jdiaz/scratch/stage5_control_resource_transition_candidate_20261011/launch_domain_2to4_CONDITIONAL.sh /home/jdiaz/scratch/stage5_guard_protocol_actual_four_roles_prepare_20261011/protocol_nominal_temporal_domain_mesh_actual_REVIEW_ONLY.json /home/jdiaz/scratch/stage5_control_resource_transition_candidate_20261011/root_actual_domain_resource_prelaunch_review.json
```


## 20261011_SPECTRAL_GUARD2B_ACTUAL_22325 — continuación espectral real

Registro UTC: 2026-10-11T22:25:47.223584+00:00

Propósito: continuar el control espectral775nm desde22.32499999999994ps con guarda local2.5e-6, conservando fuente física v6, reloj postcascada, E0, trabajo acumulado, estados e historiales. Revisión root a42978559035dc798177893995766a5220f768fcfa029127ff937b2fa3bb1f14. No acredita todavía convergencia ni cierra etapa5.

Recursos: primer lanzamiento NumBa2, CPU6/7, BLAS/OMP1; tres otras corridas ocupan10núcleos, total12/16(75%). CPU4/5 libres tras1550; posible transición posterior2→4 sólo con transferencia revisada. Scratch para datos/cachés. Comando foreground observable y checkpointed; horizonte de gatillo0.1mV,confirmación0.5ps,+5ps o censura200ps, sinreset completo.

```bash
bash /home/jdiaz/scratch/stage5_guard2_contingency_plan_20261011/four_controls_plan/launch_spectral_guard2b.sh
```

Salidas: /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_trigger_v6_guard2b_20261011/{manifest.json,checkpoint_latest.npz,trajectory.npz,snapshots.npz,execution.console.log}. El estado preparado y el padre detenido se preservan en provenance y anchor. Tiempo esperado: horas si se mantiene hasta su horizonte; la reasignación de recursos no reinicia la física.


## 20261011_SPECTRAL_RESOURCES4_ACTUAL_22795

Reanudación autorizada del control espectral desde 22.794999999999646 ps, cambiando sólo NumBa 2→4 y CPU6,7→4,5,6,7. Mismo estado, ecuaciones, origen energético, guarda2.5e-6 y controlador: umbral positivo0.1mV, confirmación0.5ps, margen5ps, techo200ps sin gatillo. Recibo real29e90285ebb58ebde94694c88b4ea35e869cb93d21f9ec41b4948bfffc6d4cd4; aprobación raíz aae3d46ba2436fc315402bd688e8472f19a79239537a302741b08ccb9e3b0f38. Prerrequisitos1550 y ancho amplio terminados. Preservada preparación incompleta anterior; corregida copia de malla archivada sin cambiar física.

```bash
bash /home/jdiaz/scratch/stage5_control_resource_transition_candidate_20261011/launch_spectral_2to4_CONDITIONAL.sh /home/jdiaz/scratch/stage5_guard_protocol_actual_five_roles_prepare_20261011/protocol_nominal_temporal_domain_mesh_spectral_actual_REVIEW_ONLY.json /home/jdiaz/scratch/stage5_spectral_resources4_stop_prepare_plan_20261011/root_actual_spectral_resource_prelaunch_review.json
```

Recursos: 4 trabajadores espectrales, BLAS/OMP1; total de esta campaña14/16 núcleos físicos, reservados14,15 completos; RAM disponible~111GiB. Duración prevista: horas, dependiente del paso aceptado real. Resultados en /home/jdiaz/scratch/stage5_effective_L320_h5_775_spectral_trigger_v6_guard2b_resources4_20261011: manifest, checkpoint, trayectoria, snapshots, controlador y log; detener sólo por su gatillo+margen o censura. Lanzamiento en primer plano por Codex autorizado por el usuario.


## 20261012_SSH_TRANSPORT_RECONNECT_READ_ONLY_OBSERVER

Usado el 2026-10-12 para recuperar observación tras perder simultáneamente los clientes SSH originales. Los PID numéricos 279800 y 296286 continuaron vivos y avanzando. El exit 1 pertenece al transporte, no al solver. El observador abre sólo lectura de stdout de los tees existentes y lee manifiestos; no reinicia simulaciones, no envía señales, no cambia estados ni afinidades. CPU 8, un hilo; duración ligada a las dos corridas existentes. Programa SHA256 ab27cbae411ee4cc2ee3319f805d33740c12e620d63afe07d8c09874979895d5; recibo de adjunción SHA256 3a7ec9d092b681d0443525c9503f489b95bde5eb641f72f46f6f6a7bad74467c. Entrada histórica usada, no un comando pendiente que deba repetirse.

```bash
ssh -T -o BatchMode=yes geminga 'env PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 taskset -c 8 /home/jdiaz/.conda/envs/snspd/bin/python -u /home/jdiaz/scratch/stage5_ssh_observation_reconnect_20261012/observe_existing_foreground_stdout.py'
```

## 20261012_STAGE5_COMPLETE_FINAL_PUBLICATION

Estado: COMPLETED. Etapas5.1,5.2,5.3 cerradas como modelo neutral efectivo exploratorio. No quedan simulaciones obligatorias ni procesos activos de esta campaña. El cuaderno operativo actualizado es /home/jdiaz/GEMINGA_COMMANDS_PRESENT.md; las órdenes largas anteriores, intentos fallidos, políticas numéricas y recibos permanecen arriba sin cambios.

Admisión: /home/jdiaz/scratch/stage5_final_evidence_completed_20261012/admission.json
SHA256: ed30c8525ecb3952e4042eca3b84e9ea0c2c78f41d7e0f17d880779b68747c4a
Cobertura: /home/jdiaz/scratch/stage5_final_evidence_completed_20261012/finite_control_coverage_review.json
SHA256: 2c404823c78594ae22fad71f1fee50e62b5d9ad3ee7beab89436b9bfd40252ac

Los cuatro controles775nm pasan en0–36.6ps; extremo crítico33ps. Todas las colas posteriores permanecen fuera de convergencia de cuatro ejes, con auditorías originales individuales hasta cada cruce+5ps. Nominal775: cruce31.94131686165ps/final36.94131686165ps. Nominal1550: censura200ps sin cruce. Señal positiva0.1mV; confirmación0.5ps; margen5ps desdecruce. La revisión espectral Vout5% fue posterior a los datos; el FAIL previo2% sigue visible.54 estados es referencia exploratoria;80 controlopcional.

Informe final: /home/jdiaz/pysnspd/output/pdf/implementation/Informe_final_etapa5_modelo_efectivo_20261011.pdf
SHA256: bdad8b47bb4da86c5e17f2d31771378f8ed0f49c83755386d42b36ea9349fd6b
14páginas, todas revisadas; build/renders/recibo numérico completo en /home/jdiaz/scratch/stage5_final_report_completed_20261012_final.

Documentos matemáticos finales: /home/jdiaz/pysnspd/Matematica/Principal.pdf (2p), Auxiliar.pdf (8p), Precomputo.pdf (5p). Compilados y revisados15/15 páginas, fuentes y procedencia en /home/jdiaz/scratch/stage5_math_documents_completed_20261012_final. Dependencia documental etoolbox instalada sólo en scratch; sistema TeX sin cambios globales.

Pipeline visual05 comprobada sobre ambos originales: /home/jdiaz/scratch/stage5_nominal775_saved_visual_final_v2_20261012 y /home/jdiaz/scratch/stage5_nominal1550_saved_visual_final_v2_20261012. PNG/PDF físicos, con detalle de lectura y normalizaciones energéticas explícitas. No se lanzó evolución adicional para generar documentos.

Este cierre admite una hipótesis BCS/KWT útil; calibradores τee(Tc),τep(Tc),η,σ,interfaz/material y comparación experimental son estudios posteriores. No se atribuye jitterestadístico, retraso óptico absoluto ni reset total a estas ventanas.


## 20261012_STAGE5_PARAMETER_SENSITIVITY_PREPARED

## Campaña de sensibilidad eléctrica — preparada el 12 de octubre de 2026

**No ejecutada.** Hasta16 casos distintos, en grupos que permiten detenerse tras
obtener información útil. No hace falta correrlos todos para utilizar el modelo
cerrado. Cada caso prepara su propio DC y mantiene ecuaciones, geometría,
circuito y perfil numérico fresh; cambia un solo input respecto de la referencia
de su misma longitud de onda. El sesgo se mantiene positivo.

| Grupo | Casos nuevos | Input físico |
|---|---:|---|
| bias — empezar aquí | 5 | Nominal775; 775 a0,65/0,75Idep; nominal1550; 1550 a0,75Idep |
| rates — calibradores | 6 | 775nm: τee(Tc)=2.5/10ps, τep(Tc)=12.35/49.4ps, τesc=7.5/30ps |
| source — deposición | 5 | 775nm: ηret=0.5/0.8, σ=5.338539125824206/10.677078251648414nm; 1550nm: ηret=0.8 |

rates incluye nominal775, source ambas referencias; si ya existen completas
en ESTA campaña, verifica identidad, todos sus SHA terminales y una auditoría
nueva antes de reutilizarlas. **No reutiliza los originales históricos como
si fueran fresh.** La unión de los tres grupos contiene16 simulaciones,
no19. El grupo all es opcional y solicita esa unión completa.

Los puntos son **sondas exploratorias escogidas para medir sensibilidad**;
no son intervalos experimentales, límites de calibración ni sus centros.
Los nominales libres siguen Allmaras: τee(Tc)=5ps y τep(Tc)=24.7ps.
El nominal de escape es15ps, retención0.667, σ=7.549834435nm y
I=21.29720614489798µA. La escala de sesgo usa la predicción bulk efectiva
Idep=30.36744305465µA de esta película de7nm:
I065=19.7388379855225µA e I075=22.7755822909875µA.
No se identifica con una corriente crítica experimental medida.
Un nuevo DC sin rama superconductora admisible detiene el caso y conserva
su evidencia.

### Órdenes normales

Ejecutar en geminga con la preparación común de este cuaderno. El directorio
sensitivity_20261012_01 debe ser nuevo al iniciar bias; los otros grupos usan
exactamente el mismo destino y --resume. Cambiar el sufijo en todas las órdenes
si ese destino ya pertenece a otra campaña.

~~~bash
cd /home/jdiaz/pysnspd
set -o pipefail
mkdir -p /home/jdiaz/scratch/photon_runs/logs

/home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/06_effective_sensitivity_campaign.py \
  --output-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01 \
  --group bias --execute \
  2>&1 | tee -a /home/jdiaz/scratch/photon_runs/logs/sensitivity_20261012_01_bias.log
~~~

Tras completar y auditar bias, continuar por separado si se desea medir
los calibradores:

~~~bash
/home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/06_effective_sensitivity_campaign.py \
  --output-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01 \
  --group rates --execute --resume \
  2>&1 | tee -a /home/jdiaz/scratch/photon_runs/logs/sensitivity_20261012_01_rates.log
~~~

Para la fuente inicial:

~~~bash
/home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/06_effective_sensitivity_campaign.py \
  --output-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01 \
  --group source --execute --resume \
  2>&1 | tee -a /home/jdiaz/scratch/photon_runs/logs/sensitivity_20261012_01_source.log
~~~

Ante una interrupción, repetir el mismo grupo/comando, añadiendo --resume
también a bias. El fotón no se deposita de nuevo en un checkpoint.
Un lock exclusivo impide lanzar dos escritores/grupos simultáneamente
sobre esta campaña. No editar fuentes físicas ni el orquestador mientras
se usa su plan: sus SHA forman parte de la identidad. Si cambian, usar
un destino nuevo.

Perfil explícito: L=320nm, hmesh=5nm, película7nm, Tb=0.9K,
cuadratura54 estados/96 fonones/48 transporte; cap0.02ps, NL2e-7,
guarda local2.5e-6 desde cero, máximo40 iteraciones.
**No son los defaults genéricos del parser.** La campaña histórica admitida
conserva sus guardas por segmentos y no se declara idéntica a este perfil fresh.

Evento: Vout−Vout,DC positivo, umbral0.1mV, confirmación0.5ps,
histéresis5% y final en cruce+5ps; sin detección, techo200ps y
latencia censurada. No espera pico, FWHM, recuperación ni reset total.
El cero de tiempo es la transferencia efectiva postcascada, no la absorción
óptica absoluta. Mantiene el circuito completo de tres estados y la
inductancia exterior congelada desde cada DC.

Recursos: **un caso en primer plano, cuatro núcleos físicos, cuatro hilos
Numba, BLAS/OpenMP=1**. El runner repite inventario/afinidad/cuota/RAM y aplica
techo90% para cada caso. Auditoría y visual son procesos únicos.
Datos, cachés, checkpoints, logs y figuras van a scratch. Reservar núcleos
completos libres y comprobar carga ajena según la política del servidor.
Planificar varias horas por caso; algunos casos censurados a200ps con este
perfil estricto pueden requerir un día o más. No hay ETA calibrada ni plazo
garantizado; los16 pueden requerir varios días. El grupo bias reduce el
compromiso inicial a cinco casos. Progreso y ETA al objetivo vigente salen
por pantalla, con checkpoints cada1ps o60s. **No se ejecutaron estas simulaciones
al preparar el cuaderno.**

Productos: campaign_plan.json fija inputs/comandos/SHA;
runs/<caso>/ conserva material, DC, trayectoria, snapshots, fuente y
checkpoint completos. audits/<caso>_<fecha>/audit.json comprueba balance,
soporte, corriente, deposición, guardas y evento.
campaign_execution.json y aggregate_<grupo>.json autentican terminales y
auditorías. Un fallo detiene la secuencia y preserva lo anterior; no se presenta
como ausencia física de detección. Los nuevos inputs no heredan los cuatro
controles de convergencia de la referencia775: es exploración de sensibilidad,
con conservación individual y necesidad de contrastar precisión donde
cambie el régimen.

### Órdenes visuales comparativas

Ejecutar tras completar el grupo correspondiente. Leen estados guardados,
sin evaluar RHS ni lanzar simulaciones. --individual añade figuras físicas
y snapshots de visual05 para cada caso. Cada destino debe ser nuevo;
para repetir, cambiar visual_..._01 por un sufijo nuevo, preservando anteriores.

~~~bash
env NUMBA_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/06_plot_effective_sensitivity_campaign.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01 \
  --output-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/visual_bias_01 \
  --group bias --individual
~~~

~~~bash
env NUMBA_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/06_plot_effective_sensitivity_campaign.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01 \
  --output-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/visual_rates_01 \
  --group rates --individual
~~~

~~~bash
env NUMBA_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/06_plot_effective_sensitivity_campaign.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01 \
  --output-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/visual_source_01 \
  --group source --individual
~~~

Omitir --individual si sólo se desean comparaciones eléctricas rápidas.
El análisis entrega PNG/PDF de Vout menos basal[mV], con recuadro[µV]
para señales pequeñas; cada familia se corta en la ventana temporal común
de sus miembros. Las líneas de cruce sólo muestran eventos confirmados
guardados. metrics.csv y comparison_manifest.json incluyen:
latencia relativa[ps], diferencia respecto de referencia sólo si ambas
detectan, máxima diferencia Vout en ventana común[mV], mínimo Is[µA]
y pendiente local[mV/ps], ajustada sobre muestras aceptadas en cruce±0.25ps.
No asigna pendiente ni latencia ficticia a censurados.
Esa pendiente no incluye filtro, ganancia, ruido ni ancho de banda del setup;
no constituye una predicción de jitter.

Los máximos Vout se llaman **máximos observados en la ventana guardada**:
una corrida que termina en cruce+5ps puede seguir subiendo. No se usan como
pico total del pulso ni para comparar amplitudes a ventanas diferentes.
Los cruces representan un discriminador ideal en la carga50Ω; la comparación
absoluta con voltajes medidos requiere calibrar el readout real.

Interpretación crítica: variar τee(Tc) cambia conjuntamente BGK y movilidad
KWT; variar τep(Tc) cambia conjuntamente normalización acústica de colisiones
y movilidad KWT. Variar σ o ηret con inicialización FD/BE común también cambia
el reparto energético electrón/fonón; σ no es una sensibilidad geométrica pura.
La extrapolación efectiva KWT/rigidez en la zona caliente por encima de Tc
sigue siendo una hipótesis calibrable, identificada en Parametros.pdf.

--plan-only inspecciona un destino NUEVO sin lanzar hijos numéricos.
Si luego se ejecuta ese mismo destino, añadir --resume.
Las pruebas rápidas cubren los16 comandos con parser real05, censuras,
ventanas comunes, checksums/identidad, lock, memoria descomprimida y
coordinación de receipts con fixtures, sin evoluciones nuevas.


## 2026-10-13 — Case report individual: pipeline visual07

# Case report individual del modelo efectivo

Pipeline: /home/jdiaz/pysnspd/plot_pipelines/07_case_report.py

Produce **un PDF autónomo por run** del formato efectivo de pipeline05
(pysnspd.effective_photon_run.v1). Acepta varios directorios como comodidad de
ejecución: nunca superpone, resta ni compara casos. Los holds sin fotón y los
formatos históricos no se presentan con el reloj de una corrida fotónica.

## Contenido estándar

1. Parámetros exclusivos del caso: geometría, espesor y resolución separados;
   material, corriente DC, absorción/deposición, energía fotónica y retenida,
   calibradores Allmaras, circuito, discriminador y precisión numérica.
2. Observables de laboratorio: Vout menos basal DC, corriente de lectura
   Ib−Is, cruce, confirmación, latencia relativa, final guardado y pendiente
   local. También Ib, Is, Vdev, vc y potencia en la carga.
3. Estado interno: curvas de gap, energías y temperaturas equivalentes;
   snapshots del gap completo y de la región central; fase relativa al DC;
   Te, Tph, exceso energético, Js,x, Jn,x, |Jtotal| y potencial.
4. Balance independiente, precisión, procedencia, hashes y límites.

Con seis snapshots automáticos son **16 páginas A4 apaisadas**. Los mapas
comparten escala de color entre los tiempos del mismo campo. Se muestran
unidades físicas: nm, ps, meV, K, kJ/m³, A/m², V y A con prefijos legibles.

Vout no incluye filtro, amplificador, ruido ni ancho de banda experimental.
Su máximo es el máximo **observado**, no el pico de un pulso completo.
La pendiente no predice jitter. El reloj empieza en la transferencia efectiva
postcascada, no en la absorción óptica absoluta. Sólo un cruce confirmado
produce latencia; un caso censurado queda sin latencia.

## Comando para cuatro casos terminados

Ejecutar en geminga. Cada repetición requiere un destino nuevo.

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir \
    /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs/nominal775 \
    /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs/bias065_775 \
    /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs/bias075_775 \
    /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs/nominal1550 \
  --output-dir /home/jdiaz/scratch/case_reports/bias_20261013_01 \
  --workers 2
~~~

Para otro caso, sustituir la lista de --run-dir; uno solo también funciona.
--workers 2 paraleliza dos PDFs independientes, con un núcleo físico y
bibliotecas de un hilo por trabajador. Revisa afinidad, cuotas y memoria,
reserva núcleos completos y respeta el techo 90 %. Los dos ejemplos nominales
tardaron unos seis segundos por caso; no es una garantía para archivos mayores.

## Caso aún en ejecución

Añadir --allow-incomplete con un destino nuevo. Ejemplo:

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs/bias075_1550 \
  --output-dir /home/jdiaz/scratch/case_reports/bias075_1550_partial_20261013_01 \
  --allow-incomplete --workers 1
~~~

Captura una copia coherente de un checkpoint atómico; no detiene ni escribe
en la simulación. Todas las páginas quedan PROVISIONALES si no hay terminal
completo. El final de trayectoria y el último snapshot se indican por separado.
No convierte ausencia de cruce hasta el checkpoint en censura a 200 ps.
Si el caso ya terminó, usa y autentica sus productos terminales.

## Elegir instantes

--snapshots-ps 0 1 5 10 15 20 solicita las muestras guardadas más cercanas;
se incluyen siempre inicial y final. Máximo 24 instantes pedidos, dentro del
intervalo de snapshots existente. No se interpola dinámica en tiempo.

Sin esa opción, una detección muestra inicial, 1 ps, primera hotbelt guardada,
muestras próximas a cruce/confirmación y final. Sin detección: inicial,
1 ps, 10 ps, cuarto, mitad y final de la ventana, eliminando duplicados.
La hotbelt es conectividad del grafo nodal con |Delta|≤0.2|DeltaDC|
entre bordes transversales; no es el gatillo eléctrico ni prueba continuidad
entre muestras.

## Datos originales y reconstrucciones

Las curvas y mapas de gap/energía/temperatura se leen de trajectory.npz y
snapshots.npz, autenticados por manifest.json, junto con el DC y la malla
original. Se rechaza corrupción o una versión distinta de las tres fuentes
utilizadas para reconstruir corrientes.

El gap se interpola como campo complejo P1 antes de tomar su magnitud.
La fase se muestra en nodos, sin interpolar ángulos y ocultando amplitudes
menores a 0.2 del DC. Un salto de fase envuelta no prueba un vórtice.

Los flujos superconductores se recuperan de la derivada gauge exacta del
gradiente energético, independiente de las ocupaciones no guardadas.
El potencial y flujos normales se recuperan de la misma neutralidad algebraica,
con Is del estado aceptado exactamente coincidente con cada snapshot.
Se contrasta Vdev y la conservación de corriente con los valores guardados.

J por triángulo es una **aproximación de visualización**: mínimos cuadrados
ponderados de proyecciones sobre caras duales conservativas, divididas por
espesor físico y longitud dual. El recibo conserva residuos de ese ajuste.
Las flechas indican dirección con longitud visual fija de 6 nm; se omiten
magnitudes menores al 0.01 % del máximo entre los snapshots seleccionados.
El color aporta la magnitud o componente firmada. No se interpretan cambios
locales de signo como inversión de corriente global.

Te y Tph son observadores de energía equivalente, no estados térmicos
impuestos. La energía electrónica guardada representa cuasipartículas y
excluye condensación/gradiente; el balance total sí incluye esos términos.
No se inventan distribuciones completas ni mapas de colisiones o calor KWT
que los snapshots no guardaron. No se evalúa el RHS cinético ni un nuevo DC,
no se evoluciona el sistema y no se aplica un ajuste experimental.

## Productos y validación

En el destino scratch: case_report_<caso>_<identidad>.pdf, un recibo JSON
por PDF, case_reports_manifest.json como inventario y .build para capturas.
Nunca reemplaza destinos existentes. No modifica datos de entrada.

Validación del 13 de octubre de 2026: 20 pruebas rápidas; ejemplos nominales
775 y 1550 nm, detección y censura respectivamente, 16 páginas cada uno.
Revisión visual de todas las páginas, escalas y etiquetas; cotejo de hashes
y reconstrucción de Vdev con error máximo aproximadamente 1.4e-17 V.
Los 13 archivos físicos del checkpoint operativo permanecen sin cambios.
Los ejemplos están en /home/jdiaz/scratch/case_report_validation_20261013_release.


## 2026-10-13 - Case report refresh: English, fixed six snapshots, winding and model summary

Pipeline: /home/jdiaz/pysnspd/plot_pipelines/07_case_report.py
Guide: /home/jdiaz/pysnspd/plot_pipelines/CASE_REPORT.md

The current report is English, 19 pages, with larger chart labels, one PDF per
case, quantity-specific palettes, phase-winding vorticity and section6 containing
the coupled model equations and complete three-state circuit schematic.
This entry supersedes the earlier snapshot-selection behavior.

--snapshots-ps takes the first six requested times [ps], repeats the last if fewer
are supplied, preserves order/repetitions and selects nearest saved states.
No initial/final snapshots are silently added. Every panel labels requested and
actual time. Out-of-range requests select the nearest stored endpoint.
Omit the flag for six automatic targets. No new simulation or full reset.

~~~bash
cd /home/jdiaz/pysnspd
CASE_RUNS=/home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir "$CASE_RUNS/nominal775" "$CASE_RUNS/bias065_775" \
    "$CASE_RUNS/bias075_775" "$CASE_RUNS/nominal1550" \
  --snapshots-ps 0 1 5 10 20 35 \
  --output-dir /home/jdiaz/scratch/case_reports/english_20261013_01 \
  --workers 2
~~~

Choose a fresh output directory for every invocation. --allow-incomplete reads
one coherent checkpoint and labels the PDF provisional. It does not interrupt
the solver or infer terminal censoring from a partial window.
Resources: two independent foreground workers, one physical core each, BLAS/
OpenMP/Numba=1, shared ceiling90%, all outputs/caches/captures under scratch.
Two nominal reports required about six seconds per case. Larger archives can differ.

Vorticity is signed condensate phase winding, with P1 zero localization.
Gray marks undefined/insufficient-amplitude contours, with a diagnostic floor
1e-6 of local DC amplitude. It is not curl(J), reconstructed magnetic flux or
continuous vortex tracking. Primary method:
https://www.mcs.anl.gov/papers/P5324-0315.pdf, Eq.(1).
Receipt contains actual/requested times, repetitions, cadence, winding checks
and a follow-up acquisition queue when saved time resolution is insufficient.
No physical solver or save-rate code was modified. Vout still lacks the
experimental filter/amplifier/noise response.
Validation:31 quick tests, two19-page examples and an actual two-time/padded CLI run.

## 2026-10-14 - Paired base cases: eta, relaxation, preparation and domain

Prepared only; the physical photon campaign has NOT been executed.
Pipeline: /home/jdiaz/pysnspd/pipelines/08_effective_basecase_sweep.py
Plan: /home/jdiaz/scratch/photon_runs/basecases_eta_20261014_01/campaign_plan.json

Purpose: compare 775/1550 nm at the same 21.5 uA and tau_escape=20 ps,
without tuning parameters separately for each wavelength. Rs=608 ohm/square
is retained: this closure has Idep proportional to 1/Rs at fixed D,Tc,Tb,w,
so replacing608 by425 would increase Idep30.3674 to43.4433 uA.
425 is a conditional340/RRR conversion, not a cryogenic measurement.
No RHS, solver, physical source pin or default material parameter was changed.

20 cases /10 wavelength pairs:
- eta=.667 reference and requested .7,.8,.9,1.0 (10 cases).
- at eta=.8, tau_ee(Tc)=2.5 and10 ps (4 cases).
- at eta=.8, sigma=7.549834435/sqrt(2) and7.549834435*sqrt(2) nm (4 cases).
- at eta=.8, length=640 nm, with the same5 um device and exterior-inductance
  accounting (2 cases; reference length320 nm is already in the eta matrix).
All remaining parameters use the admitted exploratory54/96/48 quadratures,
h_mesh5 nm, tau_ep(Tc)24.7 ps, Tb.9 K, Tc8.65 K and D5e-5 m2/s.
The source is common FD/BE at transfer with fresh DC gap unchanged.
Eta is TOTAL film energy retention at transfer, conditional on photon absorption.
It is not optical absorption probability or electronic-only conversion efficiency.
Eta1 is an upper-bound source test; later20 ps phonon escape stays active.
Sigma changes both spatial preparation and the initial electron/phonon partition.
Tau_ee changes BGK and KWT mobility together. These points are hypotheses,
not measured uncertainty intervals or a new experimental calibration.

Foreground, up to13 independent runs, one Numba thread/distinct physical core
per run; BLAS/OpenMP1. Inventory checked2026-10-14:16physical/32logical,
no CPU quota, about118 GiB available RAM and945 GiB free scratch.
Workers occupy13/16=81.25% physical cores; one admitted core is reserved for
coordination, two complete cores stay outside the87.5% allocated core set.
Conservative scheduling reservation55 GiB including coordinator; rechecked
against90% available memory. This is allocation, not a utilization guarantee:
I/O, serial setup, and the final incomplete wave can lower instantaneous usage.
SMT logical-CPU percentages can be roughly half the physical-core percentage.
The script rejects requests below50% of visible physical cores, over-budget
requests and a pre-existing effective photon solver under jdiaz.
Do not run simultaneously with another heavy campaign without a shared budget.
Expected duration: multiple hours; no measured global ETA before executing
this new concurrent matrix. Cases ending at200 ps can dominate the duration.
Each running case logs its measured progress and current-horizon ETA.

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  TMPDIR=/home/jdiaz/scratch/tmp/pysnspd \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u \
  pipelines/08_effective_basecase_sweep.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/basecases_eta_20261014_01 \
  --workers 13
~~~

Observation: positive baseline-subtracted Vout0.1 mV, confirmation at crossing+.5 ps,
end at crossing+5 ps. No crossing: censor at200 ps, latency unassigned.
Never requests full reset. No PDF/report generation.
Outputs in scratch: runs/<case>/{manifest.json,trajectory.npz,snapshots.npz,
checkpoint_latest.npz,dc/,material.*}; logs/; audits/; campaign_execution.json;
latency_pairs.csv (accepted-controller crossings and differences only when both
members detect; censored times are blank). Individual audits do not establish
new convergence or physical calibration. Failed cases remain marked failed;
independent cases continue and the campaign exits nonzero.

Resume the SAME command with --resume appended. Source/parameter/hash changes
are rejected; completed cases are authenticated and not simulated again.
SIGINT/SIGTERM stops new work and requests accepted-state checkpoints from
active physical solvers. No forced kill, background service or polling daemon.

Validation before handoff: all20 generated commands parse into the exact
expected runner parameter identities; all10 pairs differ only in wavelength;
50% minimum and global ceiling checked; synthetic foreground workers exercise
bounded concurrent slots, nonoverlap, exact completed-case resume and rejection
of changed source pins. Fixtures only under scratch/tmp, not scientific runs.
The13 physical runner pins equal the existing October12 campaign.
Static bulk derivation matches30.3674430546503 uA to floating precision;
low-temperature closed form30.38655886 uA differs by0.063%.

## 2026-10-15 - Individual reports for13 completed basecase runs

Read-only inspection found13 COMPLETE_AUDITED cases,5 RUNNING and2 FAILED before photon preparation (length640 DC Newton limit). Audit receipt hashes and terminal manifest hashes were checked. This command has NOT been executed. Produces one independent PDF per completed case, with six automatically selected snapshots adapted to the case duration. It excludes live and failed cases, preserves originals and uses a fresh scratch destination. Two plotting workers use cores12 and11 under the current full-affinity topology; active photon workers were on2,4,6,8,10. Recheck global allocation if launching other jobs. No new photon simulation, model change or full reset. Expected plotting cost: order minutes for13 cases, not a measured ETA.

~~~bash
cd /home/jdiaz/pysnspd
RUNS=/home/jdiaz/scratch/photon_runs/basecases_eta_20261014_01/runs
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir \
    "$RUNS/eta0667_775" "$RUNS/eta0667_1550" \
    "$RUNS/eta0700_775" "$RUNS/eta0700_1550" \
    "$RUNS/eta0800_775" "$RUNS/eta0800_1550" \
    "$RUNS/eta0900_775" "$RUNS/eta0900_1550" \
    "$RUNS/eta1000_775" "$RUNS/eta1000_1550" \
    "$RUNS/ee_half_775" "$RUNS/ee_half_1550" \
    "$RUNS/ee_double_775" \
  --output-dir /home/jdiaz/scratch/case_reports/basecases_completed_20261015_01 \
  --workers 2
~~~


## 2026-10-15 — 1.1.0: rigidez térmica candidata y difusión calibrada

Estado: v1.0.0 publicada en 7992f6a; thesis-v1.0.0 preservada. Revisión 1.1.0 preparada, sin admisión fotónica nueva. Dos pruebas de 0.1 ps alcanzaron su límite de 240 s: 775 nm llegó a 0.051875 ps y 1550 nm a 0.046875 ps. Conservaron soporte, capacidad positiva y balance. No completaron la duración ni aportan latencias. No se extrapola su tasa inicial a 200 ps.

Caso principal: 21.5 uA, ancho 80 nm, espesor 7 nm, longitud resuelta 320 nm, malla 5 nm, baño 0.9 K, Tc 8.65 K, escape fonónico 20 ps y sigma 7.549834435 nm independiente de D. Rs=425 ohm/square principal y 608 de control; Dcal=1.4103757147386595 y 0.689136928331325 cm2/s respectivamente, con objetivo Idep=25.8666666667 uA. N0 y ell0 se recalculan por Einstein. El cierre usa G(Delta,Te) libre, U=UBCS+G-sum(Te*GT), calor Q+PJ+Te*dotGT y Ceff>0 sin piso. Eta=0.8 es exploratoria, no ajustada al retraso experimental. Circuito completo heredado; Vout intrínseco sin filtro/amplificador.

Plan: /home/jdiaz/scratch/photon_runs/thermal_basecases_20261015_01/campaign_plan.json. Hay 12 casos, seis pares 775/1550 nm: anclajes; rigidez a Tb; Rs=608 con D calibrada; eta=0.667; eta=0.9; tauEE(Tc)=2.5 ps. Los diez controles requieren ambos anclajes completos y auditados. Una censura a 200 ps no recibe latencia ficticia; un fallo numérico o de capacidad bloquea los controles. Completar el catálogo no establece convergencia ni validación experimental.

Recursos: el ejecutor revisa topología, afinidad, cuotas y memoria. Usa hasta 12 núcleos físicos para workers, un coordinador y el límite global del 90%, contando procesos 05/08/09/10 activos y reservando dos núcleos completos. Las dos anclas usan entre 1 y 6 hilos por caso según disponibilidad; los controles usan hasta diez procesos de un hilo en cola. BLAS/OpenMP queda en 1. Datos, cachés y logs permanecen en scratch. Las corridas pueden tardar horas o días; cada caso informa progreso y ETA de su sesión actual. Checkpoints cada 60 s de pared o 1 ps. Ningún comando simula el reset eléctrico completo.

### A. Anclajes: ejecución larga pendiente

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  pipelines/10_thermal_stiffness_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/thermal_basecases_20261015_01 \
  --stage anchors --max-worker-cores 12
~~~

### B. Controles: después de A auditada

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  pipelines/10_thermal_stiffness_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/thermal_basecases_20261015_01 \
  --stage controls --max-worker-cores 12
~~~

Para reanudar el mismo stage, añadir --resume; deben conservarse fuentes y parámetros. Un fallo de integración/capacidad requiere diagnóstico, no reinicio ciego. Si cambian las fuentes identificadas por hash, usar un directorio nuevo y preservar los resultados anteriores.

### C. Reportes individuales de los anclajes

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir \
    /home/jdiaz/scratch/photon_runs/thermal_basecases_20261015_01/runs/anchor425_eta0800_775 \
    /home/jdiaz/scratch/photon_runs/thermal_basecases_20261015_01/runs/anchor425_eta0800_1550 \
  --snapshots-ps 0 1 5 10 20 35 \
  --output-dir /home/jdiaz/scratch/case_reports/thermal_basecases_20261015_01 \
  --workers 2
~~~

El reporte 07 reconstruye las corrientes con el cierre y la Te guardada de cada caso. Produce un PDF por corrida y elige los seis tiempos guardados más próximos a los solicitados.


## 2026-10-15 — v1.1.1: corrección Pauli y nuevos anclajes

Estado: bloqueo continuo de 1.1.0 identificado y corregido. Sus dos anclajes se detuvieron con aprobación expresa y conservaron checkpoints; no reanudar con ecuaciones distintas. Las pruebas frescas a malla nominal completaron 0.2 ps: 775 nm, 50 pasos/175.14 s; 1550 nm, 48 pasos/169.81 s. Defectos energéticos máximos 4.99e-25 y 2.78e-25 J. No son latencias ni validación de toda la ventana. Guardias, tolerancias, circuito y difusión no se relajaron.

El calor electrónico usa actividades f(1-p) al calentar y p(1-f) al enfriar, con normalización exacta de energía y producción entrópica adicional no negativa. Rigidez Te, Dcal y Ceff>0 permanecen. Cinco documentos matemáticos actualizados en sitio. Datos diagnósticos: /home/jdiaz/scratch/v1_1_1_debug_20261015.

Plan fresco ya preparado: /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_1_20261015_02/campaign_plan.json. Doce casos en dos etapas: anclajes 775/1550 nm a 21.5 uA, tau_escape=20 ps, Rs425, eta0.8, sigma7.549834435 nm; controles de rigidez Tb, Rs608 con Dcal, eta0.667/0.9 y tauEE2.5 ps. Primero terminar y auditar ambos anclajes. Cruce positivo 0.1 mV confirmado durante 0.5 ps, extra5 ps; censura sin cruce200 ps. Sin reset total.

Recursos: hasta12 núcleos físicos de workers, un coordinador y dos núcleos completos reservados; presupuesto90% global. El ejecutor vuelve a contar otras campañas. Al preparar quedaban activos el coordinador08 en CPU13 y un worker05 en CPU10; este inventario puede cambiar. Hasta6 hilos por anclaje cuando haya espacio, sin prometer escalamiento. Datos y cachés en scratch. El costo completo puede ser de horas/días; el tiempo de0.2 ps no es una ETA válida hasta detección. Se registran ETA por sesión y checkpoints. Estos comandos largos NO se lanzaron aquí.

### A. Anclajes frescos v1.1.1

~~~bash
cd /home/jdiaz/pysnspd
env TMPDIR=/home/jdiaz/scratch/tmp/pysnspd \
  PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  pipelines/10_thermal_stiffness_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_1_20261015_02 \
  --stage anchors --max-worker-cores 12
~~~

### B. Controles, sólo tras A completa y auditada

~~~bash
cd /home/jdiaz/pysnspd
env TMPDIR=/home/jdiaz/scratch/tmp/pysnspd \
  PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  pipelines/10_thermal_stiffness_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_1_20261015_02 \
  --stage controls --max-worker-cores 12
~~~

Añadir --resume sólo para la misma campaña1.1.1 interrumpida, con fuentes y parámetros idénticos. La antigua campaña1.1.0 queda histórica.

### C. Reportes individuales cuando terminen ambos anclajes

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir \
    /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_1_20261015_02/runs/anchor425_eta0800_775 \
    /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_1_20261015_02/runs/anchor425_eta0800_1550 \
  --snapshots-ps 0 1 5 10 20 35 \
  --output-dir /home/jdiaz/scratch/case_reports/thermal_basecases_v1_1_1_20261015_02 \
  --workers 2
~~~

Diagnóstico ligero reproducible: con el mismo entorno ejecutar python -m pytest -q tests/test_thermal_pauli_heat.py tests/test_thermal_stiffness_coupled_candidate.py; sin evolución larga ni reportes adicionales.


## 2026-10-16 — v1.1.2: reloj físico, convergencia y techo100 ps

Los anclajes1.1.1 no progresaban de manera fiable:775 nm cerca de3.02 ps y1550 nm cerca de4.51 ps. Se detuvo sólo su coordinador conSIGTERM y se preservaron checkpoints. Reproducción directa: h=1.9073486328125e-8 ps, RHS máximo de población1048.6757/ps, cero iteraciones y cero cambio en todos los campos. El residuo absoluto estaba debajo de2e-7; el reloj avanzaba sobre un predictor nulo. En el estado1550, Ceff/Cfd=5.41e-7 seguía positivo y su dirección continua incrementaba Ceff; no se concluye que ese checkpoint demuestre una capacidad física negativa.

1.1.2 exige al menos una corrección, usa epsilonNL(h)=2e-7*min(1,h/0.02ps), comprueba la respuesta térmica completa dentro de la convergencia y retrocede las correcciones internas para mantener Ceff>0. No cambia las ecuaciones físicas, Dcal, rigidezTe, fuente, calorPauli, circuito ni la guarda independiente de energía. El presupuesto no lineal se hace más estricto al bajar h. El límite adaptativo1e-6 ps produce parada diagnóstica con checkpoint; no se impone piso a h/Ceff ni se cuenta como detección. Los estados1.1.1 quedan históricos; no reanudarlos como1.1.2.

Caso nuevo:21.5 uA, ancho80 nm, espesor7 nm, longitud320 nm, malla5 nm, Rs425, D calibrada, escape20 ps, eta0.8 y sigma7.549834435 nm. Techo sin cruce100 ps; cruce positivo0.1 mV, confirmación0.5 ps, margen5 ps. Sin reset total. El plan12 casos ya está preparado; controles sólo tras ambos anclajes auditados. Resolver el bug no acredita convergencia de latencia ni validez global del cierre térmico.

Recursos: hasta12 núcleos físicos para workers, un coordinador y al menos2 núcleos completos reservados; presupuesto90% compartido. El ejecutor cuenta otros procesos antes de asignar1–6 hilos por anclaje. BLAS/OpenMP1, datos/cachés/logs en scratch. El costo puede ser de horas o días; no hay ETA completa validada. Las pruebas internas de6 ps usan dos procesos de3 hilos, primer plano, checkpoints y límite de30 min por proceso; su estado final se registra en docs/releases/v1.1.2_checks.json. No extrapolar el rendimiento inicial a100 ps.

### A. Nuevos anclajes — comando largo para el usuario

~~~bash
cd /home/jdiaz/pysnspd
env TMPDIR=/home/jdiaz/scratch/tmp/pysnspd \
  PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  pipelines/10_thermal_stiffness_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_2_20261016_01 \
  --stage anchors --ceiling-ps 100 --max-worker-cores 12
~~~

### B. Controles — sólo después de A completa y auditada

~~~bash
cd /home/jdiaz/pysnspd
env TMPDIR=/home/jdiaz/scratch/tmp/pysnspd \
  PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  pipelines/10_thermal_stiffness_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_2_20261016_01 \
  --stage controls --ceiling-ps 100 --max-worker-cores 12
~~~

Para reanudar sólo una interrupción1.1.2 con fuentes, parámetros, horizonte y hilos idénticos, añadir --resume. Una parada NUMERICAL_RESOLUTION_LIMIT requiere revisión del solver/cierre, no un reinicio ciego ni cambiar el límite para ocultar el fallo.

### C. Reportes individuales de los nuevos anclajes completos

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir \
    /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_2_20261016_01/runs/anchor425_eta0800_775 \
    /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_2_20261016_01/runs/anchor425_eta0800_1550 \
  --snapshots-ps 0 1 5 10 20 35 \
  --output-dir /home/jdiaz/scratch/case_reports/thermal_basecases_v1_1_2_20261016_01 \
  --workers 2
~~~

Diagnóstico ligero: con el entorno anterior ejecutar python -m pytest -q tests/test_midpoint_clock_progress.py tests/test_thermal_stiffness_coupled_candidate.py. Evidencia causal y checkpoints internos: /home/jdiaz/scratch/v1_1_2_debug_20261016.
