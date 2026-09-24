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

El [informe vigente](../output/pdf/implementation/Informe_etapa_4_revision_practica.pdf)
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
y el [informe](../output/pdf/implementation/Informe_etapa_4_respuesta_de_carga.pdf).

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
