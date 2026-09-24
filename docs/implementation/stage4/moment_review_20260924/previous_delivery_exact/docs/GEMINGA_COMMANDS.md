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
