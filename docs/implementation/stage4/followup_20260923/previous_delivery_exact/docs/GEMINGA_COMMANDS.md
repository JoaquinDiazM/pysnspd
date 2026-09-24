# Geminga: comandos vigentes

Actualizado el 23 de septiembre de 2026. Cuenta jdiaz, sin administrador.

## Pendiente: dos controles enfocados de etapa 4A

Los 40 controles locales y seis estados 2D anteriores terminaron; no repetirlos.
El [informe de avance](implementation/stage4/review_20260923/Informe_avance_etapa_4A.md)
explica el efecto de las reacciones de borde y la resolución todavía insuficiente
alrededor del centro suprimido. Las anclas espectrales ya son suficientes.

Propósito: comparar δ=0,05 en 561 nodos y δ=0,1 en 2145 nodos, con el mismo perfil
y población. Se mantienen fijos los cuatro bordes mediante cargas conjugadas,
se compara la derivada discreta con la exacta del perfil y se calcula corriente
integrada por sección. No se evoluciona el tiempo ni se incorpora fotón.

Estimación: **40-65 minutos**, un proceso y un hilo; reservar **2 GiB RAM** y
menos de 100 MB de salida. Se basa en los 518 s medidos para 561 nodos. El coste
puede cambiar con la malla y las consultas espectrales; la ETA es orientativa.

```bash
cd /home/jdiaz/pysnspd
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/run_campaign.py \
  --plan docs/implementation/stage4/review_20260923/next_campaign_plan.json \
  --phase spatial --output-root tmp/stage4A_core_followup_20260923 --execute
```

No requiere activar un entorno ni crear sesiones. Puede ejecutarse desde la
terminal que ya se usa. Las barras muestran progreso, tiempo transcurrido y ETA
por sección. Se guardan identity.json, progress.jsonl, result.json y fields.npz
por caso, y summary.json al finalizar. Un signo negativo queda registrado como
diagnóstico; una excepción detiene el lote sin reintentar ni aplicar recortes.

El directorio de salida debe ser nuevo; no se sobrescriben resultados. Al terminar,
basta avisar en el chat: los resultados se recuperarán de la carpeta indicada.

## Comandos ligeros opcionales

Verificar archivos actuales, resultados y entregas históricas (segundos, sin física):

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage4_core/verify_review.py
```

Reproducir sólo las restricciones algebraicas de los seis estados archivados
(segundos, cero consultas espectrales; elegir otro nombre si ya existe):

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage4_core/postprocess_completed.py \
  --output tmp/stage4_boundary_reproduction.json
```

Recalcular la referencia térmica Usadel de C.29-C.31 con los parámetros Korzh
(segundos, sin malla ni transiente; no es la respuesta a ocupaciones congeladas):

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage4_core/linear_usadel_reference.py \
  --output tmp/stage4_linear_usadel_reproduction.json
```

## Historial íntegro

Las entradas anteriores se conservan sin alteración en
[la libreta anterior](implementation/stage4/review_20260923/previous_delivery_exact/docs/GEMINGA_COMMANDS.md).
Desde ella se accede a todas las libretas previas. Los dos lotes iniciales de
etapa 4 ya terminaron y sus comandos son históricos.

## Política de cómputo

No ejecutar desde el agente cálculos razonablemente previstos de más de cinco
minutos. Registrar propósito, comando, salidas y recursos, entregar el comando
en el chat y esperar a que el usuario ejecute y avise. No usar trabajos de fondo,
sondeos o nuevas sesiones para eludir esta entrega. Un chequeo de duración
incierta usa timeout de 240 s: si vence, queda incompleto y se entrega para
corrida manual, sin reinicios sucesivos ni subdivisiones para eludir el límite.
