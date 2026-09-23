# Geminga: comandos vigentes

Actualizado el23 de septiembre de2026. Cuenta jdiaz, sin administrador.

## Etapa4 iniciada: dos lotes manuales sin fotón

La [preparación](implementation/stage4/start_20260923/README.md) incluye177 pruebas
y79 subpruebas aprobadas en Geminga, y un piloto local de145 s. Las etapas2–3
y la investigación3.5 permanecen cerradas dentro de sus alcances; no repetirlas.
Se investigan núcleo, estabilidad y disipación estática/instantánea, todavía
sin trayectorias, fotón, reservorio ni circuito en este lote.

### 1. Controles locales:40 casos

Propósito: variarδ, amplitud, gradiente y población; comparar D.36, fuerza,
corriente conjugada y tres movilidades reutilizando las mismas consultas.
Estimación:5–20 minutos;1proceso/1hilo; reservar1 GB RAM y <100 MB de salida.

```bash
cd /home/jdiaz/pysnspd
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/run_campaign.py \
  --plan docs/implementation/stage4/start_20260923/campaign_plan.json \
  --phase local --output-root tmp/stage4A_local_20260923 --execute
```

### 2. Estados espaciales2D:6 casos

Propósito: perfiles impuestos suave/suprimido en dos mallas, tresδ y tres
movilidades; obtener mapas y balances instantáneos. Los seis estados no son
seis transientes. Estimación45 minutos–3 horas;1proceso/1hilo; reservar1 GB RAM
y <100 MB de salida. La estimación usa el piloto local, no una corrida larga
en Geminga. No se requiere1D ni una transferencia fotónica.

```bash
cd /home/jdiaz/pysnspd
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/run_campaign.py \
  --plan docs/implementation/stage4/start_20260923/campaign_plan.json \
  --phase spatial --output-root tmp/stage4A_spatial_20260923 --execute
```

Los lotes son independientes. Pueden ejecutarse en orden o en terminales
separadas. No requieren activar un entorno: se usa el Python explícito.
Los dos muestran barras, tiempo transcurrido y ETA por sección; el promedio
de tareas es orientativo porque sus costes difieren.

Salidas: identidad y hashes, progresoJSONL, un resultado por caso, resumen;
el lote2D además guarda `fields.npz`. Un símbolo negativo queda como resultado,
sin intentar evolucionarlo. Una excepción conserva la evidencia y detiene
el lote sin reintentar. Los directorios existentes se rechazan: no se sobrescribe.
Al terminar basta avisar; los datos se recuperan de estas carpetas sin pegar logs.

## Comprobación ligera opcional

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage4_core/verify_delivery.py
```

Sólo lee archivos y verifica la entrega y su historial. Tarda segundos y
no ejecuta espectros ni transientes. Añadir `--execute` al corredor anterior
es la acción deliberada que inicia física; sin él sólo describe el plan.

## Historial preservado

La libreta de cierre3.5 se conserva íntegra en
[start_20260923/previous_delivery_exact/docs/GEMINGA_COMMANDS.md](implementation/stage4/start_20260923/previous_delivery_exact/docs/GEMINGA_COMMANDS.md).
Desde ella se accede a todas las entradas previas. Sus estados de no-pendiente
o de etapas abiertas son históricos; los dos comandos anteriores son la cola actual.

## Política para futuros cálculos

Un cálculo conocido o razonablemente previsto de más de cinco minutos se
prepara para ejecución del usuario. Su entrada incluirá propósito, comando
exacto, directorio nuevo, salidas, duración y memoria estimadas. El mismo comando
se incluirá explícitamente en un bloque copiable del chat; no se remitirá al
usuario únicamente a esta libreta.

Los corredores mostrarán progreso por tareas o pasos útiles, tiempo transcurrido
y ETA cuando exista una estimación. Una ETA no es una garantía de duración.
Después de entregar un cálculo largo se espera a que el usuario lo ejecute y
avise «reinicia» o equivalente. No se usan trabajos de fondo, sondeos ni sesiones
para eludir esa entrega.

Los chequeos ligeros de duración incierta usan un timeout acotado de 240 s.
Si vence, el intento queda incompleto: se conserva su evidencia y se entrega
el comando sin límite para ejecución del usuario. No se reinicia repetidamente
ni se divide un cálculo largo para eludir el límite.
