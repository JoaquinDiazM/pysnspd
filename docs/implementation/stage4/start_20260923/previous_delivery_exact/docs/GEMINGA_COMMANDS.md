# Geminga: comandos vigentes

Actualizado el 23 de septiembre de 2026. Cuenta `jdiaz`, sin administrador.

## Estado: investigación 3.5 cerrada; etapa 4 preparada sin fotón

La decisión recibida fue **«Cerrar investigación 3.5 y preparar etapa 4 sin fotón»**.
El [informe final](implementation/stage3_5/assessment_r2_20260923/Informe_cierre_investigacion_etapa_3_5.md)
y el [contrato de etapa 4](implementation/stage4/README.md) son la entrada vigente.
El ancho, reparto y reloj de transferencia de Korzh siguen abiertos. La etapa 4
comenzará con controles del núcleo, estabilidad y disipación sin fotón.

**No hay cálculo largo pendiente.** No repetir lotes de etapas 2–3. Las tareas
de implementación de etapa 4 están definidas pero aún no tienen un corredor
físico preparado. Por eso no se ofrece un comando de simulación ficticio.

## Verificación ligera de la entrega actual

```bash
cd /home/jdiaz/pysnspd
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  sandbox/stage3_5/assessment_r2_20260923/verify_delivery.py
```

Duración prevista: segundos; 1 proceso, sin RHS ni transientes. Comprueba hashes,
127 entradas, decisión recibida, preparación sin fotón y cadena histórica.
La salida es JSON en la terminal. Es opcional; la entrega ya incluye su ejecución.

## Reproducción opcional de las figuras y cuentas

```bash
cd /home/jdiaz/pysnspd
env PYTHONPATH=/home/jdiaz/pysnspd/tmp/stage3_5_document_deps \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  sandbox/stage3_5/assessment_r2_20260923/cascade_extract_partition.py
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_5/assessment_r2_20260923/cascade_plot_summary.py
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_5/assessment_r2_20260923/planning_scales.py
```

Extrae las figuras originales 2.6(a)/2.8(a) de Allmaras desde el PDF preservado
en `tmp/pdfs/modelo_v0_3/sources/Allmaras_Thesis_Final.pdf`, comprueba su hash y
genera las cotas y PNG bajo `docs/implementation/stage3_5/assessment_r2_20260923/`.
El último comando evalúa 36 estados algebraicos de movilidad y 16 combinaciones
auxiliares de longitud/ventana; no predice trayectorias. Duración: segundos por
comando; un proceso, memoria dominada por el PDF (~100 MB) y su extracción.
Las dependencias documentales están aisladas de las bibliotecas del solver.
Regenerar PNG con otras versiones puede cambiar bytes: conservar la entrega
publicada y revisar el manifiesto antes de sustituirla.

## Historial preservado

La libreta que acompañó r1 se conserva íntegra en
[assessment_r2_20260923/previous_delivery_exact/docs/GEMINGA_COMMANDS.md](implementation/stage3_5/assessment_r2_20260923/previous_delivery_exact/docs/GEMINGA_COMMANDS.md).
Incluye la reproducción anterior y los enlaces a todas las entradas previas.
Los estados de abierto o pendiente en ese historial no son una cola vigente.
La producción y `v1.0.0` se conservan.

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
