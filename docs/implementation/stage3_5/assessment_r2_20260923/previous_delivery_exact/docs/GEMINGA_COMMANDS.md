# Geminga: comandos vigentes

Actualizado el 23 de septiembre de 2026. Cuenta `jdiaz`, sin administrador.

## Estado actual: investigación 3.5; sin cálculo largo pendiente

La [revisión vigente](implementation/stage3_5/CURRENT.md) investiga las127
variables del inventario. El usuario eligió formación del hotbelt y latencia
relativa775/1550nm en el hilo80nm de Korzh; tiempo desde transferencia; referencia
material ajustada D=0,5cm²/s y608Ω por cuadrado. Después de revisar Allmaras y
Zotova–Vodolazov, eligió **caracterizar primero la cascada y dejar el ancho
gaussiano abierto**. No se adoptan5–20nm ni la equivalencia histórica1,4–1,9nm.

El [informe de investigación](implementation/stage3_5/research_20260923/Informe_investigacion_etapa_3_5_r1.md)
contiene resultados documentales y figuras analíticas. El
[plan de cascada](implementation/stage3_5/research_20260923/cascade_characterization.md)
identifica los datos de transferencia que faltan. No se ejecutaron transientes
nuevos, no se modificaron parámetros físicos y no hay un lote largo por lanzar.
No repetir las corridas de etapas2–3. El desarrollo3 permanece cerrado dentro
del alcance acordado, con requisitos dinámicos aún pendientes antes de4–5.

## Verificación ligera reproducible

```bash
cd /home/jdiaz/pysnspd
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  sandbox/stage3_5/research_20260923/verify_delivery.py
```

Comprueba integridad,127 entradas, decisiones recibidas y entregas históricas.
Lee archivos, sin resolver espectros ni trayectorias. Duración prevista:
segundos;1proceso/1hilo; memoria pequeña frente a una simulación. La salida es
JSON en terminal, con dictamen y conteos. El verificador anterior aislado usa
su libreta histórica; el actual aplica esa copia preservada al revisar la cadena.

Para regenerar las figuras analíticas y el registro, desde la raíz del repo:

```bash
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_5/research_20260923/geometry_scales.py
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_5/research_20260923/profile_comparison.py
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_5/research_20260923/build_range_ledger.py
```

Cada comando tarda segundos, usa1proceso y genera los JSON/PNG documentados en
`docs/implementation/stage3_5/research_20260923/`. Son fórmulas y contabilidad;
no caracterizan por sí solos una cascada. Regenerar PNG con otra versión de
Matplotlib puede cambiar sus bytes: conservar la entrega publicada y generar
una revisión antes de actualizar su manifiesto. Estos comandos son opcionales,
no una ejecución requerida al usuario.

La libreta que acompañó el cierre de etapa3 se conserva exactamente en
[research_20260923/previous_delivery_exact/docs/GEMINGA_COMMANDS.md](implementation/stage3_5/research_20260923/previous_delivery_exact/docs/GEMINGA_COMMANDS.md),
incluidos los resultados y comandos originales. La cadena sigue verificable.

## Lectura reproducible de la cascada histórica de Allmaras

```bash
cd /home/jdiaz/pysnspd
env PYTHONPATH=/home/jdiaz/pysnspd/tmp/stage3_5_document_deps \
  /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_5/research_20260923/digitize_a20_radial.py
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_5/research_20260923/plot_cascade_shape.py
```

Entrada: PDF de Allmaras2020 en `tmp/pdfs/modelo_v0_3/sources/Allmaras_Thesis_Final.pdf`,
con SHA verificado por el script. Extrae4 radios visibles de la figura2.7(a),
registra el margen de lectura y genera `cascade_digitization.json` y
`figures/04_cascade_shape.png` en la revisión. No usa zonas tapadas por la leyenda.
La dependencia pdfplumber0.11.8 está instalada sólo en
`tmp/stage3_5_document_deps`; el comando la activa sin modificar el entorno
físico. Su instalación aislada tardó unos3s. La extracción reproducida en
Geminga coincidió exactamente con el JSON local.
Tiempo previsto: segundos,1 proceso; sin simulación de cascada. Son comandos
opcionales. Regenerar figuras puede cambiar bytes según bibliotecas; revisar
el manifiesto antes de publicar otra versión.

## Historial preservado

La libreta anterior, incluido el comando del lote que ya terminó, se conserva
íntegra en
[implementation/stage3/closure_20260923/GEMINGA_COMMANDS_before_closure.md](implementation/stage3/closure_20260923/GEMINGA_COMMANDS_before_closure.md).
Sus enlaces mantienen accesibles las entradas previas. Las instrucciones de
«pendiente» o «siguiente cálculo» de esos archivos son históricas, no una cola
actual de trabajo. El solver de producción y `v1.0.0` no cambian.

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
