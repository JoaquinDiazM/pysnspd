# Geminga: comandos vigentes

Actualizado el 23 de septiembre de 2026. Cuenta `jdiaz`, sin administrador.

## Estado actual: sin cálculo largo pendiente

El usuario autorizó el cierre de desarrollo de la etapa 3 y la apertura de
investigación 3.5. El lote de seis instantáneas mixtas con carga de reservorio
ya terminó en **355,651 s**. **No repetir ese lote ni los lotes históricos.**
La regresión focalizada aprobó 243 pruebas y 46 subpruebas en 10,29 s.

La [guía de cierre](implementation/stage3/closure_20260923/README.md) y el
[informe final](implementation/stage3/closure_20260923/Informe_cierre_desarrollo_etapa_3_y_apertura_3_5.md)
([PDF](../output/pdf/implementation/Informe_cierre_desarrollo_etapa_3_y_apertura_3_5.pdf))
registran el alcance: energía y balances instantáneos, con sensibilidad de malla.
No se ha completado el contrato temporal original, toda D.27, la interfaz
cinética ni una admisión de producción. Los datos originales no se reclasifican.

La [secuencia vigente](implementation/SECUENCIA_VIGENTE.md) mantiene el transiente
débil y la interfaz cinética como requisitos antes de etapas 4–5. La
[etapa 3.5](implementation/stage3_5/README.md) investiga parámetros, fuentes
Korzh/Allmaras y límites de aplicación; no ejecuta barridos ni adopta rangos.
El inventario tiene 127 entradas en 15 familias. L2D/W = 1,5–6 es sólo un
ejemplo, y la reducción transversal 1D sigue siendo condicional.

## Verificación ligera reproducible

Desde la copia sincronizada de `main`:

```bash
cd /home/jdiaz/pysnspd
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  sandbox/stage3_spatial/closure_20260923/verify_delivery.py
```

Propósito: verificar la integridad y coherencia de la entrega y sus referencias.
Lee los resultados guardados; no construye espectros, no evalúa una trayectoria
y no repite el lote físico. Recursos previstos: un proceso, un hilo, memoria
pequeña frente a una simulación; duración de segundos, condicionada a la lectura
de archivos. La salida de terminal informa el dictamen y las discrepancias.

La auditoría aritmética de las seis instantáneas ya está conservada en
`docs/implementation/stage3/closure_20260923/full_reservoir_review.json`
y su explicación en `full_reservoir_review.md`. Su fuente es
`sandbox/stage3_spatial/closure_20260923/audit_full_reservoir.py`.
Ese script histórico escribe con exclusión de archivos existentes: **no es un
comando de repetición sobre esta entrega**. Para la revisión cotidiana se usa
el verificador anterior, sin sobrescribir la auditoría.

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
