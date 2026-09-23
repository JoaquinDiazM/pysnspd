# Geminga: comandos vigentes

Actualizado el 23 de septiembre de 2026. Cuenta `jdiaz`, sin administrador.

## Pendiente: etapa 3, malla mixta con carga de reservorio

**Ejecución manual; el agente no ha lanzado este lote largo.**

Los seis casos abiertos anteriores terminaron en 4,52 minutos. No repetirlos.
El piloto mixto de dos perfiles terminó en 78,34 s: fuerza, corriente, potencial,
KWT, depósito de calor y circuito de la memoria tienen identidades consistentes.
El cálculo espectral conserva el solver causal y usa una semilla para acelerarlo.

El piloto mostró que los extremos libres dominan el calor de la hélice. El
siguiente lote añade un diagnóstico con carga prescrita del reservorio,
amplitud terminal fija y su trabajo energético explícito. La corriente inyectada
y los estados del circuito se hacen coherentes con esa carga. Son instantáneas:
todavía no se integra un transiente ni se acredita la interfaz cinética completa.

### Comando que falta ejecutar

```bash
cd /home/jdiaz/pysnspd
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage3_spatial/coupled_20260923/run_reservoir_batch.py \
  --output-root tmp/stage3_reservoir_full_20260923 \
  --execute
```

Tiempo previsto: **8-15 minutos**, sujeto a la caché y a la carga de la máquina.
Estimación a partir del piloto: 55 cuadraturas por perfil frente a 55/171/595 en
las tres mallas; el contraste por diferencias finitas se hace sólo en la malla
pequeña perturbada. Recursos previstos: un proceso CPU, un hilo, menos de 2 GB
RAM y menos de 100 MB de salida. No se lanza automáticamente por superar 5 min.

La primera fase calcula las seis instantáneas del material con la misma fuente
verificada. Imprime barras, tiempo transcurrido y ETA por tareas y cuadraturas;
la ETA inicial es aproximada y puede quedar desconocida si se supera la duración
estimada. La segunda fase aplica la carga del reservorio a esos datos guardados,
resuelve de nuevo el potencial y contabiliza el trabajo del borde. No repite
las consultas espectrales. El estado de extremos libres queda como control.

La finalidad es medir la sensibilidad espacial del calor y la respuesta con
carga prescrita, y la corriente normal en las aristas interiores próximas a los
extremos. Estas aristas no son una discretización certificada del flujo normal
externo. El calor del control libre no se usa como objetivo de convergencia.

Salidas en `tmp/stage3_reservoir_full_20260923/`:

- `mixed_control/`: campos, fuerzas, poblaciones, pruebas y registro de progreso.
- `reservoir/`: respuesta con carga prescrita, potencial, calor y balance con su trabajo.
- `manifest.json`: fuentes del corredor y registros de alcance.
- `summary.json`: éxito del lote completo y alcance, cuando ambas fases terminan.
- `failure.json`: fallo preservado si se interrumpe la secuencia.

El directorio debe ser nuevo. No hay reintento, cambio automático de método ni
sobrescritura. Cuando termine, basta indicar «terminó» o «reinicia»; el agente
recogerá directamente los archivos, sin pedir que se copie toda la terminal.

## Comandos ligeros opcionales

Verificar la entrega y describir el lote sin calcular nuevos estados:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_spatial/coupled_20260923/verify_delivery.py
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage3_spatial/coupled_20260923/run_reservoir_batch.py
```

Últimas líneas de progreso guardadas por la primera fase:

```bash
tail -n 6 /home/jdiaz/pysnspd/tmp/stage3_reservoir_full_20260923/mixed_control/progress.jsonl
```

## Historial y siguiente trabajo

La libreta anterior se conserva íntegra en
`docs/implementation/stage3/coupled_20260923/GEMINGA_COMMANDS_before_coupled.md`.
Sus enlaces mantienen accesible el historial previo; no repetir sus lotes.

Informe actual:
`output/pdf/implementation/Informe_empalme_y_acoplamiento_etapa_3_20260923.pdf`.
La etapa 3 sigue abierta. Tras revisar este lote sigue la evolución débil con
reservorios dependientes de la corriente, transporte cinético en la unión y
circuito de la memoria. Producción y `v1.0.0` permanecen sin cambios.

Todo cálculo previsto de más de cinco minutos se entrega al usuario. No se usan
trabajos de fondo ni sondeos automáticos para eludir esa entrega. Un chequeo breve
incierto usa timeout de 240 s; si vence, queda incompleto y se entrega su comando
sin reiniciarlo repetidamente ni dividirlo para eludir el límite.
