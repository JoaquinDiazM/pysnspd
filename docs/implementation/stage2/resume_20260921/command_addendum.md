
## 2026-09-21 - Etapa 2: referencia recuperada y lote de convergencia RK4

La referencia manual de code_000 terminó: 1164,55 s y 6843 RHS. Sus archivos
originales siguen en tmp/stage2_user_20260921_193252 y tienen copia verificada
en docs/implementation/stage2/resume_20260921/manual_reference. No repetirla.
RK4 con 10/20/40 pasos pasó la comparación corta: error máximo fino 1,57e-6,
reducciones 3,70/4,15 y 24,94 s frente a 1164,55 s. Los kernels y umbrales
siguen iguales. La etapa completa aún está pendiente.

### Ejecutar manualmente: lote único de validaciones pendientes

El agente NO ejecutó este lote. Está preregistrado en
docs/implementation/stage2/resume_20260921/manual_validation_plan.json.
Primero compara tres resoluciones temporales con una cuarta trayectoria de
paso reducido a la mitad, para una y dos celdas hasta t=2. Si pasan, evalúa
las mallas electrónicas 630/1260/2520 y fonónicas 1025/2049/4097 en el caso de
dos celdas, más equilibrio, fronteras de ocupación, campos y soporte.

Recursos estimados: un hilo de CPU, reservar 8 GiB de RAM, unas 2-3 horas
(aproximadamente 150 minutos; puede excederse). La red de 2520 estados domina
el coste. Los datos de cada trayectoria se guardan antes de continuar. No es
un transiente espacial de SNSPD ni una calibración de NbN.

Ejecutar en la terminal de screen, como jdiaz, desde el repositorio:

```bash
cd /home/jdiaz/pysnspd
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
STAGE2_BATCH=/home/jdiaz/pysnspd/tmp/stage2_validation_$(date +%Y%m%d_%H%M%S)
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage2_cells/resume_20260921/run_acceptance_batch.py --plan docs/implementation/stage2/resume_20260921/manual_validation_plan.json --output-root "$STAGE2_BATCH" --execute
```

El proceso permanece en primer plano; screen permite desconectar la terminal.
Hay progreso cada 100 evaluaciones y mensajes por tarea. Reutiliza el piloto
TWO40 existente sólo tras verificar parámetros, fuentes, versiones y hashes.
Se detiene ante el primer fallo; no reintenta ni sobrescribe resultados. Una
interrupción deja evidencia incompleta y requiere diagnóstico antes de seguir.

Salidas: batch_manifest.json, batch_events.jsonl, recibos de tareas, JSON,
NPZ y registros .console.log/.progress.jsonl en STAGE2_BATCH. Al terminar o
detenerse, escribir «reinicia» e indicar la screen si cambió. El agente puede
leer los archivos y la screen directamente, sin pegar el output al chat.
BATCH_COMPLETE sólo significa que las tareas listadas pasaron: todavía se
revisan actividad, estacionariedad y el dictamen final antes de cerrar etapa 2.

### Vista previa sin cálculos

No crea carpetas ni lanza subprocesos/RHS; valida hashes y muestra las tareas:

```bash
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage2_cells/resume_20260921/run_acceptance_batch.py --plan docs/implementation/stage2/resume_20260921/manual_validation_plan.json --output-root tmp/stage2_preview --dry-run
```

### Diagnósticos ligeros reutilizables

Suavidad del RHS: 17 evaluaciones estáticas, 2,65 s medidos, menos de 1 GiB.
No integra una trayectoria. Elegir una salida nueva:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout --signal=TERM --kill-after=5 240 /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage2_cells/resume_20260921/diagnose_rhs_smoothness.py --initial tmp/stage2_user_20260921_193252/one_reference_probe_initial.npz --output tmp/rhs_smoothness_recheck.json
```

Comparación continua seleccionada: 3,65 s medidos; sólo integra las referencias,
sin reconstruir redes ni trayectorias. La reproducción reemplaza su JSON de
evidencia y cambia su hash por los tiempos medidos; conservar antes una copia.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout --signal=TERM --kill-after=5 240 /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage2_cells/resume_20260921/close_selected_continuum.py
```
