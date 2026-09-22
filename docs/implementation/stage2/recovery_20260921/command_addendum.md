
## 2026-09-21 — recuperación autorizada en screen code_000

El lote anterior se detuvo en `two_electron_4`: RK4 produjo una ocupación
negativa en la cola electrónica. Las trayectorias completadas y la evidencia
fallida se conservan. La bisección temporal local tampoco avanzó dentro del
presupuesto del diagnóstico. Esta continuación prueba SSPRK3 con un factor
común por evento y presupuesto compartido de partículas/huecos/fonones.
No recorta poblaciones ni modifica los kernels físicos. El limitador introduce
un error temporal que se registra y debe superar nuevos controles.

**Autorización específica:** el usuario pidió revisar, relanzar y dejar
desacoplada la screen `code_000`, sin esperar consumiendo uso de Codex.
Este relanzamiento es una excepción explícita a la entrega manual habitual
de cálculos mayores de cinco minutos. No implica permisos de administrador.

**Recursos estimados:** un hilo, hasta 12 GiB de memoria, aproximadamente 3–6 h;
la estimación temporal es incierta y depende de las mallas finas. No hay
reintentos ni continuación automática después de un fallo.

**Salidas nuevas:** `/home/jdiaz/pysnspd/tmp/stage2_recovery_20260921_ssp/`
contendrá trayectorias, condiciones iniciales, hashes, evaluaciones de criterios
y `batch_events.jsonl`. La consola completa estará en
`/home/jdiaz/pysnspd/tmp/stage2_recovery_20260921_ssp.screen.log`.
La finalización del lote no emite por sí sola la admisión de la etapa 2.

Comprobar el plan sin ejecutar física:

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage2_cells/recovery_20260921/run_limited_acceptance_batch.py --plan docs/implementation/stage2/recovery_20260921/limited_acceptance_plan.json --output-root tmp/stage2_recovery_20260921_ssp --dry-run
```

Lanzamiento autorizado, una sola vez, en una ventana nueva de la misma screen:

```bash
screen -S 3040582.code_000 -X screen -t stage2-recovery bash /home/jdiaz/pysnspd/sandbox/stage2_cells/recovery_20260921/launch_screen_batch.sh
screen -S 3040582.code_000 -X select stage2-recovery
screen -S 3040582.code_000 -X detach
```

Acceso manual y lectura del registro, sin copiarlo al chat:

```bash
screen -r code_000
# Ctrl-A D desacopla la screen sin detener el cálculo.
tail -n 35 /home/jdiaz/pysnspd/tmp/stage2_recovery_20260921_ssp.screen.log
screen -S 3040582.code_000 -p stage2-recovery -X hardcopy -h /home/jdiaz/pysnspd/tmp/stage2_recovery_screen_snapshot.txt
```

Diagnósticos ligeros reproducibles (la verificación temporal corta conserva sus
salidas y rechaza repetirlas sobre los mismos archivos):

```bash
cd /home/jdiaz/pysnspd
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 timeout --signal=TERM --kill-after=5 240 /home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage2_cells/recovery_20260921/run_short_checks.py
/home/jdiaz/.conda/envs/snspd/bin/python -m pytest -q sandbox/stage2_cells/recovery_20260921/test_limited_ssp.py sandbox/stage2_cells/recovery_20260921/test_limited_assessors.py
```

Si se detiene: conservar íntegra la carpeta y revisar el primer fallo; no
repetir a ciegas ni borrar resultados. Al terminar, indicar «reinicia» permite
recuperar los archivos directamente de Geminga y continuar la revisión.
