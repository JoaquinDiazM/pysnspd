## 2026-09-24 — Unión dinámica polarizada: ejecución pendiente

Este es el comando nuevo. No repetir Euler dual ni las campañas anteriores.
Construye dos referencias con corriente y luego cinco respuestas acopladas de
amplitud, fase, poblaciones y potencial sobre la malla dual, con el circuito
completo de la memoria. La comparación incluye un refinamiento independiente.

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/run_final_coupling.py \
  --plan docs/implementation/stage4/coupled_closure_20260924/final_coupling_plan.json \
  --output-root /home/jdiaz/scratch/stage4_final_coupling_20260924 --execute
```

**Tiempo orientativo:** una a varias horas; no es un plazo garantizado. Una
consulta de espectro/cinética con 13 perturbaciones tardó 1,92 s en la malla
completa sobre fondo uniforme. La corriente, las proximidades del gap y el
refinamiento aumentan ese coste. La ETA usa el avance de la fase activa.

**Recursos:** hasta 28 hilos de los 32 detectados, con dos núcleos físicos
reservados. Dos casos simultáneos comparten ese límite; nunca se asigna 90 % a
cada caso por separado. Un hilo BLAS/OpenMP por proceso. Reserva de RAM
aproximada de 34 GiB para dos casos, revisada contra la disponibilidad real;
disco esperado del orden de unos GB, dependiente de los registros espectrales.

**Salidas:** `references/cases/*/reference.npz`,
`responses/cases/*/results.json`, `coupled_fields.npz`,
`integrated_moments.npz`, registros por energía y `workflow_result.json`.
Todo queda bajo el destino del comando. Barras/ETA aparecen en la terminal,
y la salida detallada permanece en los archivos `logs/*.console.log` de cada
campaña. No hace falta copiar el registro al chat: basta avisar que terminó.

El destino debe ser nuevo. Si falla un caso, se guardan sus resultados ya
obtenidos, se cancelan sus propios trabajadores y los casos independientes
siguen. No hay reintento ni reinicio automático a partir de sumas parciales.
Una referencia incompleta sólo bloquea sus respuestas dependientes.

La integral de energía incluye las colas semiinfinitas mediante `x=C/E`;
el valor C separa las dos cuadraturas y no elimina estados de alta energía.
Se mantiene el margen práctico del 2 % en observables relevantes. El estado
de ejecución `COMPLETED` no declara por sí solo aprobado el cierre físico.

**Preflight opcional, segundos y sin cálculo físico:**

```bash
cd /home/jdiaz/pysnspd
/home/jdiaz/.conda/envs/snspd/bin/python sandbox/stage4_core/run_final_coupling.py \
  --plan docs/implementation/stage4/coupled_closure_20260924/final_coupling_plan.json \
  --output-root /home/jdiaz/scratch/stage4_final_coupling_20260924
```

**Diagnóstico ligero reproducible opcional:** una consulta sobre los 1712
nodos y una tangente térmica, sin recalcular el lote. Usar un destino nuevo.

```bash
cd /home/jdiaz/pysnspd
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 timeout 120 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u sandbox/stage4_core/pilot_biased_coupling.py \
  --output /home/jdiaz/scratch/stage4_coupling_pilot_reproduction
```

El piloto original ya terminó en 2,21 s y no necesita repetirse. El modelo
de producción y v1.0.0 no cambian; aún no se simula un fotón. El informe
de alcance y criterios está en
`docs/implementation/stage4/coupled_closure_20260924/README.md`.

También se completó el flujo real sobre una malla pequeña en 31,92 s, con
referencia aceptada y 1112 energías. Su comando y los planes exactos están en
`docs/implementation/stage4/coupled_closure_20260924/pipeline_smoke/integration_receipt.json`.
Es un diagnóstico de integración ya terminado; no es necesario repetirlo
antes del lote principal. Para una reproducción opcional debe elegirse otro
`--output-root` y conservar la afinidad indicada en ese recibo.
