# Geminga: comandos presentes de pysnspd 1.1.3

Sólo los comandos de esta versión son operativos. Conservar el historial anterior en /home/jdiaz/GEMINGA_COMMANDS.md y docs/GEMINGA_COMMANDS.md. No reanudar checkpoints1.1.2 bajo física1.1.3. Todos los destinos deben ser nuevos; si existen, cambiar el sufijo de forma coherente y conservar los anteriores.

## Preparación

Usar jdiaz, verificar que finalizaron las continuaciones diagnósticas y comprobar la carga de todos los usuarios. El asignador cuenta los runners registrados05/08/09/10; su detección no garantiza que trabajos ajenos o scripts diagnósticos estén libres. Nunca ejecutar otra campaña sobre sus núcleos activos.

~~~bash
cd /home/jdiaz/pysnspd
set -o pipefail
export PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache
export NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache
export MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib
export TMPDIR=/home/jdiaz/scratch/tmp/pysnspd
export NUMBA_NUM_THREADS=6
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
mkdir -p "$PYTHONPYCACHEPREFIX" "$NUMBA_CACHE_DIR" "$MPLCONFIGDIR" "$TMPDIR"
~~~

Normal vigente: electron temperature y thermal-depairing; quadratic y bath son controles identificados. El coste no está garantizado:0–8 ps proyectó11,18/12,43 h por100 ps incluyendo preparación/JIT; pasos posteriores pueden ser más costosos. No hay promesa de100 ps en diez horas.

## Campaña normal: anclajes primero

Plan nuevo de12 casos, sin evolucionar todavía:

~~~bash
/home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/10_thermal_stiffness_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_3_20261017_01 \
  --ceiling-ps 100 --max-worker-cores 12 --prepare-only
~~~

Dos anclajes simultáneos775/1550 nm,21,5 µA, ancho80 nm, escape20 ps, eta0,8:

~~~bash
/home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/10_thermal_stiffness_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_3_20261017_01 \
  --ceiling-ps 100 --max-worker-cores 12 --stage anchors
~~~

El asignador conserva la asignación verificada de hasta seis hilos por anclaje,12 núcleos de trabajo más coordinador y al menos dos núcleos físicos completos reservados. No se promovió el cambio de afinidad L3. Parámetros, fuentes y asignación efectiva quedan en el plan/recibo.

Los controles sólo se habilitan después de ambos anclajes completos y auditados:

~~~bash
/home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/10_thermal_stiffness_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_3_20261017_01 \
  --ceiling-ps 100 --max-worker-cores 12 --stage controls
~~~

## Un caso explícito, alternativa a la campaña

Ejecutar este caso sólo cuando haya terminado la campaña o no esté usando sus recursos. Todos los parámetros esenciales se hacen explícitos; length160 nm y eta0,667 siguen siendo defaults heredados de pipeline09 hasta cualquier decisión final diferente.

~~~bash
/home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/09_thermal_stiffness_photon_run.py \
  --output-dir /home/jdiaz/scratch/photon_runs/thermal_v1_1_3_775_eta0800_20261017_01 \
  --wavelength-nm 775 --length-nm 320 --spacing-nm 5 \
  --gradient-closure thermal-depairing --stiffness-temperature electron \
  --current-uA 21.5 --sheet-resistance-ohm 425 --target-depairing-uA 25.866666666666667 \
  --bath-K 0.9 --tau-escape-ps 20 \
  --tau-ee-Tc-ps 5 --tau-ep-Tc-ps 24.7 \
  --sigma-nm 7.549834435 --retained-fraction 0.8 \
  --xi-cutoff 32 --xi-order 6 --phonon-nodes 96 --transport-order 48 \
  --dt-ps 0.02 --step-tolerance 2e-7 --local-energy-tolerance 2.5e-6 \
  --minimum-adaptive-step-ps 1e-6 --max-iterations 60 --max-rejections 50 \
  --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5 --ceiling-ps 100 \
  --snapshot-every-ps 0.1 --checkpoint-every-ps 1 --checkpoint-wall-s 60 \
  --threads 6 --cpu-start 0
~~~

Para1550 nm, cambiar únicamente wavelength a1550 y el nombre del directorio a1550. Si se quieren ejecutar ambos simultáneamente, usar la campaña10 para compartir presupuesto; no lanzar dos procesos independientes con el mismo cpu-start.

La parada es primer cruce positivo0,1 mV, confirmación0,5 ps y cruce+5 ps; sin cruce queda censurado a100 ps. Si exige un paso inferior a1e−6 ps, termina por resolución insuficiente, sin asignar una latencia. El reloj comienza postcascada. Ni pico completo ni reset eléctrico total forman parte de esta orden.

## Visualización: un PDF por caso

Sólo con fuentes de reconstrucción compatibles con los manifiestos. Un hash diferente exige el productor archivado correcto; no desactivar la verificación. Este comando produce reportes individuales, no compara corridas:

~~~bash
env NUMBA_NUM_THREADS=1 taskset -c 0,1,2,3,16,17,18,19 \
  /home/jdiaz/.conda/envs/snspd/bin/python plot_pipelines/07_case_report.py \
  --run-dir \
    /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_3_20261017_01/runs/anchor425_eta0800_775 \
    /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_3_20261017_01/runs/anchor425_eta0800_1550 \
  --snapshots-ps 1 3 5 10 20 40 \
  --output-dir /home/jdiaz/scratch/case_reports/thermal_v1_1_3_20261017_01 \
  --workers 2
~~~

Exactamente seis slots: tiempos guardados más cercanos; se truncan extras y se repite el último solicitado si faltan. Para un caso todavía RUNNING añadir --allow-incomplete; el PDF será provisional y capturará un checkpoint atómico. Los terminales diagnósticos especiales deben tener un productor compatible, no reclasificarse como una detección normal. Vout se muestra como señal intrínseca sin filtro/amplificador.

## Reanudación

Tras una interrupción controlada, repetir el mismo comando, directorio y parámetros añadiendo --resume. Exige fuentes, DC, geometría y parámetros idénticos; no redeposita el fotón. No editar fuentes físicas de procesos activos, no reabrir casos completos y no mezclar versiones.
