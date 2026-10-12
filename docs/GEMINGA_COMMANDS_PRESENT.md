# Geminga: photon-runs del modelo efectivo 1.1.0

Etapa 5 completada el 12 de octubre de 2026: sistema acoplado, balance independiente y cuatro controles numéricos admitidos. No quedan simulaciones activas de esta campaña. Este cuaderno operativo contiene órdenes futuras; el historial completo, incluidos intentos fallidos y cambios de precisión, permanece en /home/jdiaz/GEMINGA_COMMANDS.md.

Admisión: /home/jdiaz/scratch/stage5_final_evidence_completed_20261012/admission.json; SHA256 ed30c8525ecb3952e4042eca3b84e9ea0c2c78f41d7e0f17d880779b68747c4a.
Cobertura: /home/jdiaz/scratch/stage5_final_evidence_completed_20261012/finite_control_coverage_review.json; SHA256 2c404823c78594ae22fad71f1fee50e62b5d9ad3ee7beab89436b9bfd40252ac.

Los cuatro ejes pasan para 775 nm en 0–36,6 ps, incluyendo formación y cruces/confirmaciones hasta 33 ps. Todas las colas posteriores a 36,6 ps quedan fuera de la convergencia de cuatro ejes, aunque cada trayectoria original completa su propio cruce+5 ps y pasa sus auditorías individuales. El nominal 775 nm cruzó a 31,94131686165 ps y terminó a 36,94131686165 ps. El nominal 1550 nm completó 200 ps sin cruce: está censurado, con conservación y un contraste de paso tardío propio, sin atribuirle los cuatro refinamientos de 775 nm.

La referencia de producción exploratoria usa 54 estados electrónicos; 80 es un control espectral opcional. El modelo es neutral efectivo BCS/KWT, con material acústico efectivo y fuente postcascada; la admisión no es calibración experimental absoluta, jitter estadístico ni reset eléctrico completo.

## Preparación común

Ejecutar esta preparación en la terminal de geminga antes de las órdenes siguientes. Usar sólo jdiaz. Todos los destinos, temporales, cachés y plots permanecen en scratch.

~~~bash
cd /home/jdiaz/pysnspd
set -o pipefail
export PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache
export NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache
export MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib
export TMPDIR=/home/jdiaz/scratch/tmp/pysnspd
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=4
mkdir -p /home/jdiaz/scratch/photon_runs/logs
mkdir -p /home/jdiaz/scratch/tmp/pysnspd/pycache /home/jdiaz/scratch/tmp/pysnspd/numba_cache /home/jdiaz/scratch/tmp/pysnspd/matplotlib
~~~

Las simulaciones se muestran en secuencia: un proceso y cuatro hilos Numba sobre cuatro núcleos físicos distintos, BLAS/OpenMP a uno. El runner comprueba el presupuesto del 90 % y aplica su afinidad; no envolverlo en un taskset estrecho antes de esa comprobación. cpu-start es el índice inicial del grupo físico admitido. Comprobar además la carga de todos los usuarios: la asignación no reserva los núcleos frente a trabajos ajenos. En el anfitrión de 16 núcleos físicos se conservan CPU14/15 y sus hermanos libres; nunca superar 14 núcleos físicos entre todos los trabajos.

Los directorios de salida deben ser nuevos. Si existe un sufijo _01, cambiarlo de forma coherente en corrida, log, auditoría y ploteo; no borrar resultados. tee -a conserva cualquier log anterior. Las órdenes quedan en primer plano con checkpoints atómicos cada 1 ps o 60 s.

## Pipeline normal

775 nm, referencia exploratoria de 54 estados:

~~~bash
/home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/05_effective_photon_run.py \
  --output-dir /home/jdiaz/scratch/photon_runs/effective_775_20261012_01 \
  --wavelength-nm 775 --length-nm 320 --spacing-nm 5 \
  --xi-cutoff 32 --xi-order 6 --phonon-nodes 96 --transport-order 48 \
  --sigma-nm 7.549834435 --retained-fraction 0.667 \
  --current-uA 21.29720614489798 --bath-K 0.9 \
  --tau-ee-Tc-ps 5 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15 \
  --dt-ps 0.02 --step-tolerance 2e-7 --local-energy-tolerance 2.5e-6 \
  --max-iterations 40 --max-rejections 50 \
  --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5 --ceiling-ps 200 \
  --snapshot-every-ps 0.1 --checkpoint-every-ps 1 --checkpoint-wall-s 60 \
  --threads 4 --cpu-start 0 \
  2>&1 | tee -a /home/jdiaz/scratch/photon_runs/logs/effective_775_20261012_01.console.log
~~~

1550 nm, una vez finalizado el caso anterior:

~~~bash
/home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/05_effective_photon_run.py \
  --output-dir /home/jdiaz/scratch/photon_runs/effective_1550_20261012_01 \
  --wavelength-nm 1550 --length-nm 320 --spacing-nm 5 \
  --xi-cutoff 32 --xi-order 6 --phonon-nodes 96 --transport-order 48 \
  --sigma-nm 7.549834435 --retained-fraction 0.667 \
  --current-uA 21.29720614489798 --bath-K 0.9 \
  --tau-ee-Tc-ps 5 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15 \
  --dt-ps 0.02 --step-tolerance 2e-7 --local-energy-tolerance 2.5e-6 \
  --max-iterations 40 --max-rejections 50 \
  --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5 --ceiling-ps 200 \
  --snapshot-every-ps 0.1 --checkpoint-every-ps 1 --checkpoint-wall-s 60 \
  --threads 4 --cpu-start 0 \
  2>&1 | tee -a /home/jdiaz/scratch/photon_runs/logs/effective_1550_20261012_01.console.log
~~~

Estas órdenes fresh usan cap de paso 0,02 ps, tolerancia no lineal 2e-7 y guarda local 2,5e-6 desde t=0. La campaña admitida tuvo políticas de precisión por segmentos y guardas prospectivas; no se afirma identidad exacta con estas trayectorias nuevas ni se atribuye esa guarda retrospectivamente. Los rechazos adaptan el paso; se conserva por separado el criterio energético global del 1 %.

La observación usa Vout positivo, umbral 0,1 mV, confirmación 0,5 ps y final a cruce+5 ps; si no cruza, termina censurada a 200 ps. No se pide pico ni reset total y no se cambia la dirección global de corriente. t=0 es el estado efectivo posterior a la cascada, no el instante óptico absoluto.

τ_ee(T_c)=5 ps y τ_ep(T_c)=24,7 ps siguen la nomenclatura y los nominales libres de Allmaras (2020): son calibradores, no constantes materiales fijadas ni centros obligatorios de un rango. σ=7,549834435 nm, retención 0,667 y τ_escape=15 ps son hipótesis exploratorias documentadas.

--spacing-nm 5 es h_mesh, resolución espacial; el espesor de película de este caso es d_film=7 nm. Reducir h_mesh no adelgaza la cinta. En un barrido de espesor a N0 y D volumétricos fijos debe mantenerse Rsquare·d_film constante, pues N0 se infiere de esa combinación; mantener Rsquare mientras cambia el espesor cambia el material inferido.

Coste: planificar varias horas y seguir el progreso de la sesión, sin plazo garantizado. Como evidencia de escala, el tramo nominal 775 nm desde 22,435 hasta aproximadamente 36,92 ps empleó 4,12 h con cuatro hilos; el tramo 1550 nm 47,3–200 ps empleó 6,59 h con dos hilos y cap tardío 0,1 ps. Son segmentos de la campaña histórica, no tiempos totales ni estimaciones de estas órdenes fresh con guarda desde cero y cap 0,02 ps.

## Reanudación exacta

Tras una interrupción controlada, repetir el mismo comando completo, el mismo destino y todos sus parámetros, añadiendo --resume. El fotón no se deposita de nuevo. Se exige identidad de fuentes, parámetros, DC y malla. No reabrir un resultado completado ni reanudar después de editar el modelo.

Si se ejecutan casos simultáneos, asignar grupos físicos disjuntos, por ejemplo cpu-start 0 y 4, sólo tras comprobar el presupuesto global. Los comandos secuenciales anteriores no necesitan cambios de recursos.

## Control opcional de 80 estados a 775 nm

Este paquete aumenta simultáneamente corte/orden electrónicos, nodos fonónicos y cuadratura de transporte: 32/6/96/48 → 64/8/144/96. No es una verdad continua ni una nueva calibración; usa salida distinta. Ejecutar después de las corridas anteriores y repetir su auditoría y ploteo con esa ruta.

~~~bash
/home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/05_effective_photon_run.py \
  --output-dir /home/jdiaz/scratch/photon_runs/effective_775_spectral80_20261012_01 \
  --wavelength-nm 775 --length-nm 320 --spacing-nm 5 \
  --xi-cutoff 64 --xi-order 8 --phonon-nodes 144 --transport-order 96 \
  --sigma-nm 7.549834435 --retained-fraction 0.667 \
  --current-uA 21.29720614489798 --bath-K 0.9 \
  --tau-ee-Tc-ps 5 --tau-ep-Tc-ps 24.7 --tau-escape-ps 15 \
  --dt-ps 0.02 --step-tolerance 2e-7 --local-energy-tolerance 2.5e-6 \
  --max-iterations 40 --max-rejections 50 \
  --threshold-mV 0.1 --confirmation-ps 0.5 --extra-time-ps 5 --ceiling-ps 200 \
  --snapshot-every-ps 0.1 --checkpoint-every-ps 1 --checkpoint-wall-s 60 \
  --threads 4 --cpu-start 0 \
  2>&1 | tee -a /home/jdiaz/scratch/photon_runs/logs/effective_775_spectral80_20261012_01.console.log
~~~

## Pipeline visual

Ejecutar sólo sobre datos guardados y con presupuesto disponible. No evoluciona el modelo. La envoltura deja al planificador dos núcleos completos; éste selecciona [4,20], reserva [5,21] y fija el único proceso de render a CPU20, hermano lógico de CPU4. No añadir --cpu, opción inexistente. Todos los hilos de librerías siguen limitados a uno; NUMBA_NUM_THREADS=1 queda explícito para este proceso.

775 nm:

~~~bash
env NUMBA_NUM_THREADS=1 taskset -c 4,5,20,21 \
  /home/jdiaz/.conda/envs/snspd/bin/python plot_pipelines/05_plot_effective_photon_run.py \
  --run-dir /home/jdiaz/scratch/photon_runs/effective_775_20261012_01 \
  --output-dir /home/jdiaz/scratch/photon_runs/effective_775_20261012_01/figures_v1
~~~

1550 nm:

~~~bash
env NUMBA_NUM_THREADS=1 taskset -c 4,5,20,21 \
  /home/jdiaz/.conda/envs/snspd/bin/python plot_pipelines/05_plot_effective_photon_run.py \
  --run-dir /home/jdiaz/scratch/photon_runs/effective_1550_20261012_01 \
  --output-dir /home/jdiaz/scratch/photon_runs/effective_1550_20261012_01/figures_v1
~~~

La selección automática incluye inicio/final y tiempos realmente guardados. Si se añade --snapshots-ps, cada tiempo pedido debe pertenecer a la ventana real; pedir 200 ps a una corrida que terminó antes impide el ploteo. El destino figures_v1 debe ser nuevo.

Productos PNG y PDF: Vout/Vdev [mV], corrientes [µA], incrementos de energía [eV], defecto firmado/Eretenida [%] y defecto global absoluto/(Eretenida+Σ|ΔWext|) [%]; mapas del gap [meV], temperaturas equivalentes [K] y densidades energéticas [J/m³]. Si Vdev o el umbral aplana Vout, un recuadro conserva la lectura menos basal con escala propia (mV o µV); 331 nV equivalen a 0,331 µV y 0,000331 mV. El gap se reconstruye primero como interpolante P1 complejo espacial y después se toma su módulo. Las temperaturas equivalentes son observadores de energía de las distribuciones. La visual05 exige una fuente fotónica con energía retenida positiva; no es el ploteador de un DC-only.

## Auditoría de las nuevas corridas guardadas

~~~bash
env NUMBA_NUM_THREADS=1 taskset -c 4 /home/jdiaz/.conda/envs/snspd/bin/python \
  sandbox/stage5_effective_validation/audit_runs.py \
  --run photon775=/home/jdiaz/scratch/photon_runs/effective_775_20261012_01 \
  --run photon1550=/home/jdiaz/scratch/photon_runs/effective_1550_20261012_01 \
  --output-dir /home/jdiaz/scratch/photon_runs/audit_775_1550_20261012_01
~~~

Este auditor comprueba por corrida potencia instantánea, energía global, guarda local reconstruida, deposición, soporte, signo de corriente, eventos y checksums. Una corrida nueva no hereda automáticamente la admisión de cuatro ejes: para variaciones de geometría, material, sesgo o calibradores se revisa la sensibilidad pertinente.

La campaña cerrada conserva los resultados originales 2 % y 3/2/2. El temporal aplica 2 % (piso absoluto de 0,1 ps sólo para su error de latencia); malla/dominio aplican 5 % al campo L2/DC y déficit/pico real, con mínimo nodal diagnóstico. Para el par espectral 775 nm 54→80 se admitieron campo/déficit 5 %, Vout 5 % y latencia 2 %; el presupuesto de señal se decidió después de los datos y su fallo anterior al 2 % permanece visible. Energía global 1 %, soporte, identidad, guardas y cobertura no se relajaron. La reconstrucción temporal interpola la magnitud física por snapshot, sin cuerda compleja temporal ni eventos inventados.

Cada corrida conserva manifest.json, trajectory.npz, snapshots.npz, source_initial_state.npz, checkpoint_latest.npz, material y DC. Las distribuciones completas quedan en la fuente inicial y checkpoint; los hashes terminales autentican los productos. Ninguna orden de este cuaderno solicita jitter estadístico, latencia óptica absoluta ni reset completo.

## Fuentes de esta entrega

- Runner: pipelines/05_effective_photon_run.py; SHA256 70b7c5b5d39c74f19a35ad816f413ecc47451c46651817708549b8e9fc8e97aa.
- Dispatcher de auditoría: sandbox/stage5_effective_validation/guard_protocol_dispatcher.py; SHA256 13663d94b35184452a274a9bad73b395154af0e0dca224b370074239dd1ce50c.

El historial /home/jdiaz/GEMINGA_COMMANDS.md permanece íntegro. Este PRESENT y docs/GEMINGA_COMMANDS_PRESENT.md deben publicarse con contenido idéntico; los PRESENT históricos de stage5 conservan su alcance fechado.
