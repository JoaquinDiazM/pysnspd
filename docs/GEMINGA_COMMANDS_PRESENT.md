Vigente desde2026-10-16: **v1.1.2**, sección final de este cuaderno. Techo sin cruce100 ps. Las entradas previas conservan historia y no son comandos actuales.

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

## Campaña de sensibilidad eléctrica — preparada el 12 de octubre de 2026

**No ejecutada.** Hasta16 casos distintos, en grupos que permiten detenerse tras
obtener información útil. No hace falta correrlos todos para utilizar el modelo
cerrado. Cada caso prepara su propio DC y mantiene ecuaciones, geometría,
circuito y perfil numérico fresh; cambia un solo input respecto de la referencia
de su misma longitud de onda. El sesgo se mantiene positivo.

| Grupo | Casos nuevos | Input físico |
|---|---:|---|
| bias — empezar aquí | 5 | Nominal775; 775 a0,65/0,75Idep; nominal1550; 1550 a0,75Idep |
| rates — calibradores | 6 | 775nm: τee(Tc)=2.5/10ps, τep(Tc)=12.35/49.4ps, τesc=7.5/30ps |
| source — deposición | 5 | 775nm: ηret=0.5/0.8, σ=5.338539125824206/10.677078251648414nm; 1550nm: ηret=0.8 |

rates incluye nominal775, source ambas referencias; si ya existen completas
en ESTA campaña, verifica identidad, todos sus SHA terminales y una auditoría
nueva antes de reutilizarlas. **No reutiliza los originales históricos como
si fueran fresh.** La unión de los tres grupos contiene16 simulaciones,
no19. El grupo all es opcional y solicita esa unión completa.

Los puntos son **sondas exploratorias escogidas para medir sensibilidad**;
no son intervalos experimentales, límites de calibración ni sus centros.
Los nominales libres siguen Allmaras: τee(Tc)=5ps y τep(Tc)=24.7ps.
El nominal de escape es15ps, retención0.667, σ=7.549834435nm y
I=21.29720614489798µA. La escala de sesgo usa la predicción bulk efectiva
Idep=30.36744305465µA de esta película de7nm:
I065=19.7388379855225µA e I075=22.7755822909875µA.
No se identifica con una corriente crítica experimental medida.
Un nuevo DC sin rama superconductora admisible detiene el caso y conserva
su evidencia.

### Órdenes normales

Ejecutar en geminga con la preparación común de este cuaderno. El directorio
sensitivity_20261012_01 debe ser nuevo al iniciar bias; los otros grupos usan
exactamente el mismo destino y --resume. Cambiar el sufijo en todas las órdenes
si ese destino ya pertenece a otra campaña.

~~~bash
cd /home/jdiaz/pysnspd
set -o pipefail
mkdir -p /home/jdiaz/scratch/photon_runs/logs

/home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/06_effective_sensitivity_campaign.py \
  --output-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01 \
  --group bias --execute \
  2>&1 | tee -a /home/jdiaz/scratch/photon_runs/logs/sensitivity_20261012_01_bias.log
~~~

Tras completar y auditar bias, continuar por separado si se desea medir
los calibradores:

~~~bash
/home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/06_effective_sensitivity_campaign.py \
  --output-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01 \
  --group rates --execute --resume \
  2>&1 | tee -a /home/jdiaz/scratch/photon_runs/logs/sensitivity_20261012_01_rates.log
~~~

Para la fuente inicial:

~~~bash
/home/jdiaz/.conda/envs/snspd/bin/python -u pipelines/06_effective_sensitivity_campaign.py \
  --output-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01 \
  --group source --execute --resume \
  2>&1 | tee -a /home/jdiaz/scratch/photon_runs/logs/sensitivity_20261012_01_source.log
~~~

Ante una interrupción, repetir el mismo grupo/comando, añadiendo --resume
también a bias. El fotón no se deposita de nuevo en un checkpoint.
Un lock exclusivo impide lanzar dos escritores/grupos simultáneamente
sobre esta campaña. No editar fuentes físicas ni el orquestador mientras
se usa su plan: sus SHA forman parte de la identidad. Si cambian, usar
un destino nuevo.

Perfil explícito: L=320nm, hmesh=5nm, película7nm, Tb=0.9K,
cuadratura54 estados/96 fonones/48 transporte; cap0.02ps, NL2e-7,
guarda local2.5e-6 desde cero, máximo40 iteraciones.
**No son los defaults genéricos del parser.** La campaña histórica admitida
conserva sus guardas por segmentos y no se declara idéntica a este perfil fresh.

Evento: Vout−Vout,DC positivo, umbral0.1mV, confirmación0.5ps,
histéresis5% y final en cruce+5ps; sin detección, techo200ps y
latencia censurada. No espera pico, FWHM, recuperación ni reset total.
El cero de tiempo es la transferencia efectiva postcascada, no la absorción
óptica absoluta. Mantiene el circuito completo de tres estados y la
inductancia exterior congelada desde cada DC.

Recursos: **un caso en primer plano, cuatro núcleos físicos, cuatro hilos
Numba, BLAS/OpenMP=1**. El runner repite inventario/afinidad/cuota/RAM y aplica
techo90% para cada caso. Auditoría y visual son procesos únicos.
Datos, cachés, checkpoints, logs y figuras van a scratch. Reservar núcleos
completos libres y comprobar carga ajena según la política del servidor.
Planificar varias horas por caso; algunos casos censurados a200ps con este
perfil estricto pueden requerir un día o más. No hay ETA calibrada ni plazo
garantizado; los16 pueden requerir varios días. El grupo bias reduce el
compromiso inicial a cinco casos. Progreso y ETA al objetivo vigente salen
por pantalla, con checkpoints cada1ps o60s. **No se ejecutaron estas simulaciones
al preparar el cuaderno.**

Productos: campaign_plan.json fija inputs/comandos/SHA;
runs/<caso>/ conserva material, DC, trayectoria, snapshots, fuente y
checkpoint completos. audits/<caso>_<fecha>/audit.json comprueba balance,
soporte, corriente, deposición, guardas y evento.
campaign_execution.json y aggregate_<grupo>.json autentican terminales y
auditorías. Un fallo detiene la secuencia y preserva lo anterior; no se presenta
como ausencia física de detección. Los nuevos inputs no heredan los cuatro
controles de convergencia de la referencia775: es exploración de sensibilidad,
con conservación individual y necesidad de contrastar precisión donde
cambie el régimen.

### Órdenes visuales comparativas

Ejecutar tras completar el grupo correspondiente. Leen estados guardados,
sin evaluar RHS ni lanzar simulaciones. --individual añade figuras físicas
y snapshots de visual05 para cada caso. Cada destino debe ser nuevo;
para repetir, cambiar visual_..._01 por un sufijo nuevo, preservando anteriores.

~~~bash
env NUMBA_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/06_plot_effective_sensitivity_campaign.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01 \
  --output-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/visual_bias_01 \
  --group bias --individual
~~~

~~~bash
env NUMBA_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/06_plot_effective_sensitivity_campaign.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01 \
  --output-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/visual_rates_01 \
  --group rates --individual
~~~

~~~bash
env NUMBA_NUM_THREADS=1 /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/06_plot_effective_sensitivity_campaign.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01 \
  --output-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/visual_source_01 \
  --group source --individual
~~~

Omitir --individual si sólo se desean comparaciones eléctricas rápidas.
El análisis entrega PNG/PDF de Vout menos basal[mV], con recuadro[µV]
para señales pequeñas; cada familia se corta en la ventana temporal común
de sus miembros. Las líneas de cruce sólo muestran eventos confirmados
guardados. metrics.csv y comparison_manifest.json incluyen:
latencia relativa[ps], diferencia respecto de referencia sólo si ambas
detectan, máxima diferencia Vout en ventana común[mV], mínimo Is[µA]
y pendiente local[mV/ps], ajustada sobre muestras aceptadas en cruce±0.25ps.
No asigna pendiente ni latencia ficticia a censurados.
Esa pendiente no incluye filtro, ganancia, ruido ni ancho de banda del setup;
no constituye una predicción de jitter.

Los máximos Vout se llaman **máximos observados en la ventana guardada**:
una corrida que termina en cruce+5ps puede seguir subiendo. No se usan como
pico total del pulso ni para comparar amplitudes a ventanas diferentes.
Los cruces representan un discriminador ideal en la carga50Ω; la comparación
absoluta con voltajes medidos requiere calibrar el readout real.

Interpretación crítica: variar τee(Tc) cambia conjuntamente BGK y movilidad
KWT; variar τep(Tc) cambia conjuntamente normalización acústica de colisiones
y movilidad KWT. Variar σ o ηret con inicialización FD/BE común también cambia
el reparto energético electrón/fonón; σ no es una sensibilidad geométrica pura.
La extrapolación efectiva KWT/rigidez en la zona caliente por encima de Tc
sigue siendo una hipótesis calibrable, identificada en Parametros.pdf.

--plan-only inspecciona un destino NUEVO sin lanzar hijos numéricos.
Si luego se ejecuta ese mismo destino, añadir --resume.
Las pruebas rápidas cubren los16 comandos con parser real05, censuras,
ventanas comunes, checksums/identidad, lock, memoria descomprimida y
coordinación de receipts con fixtures, sin evoluciones nuevas.


## 2026-10-13 — Case report individual: pipeline visual07

# Case report individual del modelo efectivo

Pipeline: /home/jdiaz/pysnspd/plot_pipelines/07_case_report.py

Produce **un PDF autónomo por run** del formato efectivo de pipeline05
(pysnspd.effective_photon_run.v1). Acepta varios directorios como comodidad de
ejecución: nunca superpone, resta ni compara casos. Los holds sin fotón y los
formatos históricos no se presentan con el reloj de una corrida fotónica.

## Contenido estándar

1. Parámetros exclusivos del caso: geometría, espesor y resolución separados;
   material, corriente DC, absorción/deposición, energía fotónica y retenida,
   calibradores Allmaras, circuito, discriminador y precisión numérica.
2. Observables de laboratorio: Vout menos basal DC, corriente de lectura
   Ib−Is, cruce, confirmación, latencia relativa, final guardado y pendiente
   local. También Ib, Is, Vdev, vc y potencia en la carga.
3. Estado interno: curvas de gap, energías y temperaturas equivalentes;
   snapshots del gap completo y de la región central; fase relativa al DC;
   Te, Tph, exceso energético, Js,x, Jn,x, |Jtotal| y potencial.
4. Balance independiente, precisión, procedencia, hashes y límites.

Con seis snapshots automáticos son **16 páginas A4 apaisadas**. Los mapas
comparten escala de color entre los tiempos del mismo campo. Se muestran
unidades físicas: nm, ps, meV, K, kJ/m³, A/m², V y A con prefijos legibles.

Vout no incluye filtro, amplificador, ruido ni ancho de banda experimental.
Su máximo es el máximo **observado**, no el pico de un pulso completo.
La pendiente no predice jitter. El reloj empieza en la transferencia efectiva
postcascada, no en la absorción óptica absoluta. Sólo un cruce confirmado
produce latencia; un caso censurado queda sin latencia.

## Comando para cuatro casos terminados

Ejecutar en geminga. Cada repetición requiere un destino nuevo.

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir \
    /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs/nominal775 \
    /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs/bias065_775 \
    /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs/bias075_775 \
    /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs/nominal1550 \
  --snapshots-ps 1 11 11.75 12.5 20 35 \
  --output-dir /home/jdiaz/scratch/case_reports/bias_20261013_04 \
  --workers 2
~~~

Para otro caso, sustituir la lista de --run-dir; uno solo también funciona.
--workers 2 paraleliza dos PDFs independientes, con un núcleo físico y
bibliotecas de un hilo por trabajador. Revisa afinidad, cuotas y memoria,
reserva núcleos completos y respeta el techo 90 %. Los dos ejemplos nominales
tardaron unos seis segundos por caso; no es una garantía para archivos mayores.

## Caso aún en ejecución

Añadir --allow-incomplete con un destino nuevo. Ejemplo:

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir /home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs/bias075_1550 \
  --output-dir /home/jdiaz/scratch/case_reports/bias075_1550_partial_20261013_01 \
  --allow-incomplete --workers 1
~~~

Captura una copia coherente de un checkpoint atómico; no detiene ni escribe
en la simulación. Todas las páginas quedan PROVISIONALES si no hay terminal
completo. El final de trayectoria y el último snapshot se indican por separado.
No convierte ausencia de cruce hasta el checkpoint en censura a 200 ps.
Si el caso ya terminó, usa y autentica sus productos terminales.

## Elegir instantes

--snapshots-ps 0 1 5 10 15 20 solicita las muestras guardadas más cercanas;
se incluyen siempre inicial y final. Máximo 24 instantes pedidos, dentro del
intervalo de snapshots existente. No se interpola dinámica en tiempo.

Sin esa opción, una detección muestra inicial, 1 ps, primera hotbelt guardada,
muestras próximas a cruce/confirmación y final. Sin detección: inicial,
1 ps, 10 ps, cuarto, mitad y final de la ventana, eliminando duplicados.
La hotbelt es conectividad del grafo nodal con |Delta|≤0.2|DeltaDC|
entre bordes transversales; no es el gatillo eléctrico ni prueba continuidad
entre muestras.

## Datos originales y reconstrucciones

Las curvas y mapas de gap/energía/temperatura se leen de trajectory.npz y
snapshots.npz, autenticados por manifest.json, junto con el DC y la malla
original. Se rechaza corrupción o una versión distinta de las tres fuentes
utilizadas para reconstruir corrientes.

El gap se interpola como campo complejo P1 antes de tomar su magnitud.
La fase se muestra en nodos, sin interpolar ángulos y ocultando amplitudes
menores a 0.2 del DC. Un salto de fase envuelta no prueba un vórtice.

Los flujos superconductores se recuperan de la derivada gauge exacta del
gradiente energético, independiente de las ocupaciones no guardadas.
El potencial y flujos normales se recuperan de la misma neutralidad algebraica,
con Is del estado aceptado exactamente coincidente con cada snapshot.
Se contrasta Vdev y la conservación de corriente con los valores guardados.

J por triángulo es una **aproximación de visualización**: mínimos cuadrados
ponderados de proyecciones sobre caras duales conservativas, divididas por
espesor físico y longitud dual. El recibo conserva residuos de ese ajuste.
Las flechas indican dirección con longitud visual fija de 6 nm; se omiten
magnitudes menores al 0.01 % del máximo entre los snapshots seleccionados.
El color aporta la magnitud o componente firmada. No se interpretan cambios
locales de signo como inversión de corriente global.

Te y Tph son observadores de energía equivalente, no estados térmicos
impuestos. La energía electrónica guardada representa cuasipartículas y
excluye condensación/gradiente; el balance total sí incluye esos términos.
No se inventan distribuciones completas ni mapas de colisiones o calor KWT
que los snapshots no guardaron. No se evalúa el RHS cinético ni un nuevo DC,
no se evoluciona el sistema y no se aplica un ajuste experimental.

## Productos y validación

En el destino scratch: case_report_<caso>_<identidad>.pdf, un recibo JSON
por PDF, case_reports_manifest.json como inventario y .build para capturas.
Nunca reemplaza destinos existentes. No modifica datos de entrada.

Validación del 13 de octubre de 2026: 20 pruebas rápidas; ejemplos nominales
775 y 1550 nm, detección y censura respectivamente, 16 páginas cada uno.
Revisión visual de todas las páginas, escalas y etiquetas; cotejo de hashes
y reconstrucción de Vdev con error máximo aproximadamente 1.4e-17 V.
Los 13 archivos físicos del checkpoint operativo permanecen sin cambios.
Los ejemplos están en /home/jdiaz/scratch/case_report_validation_20261013_release.


## 2026-10-13 - Case report refresh: English, fixed six snapshots, winding and model summary

Pipeline: /home/jdiaz/pysnspd/plot_pipelines/07_case_report.py
Guide: /home/jdiaz/pysnspd/plot_pipelines/CASE_REPORT.md

The current report is English, 19 pages, with larger chart labels, one PDF per
case, quantity-specific palettes, phase-winding vorticity and section6 containing
the coupled model equations and complete three-state circuit schematic.
This entry supersedes the earlier snapshot-selection behavior.

--snapshots-ps takes the first six requested times [ps], repeats the last if fewer
are supplied, preserves order/repetitions and selects nearest saved states.
No initial/final snapshots are silently added. Every panel labels requested and
actual time. Out-of-range requests select the nearest stored endpoint.
Omit the flag for six automatic targets. No new simulation or full reset.

~~~bash
cd /home/jdiaz/pysnspd
CASE_RUNS=/home/jdiaz/scratch/photon_runs/sensitivity_20261012_01/runs
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir "$CASE_RUNS/nominal775" "$CASE_RUNS/bias065_775" \
    "$CASE_RUNS/bias075_775" "$CASE_RUNS/nominal1550" \
  --snapshots-ps 0 1 5 10 20 35 \
  --output-dir /home/jdiaz/scratch/case_reports/english_20261013_01 \
  --workers 2
~~~

Choose a fresh output directory for every invocation. --allow-incomplete reads
one coherent checkpoint and labels the PDF provisional. It does not interrupt
the solver or infer terminal censoring from a partial window.
Resources: two independent foreground workers, one physical core each, BLAS/
OpenMP/Numba=1, shared ceiling90%, all outputs/caches/captures under scratch.
Two nominal reports required about six seconds per case. Larger archives can differ.

Vorticity is signed condensate phase winding, with P1 zero localization.
Gray marks undefined/insufficient-amplitude contours, with a diagnostic floor
1e-6 of local DC amplitude. It is not curl(J), reconstructed magnetic flux or
continuous vortex tracking. Primary method:
https://www.mcs.anl.gov/papers/P5324-0315.pdf, Eq.(1).
Receipt contains actual/requested times, repetitions, cadence, winding checks
and a follow-up acquisition queue when saved time resolution is insufficient.
No physical solver or save-rate code was modified. Vout still lacks the
experimental filter/amplifier/noise response.
Validation:31 quick tests, two19-page examples and an actual two-time/padded CLI run.

## 2026-10-14 - Paired base cases: eta, relaxation, preparation and domain

Prepared only; the physical photon campaign has NOT been executed.
Pipeline: /home/jdiaz/pysnspd/pipelines/08_effective_basecase_sweep.py
Plan: /home/jdiaz/scratch/photon_runs/basecases_eta_20261014_01/campaign_plan.json

Purpose: compare 775/1550 nm at the same 21.5 uA and tau_escape=20 ps,
without tuning parameters separately for each wavelength. Rs=608 ohm/square
is retained: this closure has Idep proportional to 1/Rs at fixed D,Tc,Tb,w,
so replacing608 by425 would increase Idep30.3674 to43.4433 uA.
425 is a conditional340/RRR conversion, not a cryogenic measurement.
No RHS, solver, physical source pin or default material parameter was changed.

20 cases /10 wavelength pairs:
- eta=.667 reference and requested .7,.8,.9,1.0 (10 cases).
- at eta=.8, tau_ee(Tc)=2.5 and10 ps (4 cases).
- at eta=.8, sigma=7.549834435/sqrt(2) and7.549834435*sqrt(2) nm (4 cases).
- at eta=.8, length=640 nm, with the same5 um device and exterior-inductance
  accounting (2 cases; reference length320 nm is already in the eta matrix).
All remaining parameters use the admitted exploratory54/96/48 quadratures,
h_mesh5 nm, tau_ep(Tc)24.7 ps, Tb.9 K, Tc8.65 K and D5e-5 m2/s.
The source is common FD/BE at transfer with fresh DC gap unchanged.
Eta is TOTAL film energy retention at transfer, conditional on photon absorption.
It is not optical absorption probability or electronic-only conversion efficiency.
Eta1 is an upper-bound source test; later20 ps phonon escape stays active.
Sigma changes both spatial preparation and the initial electron/phonon partition.
Tau_ee changes BGK and KWT mobility together. These points are hypotheses,
not measured uncertainty intervals or a new experimental calibration.

Foreground, up to13 independent runs, one Numba thread/distinct physical core
per run; BLAS/OpenMP1. Inventory checked2026-10-14:16physical/32logical,
no CPU quota, about118 GiB available RAM and945 GiB free scratch.
Workers occupy13/16=81.25% physical cores; one admitted core is reserved for
coordination, two complete cores stay outside the87.5% allocated core set.
Conservative scheduling reservation55 GiB including coordinator; rechecked
against90% available memory. This is allocation, not a utilization guarantee:
I/O, serial setup, and the final incomplete wave can lower instantaneous usage.
SMT logical-CPU percentages can be roughly half the physical-core percentage.
The script rejects requests below50% of visible physical cores, over-budget
requests and a pre-existing effective photon solver under jdiaz.
Do not run simultaneously with another heavy campaign without a shared budget.
Expected duration: multiple hours; no measured global ETA before executing
this new concurrent matrix. Cases ending at200 ps can dominate the duration.
Each running case logs its measured progress and current-horizon ETA.

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  TMPDIR=/home/jdiaz/scratch/tmp/pysnspd \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 BLIS_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python -u \
  pipelines/08_effective_basecase_sweep.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/basecases_eta_20261014_01 \
  --workers 13
~~~

Observation: positive baseline-subtracted Vout0.1 mV, confirmation at crossing+.5 ps,
end at crossing+5 ps. No crossing: censor at200 ps, latency unassigned.
Never requests full reset. No PDF/report generation.
Outputs in scratch: runs/<case>/{manifest.json,trajectory.npz,snapshots.npz,
checkpoint_latest.npz,dc/,material.*}; logs/; audits/; campaign_execution.json;
latency_pairs.csv (accepted-controller crossings and differences only when both
members detect; censored times are blank). Individual audits do not establish
new convergence or physical calibration. Failed cases remain marked failed;
independent cases continue and the campaign exits nonzero.

Resume the SAME command with --resume appended. Source/parameter/hash changes
are rejected; completed cases are authenticated and not simulated again.
SIGINT/SIGTERM stops new work and requests accepted-state checkpoints from
active physical solvers. No forced kill, background service or polling daemon.

Validation before handoff: all20 generated commands parse into the exact
expected runner parameter identities; all10 pairs differ only in wavelength;
50% minimum and global ceiling checked; synthetic foreground workers exercise
bounded concurrent slots, nonoverlap, exact completed-case resume and rejection
of changed source pins. Fixtures only under scratch/tmp, not scientific runs.
The13 physical runner pins equal the existing October12 campaign.
Static bulk derivation matches30.3674430546503 uA to floating precision;
low-temperature closed form30.38655886 uA differs by0.063%.

## 2026-10-15 - Individual reports for13 completed basecase runs

Read-only inspection found13 COMPLETE_AUDITED cases,5 RUNNING and2 FAILED before photon preparation (length640 DC Newton limit). Audit receipt hashes and terminal manifest hashes were checked. This command has NOT been executed. Produces one independent PDF per completed case, with six automatically selected snapshots adapted to the case duration. It excludes live and failed cases, preserves originals and uses a fresh scratch destination. Two plotting workers use cores12 and11 under the current full-affinity topology; active photon workers were on2,4,6,8,10. Recheck global allocation if launching other jobs. No new photon simulation, model change or full reset. Expected plotting cost: order minutes for13 cases, not a measured ETA.

~~~bash
cd /home/jdiaz/pysnspd
RUNS=/home/jdiaz/scratch/photon_runs/basecases_eta_20261014_01/runs
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir \
    "$RUNS/eta0667_775" "$RUNS/eta0667_1550" \
    "$RUNS/eta0700_775" "$RUNS/eta0700_1550" \
    "$RUNS/eta0800_775" "$RUNS/eta0800_1550" \
    "$RUNS/eta0900_775" "$RUNS/eta0900_1550" \
    "$RUNS/eta1000_775" "$RUNS/eta1000_1550" \
    "$RUNS/ee_half_775" "$RUNS/ee_half_1550" \
    "$RUNS/ee_double_775" \
  --output-dir /home/jdiaz/scratch/case_reports/basecases_completed_20261015_01 \
  --workers 2
~~~


## 2026-10-15 — 1.1.0: rigidez térmica candidata y difusión calibrada

Estado: v1.0.0 publicada en 7992f6a; thesis-v1.0.0 preservada. Revisión 1.1.0 preparada, sin admisión fotónica nueva. Dos pruebas de 0.1 ps alcanzaron su límite de 240 s: 775 nm llegó a 0.051875 ps y 1550 nm a 0.046875 ps. Conservaron soporte, capacidad positiva y balance. No completaron la duración ni aportan latencias. No se extrapola su tasa inicial a 200 ps.

Caso principal: 21.5 uA, ancho 80 nm, espesor 7 nm, longitud resuelta 320 nm, malla 5 nm, baño 0.9 K, Tc 8.65 K, escape fonónico 20 ps y sigma 7.549834435 nm independiente de D. Rs=425 ohm/square principal y 608 de control; Dcal=1.4103757147386595 y 0.689136928331325 cm2/s respectivamente, con objetivo Idep=25.8666666667 uA. N0 y ell0 se recalculan por Einstein. El cierre usa G(Delta,Te) libre, U=UBCS+G-sum(Te*GT), calor Q+PJ+Te*dotGT y Ceff>0 sin piso. Eta=0.8 es exploratoria, no ajustada al retraso experimental. Circuito completo heredado; Vout intrínseco sin filtro/amplificador.

Plan: /home/jdiaz/scratch/photon_runs/thermal_basecases_20261015_01/campaign_plan.json. Hay 12 casos, seis pares 775/1550 nm: anclajes; rigidez a Tb; Rs=608 con D calibrada; eta=0.667; eta=0.9; tauEE(Tc)=2.5 ps. Los diez controles requieren ambos anclajes completos y auditados. Una censura a 200 ps no recibe latencia ficticia; un fallo numérico o de capacidad bloquea los controles. Completar el catálogo no establece convergencia ni validación experimental.

Recursos: el ejecutor revisa topología, afinidad, cuotas y memoria. Usa hasta 12 núcleos físicos para workers, un coordinador y el límite global del 90%, contando procesos 05/08/09/10 activos y reservando dos núcleos completos. Las dos anclas usan entre 1 y 6 hilos por caso según disponibilidad; los controles usan hasta diez procesos de un hilo en cola. BLAS/OpenMP queda en 1. Datos, cachés y logs permanecen en scratch. Las corridas pueden tardar horas o días; cada caso informa progreso y ETA de su sesión actual. Checkpoints cada 60 s de pared o 1 ps. Ningún comando simula el reset eléctrico completo.

### A. Anclajes: ejecución larga pendiente

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  pipelines/10_thermal_stiffness_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/thermal_basecases_20261015_01 \
  --stage anchors --max-worker-cores 12
~~~

### B. Controles: después de A auditada

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  pipelines/10_thermal_stiffness_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/thermal_basecases_20261015_01 \
  --stage controls --max-worker-cores 12
~~~

Para reanudar el mismo stage, añadir --resume; deben conservarse fuentes y parámetros. Un fallo de integración/capacidad requiere diagnóstico, no reinicio ciego. Si cambian las fuentes identificadas por hash, usar un directorio nuevo y preservar los resultados anteriores.

### C. Reportes individuales de los anclajes

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir \
    /home/jdiaz/scratch/photon_runs/thermal_basecases_20261015_01/runs/anchor425_eta0800_775 \
    /home/jdiaz/scratch/photon_runs/thermal_basecases_20261015_01/runs/anchor425_eta0800_1550 \
  --snapshots-ps 0 1 5 10 20 35 \
  --output-dir /home/jdiaz/scratch/case_reports/thermal_basecases_20261015_01 \
  --workers 2
~~~

El reporte 07 reconstruye las corrientes con el cierre y la Te guardada de cada caso. Produce un PDF por corrida y elige los seis tiempos guardados más próximos a los solicitados.


## 2026-10-15 — v1.1.1: corrección Pauli y nuevos anclajes

Estado: bloqueo continuo de 1.1.0 identificado y corregido. Sus dos anclajes se detuvieron con aprobación expresa y conservaron checkpoints; no reanudar con ecuaciones distintas. Las pruebas frescas a malla nominal completaron 0.2 ps: 775 nm, 50 pasos/175.14 s; 1550 nm, 48 pasos/169.81 s. Defectos energéticos máximos 4.99e-25 y 2.78e-25 J. No son latencias ni validación de toda la ventana. Guardias, tolerancias, circuito y difusión no se relajaron.

El calor electrónico usa actividades f(1-p) al calentar y p(1-f) al enfriar, con normalización exacta de energía y producción entrópica adicional no negativa. Rigidez Te, Dcal y Ceff>0 permanecen. Cinco documentos matemáticos actualizados en sitio. Datos diagnósticos: /home/jdiaz/scratch/v1_1_1_debug_20261015.

Plan fresco ya preparado: /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_1_20261015_02/campaign_plan.json. Doce casos en dos etapas: anclajes 775/1550 nm a 21.5 uA, tau_escape=20 ps, Rs425, eta0.8, sigma7.549834435 nm; controles de rigidez Tb, Rs608 con Dcal, eta0.667/0.9 y tauEE2.5 ps. Primero terminar y auditar ambos anclajes. Cruce positivo 0.1 mV confirmado durante 0.5 ps, extra5 ps; censura sin cruce200 ps. Sin reset total.

Recursos: hasta12 núcleos físicos de workers, un coordinador y dos núcleos completos reservados; presupuesto90% global. El ejecutor vuelve a contar otras campañas. Al preparar quedaban activos el coordinador08 en CPU13 y un worker05 en CPU10; este inventario puede cambiar. Hasta6 hilos por anclaje cuando haya espacio, sin prometer escalamiento. Datos y cachés en scratch. El costo completo puede ser de horas/días; el tiempo de0.2 ps no es una ETA válida hasta detección. Se registran ETA por sesión y checkpoints. Estos comandos largos NO se lanzaron aquí.

### A. Anclajes frescos v1.1.1

~~~bash
cd /home/jdiaz/pysnspd
env TMPDIR=/home/jdiaz/scratch/tmp/pysnspd \
  PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  pipelines/10_thermal_stiffness_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_1_20261015_02 \
  --stage anchors --max-worker-cores 12
~~~

### B. Controles, sólo tras A completa y auditada

~~~bash
cd /home/jdiaz/pysnspd
env TMPDIR=/home/jdiaz/scratch/tmp/pysnspd \
  PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  pipelines/10_thermal_stiffness_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_1_20261015_02 \
  --stage controls --max-worker-cores 12
~~~

Añadir --resume sólo para la misma campaña1.1.1 interrumpida, con fuentes y parámetros idénticos. La antigua campaña1.1.0 queda histórica.

### C. Reportes individuales cuando terminen ambos anclajes

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  NUMEXPR_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir \
    /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_1_20261015_02/runs/anchor425_eta0800_775 \
    /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_1_20261015_02/runs/anchor425_eta0800_1550 \
  --snapshots-ps 0 1 5 10 20 35 \
  --output-dir /home/jdiaz/scratch/case_reports/thermal_basecases_v1_1_1_20261015_02 \
  --workers 2
~~~

Diagnóstico ligero reproducible: con el mismo entorno ejecutar python -m pytest -q tests/test_thermal_pauli_heat.py tests/test_thermal_stiffness_coupled_candidate.py; sin evolución larga ni reportes adicionales.


## 2026-10-16 — v1.1.2: reloj físico, convergencia y techo100 ps

Los anclajes1.1.1 no progresaban de manera fiable:775 nm cerca de3.02 ps y1550 nm cerca de4.51 ps. Se detuvo sólo su coordinador conSIGTERM y se preservaron checkpoints. Reproducción directa: h=1.9073486328125e-8 ps, RHS máximo de población1048.6757/ps, cero iteraciones y cero cambio en todos los campos. El residuo absoluto estaba debajo de2e-7; el reloj avanzaba sobre un predictor nulo. En el estado1550, Ceff/Cfd=5.41e-7 seguía positivo y su dirección continua incrementaba Ceff; no se concluye que ese checkpoint demuestre una capacidad física negativa.

1.1.2 exige al menos una corrección, usa epsilonNL(h)=2e-7*min(1,h/0.02ps), comprueba la respuesta térmica completa dentro de la convergencia y retrocede las correcciones internas para mantener Ceff>0. No cambia las ecuaciones físicas, Dcal, rigidezTe, fuente, calorPauli, circuito ni la guarda independiente de energía. El presupuesto no lineal se hace más estricto al bajar h. El límite adaptativo1e-6 ps produce parada diagnóstica con checkpoint; no se impone piso a h/Ceff ni se cuenta como detección. Los estados1.1.1 quedan históricos; no reanudarlos como1.1.2.

Caso nuevo:21.5 uA, ancho80 nm, espesor7 nm, longitud320 nm, malla5 nm, Rs425, D calibrada, escape20 ps, eta0.8 y sigma7.549834435 nm. Techo sin cruce100 ps; cruce positivo0.1 mV, confirmación0.5 ps, margen5 ps. Sin reset total. El plan12 casos ya está preparado; controles sólo tras ambos anclajes auditados. Resolver el bug no acredita convergencia de latencia ni validez global del cierre térmico.

Recursos: hasta12 núcleos físicos para workers, un coordinador y al menos2 núcleos completos reservados; presupuesto90% compartido. El ejecutor cuenta otros procesos antes de asignar1–6 hilos por anclaje. BLAS/OpenMP1, datos/cachés/logs en scratch. El costo puede ser de horas o días; no hay ETA completa validada. Las pruebas internas de6 ps usan dos procesos de3 hilos, primer plano, checkpoints y límite de30 min por proceso; su estado final se registra en docs/releases/v1.1.2_checks.json. No extrapolar el rendimiento inicial a100 ps.

### A. Nuevos anclajes — comando largo para el usuario

~~~bash
cd /home/jdiaz/pysnspd
env TMPDIR=/home/jdiaz/scratch/tmp/pysnspd \
  PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  pipelines/10_thermal_stiffness_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_2_20261016_01 \
  --stage anchors --ceiling-ps 100 --max-worker-cores 12
~~~

### B. Controles — sólo después de A completa y auditada

~~~bash
cd /home/jdiaz/pysnspd
env TMPDIR=/home/jdiaz/scratch/tmp/pysnspd \
  PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  NUMBA_CACHE_DIR=/home/jdiaz/scratch/tmp/pysnspd/numba_cache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  pipelines/10_thermal_stiffness_basecases.py \
  --campaign-dir /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_2_20261016_01 \
  --stage controls --ceiling-ps 100 --max-worker-cores 12
~~~

Para reanudar sólo una interrupción1.1.2 con fuentes, parámetros, horizonte y hilos idénticos, añadir --resume. Una parada NUMERICAL_RESOLUTION_LIMIT requiere revisión del solver/cierre, no un reinicio ciego ni cambiar el límite para ocultar el fallo.

### C. Reportes individuales de los nuevos anclajes completos

~~~bash
cd /home/jdiaz/pysnspd
env PYTHONPYCACHEPREFIX=/home/jdiaz/scratch/tmp/pysnspd/pycache \
  MPLCONFIGDIR=/home/jdiaz/scratch/tmp/pysnspd/matplotlib \
  NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 \
  /home/jdiaz/.conda/envs/snspd/bin/python \
  plot_pipelines/07_case_report.py \
  --run-dir \
    /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_2_20261016_01/runs/anchor425_eta0800_775 \
    /home/jdiaz/scratch/photon_runs/thermal_basecases_v1_1_2_20261016_01/runs/anchor425_eta0800_1550 \
  --snapshots-ps 0 1 5 10 20 35 \
  --output-dir /home/jdiaz/scratch/case_reports/thermal_basecases_v1_1_2_20261016_01 \
  --workers 2
~~~

Diagnóstico ligero: con el entorno anterior ejecutar python -m pytest -q tests/test_midpoint_clock_progress.py tests/test_thermal_stiffness_coupled_candidate.py. Evidencia causal y checkpoints internos: /home/jdiaz/scratch/v1_1_2_debug_20261016.
