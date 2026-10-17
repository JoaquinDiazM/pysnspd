# pySNSPD 1.1.3

**Neutral efectivo BCS/KWT con corrección térmica de depareo y referencia DC conservada.** La revisión sustituye la extrapolación cuadrática del gradiente térmico a gradientes grandes, mantiene la difusión calibrada al depareo y resuelve conjuntamente condensado, distribuciones electrónicas/fonónicas y circuito.

El bloqueo 1.1.2 tenía una componente constitutiva: Ceff podía acercarse a cero y cambiar de signo a gradiente fijo. La nueva energía libre evalúa la corrección térmica mediante Matsubara de depareo finito. Un término conservativo independiente de temperatura preserva el DC anterior y deja elasticidad residual a Te alta: es una hipótesis explícita. [Diagnóstico](docs/releases/v1.1.3.md), [modelo vigente](docs/implementation/MODELO_VIGENTE.md) y [estado](status.md).

Los anclajes nuevos de 775/1550 nm completaron **0–12 ps** con capacidad positiva, soporte físico y balance independiente. La continuación 8–12 ps conserva los estados, las fuentes y el trabajo externo, sin redepositar el fotón. Esa cobertura no certifica detección, latencia experimental ni jitter.

## Ejecución

La entrada operativa actual es [pipeline 09](pipelines/09_thermal_stiffness_photon_run.py), con [campaña 10](pipelines/10_thermal_stiffness_basecases.py) y [reporte 07](plot_pipelines/07_case_report.py). [Comandos presentes](docs/GEMINGA_COMMANDS_PRESENT.md) reúne las órdenes normales y visuales. Las clases térmicas especializadas permanecen en pysnspd.experimental; la API pysnspd.effective ofrece piezas compartidas y la referencia al baño, por lo que el cierre térmico debe seleccionarse explícitamente.

~~~bash
python -m pip install -e '.[performance]'
~~~

La referencia usa ancho80 nm, espesor7 nm, baño0,9 K, corriente21,5 µA y escape20 ps. Rs=425 Ω y D calibrada son hipótesis documentadas. El depareo25,8667 µA resulta de escalar una estimación para120 nm: no es una medida independiente de la cinta80 nm.

La observación termina con Vout positivo0,1 mV, confirmación0,5 ps y margen5 ps desde el primer cruce; sin gatillo se censura a100 ps. Vout no incluye filtro/amplificador experimental. No se pide reset completo.

Todos los datos pesados, logs, cachés JIT y renders van a /home/jdiaz/scratch; NUMBA_CACHE_DIR debe apuntar allí. Se conserva el presupuesto global90 % y los checkpoints antiguos con sus fuentes. La nueva física inicia casos fresh.

## Alcance

La cinética conserva espectros BCS locales q=0: la corrección de gradiente no introduce una DOS cinética depareada. Qeff reúne gradientes de fase y amplitud como hipótesis efectiva. El material fonónico es acústico; la tabla NbN/QNN-MIT no tiene normalización absoluta resuelta. No es Usadel microscópico exacto.

Tau_ee(Tc) y tau_ep(Tc), nominales libres5 y24,7 ps de nomenclatura Allmaras, siguen siendo calibradores. El reloj comienza postcascada; la cascada ultrarrápida no se evoluciona. Una trayectoria determinista no determina jitter estadístico ni latencia óptica absoluta.

Se mantienen las tolerancias físicas y las guardas de capacidad/rigidez, soporte y energía. Las sesiones0–8 ps proyectan11,18/12,43 h por100 ps para775/1550 nm, incluyendo preparación/JIT; no es garantía ni medida exclusiva posterior al DC. El benchmark del mismo estado a12 ps favorece seis hilos:1,53× frente a dos, con RHS y endpoints idénticos bit a bit. La afinidad L3 propuesta no se promovió sin una comparación específica.

## Historia

thesis-v1.0.0 conserva la memoria; v1.0.0 conserva la referencia efectiva admitida exploratoriamente en0–36,6 ps de775 nm. Esa admisión no se transfiere a1.1.3. [Changelog](CHANGELOG.md), [releases](docs/releases/) e [historial de comandos](docs/GEMINGA_COMMANDS.md) preservan decisiones anteriores.

[Licencia](LICENSE). Mantenedor: Joaquin Andres Diaz Monge.

La región crítica775 nm exigió pasos más pequeños: el tramo8–12 ps tomó71,7min, frente a28,2min en1550 nm. No extrapolar un intervalo suave al pulso entero. La ampliación diagnóstica a25ps quedó aplazada al cerrar esta entrega; no se ejecutó ni se declara completa.
