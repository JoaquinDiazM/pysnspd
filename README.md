# pySNSPD 1.1.1

**Versión vigente: v1.1.1, extensión térmica efectiva con difusión calibrada.** Corrige el bloqueo de los anclajes 1.1.0 causado por un portador de calor que podía enfriar estados vacíos. El reparto usa actividades FD/Pauli según el signo de la potencia: conserva energía y produce entropía adicional no negativa fuera de FD. Conserva la rigidez dependiente de Te y la guardia Ceff>0, sin clipping ni pisos.

Los casos nuevos usan [pipeline 09](pipelines/09_thermal_stiffness_photon_run.py) y [campaña 10](pipelines/10_thermal_stiffness_basecases.py). Los dos ensayos nominales completaron 0.2 ps en menos de tres minutos por caso; esto elimina el bloqueo observado cerca de 0.05 ps y no acredita todavía latencias completas. Iniciar una campaña nueva: las trayectorias 1.1.0 se preservan y no se reanudan con ecuaciones distintas. Véanse [corrección y evidencia](docs/releases/v1.1.1.md), [matemáticas](Matematica/), [supuestos](Matematica/Supuestos.pdf) y [comandos actuales](docs/GEMINGA_COMMANDS_PRESENT.md).

## Referencia histórica admitida: v1.0.0

**La etapa 5 de la referencia v1.0.0 se completó el 12 de octubre de 2026 con admisión exploratoria.** La decisión y su alcance están en [estado](status.md) y en [etapa 5](docs/implementation/stage5/README.md).

Modelo efectivo de investigación para estudiar la respuesta de una cinta superconductora tras la absorción de un fotón. Evoluciona conjuntamente el gap complejo, las ocupaciones de cuasipartículas y fonones, el potencial eléctrico cuasiestático y el circuito completo de polarización/lectura de la memoria. Conserva una energía común y registra un balance independiente de potencia y energía.

La referencia actual es una cinta de ancho 80 nm y espesor 7 nm. La separación de malla, h_mesh, es una resolución numérica y no el espesor. El método usa volúmenes Delaunay–Voronoi y punto medio implícito acoplado. Los campos y observables se guardan en unidades físicas; las coordenadas internas normalizadas se identifican explícitamente.

## Uso

Instalación desde el repositorio con Python ≥3.10:

~~~bash
python -m pip install -e '.[performance]'
~~~

Numba es una dependencia opcional de aceleración; la campaña de referencia utiliza su backend paralelo. NumPy, SciPy, Matplotlib y las dependencias de geometría están declaradas en pyproject.toml.

La API pública es pysnspd.effective:

~~~python
from pysnspd.effective import (
    EffectiveState, EffectiveCoupledSystem,
    EffectiveBCSFunctional, EffectiveBCSKinetics,
    EffectiveNeutralPorts, prepare_effective_dc,
)
~~~

EffectiveState contiene delta_bar = Δ/(k_B T_c), p, n y circuit_SI = (I_b, I_s, v_c). Las ocupaciones p y n son distribuciones que evolucionan; las temperaturas equivalentes son observadores de energía.

Para una ejecución reproducible:

- [Pipeline física 05](pipelines/05_effective_photon_run.py): prepara su propio DC efectivo, deposita una vez la energía retenida y sigue hasta el cruce positivo de 0,1 mV más 5 ps, con confirmación a cruce+0,5 ps; sin cruce, queda censurado a 200 ps.
- [Pipeline visual 05](plot_pipelines/05_plot_effective_photon_run.py): mapas del gap en meV, temperaturas equivalentes en K, energías en J/m³ y señales en mV/µA.
- [Comandos presentes](docs/GEMINGA_COMMANDS_PRESENT.md): ejemplos normales, visuales y de reanudación.
- [Contrato del modelo](docs/implementation/MODELO_VIGENTE.md) y [estado](status.md): ecuaciones, hipótesis y evidencia de admisión.

En Geminga todos los resultados, cachés y figuras intermedias van a /home/jdiaz/scratch. Cada ejecución guarda parámetros, fuentes identificadas por SHA-256, estado completo reanudable y trazas de balance. Las ejecuciones simultáneas comparten un máximo global del 90 % de recursos.

## Resultado admitido

Se completaron un único sistema de evolución (5.1), su balance independiente (5.2) y cuatro controles numéricos de 775 nm (5.3). La resolución por defecto tiene **54 estados electrónicos**; 80 estados constituyen un control opcional de sensibilidad, no una verdad continua.

| Caso nominal | Resultado guardado |
|---|---|
| 775 nm | Hotbelt transversal observada a 15,3 ps; cruce eléctrico a 31,9413 ps; final a 36,9413 ps |
| 1550 nm | Sin cruce de 0,1 mV durante 200 ps; latencia no asignada |

Los cuatro ejes de 775 nm pasan sus márgenes declarados en **0–36,6 ps**, con cobertura crítica hasta 33 ps. Todas las colas posteriores a 36,6 ps quedan fuera de esa comparación, aunque cada trayectoria completa conserva su auditoría individual. El margen espectral de Vout del 5 % se adoptó explícitamente después de los datos; el FAIL anterior al 2 % se conserva. Esta admisión permite exploración del modelo y calibración posterior, sin equivalencia experimental automática.

## Alcance físico

Es un cierre mesoscópico neutral BCS/KWT, con rigidez radial y cinética conservativa. El material fonónico es acústico efectivo. La tabla externa NbN/QNN-MIT se conserva como referencia documental; no proporciona las tasas de esta versión ni tiene su normalización absoluta resuelta.

τ_ee(T_c) y τ_ep(T_c) son calibradores libres con nomenclatura de Allmaras (2020). Sus nominales 5 ps y 24,7 ps no fijan el centro ni los límites de un rango. La anchura gaussiana 7,5498 nm y la retención 0,667 son referencias exploratorias ajustables.

El reloj comienza en la transferencia efectiva posterior a la cascada. Una trayectoria determinista permite estudiar respuesta y latencia relativa bajo esas hipótesis; no calcula por sí sola jitter estadístico, retraso óptico absoluto ni reset eléctrico total.

## Versiones e historia

La etiqueta [thesis-v1.0.0](docs/releases/thesis-v1.0.0.md) conserva el modelo de la memoria. Por instrucción explícita del usuario del15 de octubre de2026, [v1.0.0](docs/releases/v1.0.0.md) congela ahora esta implementación efectiva. La revisión1.1.1 incorpora la extensión térmica y mantiene separada su admisión de la referencia v1.0.0. El intento de evolución microscópica de primer orden Moyal falló sus controles y se preserva como antecedente rechazado; no forma parte de la física efectiva publicada.

La [documentación anterior](docs/releases/history/README_before_v1.1.0_20261011.md) y el [historial de comandos](docs/GEMINGA_COMMANDS.md) conservan las decisiones y resultados previos. Los datos pesados permanecen en scratch y se enlazan mediante sus manifiestos.

[Licencia](LICENSE). Mantenedor: Joaquin Andres Diaz Monge.
