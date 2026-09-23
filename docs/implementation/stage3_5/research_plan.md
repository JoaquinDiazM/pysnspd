# Etapa 3.5: investigación del dominio de confianza

Fecha: 23 de septiembre de 2026. Estado: **investigación documental autorizada, sin barridos de simulación**.

El usuario autorizó cerrar el desarrollo de la etapa 3 e iniciar esta investigación. Ese cierre acepta el alcance del trabajo realizado; no acredita el transiente espacial que todavía falta. Antes de las etapas 4–5 se mantienen como requisitos la dinámica débil acoplada y la interfaz cinética 2D–1D, con sus bordes y balances. Este plan no modifica contratos, tolerancias ni dictámenes históricos.

El [inventario](parameter_inventory.json) contiene 127 entradas en 15 familias. Cada entrada distingue nombre y símbolo, unidades, rol, valor utilizado y fuente, restricciones, acoplamientos e investigación pendiente. Un valor de ensayo no se convierte en una propiedad medida. Las cajas del catálogo tampoco son intervalos de validez física. Se incluye lo heredado únicamente para reconocer qué no debe trasladarse automáticamente al modelo nuevo.

## Qué significa «dominio de confianza»

Se busca describir para qué condiciones del material, dispositivo, preparación y lectura el modelo puede responder una pregunta concreta con aproximaciones justificadas. El dominio depende del observable: conservar energía no valida por sí solo un tiempo de detección, y un pequeño error de energía total puede ocultar una fuerza local incorrecta.

**No se trata de un intervalo estadístico.** Un intervalo estadístico requiere datos o un conjunto de realizaciones y un modelo probabilístico. En esta etapa se documentará evidencia, incertidumbre de entradas y límites de aproximación. No se atribuirán probabilidades de confianza a una caja de parámetros ni a una única trayectoria determinista.

Cada resultado de investigación tendrá una ficha breve: afirmación, fuente primaria específica y convención, parámetro afectado, condiciones de aplicación, incertidumbre conocida, qué sigue desconocido y acción propuesta. Si la literatura ofrece un rango, se conservará como **rango publicado o propuesto, no adoptado**. Si no hay evidencia suficiente, se escribirá «desconocido».

## Orden de investigación y entregables

| Prioridad | Pregunta a resolver sin simular | Entregable y condición para avanzar |
|---|---|---|
| 1. Escenario coherente | ¿Qué valores corresponden a la misma película, temperatura y montaje? | Tabla única de material y circuito con procedencia; conflictos explícitos, sin promediar escenarios. |
| 2. Datos y aproximaciones | ¿Qué limita las tasas absolutas, el espectro local y la movilidad efectiva? | Fichas de DOS/α²F, Einstein, acoplamiento débil, BGK y KWT; separar incertidumbre material de error numérico. |
| 3. Transferencia fotónica | ¿Qué representa el estado en el instante de transferencia t0? | Definiciones vinculadas de energía retenida, ancho, posición y tiempo; fuente de la cascada omitida y pérdidas. |
| 4. Dominio espacial | ¿Cómo escoger región 2D, continuaciones y ventana de observación sin contaminar el observable? | Argumentos de escala y propuestas condicionadas; pruebas dinámicas necesarias identificadas, no ejecutadas. |
| 5. Señal y tiempo | ¿Qué observable se compara con qué medición? | Definición de latencia, umbral, electrónica, jitter y referencia temporal; distinguir determinismo y distribución. |
| 6. Preparación numérica | ¿Qué controles faltan para un transiente débil útil y cuánto costarían? | Diseño de un próximo ensayo registrado y estimación de recursos; todavía sin lanzamiento ni barrido. |

La investigación puede continuar con entradas materiales desconocidas si las marca. Para un ensayo numérico se puede proponer un material sintético explícito; eso no resuelve la admisión de NbN ni autoriza una predicción experimental.

## Conflictos que deben resolverse antes de combinar parámetros

Los siguientes contrastes proceden de los archivos actuales, no de una nueva búsqueda de valores convenientes:

- El catálogo y D.2 usan `σn = 4.2×10⁵ S/m`, ligado a `N0` y `D` mediante Einstein. El ensayo mixto instantáneo declara `2.5×10⁵ S/m` como escenario. Su balance eléctrico puede ser correcto sin que ambas elecciones describan un único material coherente.
- El baño histórico es `0.9 K`; los pilotos usan `kBTb/Δ0 = 0.12`. No son la misma temperatura.
- La adenda circuital vigente especifica `Ccouple = 100 pF`. El YAML histórico v3 contiene `C_rf = 1 pF`; no debe alimentar por accidente el circuito nuevo.
- La movilidad experimental usa `0.50 ps` y `2.47 ps`. Los valores predeterminados de algunas rutas antiguas son `5 ps` y `24.7 ps`. El `τkin = 0.7 ps` del ensayo de dos celdas es una entrada sintética distinta de ambos tiempos de movilidad.
- Los `10 nH` históricos son una referencia **total**. Los `9.666757458681755 nH` del ensayo son la porción exterior obtenida para un dominio de referencia de `720 nm`. Cambiar el dominio o el estado de referencia exige revisar esa partición antes de integrar; no se resta una inductancia instantánea del hotspot.

Fuentes: [D, §D.2](../../modelo_v0_4/D_sintesis_plan_y_verificaciones_v0_4.md), [adenda CM](../../modelo_v0_4/actualizaciones/circuito_memoria_20260922.md), [registro mixto](../stage3/coupled_20260923/mixed_registration.json), [configuración histórica](../../../configs/geminga_local_v3.yaml). Las formulaciones antiguas y sus resultados se conservan íntegros.

## Material, cierres y espectro: qué investigar

La identificación de `Tc`, `D`, `σn`, `N0`, espesor, ancho y penetración London debe distinguir medición, inferencia por calibración y entrada de escenario. En particular, ajustar `D` a una corriente crítica no constituye una medición independiente de difusión. La relación débil entre `Δ0` y `Tc` pertenece al cierre seleccionado: cambiar sólo el gap sin reconstruir energía, fuerzas y corriente no es una corrección consistente de acoplamiento fuerte.

La DOS fonónica derivada de NbN conserva su estatus de forma condicionada. Deben investigarse la unidad del eje y su Jacobiano, la base por átomo o celda, la densidad de modos, el preprocesado exacto y la convención de tasas. El recorte autorizado y trazable de negativos no resuelve por sí solo la normalización absoluta ni los huecos de soporte del acoplamiento. La integral λ tampoco identifica la unidad de la DOS. Véase la [revisión material](../stage1_closure/material_findings.md).

Se investigarán por separado `τkin`, la movilidad KWT, `τesc` y la distribución espectral de calentamiento. Son cierres con propósitos distintos. La potencia Joule normal y la disipación del condensado se depositan una sola vez; el trabajo superconductivo y el trabajo de reservorio son cantidades signadas, no calor positivo adicional.

El regularizador `δ = 0.10Δ0` es una elección efectiva de núcleo, no un espaciado numérico. El regulador causal `η = 10⁻⁸Δ0` es numérico y no una tasa inelástica. Tampoco debe confundirse con una eventual fracción de retención fotónica. El dominio de amplitud positiva del catálogo `[0.08, 1.5]Δ0`, su punto normal exacto y su caja de Γ no acreditan todo el dominio físico entre ellos. Se investigarán las limitaciones del espectro local adiabático y del criterio D.36. Un símbolo principal positivo en estados muestreados es necesario para el bloque del condensado, pero no prueba estabilidad de todo el sistema acoplado ni validez física de un núcleo.

## Geometría, tiempo y protección de los bordes

`L2D/W = 1.5–6` se conserva únicamente como **ejemplo propuesto por el usuario para investigar**. No se adopta como rango correcto. La longitud de la región 2D no es la longitud total del dispositivo ni la distancia del impacto al reservorio.

Deben distinguirse `L2D`, las dos continuaciones y los planos eléctricos exteriores. Las distancias `d_center,L` y `d_center,R` se miden desde el centro de la deposición hasta esos planos. Si se define una región perturbada mediante un umbral explícito de campo o energía, sus distancias a los planos se registrarán por separado como `d_region,L` y `d_region,R`; no se sustituyen por las distancias desde el centro. Una gaussiana tiene soporte infinito y carece de un borde compacto. La ventana `tobs` se fija según el observable: respuesta temprana, desarrollo de una perturbación, recuperación o señal circuital no requieren necesariamente el mismo dominio. La modificación de longitud cambia también el almacenamiento ya resuelto y la partición de inductancia exterior.

Como filtro analítico **condicional**, para difusión lineal en espacio libre con coeficiente constante D y perfil gaussiano de desviación `sγ` en el instante de transferencia `t0`, se define el tiempo transcurrido `τ = t − t0 ≥ 0`. La varianza por coordenada es

\[
\sigma^2(\tau)=s_\gamma^2+2D\tau.
\]

La fracción instantánea de masa gaussiana fuera del intervalo entre los dos planos, usando las distancias desde su centro, es

\[
\epsilon_{\rm exterior}(\tau)=\tfrac12\operatorname{erfc}\!\left(\frac{d_{{\rm center},L}}{\sqrt{2s_\gamma^2+4D\tau}}\right)
+\tfrac12\operatorname{erfc}\!\left(\frac{d_{{\rm center},R}}{\sqrt{2s_\gamma^2+4D\tau}}\right).
\]

Esta masa exterior no es el escape acumulado a través de bordes absorbentes ni una probabilidad de primer paso. Es una identidad de ese problema lineal auxiliar, **no una garantía del modelo acoplado**, de su solución con bordes ni del error del voltaje. En el transporte espectral, el coeficiente efectivo depende del espectro; el potencial es una restricción elíptica global. Una PDE difusiva no posee un frente de propagación a velocidad finita: se estudia la pequeñez del efecto sobre un observable, no un tiempo de llegada exactamente nulo. No se usará `distancia = velocidad × tiempo` como prueba de protección.

El presupuesto de contaminación de borde se propondrá según el efecto físico que se pretende distinguir; no se prefija `10⁻¹⁰` por costumbre. Después será necesario comparar transientes con dominios extendidos a resolución comparable. Los controles estacionarios y las instantáneas actuales no sustituyen esa comprobación.

La reducción transversal 1D queda **pospuesta a la validación transiente**. El mapa R implementado identifica sólo la traza uniforme del condensado y potencial; no representa modos superiores que atraviesen una interfaz. Una escala como `W²/(π²D_eff)` puede orientar la investigación en un problema difusivo lineal, pero no demuestra por sí sola que esos modos hayan desaparecido en el sistema acoplado. Habrá que medir su historia, energía y flujo en las secciones propuestas, junto con la continuidad cinética a energía fija. No basta una razón de aspecto grande.

## Inyección, lectura y jitter

`Ein`, `sγ`, `r0` y `t0` son datos vinculados de la transferencia desde una cascada no resuelta. El perfil gaussiano de D.30 usa desviación estándar, no un radio de disco. Renormalizar una gaussiana truncada conserva la energía indicada, pero no calcula las pérdidas tempranas cerca del borde. El escape posterior se calcula después y no vuelve a descontarse de `Ein`. La inyección del modelo nuevo será un salto espectral fonónico, no la reutilización silenciosa de una temperatura fonónica del solver antiguo.

Las fluctuaciones Fano permanecen opcionales y desactivadas. Investigar un factor F no equivale a haber especificado una distribución de energías, una correlación con posición o una implementación estocástica.

Se distinguirán el voltaje entre puertos energéticos `Vdev`, el voltaje central histórico, `Vout` en la carga, y la señal después de una eventual cadena de amplificación medida. El circuito tiene tres estados `Ib, Is, vc`; `Is` es la corriente total del detector y `Ib` ya no es una fuente ideal de corriente fija. La electrónica requiere datos de impedancia, banda, ganancia, retardo, ruido y referencia temporal cuando se quiera comparar una medida. El ancho de banda y el intervalo de muestreo no son el paso de integración.

Una latencia determinista es un tiempo de cruce para una condición inicial concreta. El jitter es una distribución temporal de eventos. La [revisión experimental paralela](reference_experiment.md), con su [registro de fuentes](source_register.json), identifica que K20 informa jitter como FWHM de la IRF y latencia relativa como diferencia entre sus máximos ajustados, no entre medias. Deben mantenerse la referencia óptica, los ajustes de amplificador y trigger y los filtros calibrados. Las barras identificadas como intervalos del 95% se conservan como tales, no como desviaciones estándar. No se convertirá FWHM a desviación estándar sin justificar la forma de la distribución, ni se llamará jitter al error numérico o a una diferencia de latencias de dos mallas.

El inventario añade la identidad de muestra, resistencia de hoja y temperatura de medición, sustrato, tapers, foco, polarización, adquisición, parámetros ajustados y covarianzas. Los valores del experimento o de un modelo histórico se guardan como referencia no adoptada. Por ejemplo, el inductor adicional de 96 nH de K20 no reemplaza automáticamente la inductancia de la memoria. Su parámetro fonónico γ tampoco es Γ de Usadel, y las dispersiones energéticas ajustadas no son el ancho espacial de la inyección. Estas distinciones impiden combinar en un mismo supuesto detector observaciones procedentes de muestras y cadenas diferentes.

Los umbrales, histéresis, tiempos de confirmación y reglas de recuperación que existen en `analysis/timing.py` son definiciones operativas del software. Sus valores predeterminados no son parámetros medidos de una nueva electrónica. Recuperar el 90% de corriente tampoco demuestra recuperar el 90% de eficiencia de detección.

## Lo que falta para la validación débil antes de 4–5

1. Preparar un estado estacionario del mismo problema mixto con reservorios y corriente circuital consistentes. La hélice libre del piloto no es ese equilibrio.
2. Ensamblar el transporte electrónico 2D–1D y los reservorios a energía común, conservando capacidades, energía y ocupaciones. El mapa R de campos no completa esta tarea.
3. Integrar conjuntamente el condensado complejo, ocupaciones electrónicas y fonónicas, potencial instantáneo y los tres estados circuitales. El cambio de Γ y su trabajo espectral deben estar contenidos una sola vez.
4. Incorporar la rama móvil `a_b(Is(t)), q_b(Is(t))`, su respuesta radial y su trabajo. El postproceso actual acredita sólo radio constante y balance instantáneo; el pequeño flujo normal de una arista interior no es una traza exterior D.27 certificada.
5. Registrar un ensayo débil con observables, referencias, incertidumbres y presupuestos de tiempo/espacio/espectro. Comprobar soporte, causalidad, positividad, D.36, balances acumulados y conservación de corriente durante la evolución, además de sensibilidad a extensión del dominio.
6. Sólo después decidir el siguiente ensayo de depósito sintético localizado. La identificación material y la transferencia fotónica serán necesarias antes de presentar umbrales, pulsos o latencias como predicciones del detector.

Estos puntos forman una lista de requisitos pendientes, no una nueva aceptación automática ni un plan de barrido. El certificado estricto de mallas de la etapa 2 sigue incompleto; se usará su evidencia en el alcance probado y se controlará lo que efectivamente utilice el nuevo ensayo.

## Política de entrega y coste

Esta etapa produce inventario, fichas de fuentes, conflictos, propuestas de rango y diseño del siguiente ensayo. No ejecuta transientes ni barridos exploratorios. No promueve producción y no modifica `v1.0.0`.

Antes de un cálculo futuro se estimarán duración, memoria y fuentes dominantes de coste. Los cálculos que razonablemente superen cinco minutos se anotarán en `/home/jdiaz/GEMINGA_COMMANDS.md`, preservando entradas anteriores, y su comando exacto se entregará también en un bloque copiable del chat. Se esperará la ejecución del usuario y su mensaje de reinicio; no se usarán trabajos de fondo, sondeos ni subdivisiones para eludir esa política.
