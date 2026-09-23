# NbN y cinética: qué puede fijarse y qué sigue condicionado

23 de septiembre de 2026. Investigación documental, sin transientes ni cambios al modelo. La orientación acordada es el nanohilo de **80 nm de Korzh 2020**, la latencia relativa **775/1550 nm** y la formación de un hotbelt antes de estudiar jitter. Esta elección identifica la pregunta experimental; no identifica por sí sola todos los parámetros materiales.

El [registro complementario](material_kinetics_sources.json) conserva fuentes, localizadores, clases de evidencia y aritmética reproducible. La referencia experimental general permanece en [reference_experiment.md](../reference_experiment.md).

## 1. Resultado de la revisión

Se puede cerrar la conversión de unidades y evitar tres mezclas: **340 Ω/□ a temperatura ambiente no se transforma en 608 Ω/□ mediante RRR=0.8**; la DOS electrónica por espín es la mitad de la DOS que incluye ambos espines; **0.50/2.47 ps frente a 5/24.7 ps es una elección física de movilidad**, no una conversión de unidades. No aparece evidencia que admita una caja continua de parámetros para la muestra de 80 nm. Sí hay conjuntos históricos coherentes y observaciones de otras películas útiles para evaluar plausibilidad, siempre identificadas.

## 2. Resistencia, difusión y número de estados

La resistencia de hoja $R_\square$ mide una película por cuadrado; con espesor $d$, la resistividad es $\rho_n=R_\square d$ y la conductividad $\sigma_n=1/\rho_n$. La difusión $D$ describe la dispersión espacial de electrones tras muchas colisiones elásticas. La DOS normal $N_0$ cuenta estados electrónicos por energía y volumen **para un espín**; no es la DOS fonónica ni la DOS superconductora dependiente de energía.

En el cierre difusivo usado aquí,

\[
\sigma_n=2e^2N_0D,\qquad
N_0=\frac{1}{2e^2D R_\square d}.
\]

Allmaras usa esta convención por espín [A20, p. 16, PDF 28]. Sidorova et al. definen en cambio una DOS total $N_{\rm total}=2N_0$ y escriben $\sigma_n=e^2N_{\rm total}D$ [S20, tabla II y ecuación (3)/texto contiguo, PDF 4–5]. Comparar sus valores sin el factor dos produciría una discrepancia artificial.

K20 informa $R_\square(\mathrm{ambiente})=340$ Ω/□, RRR=0.8, $d=7$ nm nominal y $T_c=8.65$ K [métodos, PDF 9–10]. **Bajo la convención** $RRR=R_{300}/R_{T_n}$, se obtiene $R_\square(T_n)=340/0.8=425$ Ω/□. K20 no especifica en ese párrafo la temperatura normal $T_n$; no se debe rotular automáticamente 20 K ni extrapolar este dato al estado normal a 0.9 K. S20 sí define explícitamente $R_{300}/R_{20}$ para sus propias muestras [tabla II]. Incluso invertir la convención daría 272 Ω/□, no 608.

| Conjunto identificado | $R_\square$, Ω/□ | $D$, cm²/s | $\sigma_n$, S/m | $N_0$, eV⁻¹ cm⁻³ por espín |
|---|---:|---:|---:|---:|
| Conversión de la película K20 a temperatura ambiente, usando **D heredado sólo para la comparación** | 340 | 1.581 | $4.20168\times10^5$ | $8.29375\times10^{21}$ |
| Conversión RRR condicional anterior, conservando ese D | 425 | 1.581 | $3.36134\times10^5$ | $6.63500\times10^{21}$ |
| Simulación hotbelt K20, tabla suplementaria 1 | 608 | 0.5 | $2.34962\times10^5$ | $1.46652\times10^{22}$ |
| Conjunto A20 de pp. 15–16 y 81 | 600 | 0.5 | $2.38095\times10^5$ | $1.48607\times10^{22}$ |

Las dos últimas columnas son **deducciones algebraicas**, no nuevas mediciones; todas usan 7 nm. Los 600/608 Ω/□ y $D=0.5$ cm²/s son entradas de modelos históricos. El escenario vigente $\sigma_n=4.2\times10^5$ S/m coincide numéricamente con convertir 340 Ω/□, pero esa coincidencia no acredita su conductividad a baja temperatura. El D heredado $1.581\times10^{-4}$ m²/s equivale a **1.581 cm²/s = 158.1 nm²/ps**; el histórico equivale a **0.5 cm²/s = 50 nm²/ps**. Difieren por 3.162, no por diez.

No hay motivo para variar $D,\sigma_n,N_0$ independientemente si se conserva Einstein. Por ejemplo, la conversión RRR condicional reduce $\sigma_n$ y $N_0$ un 20% a D fijo. Cambiar D altera además difusión espacial, ruptura de pares y escalas de longitud; mantener sólo una corriente crítica ajustada no controla todos esos efectos.

## 3. Gap, Tc y espectros fonónicos

El catálogo vigente fija la relación BCS débil $\Delta_0/(k_BT_c)=\pi e^{-\gamma_E}\simeq1.76388$; para 8.65 K resulta **1.31479 meV**. No se encontró una medición de gap específica de la muestra de 80 nm que convierta esa hipótesis en dato experimental.

Simon et al. emplean una corrección de escala a $\Delta_0/(k_BT_c)=2.1$ y una película de referencia de 5 nm, $T_c=10$ K, distinta de K20 [SI25, apéndices VII.3–VII.5]. Aplicar 2.1 a 8.65 K daría **1.56534 meV** por aritmética, sin demostrar que corresponda a esa película. El $\alpha^2F$ usado en colisiones no transforma por sí mismo el funcional estático débil en Eliashberg. Cambiar sólo el gap dejaría sin verificar energía, fuerza, corriente y equilibrio; no se propone hacerlo en esta investigación.

La DOS fonónica $g_{\rm ph}(\Omega)$ cuenta modos de vibración por energía y volumen. La función $\alpha^2F(\Omega)$ pondera esos modos por su acoplamiento a electrones. El número de modos determina capacidad térmica; el acoplamiento determina tasas. Una misma integral $\lambda=2\int\alpha^2F(\Omega)d\Omega/\Omega$ no identifica ambas cosas ni fija toda la dinámica.

Para frecuencia ordinaria $\nu$, $\Omega=h\nu$: **1 THz = 4.13567 meV** y $g_\Omega d\Omega=g_\nu d\nu$. Una DOS por celda exige además la densidad de celdas; una DOS por átomo exige la densidad atómica. En la convención adimensional de $\alpha^2F$, la transformación lineal del eje conserva λ sin el mismo Jacobiano que la DOS. Por ello λ no resuelve la unidad de la tercera columna del archivo.

La auditoría ya guardada del archivo público `nbn-a2f-ph.dat` sigue siendo la evidencia utilizable: el derivado v1 es una **forma condicionada**, con unidades/base/densidad absolutas pendientes. Recalcular aquí los mismos momentos no resolvería esa procedencia. Los huecos de soporte cerca de 30.47 meV afectan el cociente acoplamiento/DOS y una banda no térmica aunque el cambio de capacidad térmica sea pequeño. Véanse [hallazgos y manifiesto previos](../../stage1_closure/material_findings.md). No se adopta una normalización por proximidad de la integral a tres.

Babu y Guo calculan espectros y acoplamiento para distintos politipos cristalinos de NbN mediante DFT [BG19, figuras 5–6 y tabla III]. Su resultado $\lambda=0.98$ para δ-NbN es un cálculo de ese sistema, no una medición de la película K20. SI25 calcula propiedades **bulk** por DFPT y cambia el ensanchamiento electrónico para estabilizar el cálculo de NbN; no demuestra la identidad de ese espectro con una película desordenada de 7 nm. A20 pp. 53–54 ya identifica esta limitación. Esto justifica estudiar formas espectrales mejores que Debye, sin declarar admitidas sus tasas absolutas.

## 4. Cuatro tiempos con trabajos diferentes

| Tiempo | Qué describe | Evidencia y límite |
|---|---|---|
| $\tau_{ee}$ | Interacciones entre electrones; una tasa de pérdida de coherencia y una de redistribución de energía no son automáticamente iguales | A20 usa 5 ps a Tc como parámetro de movilidad ajustado y explora 0/5/10 ps; K20 usa 6 ps. El cero representa un límite del modelo, no una medición de tiempo nulo. |
| $\tau_{ep}$ | Interacción electrón–fonón; distinguir tiempo de dispersión de partícula y tiempo de relajación de energía | A20 parte de 16 ps a 10 K y supone $T^{-3}$: $16(10/8.65)^3=24.7213$ ps. Es una extrapolación modelada, no una medida directa en la película a Tc. |
| $\tau_{esc}$ | Pérdida de fonones hacia el sustrato; depende de interfaz, espesor y modos | A20 obtiene 9.4 ps con su modelo de desajuste acústico NbN/SiO₂; K20 usa 20 ps en el ajuste hotbelt; versiones 2D de A20 llegan a 80 ps ajustados. No constituyen un intervalo de confianza 9.4–80 ps. |
| $\tau_{kin}$ BGK | Tiempo efectivo que acerca p a una FD de la misma energía, sin añadir otra variable térmica | Los 0.7 ps de los ensayos son sintéticos. Ninguna fuente revisada mide ese parámetro del cierre nuevo; no se identifica con $\tau_{ee}$ de KWT. |

Fuentes: [A20, pp. 15, 38, 94–99; K20, suplemento, tabla 1]. Un proceso puede cambiar la distribución conservando energía: por eso BGK puede actuar mientras $T_E$, definida sólo por el momento energético, permanece constante. Ese cierre no conserva necesariamente el número de quasipartículas y no reproduce automáticamente un integral microscópico electrón–electrón.

También hay una coincidencia de símbolos que debe evitarse: el `tau0 = 1.87 ns` del capítulo 2 de A20 normaliza el acoplamiento electrón–fonón mediante `tau_ep(Tc)*720*zeta(5)/pi²`. El `tau0 = pi*hbar/(8*kB*Tc)` de D.11 es la escala temporal GL. Son definiciones diferentes, no dos estimaciones del mismo tiempo.

Un control de plausibilidad externo útil es Sidorova et al. [S20]: cuatro películas NbN/SiO₂ de 5.0–9.5 nm tienen tiempos **de dispersión** electrón–fonón extrapolados a 10 K de 11.9–17.5 ps, exponentes ajustados 3.21–3.77 a partir de datos 14–30 K, y tiempos **de relajación energética** en sus respectivos Tc de 1.4–4.2 ps [tabla III; ecuación (1)]. Esas magnitudes son distintas aun dentro del mismo artículo. Son otras muestras, con D entre 0.339 y 0.474 cm²/s y Tc entre 8.35 y 10.94 K; no son barras de incertidumbre de K20 ni justifican extrapolar la ley hasta el baño de 0.9 K. Sus tiempos de escape se calculan con parámetros acústicos inferidos, no mediante un cronómetro directo de salida fonónica.

## 5. La diferencia de diez veces en movilidad queda identificada

D.11 y `KWTMobility.coefficients` usan

\[
\tau_\psi^{-1}(T)=\frac{T/T_c}{0.50\;\mathrm{ps}}
+\frac{(T/T_c)^3}{2.47\;\mathrm{ps}},\qquad T=\max(T_E,T_b).
\]

`configs/geminga_local_v3.yaml` selecciona esos valores; los defaults de `config.py` son 5 y 24.7 ps. A20 p. 94 define el mismo tipo de suma de tasas con tiempos de dispersión, y sus casos posteriores utilizan 5/24.7 ps. Reducir **ambos** tiempos por diez reduce $\tau_\psi$ exactamente por diez a la misma T: a 0.9 K son 4.795 frente a 47.950 ps; a 4 K, 1.036 frente a 10.364 ps. Son valores derivados de las fórmulas, no latencias simuladas.

No debe repararse esta diferencia cambiando sólo la unidad del archivo ni reescalando el eje temporal. En KWT, $R=\sqrt{1+4|\Delta|^2\tau_\psi^2/\hbar^2}$ frena la respuesta radial y modifica de otra manera la respuesta de fase. El cambio afecta la trayectoria acoplada. Tampoco es válido reinterpretar automáticamente 2.47 ps como el tiempo energético de S20: la ecuación de movilidad necesita su propia identificación y el factor diez se aplica también al término ee.

**Decisión documental:** 0.50/2.47 ps se conserva como movilidad efectiva heredada **no calibrada para el objetivo K20**, mientras 5/24.7 y el 6/24.7 del ajuste K20 son referencias históricas separadas. No se modificó ningún coeficiente. Antes de atribuir una diferencia de latencia al nuevo modelo microscópico deberá distinguirse de la sensibilidad a esta movilidad.

## 6. Qué puede decidirse ahora

Para la orientación de 80 nm conviene preparar **un escenario K20/Allmaras explícito**, conservando el escenario de la memoria como comparación. Esto no adopta sus ajustes como verdad material. La etiqueta del escenario debe incluir conjuntamente $T_c,d,R_\square,D,N_0$, cierre de gap, DOS/acoplamiento, movilidad y escape. No hace falta inventar rangos independientes para avanzar con un ensayo sintético declarado.

Quedan pendientes una caracterización de $D$/resistencia normal de la misma muestra, gap experimental, normalización absoluta del espectro y una identificación de la movilidad. Una futura comparación pequeña debe medir cuánto cambia la **latencia relativa** y el momento de homogeneización transversal al cambiar los cierres con procedencia identificada; sólo después tiene sentido invertir más cómputo en jitter. No se ha lanzado esa comparación.

**Decisión de orientación recibida:** se usará como referencia provisional el escenario Korzh/Allmaras, sin tratar sus ajustes como mediciones. Por tanto, para preparar investigación se vinculan D=0.5 cm²/s, R□=608 Ω/□, d=7 nm y N0 por Einstein. El escenario heredado se conserva como comparación. Esto no cambia parámetros en código ni admite una caja física. No hace falta volver a preguntar por esa elección.

La inductancia tampoco es una entrada material independiente si se deriva del mismo cierre: en el límite BCS difusivo, homogéneo, de baja temperatura y corriente pequeña, `L_square = hbar*R_square/(pi*Delta0)`. Con 608 Ω/□ y el gap débil anterior da aproximadamente 97 pH/□, frente a los 64 pH/□ estimados para el inductor experimental K20. La conversión RRR de 425 Ω/□ daría aproximadamente 68 pH/□. Esta comparación condicional no mide inductancia de la región perturbada; evidencia que no debe heredarse también una profundidad de penetración de London independiente y declarar todas las escalas simultáneamente calibradas.

## 7. El soporte cinético debe acompañar al material

El CSV derivado v1 llega en su eje almacenado hasta 20; bajo la interpretación condicional THz, el último nodo de soporte común positivo está en 16.3901 THz (67.78401 meV). Esto describe esos bytes y esa conversión todavía condicionada, no un corte universal de NbN. El catálogo electrónico de conteo hasta 12 tiene, por ejemplo, $E_{\max}\simeq15.832$ meV a corriente cero y amplitud $\Delta_0$. La suma de dos energías no es el mismo límite que la absorción electrónica $E\to E+\Omega$. Un fonón de alta energía puede requerir estados que la tabla no contiene. La rejilla sintética de 1025 nodos hasta $4\Delta_0\simeq5.259$ meV tampoco cubre por ello el espectro material. Antes de trasladar la inyección D.30 o las colisiones se necesita admitir soportes y colas con el nuevo material; no basta conservar el número de nodos.
