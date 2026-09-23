# Caracterizar la transferencia antes de escoger el ancho

23 de septiembre de 2026. **Decisión vigente: caracterizar primero la cascada y dejar $s_\gamma$ abierto.** Resuelve la elección planteada en [gaussian_initial_condition.md](gaussian_initial_condition.md): su equivalencia histórica de radios sigue siendo una comparación, no una condición inicial adoptada. Este plan no modifica D.30, no inicia simulaciones ni fija una fracción retenida.

La tarea tiene un final concreto: obtener, para 775 y 1550 nm, una descripción suficiente de la energía que llega al modelo después de la etapa omitida, junto con el instante de transferencia y sus limitaciones. No se busca resolver todos los detalles de la absorción óptica. El producto puede ser una fuente efectiva respaldada por datos o un cálculo cinético reducido; debe permitir decidir si la fuente puramente fonónica y gaussiana de D.30 representa lo necesario para estudiar latencia relativa y formación del hotbelt.

## 1. Lo que las fuentes ya permiten afirmar

La cinética radial de Allmaras distingue energía electrónica que se extiende y energía fonónica más localizada. Sus figuras 2.6–2.8 corresponden a un caso modelado de 1 eV, baño de 4.325 K y $D=0.5$ cm²/s; omiten ee, difusión fonónica y escape en esa comparación. No son datos de 775/1550 nm a 0.9 K. Sus perfiles cilíndricos y esféricos tampoco resuelven una película de 7 nm con ambas superficies y una profundidad de absorción definida. [A20, pp. 15–16, 20–24; PDF 27–28, 32–36](https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf).

Zotova–Vodolazov 2012 parte de un disco térmico tras una termalización omitida, mientras su estimación gaussiana representa difusión electrónica libre. Ese antecedente no identifica el perfil de una burbuja fonónica. [ZV12, sección II y ec. (13)](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.85.024509/fulltext). Simon et al. fundamenta una forma espectral de la fuente fonónica, pero presupone su volumen espacial y omite difusión en la cinética inicial: no aporta la caracterización radial que falta. [SI25, sección III y apéndice de condición inicial](https://arxiv.org/html/2501.13791v3).

Para Korzh quedan por identificar el reparto electrónico/fonónico, su extensión lateral, su variación en profundidad y las pérdidas **antes** del instante elegido. La coincidencia de un D histórico no convierte las demás hipótesis de A20 en mediciones de esa muestra.

## 2. Qué se debe entregar para cada color

Se usará tiempo desde la absorción, $\tau=t-t_{\rm abs}$, mientras se caracteriza la etapa omitida. La transferencia ocurrirá en $\tau_{\rm tr}$; el reloj del modelo posterior puede ponerse a cero allí. La diferencia $\tau_{{\rm tr},775}-\tau_{{\rm tr},1550}$ debe conservarse como resultado o incertidumbre, pues contribuye a la latencia relativa aunque se omita de un transiente posterior. Esta notación evita confundir el instante de transferencia con el $\tau_0$ cinético usado para los ejes de A20.

Todas las energías siguientes son **excesos respecto al mismo fondo de baño**, con unidades y volumen explícitos. Una fracción fonónica de equilibrio preexistente no se cuenta como energía fotónica.

| Producto mínimo | Definición o propósito |
|---|---|
| $E_e(\tau)$, $E_{\rm ph}(\tau)$ | Energía retenida en cada subsistema; guardar valores y fracciones de $E_\gamma=hc/\lambda$. |
| $E_{\rm esc}(\tau)$ y flujos externos | Pérdida acumulada al sustrato, otras superficies o límite espacial del cálculo, separando cada destino. |
| Perfil energético de cada subsistema | Densidad lateral integrada en profundidad $\mathcal U_\alpha(\mathbf r,\tau)=\int dz\,u_\alpha(\mathbf r,z,\tau)$, en J/m²; guardar también el perfil por profundidad si está resuelto. |
| $R_{50,\alpha}$, $R_{90,\alpha}$ | Radios centrados en el depósito que contienen 50%/90% de la **energía de ese subsistema**, además de los radios de energía total. |
| Centro, segundo momento y pico | Centro energético, matriz de covarianza lateral, densidad máxima y resolución a la que se midió; detectan colas, anisotropía y diferencias que un único radio oculta. |
| Contenido espectral | $p(E,\mathbf r,\tau)$ y exceso $\delta n(\Omega,\mathbf r,\tau)$, o una reducción espectral con error declarado; comprobar si espectro y perfil espacial pueden separarse. |
| Sector y profundidad | Electrónico/fonónico; profundidad inicial o distribución usada, geometría cilíndrica/esférica/película y condición de cada superficie. No promediar escenarios sin pesos físicos. |
| Registro de transferencia | $\tau_{\rm tr}$, energía retenida, reparto, perfil, pérdidas anteriores, soporte espectral y espacial, balance, incertidumbre y fuentes. |

Para energía no negativa de un subsistema $\alpha$,

$$
F_\alpha(R,\tau)=
\frac{\int_{|\mathbf r-\mathbf r_c|\le R}\mathcal U_\alpha(\mathbf r,\tau)d^2r}
{E_\alpha(\tau)},\qquad
F_\alpha(R_{q,\alpha},\tau)=q.
$$

El centro y la covarianza se normalizan por el mismo $E_\alpha$. Sólo si el perfil es isotrópico, su anchura por momento sería $s_{{\rm mom},\alpha}^2=\mathrm{tr}(\Sigma_\alpha)/2$. No se declara que el perfil sea gaussiano por calcular ese número. Si el exceso tiene zonas negativas o $E_\alpha$ es despreciable, esos percentiles dejan de ser medidas válidas de una distribución positiva: se informa el perfil firmado y el denominador, sin recortar valores para producir un radio.

**Por qué el radio total no basta.** Con $E_{\rm ret}=E_e+E_{\rm ph}$ y un centro común,

$$
F_{\rm total}(R)=
\frac{E_e}{E_{\rm ret}}F_e(R)+
\frac{E_{\rm ph}}{E_{\rm ret}}F_{\rm ph}(R).
$$

Conocer $F_{\rm total}$ y su $R_{90}$ no determina $F_{\rm ph}$. Una pequeña cola electrónica puede mover el percentil total; además, la energía fonónica y su espectro pueden seguir evolucionando. Si $f_e=E_e/E_{\rm ret}$ es conocido y $0\le F_e\le1$, se puede acotar

$$
\max\!\left(0,\frac{F_{\rm total}-f_e}{1-f_e}\right)
\le F_{\rm ph}\le
\min\!\left(1,\frac{F_{\rm total}}{1-f_e}\right),
$$

para $f_e<1$. Es una cota algebraica de información incompleta, no una corrección del perfil.

## 3. Primer trabajo: recuperar la referencia histórica sin simular

Este es el siguiente paso ligero y ejecutable como trabajo documental:

1. Archivar las páginas 27–28 y 32–36 del PDF local de A20 con su hash; extraer sus curvas vectoriales si existen. Si se digitalizan imágenes, registrar resolución, transformación de ejes, grosor de líneas y zonas ocultas por leyendas.
2. Recuperar de la figura 2.6 curvas electrónicas y fonónicas en los instantes legibles, de la 2.7 los radios de energía **total**, y de la 2.8 las fracciones dentro de radios indicados. Conservar separado cada panel cilíndrico/esférico y la normalización que dice la figura. Contrastar solapamientos entre figuras; no rellenar detalles ocultos con números supuestos.
3. Convertir el reloj con el $\tau_0$ cinético de 1.87 ns del capítulo, sin emplear la escala GL. Estimar un intervalo de lectura a partir de resolución/grosor y comprobar al menos puntos representativos con una segunda extracción.
4. Entregar CSV de puntos extraídos, inventario de lo recuperable y una comparación de los percentiles y del reparto. Si una curva acaba antes de alcanzar una meseta, conservar energía exterior desconocida: no normalizar por su último punto como si fuese la energía total.
5. Usar esa reconstrucción para decidir qué información adicional requiere la transferencia. **No** multiplicar sus radios o tiempos por cocientes de energía, D o temperatura para presentarlos como curvas de Korzh.

La digitalización no recupera automáticamente un pico, un segundo momento completo ni un espectro. Derivar una densidad a partir de una curva acumulada amplifica el error; sólo se hará si la resolución lo permite y se propagará esa incertidumbre. Para una cola no mostrada, los momentos se reportarán parciales o acotados. Si las curvas no bastan, el dato preferible es la salida numérica original del autor; no se enviarán solicitudes externas sin autorización.

**Salidas previstas**, aún no generadas: a20_radial_digitization.csv, a20_radial_digitization.json y a20_transfer_diagnostic.md. Son referencia de su propio caso histórico. La primera entrega debe poder decir también “este observable no puede recuperarse”, sin inventar un ajuste.

### 3.1. Primer diagnóstico documental ya obtenido

Se extrajeron dos pares de radios de la imagen RGB original de la figura 2.7(a), de 1750×1315 píxeles, incrustada en el PDF. No existen curvas vectoriales en esa página. El [registro de extracción](cascade_digitization.json) conserva hashes, colores, ejes y píxeles utilizados; el script reproducible es [digitize_a20_radial.py](../../../../sandbox/stage3_5/research_20260923/digitize_a20_radial.py). No se recalculó la cascada.

| Tiempo de A20, ps | $R_{50}$, nm | $R_{90}$, nm | $R_{90}/R_{50}$ | Envolvente de lectura del cociente |
|---:|---:|---:|---:|---:|
| 0.0561 | 1.603 | 3.556 | 2.218 | 2.147–2.294 |
| 0.0935 | 1.635 | 3.734 | 2.284 | 2.212–2.360 |

Una gaussiana 2D tiene necesariamente

$$
\frac{R_{90}}{R_{50}}=\sqrt{\frac{\ln 10}{\ln 2}}\simeq1.823.
$$

Ese valor queda fuera de ambas envolventes. Por tanto, una sola gaussiana no reproduce simultáneamente estos dos percentiles de **energía total** del caso cilíndrico A20 a esos instantes. Ajustar sólo $R_{90}$ ocultaría la concentración mayor del núcleo. La conclusión no identifica el perfil fonónico ni prueba que una reducción gaussiana futura sea inútil: cuantifica un desacuerdo que deberá conservarse en su presupuesto.

La envolvente combina grosor visible, ±2 píxeles de calibración de ejes y variación horizontal local. Es un presupuesto de lectura, no un intervalo estadístico ni incertidumbre física de A20. Cambiar el umbral de selección RGB entre 4, 12 y 20 no cambió los centros extraídos. El instante final $10^{-4}\tau_0$ no se reporta para $R_{50}$ ni para el cociente, debido a la leyenda y el borde: no se extrapola la curva oculta.

La [comparación de anchuras inferidas de cada percentil](figures/04_cascade_shape.png) muestra esta incompatibilidad con sus barras de lectura. Son dos conversiones del mismo gráfico histórico; no son dos nuevas simulaciones ni extremos de un rango de $s_\gamma$ aceptado.

## 4. Cuándo hace falta un cálculo reducido de cascada

La segunda vía se activa después del diagnóstico documental y la admisión material necesaria. El cálculo mínimo útil sigue la cinética espectral y el transporte radial suficiente para estimar los productos anteriores, no la cadena óptica completa.

- **Inicio energético:** identificar cómo se distribuye la energía de cada fotón entre excitaciones de alta energía, y qué incertidumbre deja omitir la primera multiplicación. Un radio electrónico impuesto en A20 es parte de su referencia, no un dato que deba copiarse a Korzh. Una fuente puntual regularizada numéricamente requeriría demostrar independencia del núcleo inicial a las escalas de transferencia; no se fija aquí una regularización física.
- **Colisiones:** conservar el intercambio electrón–fonón con su DOS/acoplamiento admitidos. Evaluar si ee altera el reparto o la extensión en la ventana de interés; el BGK ensayado no representa por sí solo una cascada de eV ni identifica sus tasas.
- **Transporte:** tratar primero los canales que llevan energía fuera del núcleo. Difusión electrónica no equivale a difusión de toda la energía. Si propagación fonónica o superficies compiten con el reloj de transferencia, introducir el cierre mínimo que permita acotarlas con velocidades, caminos libres y transmisión justificados.
- **Profundidad:** si la extensión es comparable o menor que 7 nm, el promedio en profundidad exige un control propio. Una alternativa reducida puede separar posiciones de absorción o resolver radio y profundidad; pasar directamente a un problema 3D completo sólo se justifica si ese control cambia los observables relevantes.
- **Campos:** un espectro de metal normal o campos congelados pueden servir para una etapa temprana de alta energía, pero se deben contrastar al aproximarse a las energías que suprimen el condensado. No se prolonga esa aproximación hasta una latencia de detección por defecto.

No hay todavía un programa de esta cascada radial admitido y listo para lanzar. Por eso no se escribe un comando largo ficticio. Antes de cualquier corrida se fijarán el modelo reducido, la procedencia de sus entradas, un piloto de coste y los criterios de lectura; los cálculos previstos de más de cinco minutos se entregarán al usuario según la política del repositorio.

## 5. Material y herramientas: qué falta realmente

El escenario K20 elegido para investigación mantiene $D=0.5$ cm²/s y el resto de su conjunto coherente, identificando los ajustes como tales. Eso no completa los datos de una cascada:

| Falta | Consecuencia práctica |
|---|---|
| Unidades/base y normalización absoluta del $g_{\rm ph}$ y $\alpha^2F$ disponibles | La forma condicional del archivo NbN permite comparar espectros; todavía no certifica capacidades ni tasas absolutas. |
| Interacciones ee y e–ph a las energías iniciales | Los tiempos efectivos KWT y BGK no sustituyen kernels de multiplicación y redistribución; tampoco basta una ley térmica extrapolada. |
| Escape/transmisión y transporte fonónico en la película/interfaz | Determinan simultáneamente energía retenida, profundidad y ensanchamiento; un único tiempo ajustado no los identifica por separado. |
| Fuente óptica efectiva por color y profundidad | No está identificado el perfil energético que deja la etapa de absorción omitida. El tamaño del foco óptico no es el tamaño de la cascada de un fotón. |
| Datos radiales a 0.9 K y 775/1550 nm | Las curvas A20 verificadas no llenan este hueco mediante una conversión algebraica. |

Los módulos experimentales actuales ofrecen energía electrónica común, eventos locales conservativos e intercambio entre celdas. No constituyen un solucionador admitido de multiplicación de alta energía, propagación fonónica ni cascada con superficies. El catálogo con coordenada de conteo normalizada $x_{\max}=12$ alcanza aproximadamente 15.8 meV en el estado de referencia a corriente cero; ese12 no significa literalmente $E_{\max}/\Delta_0$. La rejilla fonónica sintética hasta $4\Delta_0$ alcanza aproximadamente 5.26 meV. El soporte positivo del archivo NbN condicional llega aproximadamente a 67.8 meV, y la excitación óptica está en escala eV. Aumentar el número de nodos sin ampliar/admitir el soporte no resuelve esta diferencia. Véanse [material_kinetics.md](material_kinetics.md) y su registro de fuentes.

## 6. Criterios para decidir una transferencia útil

Son **criterios propuestos para registrar antes de futuros resultados**, no nuevas tolerancias obligatorias ni una certificación:

1. **Balance comprensible:** en un cálculo de cascada sin trabajo externo, contabilizar $E_e+E_{\rm ph}+E_{\rm esc}+E_{\rm salida\ lateral}$ respecto a $E_\gamma$. Si evolucionan condensado o corriente, añadir su variación de energía y el trabajo de puertos una sola vez. El residual debe ser pequeño frente a la energía cuya partición o pérdida se pretende distinguir.
2. **Soporte suficiente:** estimar energía y efecto cinético de las colas excluidas y de cada frontera numérica. Un flujo cero por construcción fuera de la malla no demuestra que la cola física sea despreciable.
3. **Perfil resoluble:** al refinar tiempo, radio, profundidad cuando proceda y energía, el cambio en reparto, radios, momento y pico no debe alterar la decisión de transferencia ni superar la incertidumbre de las fuentes. El pico requiere su propio control espacial.
4. **Reloj identificable:** escoger una ventana donde la reducción efectiva conserve la evolución relevante al iniciar el modelo posterior. Desplazar la transferencia dentro de esa ventana debe dar resultados compatibles al contabilizar el tiempo omitido. No se exige que electrones y fonones hayan adquirido temperaturas térmicas iguales.
5. **Compatibilidad con D.30:** comprobar si los electrones pueden comenzar en el baño, si la energía retenida puede asignarse a fonones y si el espectro permite la forma separable propuesta. Una fracción electrónica relevante exige revisar explícitamente la transferencia, no ocultarla en $\chi$ o $s_\gamma$.
6. **Forma gaussiana evaluada después:** comparar el perfil obtenido con una gaussiana conservando energía y un criterio declarado —por ejemplo, segundo momento—, y medir el desacuerdo en $R_{50}$, $R_{90}$, pico y colas. Si ese desacuerdo afecta hotbelt o latencia relativa más que el objetivo útil del ensayo, conservar el perfil tabulado o revisar la reducción. No ajustar la gaussiana para forzar la respuesta deseada.

La salida de esta caracterización será un estado de transferencia por color, con sus incertidumbres y alcance. **Hasta entonces, $s_\gamma$, $\tau_{\rm tr}$ y el reparto de entrada permanecen abiertos.** La decisión de estudiar Korzh y mantener una descripción espacial gaussiana como candidata no sustituye esos datos.
