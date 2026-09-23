# Circuito y observables para el caso de 80 nm

Investigación documental, 23 de septiembre de 2026. Se adopta la elección comunicada del caso Korzh de 80 nm y del contraste 775/1550 nm, estudiando la formación del hotbelt a partir de un depósito localizado. **Se conserva el circuito CM de la memoria**; la elección de parámetros materiales ajustados de K20 no autoriza cambiarlo por el montaje de 96 nH. Aquí se resuelve la procedencia de la lectura; no se ejecutan transientes ni se cambia el circuito implementado. El registro de fuentes y decisiones está en [circuit_sources.json](circuit_sources.json).

## Qué representa el circuito disponible

La adenda [CM de la memoria](../../../modelo_v0_4/actualizaciones/circuito_memoria_20260922.md) define tres **variables de estado circuitales**: las corrientes de las ramas de polarización y detector, $I_b$ e $I_s$, y la tensión del condensador de acoplamiento, $v_c$. $I_s$ incluye corriente normal y superconductora. Con orientación pasiva izquierda–derecha:

$$
I_{\rm RF}=I_b-I_s,\quad V_{\rm out}=R_L I_{\rm RF},\quad V_d=v_c+V_{\rm out},\quad V_{\rm dev}=\phi_L-\phi_R,
$$
$$
L_b\dot I_b=V_b-R_bI_b-v_c-R_L(I_b-I_s),\qquad
L_{\rm ext}\dot I_s=v_c+R_L(I_b-I_s)-V_{\rm dev},\qquad
C\dot v_c=I_b-I_s.
$$

$V_{\rm dev}$ se mide sobre **todo el dominio cuya energía ya se resuelve**, incluidas sus continuaciones 1D. $V_d$ está antes del condensador; $V_{\rm out}$ es la señal en la carga. La analogía útil es una contabilidad de tres depósitos: cada inductancia almacena energía asociada a una corriente y el condensador almacena energía asociada a una tensión. La señal aparece cuando las corrientes dejan de coincidir; no es directamente la caída interna del nanohilo.

| Parámetro | Significado físico | Escenario CM de la memoria | Restricción |
|---|---|---:|---|
| $R_b$ | Resistencia serie que polariza desde $V_b$ | 10 kΩ | Su disipación es externa al film. |
| $L_b$ | Inductancia de la rama de polarización | 1 μH | No es la inductancia cinética del detector. |
| $R_L$ | Carga de lectura después de $C$ | 50 Ω | No equivale a una impedancia de entrada medida en todas las frecuencias. |
| $C$ | Acoplamiento capacitivo entre $V_d$ y la carga | 100 pF | Bloquea la componente continua de esa rama. |
| $V_b$ | Tensión de la fuente del modelo | 0.300 V histórica | Se fija por el equilibrio elegido; no conservarla automáticamente al cambiar de corriente. |
| $L_{\rm total,mem}$ | Referencia inductiva total del escenario histórico | 10 nH | Se particiona entre material resuelto y exterior; no se suma íntegra otra vez. |
| $L_{\rm ext}$ | Inductancia situada fuera de la energía espacial resuelta | Derivada en la referencia | Constante positiva; no recalcularla a partir del hotspot instantáneo. |

Estos valores proceden de la memoria, tabla 4.3 y ecuaciones 4.16–4.18 (PDF 85, 114–116), y de CM. No son mediciones del montaje de Korzh. En un equilibrio superconductivo, $V_{\rm dev}=v_c=V_{\rm out}=0$ e $I_b=I_s=V_b/R_b$. Así, **si se conserva el escenario CM** de 10 kΩ, 15.5 y 21.5 μA requieren 0.155 y 0.215 V; son conversiones del modelo, no tensiones experimentales identificadas.

La identidad que mantiene esa interpretación es

$$
\frac{d}{dt}\left(\frac{L_b I_b^2}{2}+\frac{L_{\rm ext}I_s^2}{2}+\frac{Cv_c^2}{2}\right)
+I_sV_{\rm dev}+R_bI_b^2+R_L(I_b-I_s)^2=V_bI_b.
$$

El término $I_sV_{\rm dev}$ transfiere potencia al dominio espacial. Al unir su balance, el intercambio interno $Q_\Delta$ no se agrega como una segunda fuente externa. Esta comprobación es algebraica; por sí sola no demuestra que la red reproduzca un amplificador experimental.

## Correspondencia experimental y partición inductiva

La referencia seleccionada es NbN de 80 nm, longitud activa 5 μm, espesor nominal 7 nm y baño 0.9 K. A 21.5 μA, la diferencia publicada 1550−775 es 4.2 ± 0.4 ps; también existe la comparación a 15.5 μA. No corresponde al dispositivo del récord absoluto. [Korzh 2020, PDF 4–5 y 22–24](https://eprints.lancs.ac.uk/id/eprint/140252/3/Binder1.pdf).

El montaje añade un inductor de 1.5 mm × 1 μm: 96 nH **estimados**, usando 64 pH/cuadro de películas similares. Lectura: CITLF1 a 4 K, ganancia nominal 50 dB, banda nominal 1.5 GHz, corte inferior 1 MHz y temperatura de ruido nominal <7 K; osciloscopio de 40/80 GS/s, banda analógica elegida 6 GHz; subida 20–80% de unos 80 ps. [Korzh, métodos y figura 2, PDF 10–12 y 23](https://eprints.lancs.ac.uk/id/eprint/140252/3/Binder1.pdf).

La tesis de Allmaras confirma que ese inductor evita el enclavamiento resistivo (*latching*) y distingue el protocolo de latencia de las medidas de jitter. Sus otros experimentos a 1/4 K y sus parámetros de simulación no sustituyen las condiciones anteriores. [Allmaras 2020, pp. 63–65 y 95; PDF 75–77 y 107](https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf).

La correspondencia defendible es estructural: una rama detectora, polarización y lectura. **No se ha identificado una equivalencia completa entre la red CM y el bias-T resistivo del experimento.** Faltan componentes y transferencia compleja de ese montaje. En particular:

- Los 96 nH describen el inductor añadido, no toda la rama de detector, ni $L_b$. Si todo él queda fuera del dominio espacial, contribuye íntegro a $L_{\rm ext}$.
- Aplicar la misma estimación de hoja a la parte activa da $64\,\mathrm{pH}\times(5\,\mu\mathrm m/80\,\mathrm{nm})=4\,\mathrm{nH}$. Añadido más tramo activo suman aproximadamente 100 nH **antes** de tapers y conexiones, dentro de esa aproximación. No es una nueva medición ni la inductancia diferencial certificada a 21.5 μA.
- En una reducción de elementos concentrados, $L_{\rm ext}=L_{\rm añadido}+L_{\rm activo\ no\ resuelto}+L_{\rm transiciones/conexiones}$, evaluada en la referencia. Si se dispone de una inductancia total consistente, también puede usarse $L_{\rm ext}=L_{\rm total}-L_{\rm res,dif}$. Son dos contabilidades alternativas de los mismos segmentos.
- El ejemplo anterior de etapa 3, $L_{\rm res,dif}\simeq0.33324$ nH y $L_{\rm ext}\simeq9.66676$ nH, pertenece al escenario espacial de 120 nm/720 nm y referencia de 10 nH. No se transporta al nanohilo experimental de 80 nm/5 μm.

Un intervalo 10–100 nH mezclaría dispositivos y procedencias; **no sería un rango físico de incertidumbre**. La incertidumbre de 96 nH y su dependencia de corriente no están cuantificadas en la fuente consultada. Tampoco la estimación estática acredita una impedancia inductiva ideal a cualquier frecuencia. Para la primera comparación de latencia se conserva esa limitación explícita, sin exigir por adelantado un modelo electromagnético distribuido de todo el montaje.

La referencia material ajustada elegida ($D=0.5$ cm²/s y $R_\square=608$ Ω) no convierte la inductancia de hoja experimental estimada en un parámetro microscópico compatible por definición. Una resistencia de ajuste, una inductancia extraída de películas similares y la inductancia diferencial del catálogo tienen procedencias distintas. La variante exterior de 96 nH queda como investigación separada; el circuito activo de referencia conserva los 10 nH totales de CM y su partición coherente.

## Conflictos heredados que quedan resueltos documentalmente

**100 pF frente a 1 pF.** CM y los valores por defecto de `CircuitParams`, `ThesisCircuitParameters` y los dos corredores históricos SS/fotón usan 100 pF. Los YAML contienen `C_rf_F=1e-12`, pero la búsqueda del código fuente encuentra esa clave en la validación del YAML, sin conexión a `C_couple_F` en los corredores inspeccionados: estos construyen el circuito desde `--circuit-Ccouple-pF`, cuyo defecto es 100. Por tanto, no hay evidencia para afirmar que las trayectorias históricas usaron 1 pF sólo porque aparece en el YAML. La autoridad de una corrida es su `circuit.params` guardado. Se documenta 1 pF como entrada heredada sin correspondencia demostrada con el circuito ejecutado; no se promueve a nuevo valor físico ni se modifica ningún archivo histórico.

**10 nH frente a 96 nH.** Son referencias de escenarios diferentes y además significan total histórico frente a componente añadido experimental. Se conserva CM como banco de desarrollo identificado. Una parametrización experimental posterior deberá declarar los segmentos y la referencia de $L_{\rm res,dif}$; no basta cambiar el número de `Lk_ext_H`.

**100 μV frente al umbral de Korzh.** La memoria y el corredor fotónico usan un cruce de $V_{\rm out}$ relativo a la línea basal; el corredor exige 0.5 ps de confirmación por defecto. No se ha localizado el valor absoluto del umbral experimental. No convertir 100 μV en un ajuste medido ni trasladarlo sin más a la señal amplificada.

## Escalas útiles, sin confundirlas con polos del circuito completo

Las siguientes conversiones son analíticas; sirven para planificar observación y almacenamiento. Los polos reales dependen de la red completa y de la respuesta dinámica del dispositivo.

| Cálculo orientativo | Valor | Qué permite anticipar |
|---|---:|---|
| $L_b/R_b$ para CM | 100 ps | Escala aislada de su rama de polarización. |
| $R_LC$ y $1/(2\pi R_LC)$, 100 pF/50 Ω | 5 ns; 31.83 MHz | Escala de acoplamiento bajo esa terminación simplificada. |
| Las mismas expresiones, 1 pF/50 Ω | 50 ps; 3.183 GHz | El cambio heredado es un factor 100; no una corrección despreciable. |
| $L/R_L$, 10 nH/50 Ω | 200 ps | Estimación elemental de retorno de corriente. |
| $L/R_L$, 96 nH/50 Ω | 1.92 ns | Escala del inductor añadido con carga ideal supuesta; no tiempo de reset medido. |
| $1/f_s$, 40/80 GS/s | 25/12.5 ps | Espaciado de adquisición, no límite directo del ajuste de tiempos ni paso del integrador. |
| $10^{50/20}$ | 316.23 | Ganancia de tensión equivalente bajo igual referencia de impedancia; no transferencia plana demostrada. |

El corte inferior nominal de 1 MHz del amplificador **no determina $C$**: invertir $1/(2\pi R_LC)$ sin conocer la red atribuiría su filtrado a un componente distinto. Tampoco los 6 GHz del osciloscopio reemplazan los 1.5 GHz nominales del amplificador. La duración de un pulso, su subida, su latencia y su reset son observables diferentes; no debe cerrarse el reset cuando sólo se observó el comienzo de la señal.

## Observable inicial y presupuesto numérico proporcional

Se acepta como origen temporal la **transferencia de energía al modelo**, $t_{\rm tr}$, no la absorción óptica, $t_{\rm abs}$. El retardo de cascada $\tau_{\rm casc}=t_{\rm tr}-t_{\rm abs}$ queda desconocido y separado. Por tanto, una latencia óptica contendría $\tau_{\rm casc}+\tau_{\rm modelo}$, y su diferencia entre longitudes de onda contendría también $\tau_{{\rm casc},1550}-\tau_{{\rm casc},775}$. No se supone que ese término sea cero ni que quede acotado por el presupuesto numérico de 0.1 ps.

Para cada longitud de onda se registrarán por separado $t_{\rm tr}$, $V_{\rm dev}$, $V_{\rm out}$ y las tres variables circuitales. Una cadena de medida, cuando esté identificada, produciría

$$u_\lambda(t)=(h_{\rm lectura}*V_{{\rm out},\lambda})(t)+n(t),\qquad
t_{\times,\lambda}=\inf\{t\ge t_{\rm tr}:u_\lambda(t)-u_{\rm basal}=U_*\}-t_{\rm tr},
$$

donde $h_{\rm lectura}$ es su respuesta impulsional, $n$ el ruido y $U_*$ un umbral fijo en un puerto explícito. Deben acompañarlo polaridad, regla de confirmación e interpolación. Sin transferencia conocida se informa **cruce del modelo** en $V_{\rm out}$, sin fingir señal de osciloscopio. Ausencia de cruce antes del horizonte significa dato censurado, no latencia igual al horizonte.

Korzh obtiene la latencia relativa del desplazamiento de máximos de IRF ajustadas, manteniendo ajustes y verificando pulsos equivalentes. El umbral numérico exacto no está dado. [Korzh, PDF 12](https://eprints.lancs.ac.uk/id/eprint/140252/3/Binder1.pdf). Así se distinguen

$$\Delta t_\times=t_{\times,1550}-t_{\times,775},\qquad
\Delta t_{\rm IRF}=\arg\max P_{1550}(t)-\arg\max P_{775}(t).$$

Dos transientes deterministas permiten la primera cantidad. Para obtener la segunda hace falta una distribución de eventos: si $T(E)$ transforma energía retenida en tiempo y es monótona, $P_t(t)=P_E(E(t))/|T'(E(t))|$. El máximo resultante no coincide, en general, con $T(\langle E\rangle)$. El jitter FWHM requiere esa distribución; no se obtiene de la separación entre dos curvas ni de su error numérico. Esta es también la distinción pedagógica de [Allmaras, pp. 60 y 64–65 (PDF 72 y 76–77)](https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf).

Una cadena común cancela un retardo fijo, pero no garantiza cancelar diferencias de forma: una pendiente distinta produce *time walk*, es decir, desplazamiento del cruce dependiente de la amplitud o forma. La linealización $\delta t\simeq\delta u/|du/dt|_{U_*}$ muestra por qué hay que mirar pendiente y puerto. Si las curvas son traslaciones de una misma forma, el desplazamiento relativo sí sobrevive a un filtrado lineal común.

Se propone **0.1 ps como objetivo absoluto del error numérico estimado de $\Delta t_\times$**: aproximadamente 2.4% del contraste de 4.2 ps y un cuarto de la escala publicada ±0.4 ps. Es un presupuesto de ingeniería propuesto antes de nuevos resultados, no una barra estadística, ni una precisión física garantizada. **0.2 ps** puede servir para una exploración inicial, etiquetada como tal; consume la mitad de aquella escala y no sustenta afirmar que se resolvieron diferencias de 0.1 ps. Esta elección no obliga a cada trayectoria a tener error 0.05 ps, ni prescribe $dt=0.1$ ps: se comprueba el observable combinado mediante refinamiento y localización del cruce, registrando también cada tiempo para detectar cancelaciones accidentales.

Las colas espaciales al 1%/0.1%, si se usan para delimitar el hotbelt o alejar bordes, son **indicadores geométricos propuestos**, no probabilidades de detección ni garantías automáticas sobre la latencia. Su utilidad consiste en reducir el trabajo a contrastes con efecto observable; no sustituyen medir si el corte altera $\Delta t_\times$.

## Decisiones pendientes que sí afectan al panorama general

El caso de 80 nm y la prioridad de latencia relativa están elegidos. Para comenzar investigación del hotbelt no hace falta pedir valores arbitrarios de componentes. Se puede avanzar con observables internos y con el circuito CM identificado como banco de desarrollo. Antes de presentar una reproducción eléctrica de Korzh hay que resolver:

1. **Comparación de detector o reproducción del instrumento.** La primera admite una comparación cualitativa/condicional de latencia determinista; la segunda necesita $Z_{\rm entrada}(\omega)$, transferencia/ruido, umbral y trazas de la adquisición seleccionada. Son objetivos de alcance distinto. La prioridad ya indicada favorece el primero por ahora.
2. **Reducción del exterior.** Documentar dónde terminan los 5 μm activos y qué tramos quedan en $L_{\rm ext}$. Obtener geometría de transiciones o conservar su contribución como incertidumbre no cuantificada; no inventar un margen porcentual.
3. **Datos adicionales.** Si se desea igualdad cuantitativa con el máximo IRF, harán falta eventos o ajustes y una distribución de energía/posiciones justificada. Hasta entonces, no ajustar componentes, umbral o ruido para forzar 4.2 ps.

La limitación de información circuital no invalida la mejora del modelo de SNSPD; delimita qué parte de la comparación puede atribuirse a su física. El siguiente registro debe separar precisión del cálculo, incertidumbre de entradas y diferencia entre el observable determinista y el experimental.
