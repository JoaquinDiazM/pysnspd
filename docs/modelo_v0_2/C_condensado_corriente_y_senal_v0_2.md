---
title: "C. Condensado, corriente y señal"
subtitle: "Reestructuración física de pySNSPD · Documento técnico 3 de 4"
author: "Desarrollo del modelo pySNSPD"
date: "8 de septiembre de 2026 · Revisión 0.2"
lang: es
---

# C.0. Qué se conserva y qué se sustituye

Se conserva una dinámica espacial compleja del condensado, la conservación local de corriente, el régimen difusivo, la respuesta subkelvin como objetivo físico y el acoplamiento eléctrico que determina el pulso antes de amplificación. No se exige resolver un meandro completo ni interfaces metálicas que estén lejos de la absorción.

Se conserva la motivación del cierre Vodolazov–Allmaras: utilizar una corriente de superfluido más adecuada a baja temperatura y hacerla compatible con la ecuación de fase. Vodolazov explica el término correctivo que permite esa compatibilidad [V, p. 11]. La propuesta de A permite dar un paso adicional: obtener la fuerza local de amplitud y la corriente como derivadas de un mismo funcional reducido. Esto no invalida la justificación de la corriente usada en la memoria; organiza la extensión energética y espectral que se desea evaluar. Se retiene una movilidad inspirada en Allmaras/Kramer–Watts–Tobin como cierre disipativo efectivo. La parte termodinámica uniforme queda fijada microscópicamente dentro del modelo de A; la movilidad y la aproximación de gradientes tienen un estatuto distinto y necesitan pruebas propias.

No se condiciona el programa a operar a 4 K. La modificación de Vodolazov mejora el comportamiento de equilibrio a baja temperatura y reproduce bien la corriente de depairing [V]. La precisión del transitorio completo debe demostrarse en el régimen subkelvin de interés, sin convertir una reserva sobre su derivación en una prohibición de modelarlo.


**Cómo leer C.** A fija el espectro y las fuerzas; B decide qué información de las poblaciones se conserva. Aquí se conectan ambas piezas con el movimiento de $\Delta$, el trabajo eléctrico y el voltaje medido. C.1–C.3 dicen cómo se mueven amplitud, fase y corriente; C.4 comprueba dónde va la energía; C.9 muestra por qué el circuito también determina la señal. La revisión 0.2 conserva el cierre dinámico propuesto y añade derivaciones y pruebas pequeñas, sin modificar el solver de producción.

Un resorte con rozamiento ofrece una analogía limitada pero útil: la energía define hacia dónde empuja el resorte y el rozamiento decide qué tan rápido se mueve. Aquí el funcional define la fuerza del condensado; la movilidad de C.6 fija su ritmo. La energía que pierde por esa disipación debe reaparecer en las excitaciones o salir por un flujo identificado. La analogía no representa pares como objetos mecánicos ni sustituye el espectro de A.

# C.1. Funcional espacial reducido

## C.1.1. Separación entre respuesta uniforme y gradientes

El propósito de esta sección es construir una energía de la región espacial resuelta que produzca tanto fuerzas locales como flujos de trabajo. No se añade un segundo depósito de energía independiente al de A. En el régimen de electrones térmicos se propone

$$
\begin{gathered}
\mathcal F[\Delta,T_e,\mathbf A]=\int_{\mathcal V}
\left[f_e^{\rm FD}(T_e,|\Delta|,q)+K_0|\nabla |\Delta||^2\right]dV,\\
K_0=\frac{\pi N_0\hbar D}{8k_BT_c}.
\end{gathered} \tag{C.1}
$$

La dependencia de superflujo ya está dentro de $f_e^{\rm FD}$: añadir además $K_0|\Delta|^2q^2$ contaría dos veces esa contribución. En el intervalo no térmico se utiliza

$$
\mathcal U[\Delta,p]=\int_{\mathcal V}
\left[U_{\rm vac}(|\Delta|,q)+4N_0\int E(x;|\Delta|,q)p(x)dx
+K_0|\nabla |\Delta||^2\right]dV. \tag{C.2}
$$

La variación no térmica se realiza a ocupaciones $p(x)$ fijas. En la variante térmica se realiza a $T_e$ fijo. Las fuerzas coinciden cuando las ocupaciones son Fermi–Dirac.

La elección de $K_0$ conserva el término dominante de gradiente GL y permite reproducir la escala $\xi_{\rm mod}$ de Allmaras con la movilidad apropiada. **No es una derivación exacta de la rigidez espacial de Usadel a cualquier temperatura, amplitud y corriente.** Esa distinción se cuantifica en C.8.

## C.1.2. Derivadas funcionales

Sea $X_{|\Delta|}$ la fuerza uniforme de A y $\boldsymbol\Pi=\hbar\mathbf j_s/(2e)$. Para $K_0$ constante,

$$
\begin{gathered}
\mathcal Q_{|\Delta|}=X_{|\Delta|}-2K_0\nabla^2|\Delta|,\\
\mathcal Q_\theta=-\nabla\cdot\boldsymbol\Pi.
\end{gathered} \tag{C.3}
$$

La primera identidad combina la fuerza local con el coste de que la amplitud cambie entre puntos vecinos. La segunda identidad resulta de $\delta\mathbf q=\nabla\delta\theta$ y de una integración por partes. Los términos de borde no se descartan físicamente: especifican el trabajo o las condiciones naturales de la región acoplada al exterior.

En variables complejas, usando la derivada de Wirtinger,

$$
\frac{\delta\mathcal F}{\delta\Delta^*}
=\frac{e^{i\theta}}2\left(\mathcal Q_{|\Delta|}+
\frac{i\mathcal Q_\theta}{|\Delta|}\right). \tag{C.4}
$$

C.4 es una representación polar donde $|\Delta|>0$; no autoriza dividir por cero en un núcleo normal. La construcción regular del límite se trata en C.7.

# C.2. Dinámica variacional con movilidad Allmaras/KWT

## C.2.1. Amplitud y fase

Defínanse

$$
\begin{gathered}
\tau_0=\frac{\pi\hbar}{8k_BT_c},\\
A_0(T_E)=N_0\sqrt{\frac{1+T_E/T_c}{2}},\\
R=\sqrt{1+\frac{4|\Delta|^2\tau_\psi^2}{\hbar^2}},\\
w=\dot\theta+\frac{2e}{\hbar}\phi.
\end{gathered} \tag{C.5}
$$

$\tau_\psi$ representa el tiempo efectivo que entra en la movilidad del condensado. Se distingue de un tiempo de pérdida de energía y del tiempo de termalización explícito de B. La propuesta de dinámica es

$$
\begin{gathered}
\Gamma_{|\Delta|}\partial_t|\Delta|=-\mathcal Q_{|\Delta|},\\
\Gamma_\theta w=-\mathcal Q_\theta=\nabla\cdot\boldsymbol\Pi,\\
\Gamma_{|\Delta|}=2A_0\tau_0R,\\
\Gamma_\theta=\frac{2A_0\tau_0|\Delta|^2}{R}.
\end{gathered} \tag{C.6}
$$

La energía uniforme y las fuerzas se derivaron en A. Los coeficientes positivos de C.6 constituyen una **elección de dinámica disipativa**, normalizada para conectar con Allmaras. No se presenta esta movilidad como una solución exacta del problema quasiclásico no estacionario a todas las frecuencias.

Al dividir la ecuación de amplitud por $2A_0$,

$$
\begin{gathered}
\tau_0R\partial_t|\Delta|=\frac{K_0}{A_0}\nabla^2|\Delta|-\frac{X_{|\Delta|}}{2A_0},\\
\frac{K_0}{A_0}
=\frac{\pi\hbar D}{4\sqrt2 k_BT_c\sqrt{1+T_E/T_c}}
=\xi_{\rm mod}^2(T_E).
\end{gathered} \tag{C.7}
$$

El término de gradiente conserva exactamente la escala modificada antigua. La diferencia reside en la fuerza local: se usa $X_{|\Delta|}$ del funcional en lugar del polinomio interpolado de la memoria.

Para la fase,

$$
\frac{\tau_0|\Delta|}{R}w=\frac\hbar{4eA_0|\Delta|}\nabla\cdot\mathbf j_s
=\frac{\hbar eD}{\sigma_n\sqrt2\sqrt{1+T_E/T_c}}
\frac{\nabla\cdot\mathbf j_s}{|\Delta|}. \tag{C.8}
$$

Así se recupera la normalización de la ecuación de fase de Allmaras utilizada en el anexo C de la memoria. La misma estructura de fase se obtiene ahora por variación. El cambio propuesto afecta al origen común de las derivadas y al balance, no a la validez de usar la corriente de superfluido de Usadel en el cierre anterior.

## C.2.2. Forma compleja equivalente

Como $R^2-1=4|\Delta|^2\tau_\psi^2/\hbar^2$ y $\partial_t(|\Delta|^2)=2|\Delta|\partial_t|\Delta|$,

$$
\frac{A_0\tau_0}{R}\left[
\left(\partial_t+\frac{2ie}{\hbar}\phi\right)\Delta
+\frac{2\tau_\psi^2}{\hbar^2}(\partial_t(|\Delta|^2))\Delta\right]
=-\frac{\delta\mathcal F}{\delta\Delta^*}. \tag{C.9}
$$

En el régimen no térmico se reemplaza la derivada de $\mathcal F$ por la de $\mathcal U$ a $p$ fijo. Para comprobar C.9, se multiplica por $e^{-i\theta}$: la parte real del corchete vale $R^2\partial_t|\Delta|$ y la imaginaria vale $|\Delta|w$. Se obtienen las dos filas de C.6 sin introducir un segundo cierre de corriente.

El uso de $T_E$ en los coeficientes parametriza la movilidad mediante la energía local. Esta elección efectiva se mantiene en la revisión 0.2; no obliga a evolucionar una ecuación de temperatura. La respuesta a la forma de $f$ no queda reducida a $T_E$: permanece en $X_{|\Delta|}[p]$ y $\mathbf j_s[p]$.

# C.3. Corriente y potencial

En el sector cuasineutro,

$$
\begin{gathered}
\mathbf E=-\nabla\phi-\partial_t\mathbf A,\\
\mathbf j_n=\sigma_n\mathbf E,\\
\nabla\cdot(\mathbf j_s+\mathbf j_n)=0.
\end{gathered} \tag{C.10}
$$

Para $\mathbf A=0$ y $\sigma_n$ constante,

$$
\sigma_n\nabla^2\phi=\nabla\cdot\mathbf j_s. \tag{C.11}
$$

La corriente empleada en C.11 debe ser la misma derivada de energía que conduce a C.8. Conservar carga no demuestra por sí solo conservar energía. Si se emplean corrientes diferentes, debe identificarse también el trabajo de su diferencia; sin él, no se obtiene el balance C.18.

La aproximación $\mathbf A=0$ desprecia el campo propio magnético en la región local. **No prohíbe vórtices ni phase slips deterministas**: la fase de $\Delta$ puede adquirir enrollamiento y la amplitud puede anularse. Se distingue la existencia de estas soluciones de la incorporación de nucleación estocástica, que no es necesaria para el alcance determinista de esta primera publicación.

# C.4. Derivación del balance energético y del calentamiento

## C.4.1. Identidad cinemática gauge-invariante

De la definición de $\mathbf q$ y $w$,

$$
\dot{\mathbf q}=\nabla w+\frac{2e}{\hbar}\mathbf E. \tag{C.12}
$$

\Needspace{6\baselineskip}

Defínase la densidad electrónica total local, incluida la energía de gradiente,

$$
u_* =u_e+K_0|\nabla |\Delta||^2.
$$

A ocupaciones $p$ fijas, la contribución del condensado a su derivada temporal es

$$
X_{|\Delta|}\partial_t|\Delta|+\boldsymbol\Pi\cdot\dot{\mathbf q}
+2K_0\nabla |\Delta|\cdot\nabla\partial_t|\Delta|. \tag{C.13}
$$

Se usan dos integraciones por partes locales:

$$
2K_0\nabla |\Delta|\cdot\nabla\partial_t|\Delta|
=\nabla\cdot(2K_0\partial_t|\Delta|\nabla |\Delta|)-2K_0\partial_t|\Delta|\nabla^2|\Delta|,
$$

$$
\boldsymbol\Pi\cdot\nabla w
=\nabla\cdot(\boldsymbol\Pi w)-w\nabla\cdot\boldsymbol\Pi.
$$

Con C.6 y C.12 resulta

$$
\left.\dot u_*\right|_{\rm campos}
=-Q_\Delta+\mathbf j_s\cdot\mathbf E
+\nabla\cdot\mathbf J_{\rm tr}, \tag{C.14}
$$

$$
\begin{gathered}
Q_\Delta=\Gamma_{|\Delta|}(\partial_t|\Delta|)^2+\Gamma_\theta w^2\geq0,\\
\mathbf J_{\rm tr}=2K_0\partial_t|\Delta|\nabla |\Delta|+\boldsymbol\Pi w.
\end{gathered} \tag{C.15}
$$

$Q_\Delta$ es la disipación de la dinámica elegida, no una nueva inyección externa. En una celda cerrada, verla con signo negativo en C.14 y positivo en C.16 equivale a una transferencia entre dos cuentas: disminuye la energía asociada a los campos y aumenta la de las ocupaciones, sin crear un depósito adicional. $\mathbf J_{\rm tr}$ transporta trabajo asociado a los campos; no se debe confundir su signo con el flujo de calor saliente.

## C.4.2. Suma con la energía de las ocupaciones

La energía transferida a las ocupaciones por calor y disipación es

$$
4N_0\int E\dot p\,dx=-\nabla\cdot\mathbf Q_e
+\sigma_n E^2+Q_\Delta-P_{e\text{-ph}}+S_{\gamma,e}^{(E)}.
\tag{C.16}
$$

Al sumar C.14, se cancelan las dos apariciones de $Q_\Delta$:

$$
\dot u_*=-\nabla\cdot\mathbf Q_e
+(\mathbf j_s+\mathbf j_n)\cdot\mathbf E
-P_{e\text{-ph}}+\nabla\cdot\mathbf J_{\rm tr}
+S_{\gamma,e}^{(E)}. \tag{C.17}
$$

Finalmente, al añadir la fila fonónica de B,

$$
\partial_t(u_*+u_{\rm ph})=
-\nabla\cdot\mathbf Q_e+\nabla\cdot\mathbf J_{\rm tr}
+\mathbf j_{\rm tot}\cdot\mathbf E-P_{\rm esc}+S_\gamma^{(E)}.
\tag{C.18}
$$

Esta es la identidad de conservación del cierre propuesto. Si se incorporan transporte fonónico o campos magnéticos resueltos, deben añadirse sus energías y flujos correspondientes; C.18 no pretende contenerlos tácitamente.

## C.4.3. Por qué aparece una temperatura y qué se puede evolucionar

**La variable conservada puede ser la energía.** Escribir una ecuación para $T_e$ no afirma que la energía sea una elección peor ni obliga a usar dos temperaturas desde el instante de absorción. Hay dos operaciones diferentes: elegir una coordenada para la energía y aproximar toda una distribución por Fermi–Dirac. La primera es invertible bajo las condiciones de B.20–B.21; la segunda pierde información y necesita los criterios de B.7.

\Needspace{6\baselineskip}

Si se avanza $u_*$ con C.17, la operación local es

$$
u_{\rm qp}=u_*-K_0|\nabla|\Delta||^2-U_{\rm vac}(|\Delta|,q),
\qquad
T_E=\big[u_{\rm qp}^{\rm FD}(\,\cdot\,;|\Delta|,q)\big]^{-1}(u_{\rm qp}).
\tag{C.18a}
$$

La inversa mantiene **los valores instantáneos** de $|\Delta|$ y $q$. No se sustituye antes $|\Delta|$ por su valor de equilibrio a cada temperatura: eso eliminaría la amplitud como variable dinámica. Para un espectro fijo, $\partial_Tu_{\rm qp}^{\rm FD}>0$; por tanto, $u_{\rm qp}$ y $T_E$ son dos coordenadas de la misma cantidad de energía. A baja temperatura esa derivada puede ser muy pequeña, y avanzar energía puede resultar más conveniente numéricamente que avanzar $T_E$.

La analogía es la escala de un recipiente: una misma cantidad de líquido puede indicarse por volumen o por altura **si se conoce la forma del recipiente**. Aquí la forma corresponde al espectro determinado por $|\Delta|,q$. Si cambia, también cambia la traducción entre energía y temperatura. Y conocer cuánto líquido hay no informa cómo se reparte en compartimentos: del mismo modo, $u_{\rm qp}$ o $T_E$ no reconstruyen una distribución no térmica.

Para calcular la figura C.1 se usa únicamente el caso BCS $q=0$, con $E(x)=\sqrt{x^2+|\Delta|^2}$ y $dx=\rho(E)dE$:

$$
u_{\rm qp}^{\rm FD}(T_E,|\Delta|)=4N_0\int_0^\infty
\frac{\sqrt{x^2+|\Delta|^2}}{\exp[\sqrt{x^2+|\Delta|^2}/(k_BT_E)]+1}\,dx.
\tag{C.18b}
$$

La escala $\Delta_{\rm ref}$ sólo fija las unidades del gráfico. Cada curva se obtiene invirtiendo C.18b a distintas amplitudes para una energía dada; no se impone la ecuación de gap.

![Figura C.1. Curvas de igual energía de excitaciones BCS a q=0. Cada punto usa el espectro de su amplitud instantánea; las amplitudes no se impusieron autoconsistentes. Los dos puntos destacados tienen la misma energía y diferentes temperaturas equivalentes. La curva es un mapa de coordenadas, no una trayectoria del detector.](figuras/C_01_coordenadas_energia.png){width=88%}

En la etapa espectral se avanza $p(x)$, o sus momentos validados, y se calcula su energía; si además se avanza C.17 para controlar conservación, ambas representaciones deben ser compatibles. $T_E$ puede seguir usándose para presentar resultados o parametrizar la movilidad, conservando $X_{|\Delta|}[p]$, $\mathbf j_s[p]$ y las colisiones de la distribución retenida. Una energía, un número de cuasipartículas y la fuerza sobre el condensado son momentos distintos: el número por sí solo tampoco es una coordenada energética biyectiva para distribuciones arbitrarias.

## C.4.4. Derivación de la forma térmica, paso a paso

Cuando los electrones sí son térmicos, se escribe $f=f_e^{\rm FD}(T_e,|\Delta|,\mathbf q)$, con $s_e=-f_T$ y $C_e=-T_ef_{TT}$. La disipación satisface

$$
T_e\dot s_e=\nabla\cdot(\kappa_s\nabla T_e)
+\sigma_nE^2+Q_\Delta-P_{e\text{-ph}}+S_{\gamma,e}^{(E)}.
\tag{C.19}
$$

La fuente se anula después de la deposición cuando el fotón se representa sólo mediante una condición inicial. La regla de la cadena, sin mover todavía el espectro a su equilibrio, da

$$
T_e\dot s_e=C_e\dot T_e
-T_ef_{T|\Delta|}\partial_t|\Delta|
-T_e\partial_T\partial_{\mathbf q}f\cdot\dot{\mathbf q}.
\tag{C.19a}
$$

Despejando se recupera la forma de la revisión anterior, con la fuente ahora explícita:

$$
\boxed{\begin{aligned}
C_e\dot T_e={}&\nabla\cdot(\kappa_s\nabla T_e)
+\sigma_nE^2+Q_\Delta-P_{e\text{-ph}}+S_{\gamma,e}^{(E)}\\
&+T_ef_{T|\Delta|}\partial_t|\Delta|
+T_e\partial_T\partial_{\mathbf q}f\cdot\dot{\mathbf q}.
\end{aligned}}\tag{C.20}
$$

La equivalencia con avanzar energía se comprueba sin invocar una analogía. Como $u_e=f-T_ef_T$, sus derivadas a variables instantáneas fijas son

$$
\begin{aligned}
\dot u_e={}&C_e\dot T_e
+(f_{|\Delta|}-T_ef_{T|\Delta|})\partial_t|\Delta|\\
&+(\partial_{\mathbf q}f-T_e\partial_T\partial_{\mathbf q}f)
\cdot\dot{\mathbf q}.
\end{aligned}\tag{C.20a}
$$

Al insertar C.20 en C.20a, los dos términos cruzados se cancelan. Quedan la potencia que entra en ocupaciones y $f_{|\Delta|}\partial_t|\Delta|+\partial_{\mathbf q}f\cdot\dot{\mathbf q}$. Estos últimos son precisamente los términos de campos de C.13, pues A establece $f_{|\Delta|}=X_{|\Delta|}$ y $\partial_{\mathbf q}f=\boldsymbol\Pi$ en la distribución térmica. Añadir la energía de gradiente recupera C.17. Los términos cruzados traducen el cambio del espectro a un cambio de temperatura: no representan un segundo calentador.

**Regla de contabilidad energética:** almacenamiento reversible, trabajo $\mathbf j_s\cdot\mathbf E$, disipación $Q_\Delta$ y flujos de borde se cuentan una sola vez. Añadir solamente dos derivadas de almacenamiento a una fila térmica antigua no basta para demostrar la conservación del sistema acoplado.

## C.4.5. Una prueba dinámica pequeña de esa contabilidad

Para aislar la idea se usa un **modelo termodinámico ilustrativo adimensional**, no el funcional BCS ni un ajuste de NbN. Sólo en esta prueba $d=|\Delta|/\Delta_{\rm ref}$ y $\vartheta$ es una temperatura sin unidades. Se define

$$
f_{\rm il}(d,\vartheta)=\frac{(1-d^2)^2}{4}
-\frac{1-0.4d^2}{2}\vartheta^2,
\qquad
u_{\rm il}=\frac{(1-d^2)^2}{4}
+\frac{1-0.4d^2}{2}\vartheta^2.
\tag{C.20b}
$$

Se integra $\dot d=-\partial_df_{\rm il}$ con movilidad unidad, partiendo de $d=0.35$, $\vartheta=0.2$, e invirtiendo $u_{\rm il}(d,\vartheta)=u_{\rm il}(0)$ en cada evaluación. El tiempo del ejemplo carece de calibración física. En el intervalo utilizado $C_{\rm il}=(1-0.4d^2)\vartheta>0$, y se verifica $\dot u_{\rm il}=0$ y $\vartheta\dot s_{\rm il}=(\partial_df_{\rm il})^2\geq0$.

![Figura C.2. Celda ilustrativa aislada. Izquierda: recuperación de amplitud y aumento de temperatura. Centro: la energía de orden disminuye y la de excitaciones aumenta, manteniendo constante su suma. Derecha: evolución de energía si se conserva el balance y si se congela artificialmente la temperatura durante la relajación. El segundo cálculo pierde energía hacia un reservorio implícito; no describe una celda aislada.](figuras/C_02_celda_balance.png){width=100%}

El estado al tiempo ilustrativo 12 tiene $d\simeq0.87784$ y $\vartheta\simeq0.75730$. Una segunda integración avanza la ecuación térmica de la misma energía, independientemente de la inversión: las dos representaciones difieren menos de $1.2\times10^{-12}$, y la deriva absoluta de energía de esa segunda integración es menor que $6.8\times10^{-13}$. El paso temporal es 0.004; reducirlo a la mitad cambia el estado menos de $1.6\times10^{-12}$ [Q].

Esta prueba verifica la identidad algebraica y su realización numérica en un problema controlado. No predice una temperatura máxima, una escala de picosegundos ni una latencia de detección.

# C.5. Tiempos de relajación: significado y calibración admisible

Allmaras utiliza como referencia una combinación

$$
\begin{gathered}
\frac1{\tau_{sc}(T)}=\frac1{\tau_{ee}(T)}+\frac1{\tau_{ep}(T)},\\
\tau_{ee}(T)=\tau_{ee,c}\frac{T_c}{T},\\
\tau_{ep}(T)=\tau_{ep,c}\left(\frac{T_c}{T}\right)^3.
\end{gathered} \tag{C.21}
$$

En la memoria se usaron $\tau_{ee,c}=0.50$ ps y $\tau_{ep,c}=2.47$ ps; la referencia de Allmaras de 2019 emplea escalas de 5 y 24.7 ps. Recuperar esas escalas como punto de partida físico es razonable, pero no equivale a demostrar que cualquier película de NbN tiene exactamente esos valores [M, A].

La tasa de scattering de una cuasipartícula, el tiempo de relajación de la energía total y la movilidad de amplitud no son el mismo observable. Pueden provenir de distintos momentos o modos del mismo operador de colisión. Por eso la potencia calculada desde $\alpha^2F$ no debe sustituirse por $C_e(T_e-T_{\rm ph})/\tau_{sc}$, y tampoco se exige que cada tiempo efectivo sea numéricamente idéntico.

## C.5.1. No contar dos veces los modos cinéticos

KWT incorpora el efecto de variables eliminadas. Si se evolucionan explícitamente esas mismas variables de población y se conserva íntegramente su corrección KWT, puede duplicarse una memoria física. Una manera de organizar la reducción es separar variables retenidas $P$ y eliminadas $Q$. Si el operador linealizado tiene bloques,

$$
\partial_t\begin{pmatrix}y_P\\y_Q\end{pmatrix}
=\begin{pmatrix}L_{PP}&L_{PQ}\\L_{QP}&L_{QQ}\end{pmatrix}
\begin{pmatrix}y_P\\y_Q\end{pmatrix},
$$

eliminar $y_Q$ a frecuencia $\omega$ produce

$$
L_{\rm eff}(\omega)=L_{PP}
+L_{PQ}(i\omega-L_{QQ})^{-1}L_{QP}. \tag{C.22}
$$

Su expansión de baja frecuencia genera restauración, movilidad y memoria efectivas. Los modos ya incluidos en $P$ no se vuelven a integrar en ese término. C.22 es la regla matemática de reducción; obtener sus bloques a partir del operador material completo es una tarea posterior, no un cálculo ya realizado en este paquete.

Para la primera versión se propone identificar $\tau_\psi$ como movilidad efectiva de los grados no resueltos y $\tau_{ee}$ como tiempo del operador explícito de poblaciones. Si se usa literalmente C.21 para $\tau_\psi=\tau_{sc}$, esos tiempos quedan vinculados: no son tres ajustes independientes. Una modificación independiente de $\tau_\psi$ debe declararse como calibración efectiva del cierre dinámico y congelarse después en las comparaciones fuera de calibración.

## C.5.2. Lo que sí indica la comprobación de escalas

La revisión 0.1 registró, a $T_b=0.9$ K, que multiplicar ambos tiempos de referencia por diez llevaba $\tau_{sc}$ de aproximadamente 4.795 ps a 47.95 ps y el factor KWT, para el gap de referencia usado en la prueba, de 19.18 a 191.57. Esas cifras proceden del paquete anterior y no se revalidaron en la revisión 0.2; aquí se conserva sólo su función orientativa. El cambio de escala modifica la movilidad de amplitud y fase. **No prueba que la latencia ni el tiempo al pico cambien por un factor diez**, porque la energía, la corriente y el espectro evolucionan simultáneamente.

# C.6. Inicialización y fronteras superconductoras

## C.6.1. Estado inicial

Se obtiene una solución uniforme o débilmente inhomogénea con el mismo funcional y corriente que se usarán en el transitorio:

$$
\begin{gathered}
X_{|\Delta|}^{\rm FD}(T_b,|\Delta|_b,q_b)=0,\\
j_s(T_b,|\Delta|_b,q_b)=I_b/(wd),\\
T_e=T_{\rm ph}=T_b.
\end{gathered} \tag{C.23}
$$

En una región uniforme sin disipación se puede iniciar $\mathbf j_n=0$. No se requiere crear zonas de conversión normal–superconductora para recuperar un interior bulk. La estabilidad de C.23 se comprueba con perturbaciones y con la rama uniforme; el mismo voltaje de fondo no debe heredarse de una trayectoria parcialmente estacionaria y luego interpretarse como respuesta del fotón.

## C.6.2. Qué sustituye a las fronteras metálicas

En los extremos artificiales se conecta la región 2D con una continuación superconductora que porta la corriente del circuito. La amplitud y el gradiente de fase deben ser compatibles con el estado exterior; no se impone $\Delta=0$ ni que toda la corriente entre como corriente normal. En los bordes laterales físicos se conserva el flujo normal de corriente nulo y el aislamiento térmico compatible con el dispositivo.

Si la corriente cambia rápidamente, la parte superconductora exterior puede generar una tensión inductiva aunque permanezca fría. Eliminar los contactos no significa eliminar ese trabajo: se representa en la continuación 1D o en una inductancia exterior, sin doble contabilización.

## C.6.3. Menos área 2D no significa un reservorio frío arbitrariamente próximo

La escala normal de difusión da un diagnóstico sencillo. Para un núcleo gaussiano, exigir una amplitud relativa $\varepsilon$ a distancia $\ell$ conduce a

$$
\ell\simeq\sqrt{4Dt\ln(1/\varepsilon)}. \tag{C.24}
$$

Es una estimación de alcance espacial de un problema normal idealizado, no una cota rigurosa del error de toda la simulación superconductora. Con $D=1.581$ cm$^2$/s y $\varepsilon=1\%$, se obtiene aproximadamente 121 nm a 5 ps, 241 nm a 20 ps y 382 nm a 50 ps. Imponer $T_b$ a 50 nm del impacto podría evacuar energía demasiado pronto, aunque ya no existan contactos metálicos.

Una solución compatible con la reducción de coste buscada es una región 2D próxima a la absorción, más continuaciones longitudinales 1D o condiciones térmicas transparentes. Se conserva la física exterior sin mallar con la misma resolución todo el nanohilo.

Para un exterior normal lineal semiinfinito, inicialmente a $T_b$ y con $\kappa,C$ constantes, la ecuación de calor en Laplace da un flujo saliente

$$
\widehat Q_{\rm out}(s)=\sqrt{\kappa Cs}\,
\widehat{(T_{\partial}-T_b)}(s). \tag{C.25}
$$

Invirtiendo,

$$
Q_{\rm out}(t)=\sqrt{\frac{\kappa C}{\pi}}
\int_0^t\frac{\dot T_{\partial}(s)}{\sqrt{t-s}}\,ds. \tag{C.26}
$$

C.26 supone $T_\partial(0)=T_b$ sin salto inicial. Si se prescribe un salto $\delta T_0$, se añade $\sqrt{\kappa C/\pi}\,\delta T_0/\sqrt t$. Es una referencia transparente exacta para ese problema lineal, no una condición exacta para un exterior superconductivo no lineal con escape y recombinación. En ese caso una continuación 1D física suele ser una opción más directa. El ahorro debe evaluarse con la longitud 2D mínima que deje invariantes las métricas, no con una reducción de longitud elegida exclusivamente por tiempo de cómputo.

# C.7. Amplitudes pequeñas, grandes valores de q y núcleos

Cuando $|\Delta|\to0$, la fase pierde significado local. Una tabla limitada en $q$ no puede resolver el problema mediante recorte del último valor de corriente. Tampoco es correcto imponer instantáneamente $|\Delta|=0$ cuando una estimación local excede la corriente crítica: eso elimina precisamente la dinámica de nucleación o colapso que se desea estudiar.

La rama espectral debe admitir estados fuertemente depareados. A $\Gamma$ grande y amplitud pequeña, la corriente uniforme decae en vez de quedar fijada al valor del extremo del catálogo. Esa asintótica ayuda a extender el cierre físico, pero no convierte un núcleo inhomogéneo en un material uniforme.

La representación espacial debe operar con $\Delta$ complejo y variables covariantes. Si se introduce una regularización alrededor de $|\Delta|=0$, debe aplicarse al **funcional completo** y derivar de él tanto la fuerza como la corriente regularizadas. Regularizar únicamente $q$ o únicamente el cociente de fase rompe A.10 y el balance de C.18.

No se declara resuelto el problema de los núcleos mediante una fórmula universal nueva en esta propuesta. Se requiere una prueba de independencia respecto del umbral de regularización y una comparación espacial de referencia, como la de C.8. Si esas pruebas no convergen, la aproximación espectral local debe sustituirse en esa región; aumentar el rango de la tabla no basta.

# C.8. Prueba de la aproximación de gradientes

El funcional espacial quasiclásico también contiene gradientes de los ángulos espectrales. Para $q=0$ y una amplitud lentamente variable,

$$
\begin{gathered}
\Theta_n=\arctan(|\Delta|/\epsilon_n),\\
\partial_{|\Delta|}\Theta_n=\frac{\epsilon_n}{\epsilon_n^2+|\Delta|^2}.
\end{gathered}
$$

El coeficiente de $|\nabla |\Delta||^2$ obtenido del término espectral de gradiente es

$$
K_{|\Delta|}(T,|\Delta|)=\pi N_0\hbar D k_BT
\sum_{n\geq0}\frac{\epsilon_n^2}{(\epsilon_n^2+|\Delta|^2)^2}. \tag{C.27}
$$

Cerca de $T_c$, con $|\Delta|\to0$, la suma recupera $K_0$. A $T\to0$ se usa $k_BT\sum_n\to(2\pi)^{-1}\int_0^\infty d\epsilon$ y

$$
\begin{gathered}
\int_0^\infty\frac{\epsilon^2}{(\epsilon^2+|\Delta|^2)^2}d\epsilon
=\frac\pi{4|\Delta|},\\
K_{|\Delta|}(0,|\Delta|)=\frac{\pi N_0\hbar D}{8|\Delta|}.
\end{gathered} \tag{C.28}
$$

Por tanto, para $|\Delta|\simeq1.764k_BT_c$, $K_{|\Delta|}/K_0\simeq0.567$. La longitud proporcional a $\sqrt K$ sería aproximadamente 0.753 de la que usa $K_0$, o, a la inversa, el cierre con $K_0$ daría una longitud aproximadamente un 33% mayor. No es un error demostrado del tiempo de detección; muestra que la rigidez espacial merece una prueba antes de reclamar precisión del 10–20%.

Tampoco conviene sustituir sin más $K_0$ por $K_{|\Delta|}(T,|\Delta|,q)$ en la ecuación de amplitud. Aparecen nuevas derivadas variacionales y contribuciones de entropía y energía de gradiente. Para un coeficiente térmico local, la energía de gradiente contiene $(K_{|\Delta|}-T\partial_TK_{|\Delta|})|\nabla |\Delta||^2$. Ignorar esa modificación volvería a romper el balance que se acaba de construir.

La recomendación cercana es mantener C.1 como cierre explícito, confrontarlo con un problema estacionario espacial pequeño de Usadel y promover el término de gradiente solamente si su error en la barrera, el perfil o la escala de nucleación supera el presupuesto físico. No se necesita un Keldysh espacial completo de todo el evento para realizar esa prueba preliminar.

# C.9. Desde la absorción al pico antes de amplificación

## C.9.1. Dos voltajes diferentes

Se define un par de planos de medida en las continuaciones superconductoras y

$$
V_{\rm patch}(t)=\int_L^R\mathbf E\cdot d\mathbf l.
\tag{C.29}
$$

Con $\mathbf A=0$ y la orientación indicada, $V_{\rm patch}=\phi_L-\phi_R$. La parte exterior no resuelta puede contribuir

$$
\begin{gathered}
V_{\rm port}=V_{\rm patch}+L_{\rm ext}^{\rm diff}(I)\dot I,\\
L_{\rm ext}^{\rm diff}=\frac{d\Phi_{\rm ext}}{dI}.
\end{gathered} \tag{C.30}
$$

La inductancia exterior excluye la energía inductiva ya descrita por la región resuelta. La energía exterior cumple $dU_{\rm ext}/dI=I L_{\rm ext}^{\rm diff}(I)$; no se usa $U=LI^2/2$ cuando $L$ depende de corriente sin comprobar cuál inductancia se está definiendo.

Para una fuente de corriente ideal y una resistencia física $R_L$ en paralelo,

$$
\begin{gathered}
L_{\rm ext}^{\rm diff}(I)\dot I+V_{\rm patch}
=R_L(I_b-I),\\
V_{\rm port}=R_L(I_b-I).
\end{gathered} \tag{C.31}
$$

Este circuito no representa un amplificador ni un filtro. Es la carga que modifica la corriente y el crecimiento del estado resistivo. Aun antes de amplificar, la dinámica del detector depende de ella. Si se utiliza otra topología, se conserva su trabajo y sus energías almacenadas y se declara el plano de salida.

## C.9.2. Qué tiempo representa cada fenómeno

Deben registrarse por separado el primer phase slip o cruce de vórtice identificado, la formación de una región resistiva que atraviese el ancho cuando exista, el desarrollo del pulso de puerto y el tiempo de su pico. La formación de un hotbelt es una métrica física especialmente relevante para el régimen físico considerado; no se impone como condición universal para toda señal resistiva posible.

Los tiempos de cruce de un umbral de tensión son observables instrumentales o convencionales distintos. En la memoria el mismo umbral numérico de 100 µV se aplicó a una señal interna y a otra atenuada y transducida por el circuito. Eso **no permite deducir** que el estado resistivo intrínseco se formó en 0.06–0.10 ps ni que la latencia restante fuese únicamente un retardo externo. Esa interpretación temporal no se adopta en este documento.

Para el alcance solicitado se entregarán $V_{\rm patch}(t)$, $V_{\rm port}(t)$, $I(t)$, evolución energética, topología y extensión del estado resistivo, $t_{\rm HB}$ cuando sea identificable, amplitud del pico y tiempo hasta ese pico. Un pico asociado a oscilaciones rápidas de phase slip se distingue del máximo del pulso macroscópico declarado; el criterio se fija antes de comparar modelos o datos.

No se infiere eficiencia, probabilidad de absorción o jitter estadístico a partir de estas trayectorias deterministas. Tampoco se elimina la comparación experimental: trazas preamplificador referidas al mismo puerto, latencias relativas bien definidas y escalas de crecimiento ofrecen objetivos verificables.


## C.9.3. Resultado ligero: el mismo evento prescrito produce puertos distintos

Para visualizar C.31 se prescribe únicamente $V_{\rm patch}=R_{\rm patch}(t)I$. Se utiliza la resistencia positiva ilustrativa

$$
\begin{gathered}
\zeta=tR_L/L_{\rm ref},\quad i=I/I_b,\quad \lambda=L_{\rm ext}/L_{\rm ref},\\
r(\zeta)=R_{\rm patch}/R_L
=4(1-e^{-\zeta/0.08})e^{-\zeta/0.5},\\
\lambda\,\frac{di}{d\zeta}=1-(1+r)i,\qquad
v_{\rm patch}=ri,\quad v_{\rm port}=1-i.
\end{gathered}\tag{C.31a}
$$

Ambos voltajes se normalizan por $I_bR_L$. La condición inicial es $i(0)=1$ y se comparan $\lambda=0.25,1,4$. La función $r$ es una entrada elegida para esta prueba; no se calcula a partir de un hotbelt ni se afirma que toda región resuelta se reduzca a una resistencia instantánea.

![Figura C.3. Respuesta calculada del circuito a una resistencia prescrita. Arriba se mantiene idéntica la entrada. Abajo, la inductancia cambia la descarga de corriente y el voltaje de puerto. Las curvas discontinuas muestran el voltaje de la región para lambda=1 y permiten compararlo con su puerto. Las marcas indican máximos del puerto, sin asignarles una latencia microscópica.](figuras/C_03_circuito_prescrito.png){width=95%}

Al inicio $i=1$, de modo que el puerto aún tiene tensión nula aunque comience a crecer la caída en la región: la tensión inductiva compensa esa caída. Durante la descarga, la energía almacenada en la inductancia disminuye. Más adelante la corriente se desvía hacia $R_L$, y ese desvío genera $V_{\rm port}$. Los máximos de ambas tensiones no tienen por qué coincidir. En el detector acoplado, además, la corriente vuelve a modificar la resistencia y el condensado; ese retorno de información no se simuló aquí.

| $L_{\rm ext}/L_{\rm ref}$ | Pico del puerto: tiempo | Pico del puerto: voltaje | Pico de la región: tiempo |
|---:|---:|---:|---:|
| 0.25 | 0.308 | 0.6791 | 0.072 |
| 1 | 0.644 | 0.5242 | 0.110 |
| 4 | 1.198 | 0.2673 | 0.140 |

Son valores adimensionales de C.31a sobre una malla temporal de paso 0.002; los tiempos de máximo tienen esa resolución. La tabla permite comparar efectos del circuito sin atribuirlos a la cinética del detector.

La comprobación independiente del circuito es su balance de potencia,

$$
i\lambda\frac{di}{d\zeta}+ri^2+(1-i)^2=1-i.
\tag{C.31b}
$$

A la izquierda aparecen la variación de energía inductiva, la potencia en la región y la potencia en la carga; a la derecha, la potencia de la fuente ideal. La prueba integra también este balance hasta $\zeta=6$: el residuo absoluto máximo es $1.7\times10^{-5}$, dominado por la cuadratura trapezoidal de potencia. Al repetir la evolución con medio paso, las corrientes cambian menos de $6.5\times10^{-10}$. Confirma la lectura de los dos voltajes dentro del circuito prescrito; la predicción del pulso completo requiere C.1–C.18 y B.


# C.10. Qué demuestra este documento

C.6–C.20 demuestran que la movilidad positiva elegida puede acoplarse a la termodinámica uniforme y a una energía térmica invertible sin perder ni duplicar trabajo eléctrico. La derivación común de amplitud y corriente hace explícita esa contabilidad. Su utilidad es permitir nuevas comparaciones de espectros, energía y señales; no se exige convertir todo el detector en una formulación microscópica exacta.

No demuestra por sí sola que KWT describa exactamente toda la memoria no térmica, que el gradiente constante sea suficiente dentro de un núcleo, ni que una trayectoria convergida coincida con un dispositivo real al 10–20%. Esas preguntas tienen pruebas previas concretas: reducción cinética, referencia espacial y caracterización material. D las organiza sin exigir resolver fenómenos fuera del alcance de la primera publicación.

# C.11. Estimación de coste para el dominio físico reducido

Sea $N$ el número de nodos 2D activos, $N_t=t_{\rm fin}/\delta t$ el número de pasos y $c_{\rm loc}$ el coste constitutivo por nodo. Una descomposición orientativa es

$$
W\simeq W_{\rm PRE}+N_t\,[c_{\rm loc}N+C_{\rm campo}(N)]
+W_{\rm cin,temprana}+W_{\rm exterior}. \tag{C.32}
$$

$C_{\rm campo}$ incluye las restricciones de corriente y cualquier solución global necesaria. Su escalado depende del algoritmo, la geometría y la reutilización de operadores; no se identifica universalmente con $N$ ni con $N^3$.

Para explorar escenarios, si el coste del bloque dominante se aproxima localmente por $N^p$, se obtiene

$$
\frac{W_{\rm nuevo}^{\rm bloque}}{W_{\rm antiguo}^{\rm bloque}}
\simeq\frac{t_{\rm fin,n}}{t_{\rm fin,previo}}
\frac{\delta t_{\rm previo}}{\delta t_n}
\left(\frac{N_n}{N_{\rm previo}}\right)^p r_{\rm constitutivo}. \tag{C.33}
$$

No es una predicción de wall-time hasta medir sus factores. Una comparación favorable por detenerse al pico, y no al reset, no debe presentarse como aceleración del mismo problema físico. Es un alcance menor que responde mejor a la publicación solicitada. La conservación de las continuaciones térmicas y eléctricas exteriores sigue siendo obligatoria.

En una evaluación directa de las colisiones espectrales, un orden de coste adicional es

$$
\begin{gathered}
W_{\rm cin}\sim c_{\rm col}N_{\rm cel,cin}N_{t,\rm cin}N_EN_\Omega,\\
M_{\rm ocup}\sim N_{\rm cel,cin}(N_E+N_\Omega).
\end{gathered} \tag{C.34}
$$

Aquí se cuentan solamente ocupaciones y evaluaciones de colisiones, no toda la memoria del programa. La estructura separable de los núcleos puede reducir el coste mediante reorganización de integrales; no se presupone una aceleración concreta. La duración, el número de celdas y la resolución espectral del intervalo no térmico deben medirse: ese término puede dominar aunque el dominio 2D sea menor.

# Referencias y procedencia

[M] J. A. Díaz Monge, `memoria_02.pdf`, anexo C, pp. impresas 120–125, y resultados de latencia, pp. 95–97. Fuente señalada por la revisión 0.1 para la ecuación antigua, fronteras y voltajes comparados. Los valores históricos de latencia y tiempos citados aquí no se reprodujeron en esta revisión.

[V] D. Y. Vodolazov, Physical Review Applied **7**, 034014 (2017), [DOI: 10.1103/PhysRevApplied.7.034014](https://doi.org/10.1103/PhysRevApplied.7.034014); [preprint](https://arxiv.org/abs/1611.06060). Modificación del modelo y prueba de corriente de depairing a diferentes temperaturas.

[A] J. P. Allmaras, `Allmaras_thesis.pdf`, pp. impresas 91–114. Ecuación 3.24, movilidad generalizada, balance de energía y comparación con modelos cinéticos. J. P. Allmaras et al., Physical Review Applied **11**, 034062 (2019), [DOI: 10.1103/PhysRevApplied.11.034062](https://doi.org/10.1103/PhysRevApplied.11.034062).

[F] P. Virtanen, A. Vargunin y M. Silaev, Physical Review B **101**, 094507 (2020), [DOI: 10.1103/PhysRevB.101.094507](https://doi.org/10.1103/PhysRevB.101.094507); [preprint](https://arxiv.org/abs/1909.00992). Marco variacional quasiclásico.

[R] Inspección histórica de 0.1: `JoaquinDiazM/pysnspd`, revisión `f3c26b95ff4e4a93504371e78b46ad3a20e06273`; sectores de material, Allmaras, solver, circuito, térmica y tiempos. Los ejemplos nuevos se ejecutan mediante [Q] sin modificar esos módulos. D registra la revisión local de contexto.

[Q] Documentos A y B de la revisión 0.2. Las figuras C.1–C.3 se generan con `sandbox/model_v0_2/checks_c.py`; las entradas, salidas y residuos se conservan en `verificaciones/C_*.csv` y `verificaciones/C_checks.json`. C.1 usa integrales BCS a amplitud fija; C.2 utiliza exclusivamente C.20b; C.3 utiliza exclusivamente C.31a. No son transitorios del solver de producción. C.6, la aplicación de C.22 a variables retenidas y los criterios de acoplamiento espacial constituyen propuestas de reestructuración, no resultados ya validados contra transitorios experimentales.
