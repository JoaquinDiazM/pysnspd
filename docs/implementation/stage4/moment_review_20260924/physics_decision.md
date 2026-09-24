# Decisión física tras resolver los momentos

La campaña de 181 espectros, 181 respuestas cinéticas y siete proyecciones
terminó correctamente. Su resultado resuelve la pregunta que motivó este
refinamiento: **la diferencia de corriente y fase de la proyección con un
solo potencial es real en el ensayo congelado; no la explica la cuadratura
energética gruesa**. No hace falta repetir ahora otra campaña estática para
obtener esa conclusión. La etapa 4 permanece abierta porque todavía falta
el acoplamiento temporal del nuevo núcleo con poblaciones, trabajo
espectral y dispositivo.

## Qué decisión permiten los datos

La comparación de referencia usa el núcleo radial de 65×65, la malla fina
de 50 energías y η/Δref=0,01. La perturbación angular rompe la simetría de
la distribución alrededor del núcleo. Las magnitudes son normas de campos
integrados en energía; no son errores de una latencia del detector.

| Momento | Diferencia proyección/espectral, relativa al resultado espectral | Variación numérica observada, en la misma escala | Interpretación |
|---|---:|---:|---|
| Corriente | 20,55 % | 2,87 % | Diferencia resuelta; mejora respecto del 50,27 % al omitir hT |
| Torque de fase | 77,53 % | 3,64 % | Diferencia resuelta; la norma proyectada alcanza 1,717 veces la espectral |
| Fuerza compleja | 3,23 % | 1,53 % | Su norma conjunta oculta parte del error de fase |
| Fuerza de amplitud | 0,00154 % | 1,52 % | Diferencia no resuelta frente a la variación observada |
| Flujo energético, ∫E jL dE | 0,438 % | 0,196 % | Diferencia pequeña pero resuelta; mejora respecto de omitir hT incierta |

Los números y controles completos están en [analysis.json](analysis.json).
La variación es una suma de cambios observados de malla, energía y η; no
una cota matemática del error. El control radial casi nulo continúa sin
resolver esos pequeños efectos y no debe rechazarse mediante cocientes de
números próximos a cero. El defecto de la integral de χ bajó del 18,02 %
al 0,4243 %, mientras la diferencia de corriente permaneció cerca del 20 %.

La decisión es conservar la respuesta espectral completa como referencia
para la siguiente prueba y **no promover todavía la proyección de un solo
potencial como respuesta cuantitativa de corriente/fase**. Tampoco se
descarta todo el cierre mesoscópico de la memoria: aquí se mantuvieron
congelados el gap, su fase y el potencial. La fase KWT y el potencial normal
pueden responder en el sistema completo. La superioridad de una solución
estática espectral no demuestra que haga falta añadir una nueva variable
dinámica al detector.

## Siguiente cálculo: medir la memoria del modo de carga

El paso más directo aprovecha los espectros y operadores existentes, sin
resolver otra vez la autoconsistencia ni refinar otra vez la DOS. Se compara
la respuesta del modo de carga a una perturbación longitudinal oscilante
con su eliminación algebraica. Así se distingue el error de la forma
energética χ(E)v del error adicional de suponer una respuesta instantánea.

Escribiendo el operador real ya implementado por bloques,

\[
\begin{pmatrix}r_L\\r_T\end{pmatrix}
=\begin{pmatrix}L_{LL}&L_{LT}\\L_{TL}&L_{TT}\end{pmatrix}
\begin{pmatrix}h_L\\h_T\end{pmatrix},
\]

la eliminación estática completa es

\[
h_T^{(0)}=-L_{TT}^{-1}L_{TL}h_L,\qquad
L_{\rm red}=L_{LL}-L_{LT}L_{TT}^{-1}L_{TL}.
\]

Esta reducción de Schur reproduce el problema **estático** sin introducir
estados dinámicos adicionales. Su coste pequeño no acredita una separación
de tiempos. El siguiente diagnóstico debe comprobar esa separación antes
de usarla en un transiente.

### Almacenamiento y escala de tiempo del diagnóstico

En un cierre cinético adiabático con espectros fijos, se retiene el
almacenamiento diagonal proporcional a la DOS,

\[
M(E)=\operatorname{diag}\{m_i\operatorname{Re}g_i(E)\},\qquad
t_D=\frac{\ell_0^2}{D}=\frac{\hbar}{2k_BT_c}\simeq0,4415\;\mathrm{ps},
\]

y se estudia el modelo lineal

\[
t_D M\,\partial_t h_T=L_{TL}h_L+L_{TT}h_T.
\tag{P1}
\]

Aquí t y tD deben emplear la misma unidad física. Con τ=t/tD se escribe
simplemente M∂τhT=rT. El control normal fija el prefactor sin ajuste:
g=1, M=m y el operador de caras recupera ∂t h=D∇²h.

El origen del término de almacenamiento se ve en la convención matricial
ya usada: K=Rh−hA, h=hL I+hT τ3. Sus elementos diagonales dan

\[
\frac14\operatorname{Tr}\{\tau_3,K\}=2\operatorname{Re}g\,h_L,
\qquad
\frac14\operatorname{Tr}\tau_3\{\tau_3,K\}
=2\operatorname{Re}g\,h_T.
\]

El término temporal de Usadel aporta ℏ/(4kBTc) a esta anticomutación;
con la masa geométrica resulta el tiempo tD anterior. **P1 es el diagnóstico
con almacenamiento DOS y coeficientes congelados de una expansión
adiabática**, no una afirmación de que esos términos agoten la respuesta
microscópica AC a cualquier frecuencia. Los productos temporales completos
son convoluciones; las correcciones espectrales de orden frecuencia y el
movimiento del gap requieren tratamiento adicional. La formulación de dos
tiempos y su expansión lenta están expuestas en Larkin y Ovchinnikov,
ecuaciones 5 y 10–11 de
[Nonlinear conductivity of superconductors in the mixed state](https://www.jetp.ras.ru/cgi-bin/dn/e_041_05_0960.pdf).

Para la convención exp(−iωt), el cálculo concreto es

\[
[-L_{TT}-i\omega t_D M]h_T(\omega)=L_{TL}h_L.
\tag{P2}
\]

Se usan los mismos contactos, el mismo perfil infinitesimal y las mismas
proyecciones de corriente, fuerza cartesiana y flujo de energía. ω=0 debe
recuperar la respuesta ya guardada. No se divide por una DOS mínima ni se
añade el escape artificial η como relajación. Si una masa tiende a cero,
ese límite debe conservar su carácter algebraico; no se fabrica capacidad
de almacenamiento en un intervalo sin estados.

Para reconstruir observables con fasores complejos se evalúan por separado
las partes real e imaginaria de h. Las dos componentes cartesianas de la
fuerza del gap son dos fasores: no se puede comprimir esa información
mezclando la unidad imaginaria del parámetro de orden con la del fasor.

### Banda y criterio de interpretación

Una selección económica es ν=ωtD en {0,10⁻⁴,10⁻³,10⁻²,10⁻¹,1}.
Los primeros puntos anclan el límite lento; los mayores exploran dónde el
cierre local en tiempo deja de parecer instantáneo. Para cada punto se
declaran también

\[
\frac{\hbar\omega}{\Delta_{\rm ref}}=\frac{2\nu}{1,76387694},\quad
\frac{\hbar\omega}{k_BT_b}=\frac{2\nu}{0,10404624},\quad
\frac{\hbar\omega}{\eta}=\frac{2\nu}{1,76387694\eta_{\rm rel}}.
\]

La última razón describe la resolución del contorno numérico, no un tiempo
de colisión físico. Los rasgos espectrales estrechos impiden deducir una
banda adiabática uniforme sólo de que ℏω sea pequeño frente al gap global.
El resultado a frecuencias mayores se etiqueta como estrés del cierre P1;
no acredita por sí mismo la respuesta subpicosegundo del dispositivo.

Se comparan la amplitud y el desfase de los momentos respecto del límite
ω=0, junto con la variación disponible en η. Si la memoria ya modifica los
momentos dentro de la banda lenta, la eliminación instantánea necesita una
corrección o un estado reducido. Si sólo aparecen diferencias fuera de
esa banda, el resultado no decide todavía la dinámica ultrarrápida. No se
añade hT dinámico por defecto ni se convierte una respuesta casi estática
a frecuencias bajas en validación del transiente de Korzh.

## El balance que falta al liberar el condensado

Con gap fijo, el modelo cinético tiene un balance útil y verificable:

\[
\delta U_e/E_*=-2\sum_i m_i\int_0^\infty
\bar E\operatorname{Re}g_i\,\delta h_{L,i}\,d\bar E.
\]

La traza del conmutador de conversión es cero. Por eso su derivada temporal
se contabiliza mediante flujos de caras y reservorios; la conversión de
carga no crea energía longitudinal. El balance algebraico de torque
continúa siendo

\[
B\delta I+\operatorname{Im}(d^*\delta G)
=2\int r_T\,d\bar E.
\]

En P1 su lado derecho es un término de almacenamiento del desequilibrio
de ramas. No equivale automáticamente a una capacitancia ni a la carga
eléctrica neta: en el dispositivo, neutralidad, fase y potencial deben
cerrarse conjuntamente. La distinción entre distribución y potencial
invariante de calibre aparece explícitamente en las ecuaciones 9–10 de
[Golub, Dynamic properties of short superconducting filaments](https://www.jetp.ras.ru/cgi-bin/dn/e_044_01_0178.pdf).
Sus aproximaciones próximas a Tc no se importan como calibración a 0,9 K.

Cuando Δ cambia, no basta multiplicar ∂t hL por la DOS instantánea ni
depositar QΔ: hay trabajo por el movimiento del espectro. Las ecuaciones
9–13 de
[Kozorezov y colaboradores, PRB 92, 064504](https://link.aps.org/accepted/10.1103/PhysRevB.92.064504)
contienen separadamente el almacenamiento proporcional a Re cosθ y el
término que acompaña a ∂tΔ. Ese término, sus flujos y la transferencia de
QΔ deben formar un mismo balance cuando se libere el gap. Una identidad
de torque no demuestra la existencia de una energía escalar para toda
fuerza Keldysh no térmica.

Hay además un control térmico implementable: relajar una perturbación débil
de fase con la acción Matsubara, KWT heredado y potencial normal. Con
v=2eφ/(kBTc), el problema B(c Bᵀv)=−B Ibar y la derivada covariante
∂τd+i(v/2)d dan una pérdida de energía libre normal igual a
(1/2)Σc(Δv)², junto con la pérdida KWT positiva. Ese control prueba la
dinámica térmica de la nueva fuerza espacial, pero no reemplaza el balance
de energía interna con poblaciones no térmicas.

## Control siguiente con gap móvil: límite térmico débil

Es útil avanzar en paralelo hacia un transiente térmico débil de la fuerza
espacial aceptada. No necesita una nueva hipótesis de cascada ni una nueva
movilidad: conserva la acción térmica de Matsubara y KWT heredado. Su
alcance es una respuesta condicionada a T fijo; no se presupone que la
energía interna del sistema no térmico permanezca en equilibrio durante
la relajación. La temperatura fija y el libro de energía libre deben
declararse como condiciones de este control.

No conviene volver a minimizar 256 problemas espectrales en cada paso
temporal diminuto. Alrededor del campo de referencia d0 puede obtenerse
el Hessiano reducido mediante diferenciación implícita de las ecuaciones
espectrales. Con s el conjunto de variables espectrales libres,

\[
H=F_{dd}-F_{ds}F_{ss}^{-1}F_{sd}.
\]

Es el Hessiano de la misma energía después de resolver Usadel, no una
rigidez GL añadida. La forma directamente implementable en las coordenadas
u del módulo térmico es

\[
J_n\delta u_n=m\,\delta d,\qquad
\delta f_n=g_n\delta u_n-g_n^3u_n
\operatorname{Re}(u_n^*\delta u_n),
\]

\[
H\delta d=2m\ln(T/T_c)\delta d+
4\pi(T/T_c)m\sum_n
\left(\frac{\delta d}{\epsilon_n}-\delta f_n\right).
\tag{P3}
\]

El Jacobiano Jn ya está implementado; δu es cero en los contactos
espectrales prescritos. Se factoriza una vez cada Jn y se reutilizan sus
factores para todas las direcciones. Se paralelizan frecuencias
independientes dentro del mismo presupuesto. Si no se conservaron todos
los espectros de referencia, basta una preparación de esos espectros,
no una relajación nueva del gap.

Para el acoplamiento normal, sea L=BcBᵀ el Laplaciano de potencial con
sus valores fijos de contorno, y S la aplicación v↦i d0 v escrita en
coordenadas reales. Sea M la matriz temporal KWT positiva, incluidos los
pesos de volumen. En un equilibrio exacto, eliminar el potencial da

\[
\partial_\tau\delta d=-A H\delta d,\qquad
A=M^{-1}+\tfrac12 S L^{-1}S^T.
\tag{P4}
\]

La segunda contribución es la disipación normal que sigue del balance
de corriente; no modifica la física del circuito. Para una reducción de
Galerkin que preserve la disipación cuadrática, la métrica correspondiente
es

\[
Q=A^{-1}=M-MS(2L+S^TMS)^{-1}S^TM,
\qquad
(V^TQV)\dot a=-(V^THV)a.
\tag{P5}
\]

Esta fórmula permite una base temporal pequeña, pero **una base arbitraria
de pocos perfiles no contiene todas las variaciones locales de fase**.
No debe declararse continuidad eléctrica nodo a nodo sólo porque P5
conserve la energía reducida. La propuesta práctica es construir una base
Krylov/Galerkin para aproximar el operador de todos los nodos, informar
su residuo en el sistema espacial completo y ampliarla cuando ese residuo
impida interpretar la respuesta. Las corrientes se reconstruyen con la
derivada espectral P3; cualquier defecto de continuidad asociado a la
reducción se muestra y no se cancela mediante cargas de restricción
inventadas. Resolver la exponencial de la pequeña matriz permite muestrear
tiempo físico sin imponer el límite de estabilidad del Euler explícito.

El núcleo aceptado tiene un residuo G0 pequeño pero distinto de cero.
No se lo redefine como equilibrio exacto ni se repite la relajación con
una tolerancia más exigente por ese motivo. Para el Jacobiano exacto del
RHS térmico F(d)=−A(d)G(d), la aproximación afín es

\[
\dot{\delta d}=-A_0G_0+
\left[-A_0H\delta d-(DA[\delta d])G_0\right].
\tag{P6}
\]

Se calculan una evolución base y otra inicialmente perturbada y se compara
su diferencia. Así el pequeño movimiento de la referencia permanece
contabilizado. La derivada DA usa los coeficientes KWT y S; no exige
espectros adicionales. Congelar A en P5 descarta un término de orden G0δd:
puede ser una aproximación útil, pero su tamaño se registra antes de
atribuirle el Jacobiano exacto.

Una perturbación inicial concreta es
d=d0 exp[iε(x/4)(1−r²/16)³_+], con ε=0,001 rad, complementada por una
variación de amplitud del mismo soporte. Ambas dejan intacto el contorno.
El horizonte se elige a partir de las tasas obtenidas del propio operador
KWT, no acortando tiempos físicos ni suponiendo que una escala de difusión
representa también la relajación de fase. Basta contrastar el límite
lineal con unas pocas evaluaciones no lineales de fuerza en estados de la
trayectoria y una perturbación de amplitud menor; no se requiere un
transiente fotónico ni una recuperación de nanosegundos para ese control.

Este trabajo es un hito legítimo de etapa 4: verifica que el núcleo térmico
aceptado tenga una evolución mesoscópica coherente con KWT, potencial y
disipación. No obliga todavía a una decisión física del usuario. Tampoco
sustituye el trabajo espectral y el acoplamiento con poblaciones que siguen
pendientes para cerrar el sistema no térmico y avanzar al detector.

La preparación de P3/P6 ya se ejecutó: 256 frecuencias y 2304 raíces
espectrales independientes terminaron en 32,94 s. Las diferencias finitas
con pasos 0,02 y 0,01 verifican el Hessiano, la corriente y el RHS completo;
el máximo error relativo del RHS a paso 0,01 es 4,86×10⁻⁵ y disminuye
aproximadamente por cuatro al dividir el paso por dos. La corrección por
G0 y el potencial base representa 1,10–1,13 % del RHS de las direcciones
ensayadas: conservarla evita una aproximación visible e innecesaria.
La continuidad nodal queda a 1,67×10⁻¹⁶ y el residuo de disipación de
energía libre por debajo de 1,26×10⁻¹². Estos números verifican el operador
del sistema térmico discretizado; todavía no son una trayectoria ni una
precisión atribuida al detector. Los comprobantes y sus hashes están en
[thermal_weak/verification_receipt.json](thermal_weak/verification_receipt.json).

El siguiente control reutiliza las matrices guardadas en todos los nodos
para evolucionar el sistema afín P6. Una ventana de hasta 1 ps y amplitud
inicial 10⁻³ permiten observar el inicio de la relajación sin modificar
ningún tiempo material. La referencia y la trayectoria perturbada deben
permanecer dentro de la vecindad declarada, inicialmente 2 % de Δref.
El balance no lineal de disipación se verificó en los estados de prueba
anteriores; no se atribuye automáticamente como identidad exacta a la
trayectoria del Jacobiano afín.

## Contratos del dispositivo que permanecen vigentes

Los ensayos anteriores son problemas de núcleo con contorno prescrito, sin
puerto eléctrico del nanohilo. No generan Vout ni permiten conectar un
circuito ficticio al borde circular para atribuirles una señal de detector.
Al volver al ensayo espacial con terminales se conserva la
[adenda del circuito completo de la memoria](../../../modelo_v0_4/actualizaciones/circuito_memoria_20260922.md):
los estados Ib, Is y vc, la fuente, resistencias, inductancias y capacitancia
declaradas, y su balance de potencia. La energía de superflujo de la región
resuelta no se vuelve a añadir mediante una inductancia exterior.

No se cambian constantes circuitales ni parámetros materiales para acelerar
la respuesta. El horizonte posterior del fotón seguirá siendo el cruce de
Vout acordado más su margen; acortar ese horizonte no modifica el sistema.
La presente decisión no necesita una nueva elección física del usuario:
autoriza avanzar con un diagnóstico temporal acotado y mantiene explícito
qué falta antes del transiente acoplado.
