# Cinética con espectro móvil: una construcción matricial comprobable

24 de septiembre de 2026. Este documento fija convenciones y una construcción
ejecutable de los términos temporales; no presenta un transiente de detector
ni declara terminada la unión física. La implementación independiente está en
`pysnspd/experimental/adiabatic_matrix_usadel.py`.

La ecuación de Usadel y su normalización contienen el mismo producto temporal.
Por ello no basta añadir una masa a la cinética estacionaria: también cambia
la relación entre espectro, distribución y fuerza del gap. El punto de partida
publicado es la convolución de [Sauls, ecuaciones 7–12 y 47–50](https://arxiv.org/pdf/2202.02260).
Las ecuaciones estacionarias y la parametrización de distribución proceden de
[Belzig et al., sección 2.4](https://arxiv.org/html/cond-mat/9812297v2).
La traducción al grafo y las fórmulas siguientes se derivan aquí con las
convenciones del repositorio; no son una transcripción de otra convención Nambu.

## Variables, signos y unidades

Se usa `E = energía/(kB Tc)` y tiempo físico `t`. Definimos
`κ = ħ/(2 kB Tc)`; si el tiempo se expresa en unidades de κ, el código recibe
`kappa=1`. La matriz del gap es

\[
D=\begin{pmatrix}0&d\\d^*&0\end{pmatrix},\quad
d=\Delta/(k_BT_c),\quad B=D-iE\tau_3+i v I,
\quad v=e\phi/(k_BT_c),\quad e>0.
\tag{K.1}
\]

`φ` es potencial eléctrico, `I` la identidad 2×2, y `τ3=diag(1,−1)`.
El término `η τ3` se permite exclusivamente como regulador del problema
retardado/avanzado. No se lo convierte en colisión, baño o tasa de relajación.
La transformación de calibre es `d → exp(iχ)d`,
`v → v−κ ∂tχ`. Para un cambio lento de calibre, el último desplazamiento es
de primer orden temporal. No se debe esconderlo dentro del campo estacionario
y después volver a añadirlo como corrección.

## Producto de primer orden, aplicado una sola vez

Escribimos `X=X0+X1`, donde `X1` es de primer orden temporal. El producto
retenido es

\[
(X\circ Y)_0=X_0Y_0,
\quad
(X\circ Y)_1=X_1Y_0+X_0Y_1
 +i\kappa(\partial_EX_0\partial_tY_0
             -\partial_tX_0\partial_EY_0).
\tag{K.2}
\]

Se descartan productos `X1 Y1` y derivadas temporales de correcciones de
primer orden. `Jet` almacena precisamente el valor de orden cero, sus dos
derivadas y la corrección independiente. La operación es asociativa hasta
el orden retenido; esto está comprobado con matrices no conmutativas.

En el enlace `i→j`, `Pij=diag(exp(−iαij/2),exp(iαij/2))` y
`X̄j=Pij∘Xj∘Pij†`. Si cambia el enlace, también se deriva Pij. Usar sólo
`Pij Xj Pij†` durante una transformación temporal general perdería términos.

## Espectro: la nueva corrección no es siempre traceless

Sean `mi` las áreas duales y `cij` las conductancias del grafo. El residuo
retardado y su normalización son

\[
Q^R_i=\sum_j\frac{c_{ij}}2[R_i,\bar R_j]_\circ
       -\frac{m_i}{2}[B_i,R_i]_\circ=0,
\qquad R_i\circ R_i=I.
\tag{K.3}
\]

El oráculo existente proporciona R0. La corrección R1 satisface el sistema
lineal formado por el Jacobiano de K.3 y

\[
\{R_0,R_1\}=-i\kappa[\partial_ER_0,\partial_tR_0].
\tag{K.4}
\]

Un particular es
`R1,p=−(iκ/2) R0 [∂ER0,∂tR0]`; la parte restante pertenece al espacio
tangente ` {R0,R1,t}=0`. En general el particular tiene componente identidad.
Por eso **no basta reutilizar sólo las tres componentes traceless `(g,f,f̃)`
para R1**. El nuevo solucionador conserva las cuatro entradas, apila residuo
espectral y normalización y resuelve con LSMR disperso. Entrega ambos residuos
por separado: un código de terminación de mínimos cuadrados no es aceptación.

Este es un problema lineal por energía, paralelizable; no es un solver de
dos tiempos. Las derivadas de R0 deben provenir del mismo oráculo, usando
sus tangentes o diferencias verificadas. En los contactos hay que prescribir
la corrección del reservorio correspondiente, no imponer cero durante un
cambio de calibre.

## Distribuciones dinámicas y vértice de trabajo

\[
A=-\tau_3R^\dagger\tau_3,
\quad h=h_LI+h_T\tau_3,
\quad K=R\circ h-h\circ A.
\tag{K.5}
\]

Aquí `hL` es modo de energía y `hT` desequilibrio entre ramas. En particular,

\[
K_1=R_1h-hA_1+i\kappa
 (R_Eh_t-R_th_E-h_EA_t+h_tA_E).
\tag{K.6}
\]

El flujo cinético y el residuo son la extensión directa del operador congelado:

\[
J^K_{ij}=\frac{c_{ij}}2
 (R_i\circ\bar K_j+K_i\circ\bar A_j
  -\bar R_j\circ K_i-\bar K_j\circ A_i),
\quad
Q^K_i=\sum_jJ^K_{ij}-\frac{m_i}2[B_i,K_i]_\circ.
\tag{K.7}
\]

En K.7 B excluye el η numérico. Los operadores de colisión físicos deben
añadirse con sus momentos conservados y la misma autoenergía cuando proceda;
el módulo no inventa ese cierre. Las proyecciones son `Tr(Q)/4` y
`Tr(τ3 Q)/4`, tomando la parte real para campos físicos instantáneos.

Dos reducciones analíticas comprueban signo y normalización:

* En el metal normal, `R=τ3, A=−τ3`, se obtiene exactamente
  `rL/T = Σj cij(hj−hi) − mi κ ∂t hL/T`. La masa no se deduce de un
  ajuste ni de la parte imaginaria del regulador.
* Para un gap real uniforme móvil, con `ρ=Re g` y `R2=Im f`, la proyección
  temporal longitudinal es
  `rL,1 = −mi κ [ρ ∂t hL + R2 ∂t d ∂E hL]`.
  Se cancela la aparente fuente `hL ∂tρ` mediante la identidad espectral
  `∂d Re g=−∂E Im f`. Se recupera así el trabajo publicado por
  [Vodolazov 2017, ecuación 1](https://arxiv.org/pdf/1611.06060), con la
  conversión `hL=1−2p`. No se añadió calor para forzar esa igualdad.

El segundo resultado exige derivar también K. Ignorar K1 o usar una DOS
móvil multiplicada por una masa arbitraria no recupera esa cancelación.

## Potencial, neutralidad y fuerza del gap

Con los signos K.1, el reservorio normal de pequeño potencial tiene
`hT=−v ∂E tanh[E/(2T)]`. El momento de neutralidad se escribe

\[
v+\int_0^\infty\frac14\operatorname{Re}\operatorname{Tr}K(E)\,dE=0.
\tag{K.8}
\]

La unidad común de carga se cancela en esa ecuación. Deben conservarse K1,
la sustracción del fondo y las condiciones de cola. En el límite normal
K.8 recupera exactamente el reservorio anterior. Bajo un calibre temporal,
el cambio del momento integral compensa `−κ ∂tχ`: la neutralidad no impone
`hT=0`. Este mecanismo es la traducción del término de carga local de
[Sauls, ecuación 5](https://arxiv.org/pdf/2202.02260); su prefactor depende de
la convención de densidad normal, pero el signo relativo aquí queda fijado
por K.1 y el límite normal.

La identidad estacionaria entre torque y divergencia de corriente se amplía
con almacenamiento y términos temporales. El criterio correcto es obtenerla
de K.7 y de la **misma** ecuación de gap y neutralidad. No se añade en paralelo
la vieja ecuación de Poisson óhmica: la corriente normal ya está dentro del
momento Keldysh. El circuito de la memoria se conecta por su puerto de tensión
y corriente, conservando sus tres estados físicos.

## KWT: dos cierres que no deben sumarse sin explicar el reparto

Un cierre puramente Usadel impone autoconsistencia del momento anómalo,
incluido K1. Ese momento ya determina una respuesta temporal del gap; añadir
encima la misma respuesta como fricción KWT no está justificado. Esto no
significa que cualquier modelo efectivo que combine cinética y movilidad
fenomenológica esté prohibido. Significa que debe especificar qué parte queda
resuelta y qué disipación residual se añade.

Una opción híbrida conservativa tiene ecuación de fase `Im G=0`, ecuación
radial `G_rad+Γ_res ∂t|d|=0` y conversión a calor
`Γ_res(∂t|d|)^2`, donde **G incluye el momento temporal retenido**.
Para reproducir la movilidad térmica KWT admitida se requiere una comparación
de respuesta que determine el residuo Γ_res y compruebe su pasividad. Usar
G0 en la ecuación radial pero G0+G1 en el balance perdería
`G1·∂td`, del mismo orden cuadrático que ese calor. Tampoco se puede asumir
que Γ_res es positivo o igual al Γ_KWT heredado sin calcular la respuesta.

No se ha adoptado ni calibrado esa opción en este módulo. El aviso anterior
es una condición de ensamblaje, no otra tolerancia exigida al benchmark Euler.
Un problema collisionless gapped puede tener una matriz de movilidad de
primer orden singular; rellenar ese rango con η equivaldría a inventar un baño.

## Ensayo acoplado a frecuencia finita: alternativa concreta

Para un ensayo débil sobre una referencia suave polarizada, se puede evitar
la expansión temporal y conservar exactamente los desplazamientos energéticos.
Con perturbación `δX exp(−iωt)`, `E±=E±κω`,

\[
\delta(X\circ Y)=X_+(E)\delta Y(E)+\delta X(E)Y_-(E).
\tag{K.9}
\]

Se resuelven, por energía, la linealización de K.3 y
`R+ δR+δR R−=0`, y por separado la ecuación avanzada. La componente de
distribución es

\[
\delta K=R_+\delta h-\delta h A_-+
          \delta R h_- -h_+\delta A.
\tag{K.10}
\]

La proyección de la ecuación cinética, el momento del gap, la neutralidad y
el circuito forman entonces un único sistema complejo. No es una suma de
ensayos desacoplados. Permite observar amplitud, fase, ambos modos cinéticos,
corriente y Vout, además de la potencia media del puerto; no pretende
reemplazar un transiente fuerte con fotón.

Tres detalles de implementación son esenciales:

1. No usar `.real` al proyectar amplitudes de Fourier. Son complejas; sólo la
   reconstrucción del campo físico toma la parte real.
2. Usar las dos amplitudes de cuadratura `x,y` del gap:
   `δD12=x+i y`, `δD21=x−i y`. No conjugar `δD12` para construir la segunda,
   porque conjugaría también la frecuencia. Los momentos correspondientes son
   `Gx=i m(δK12+δK21)/2` y `Gy=m(δK12−δK21)/2`.
3. Usar un intervalo de energías firmado o demostrar explícitamente el
   plegado; E+ y E− pueden cruzar cero. La corrección dinámica del momento
   se resta a su respuesta estática de equilibrio y se suma al Hessiano
   Matsubara existente, para conservar su referencia y evitar cancelaciones
   ultravioletas. Se debe comprobar el límite ω→0.

La existencia de esta formulación no sustituye la preparación numérica del
ensayo. Sí da un camino de implementación publicado y verificable, sin crear
un potencial ficticio ni una energía escalar a partir del resultado deseado.

## Contrato concreto de un híbrido radial, si se admite tras la comparación

Esta alternativa queda completamente especificada sólo si se satisfacen
conjuntamente las siguientes condiciones; no se sustituye por una suma de
fuentes independientes:

\[
\begin{aligned}
 &Q^R_0=0,\quad Q^R_1=0,\quad R\circ R=I+O(\partial_t^2),\\
 &\Pi_{L,T}Q^K=\mathcal C_{L,T}+\mathcal S_{L,T}^{\rm calor},\\
 &G_r+\Gamma_{\rm res}\,\partial_t|d|=0,\qquad G_\vartheta=0,\\
 &v+\int_0^\infty\operatorname{ReTr}K/4\,dE=0,\\
 &P_{\rm radial}=\Gamma_{\rm res}(\partial_t|d|)^2\ge0.
\end{aligned}
\tag{K.11}
\]

`G_r=Re(exp(−iϑ)G)` y `Gϑ=Im(exp(−iϑ)G)` son las componentes de la
fuerza anómala **incluyendo K1**. El ángulo ϑ es la fase del gap; no se
añade otra relajación KWT de fase. La dinámica de esa fase procede de las
ecuaciones cinéticas, neutralidad y la restricción imaginaria. El calor
radial se reparte con el cierre B.41 ya admitido y se incorpora una vez a
las poblaciones. El término de Joule no se duplica si el momento de la
cinética ya contiene el trabajo electromagnético. La referencia de energía
necesaria es energía interna fuera de la autoconsistencia, no sólo la energía
libre del benchmark térmico.

En un reservorio se prescriben el espectro, la distribución y su potencial
del **mismo** reservorio. Una perturbación de tensión en una superconductor
también cambia su fase según la convención de calibre. En los bordes aislantes
se anulan los flujos de los mismos enlaces. El puerto utiliza la corriente
integrada Keldysh y la diferencia de potencial eléctrico; el circuito conserva
sus estados, resistencias e inductancias. Ningún contacto ideal absorbe energía
sin que aparezca su potencia/flujo en el balance.

La integral sobre energías positivas supone el plegado de simetría usual para
campos reales instantáneos. Para el ensayo armónico se implementa primero
`E ∈ [−Emax,Emax]`: el momento equivalente de carga es
`v+∫ReTr K/8 dE`. Los momentos anómalos y de corriente requieren el mismo
factor de plegado y las dos cuadraturas. Es obligatorio medir el residuo de
cola y comparar las dos cuadraturas, porque un error de factor dos podría
aparentar una compensación de calor o carga.

### Cómo determinar la movilidad residual sin confundir coordenadas y disipación

La comparación debe resolver la respuesta de distribución provocada por el
gap móvil **antes** de extraer el coeficiente del momento anómalo. Mantener
h fijo no reproduce la evolución que se quiere calibrar. En una referencia
con el mismo material, temperatura y polarización, se elimina el bloque
cinético/neutralidad de la respuesta lineal y se obtiene el coeficiente radial
efectivo de baja frecuencia. Su parte disipativa se compara con Γ_KWT en las
mismas unidades. Sólo entonces puede usarse
`Γ_res=Γ_KWT−Γ_resuelto`.

Si la parte resuelta no es local, no puede reemplazarse sin explicación por
un número nodal. Si Γ_res tiene un autovalor disipativo negativo que supera
la incertidumbre numérica relevante, el híbrido no tiene una interpretación
de baño pasivo: debe mantenerse el cierre microscópico o revisar la movilidad
efectiva. Esto no se arregla recortando autovalores. Una corrección menor que la
escala de error del observable puede admitirse con la política práctica, pero
no se silencia una inestabilidad o una identidad de conservación incumplida.

### Resultado del piloto uniforme ejecutado

El script `sandbox/stage4_core/probe_dynamic_gap_mobility.py` realiza una
cuadratura analítica, de 0,016 s, con `d=1,76392593`, `Tb/Tc=0,9/8,65` y gap
real uniforme. El coeficiente mostrado es `G1/(m κ ∂td)`, sin unidades; no es
un tiempo medido. Para una velocidad radial unitaria resuelve

\[
\partial_t h_L=-\frac{\operatorname{Im}f}{\operatorname{Re}g}\,
                \partial_Eh_0,
\quad
\frac{G_{r,1}}{m\kappa\partial_td}
=2\int_0^\infty\left[
 \operatorname{Re}(\partial_df)\,\partial_Eh_0
 -\operatorname{Re}(\partial_Ef)\,\partial_th_L\right]dE.
\tag{K.12}
\]

| η/|d| | Manteniendo h fijo | Evolucionando h con K.12 |
|---:|---:|---:|
| 0,02 | −0,0135781 | 0,000473890 |
| 0,01 | −0,0139617 | 0,000118702 |
| 0,005 | −0,0140580 | 0,000029812 |
| 0,002 | −0,0140852 | 0,000004875 |

El coeficiente aparentemente negativo con h fijo desaparece al conservar la
respuesta cinética, y el restante disminuye aproximadamente como η² en esta
secuencia. El residuo cinético máximo es menor que 4×10⁻¹⁸. Las cuadraturas
Gauss de 32 y 64 puntos por panel están guardadas en
[uniform_mobility_probe.json](uniform_mobility_probe.json).

En el límite de gap ideal la masa subgap es cero: la distribución allí no
es una población observable independiente. El regulador permite evaluar
coherentemente ese límite, pero su pequeña contribución no debe identificarse
como amortiguamiento físico. Este piloto respalda la ausencia de una gran
fricción radial resuelta en el caso uniforme collisionless; **no calibra** el
caso polarizado espacial ni autoriza extrapolar una Γ_res nodal a una hotbelt.
Una sola comparación acoplada alrededor de la referencia espacial prevista
es el siguiente cálculo útil; no hace falta repetir los benchmarks ya aprobados.

### Tratamiento consistente de η y de las colisiones

Excluir η del residuo cinético evita la fuga artificial ya conocida, pero
no vuelve exacta a η finito una respuesta cuyo espectro sí usa η. Los residuos
de las identidades dinámicas y la variación con η deben acompañar el resultado.
No se extrapola un dato a η=0 inventando una ley; dos o más valores permiten
medir si el cambio afecta los momentos observables. Las identidades de calibre
usan el mismo η retardado/avanzado y no lo reinterpretan como un baño.

Un tiempo inelástico físico, en cambio, necesita su autoenergía y componente
Keldysh o una aproximación de colisión que conserve los momentos declarados.
No basta sustituir η por `ħ/(2τ)` sólo en el propagador. Esa sustitución haría
desaparecer energía/carga sin registrar su destino y no definiría una mejora
numérica del mismo sistema.

## Evidencia ejecutada

`tests/test_adiabatic_matrix_usadel.py`: seis pruebas pasaron en 0,016 s
en el entorno Python local. Comprueban difusión y almacenamiento normales,
trabajo radial publicado, invariancia bajo calibre temporal con componente
identidad no nula, reconstrucción dispersa de esa corrección y asociatividad
de la expansión, incluyendo un calibre que depende de espacio y tiempo con
derivadas de los enlaces. Son controles analíticos de las ecuaciones nuevas; no se
presentan como aceptación de latencia, material NbN o transiente acoplado.

El bloque de frecuencia finita ya está implementado en
`pysnspd/experimental/harmonic_kinetic_usadel.py`. Su matriz dispersa tiene
dos amplitudes complejas de distribución por nodo y conserva las fuentes de
gap y potencial incluidas en K.10. No aplica `.real` a esas amplitudes.
Seis pruebas de `tests/test_harmonic_kinetic_usadel.py` pasaron en 0,024 s:
difusión/almacenamiento normales, calibre temporal con resolución de contactos,
límite estacionario del operador congelado, conservación de ambas cuadraturas
del gap, identidad dinámica de carga en cada energía y equivalencia entre un
modo uniforme reducido y el grafo completo. El control de calibre utiliza
también las respuestas R/A obtenidas del solucionador espectral, conectando
las implementaciones efectivas en vez de comprobar sólo fórmulas aisladas.

Para un espectro de referencia uniforme y un modo `L u = λ M u` del laplaciano
positivo, la divergencia matricial por masa y amplitud de modo es

\[
\delta\operatorname{div}J^K/M=
\frac\lambda2(\delta R K_-+\delta K A_- -R_+\delta K-K_+\delta A).
\tag{K.13}
\]

El helper `modal_coefficients` conserva K.13 y el término local completo,
devolviendo la matriz 2×2, la fuente y los coeficientes para reconstruir δK.
Su equivalencia con el grafo se verificó a 3×10⁻¹⁴, incluyendo fuentes
espectrales y de potencial. Esta reducción sólo vale cuando el espectro de
referencia es uniforme: no se usan promedios de coeficientes para simular un
fondo espacial polarizado como si sus modos estuvieran desacoplados.

El ensamblaje vectorizado del ensayo se revisó mediante
`tests/test_coupled_response.py`: cinco pruebas pasaron en 0,630 s. La revisión
comprueba seis energías firmadas y tres columnas de forzamiento frente al
operador modal independiente, conductancia normal positiva `Ggrafo/R□`,
identidad dinámica integrada de carga, desplazamiento puro de calibre y
límites de rigidez Matsubara. La función `film_port_coefficients` permite
comprobar el metal normal antes de asignar inductancias externas; no se
inventa una inductancia superconductora para un gap exactamente cero.

## Consulta ejecutable sobre una referencia espacial polarizada

`sandbox/stage4_core/biased_harmonic_query.py` resuelve ahora cada consulta
energética sobre **todos los nodos** de una referencia polarizada. Sigue un
camino recto declarado desde el ancla Matsubara hasta E+, E− y E, con
contactos BCS, comprobación de DOS causal y normalización. Factoriza las
ecuaciones espectrales retardada, avanzada y estática una vez por energía,
y reutiliza esas factorizaciones y la cinética para todas las columnas de
gap y potencial.

La base espacial únicamente limita el número de coordenadas del gap y del
potencial para el ensayo exterior. No reemplaza por promedios las respuestas
espectrales ni las dos distribuciones: ambas se resuelven en todo el grafo.
Las columnas interiores desaparecen en los contactos. La columna del puerto
conserva la fase Josephson del reservorio y su distribución desplazada,
`δhT=−lift·(h+−h−)/Ω`. Las derivadas estáticas de contactos se incluyen
también en la sustracción de referencia.

La salida por energía distingue fuerza nodal integrada, densidad de fuerza,
densidad de carga y corriente de enlace. Las integrales firmadas utilizan el
factor 1/2 común. Se guardan las partes dinámica y estática por separado y
los residuos del sistema realmente resuelto. Una amplitud de respuesta
unitaria no se etiqueta como población finita admisible: su escalamiento
pertenece al ensayo débil exterior.

El piloto `tests/test_biased_harmonic_query.py` pasó en 0,207 s, con 15 nodos,
fase impuesta de 0,2 rad, dos modos interiores y siete columnas de forzamiento.
Incluye la unión real de referencia Matsubara, continuación, respuesta
retardada/avanzada, cinética, contactos y tangente térmica del orquestador.
Detecta una respuesta hT no nula producida por una perturbación radial en el
fondo polarizado. Este piloto valida la conexión de las implementaciones;
la aceptación cuantitativa requiere la campaña espacial y sus comparaciones
independientes, que ejecutará el usuario si supera cinco minutos.

La cola energética se revisó con la misma consulta espacial, sin ajustar una
corrección analítica al resultado. Para la referencia pequeña con fase de
0,4 rad, los momentos proyectados sumados en ±E disminuyen por un factor
7,992 entre E=128 y 256, compatible con la potencia asintótica E⁻³. La
variación de E³ por el momento es 0,181 % entre esos puntos y 0,0632 % hasta
E=12000. Al ajustar la tolerancia espectral de 10⁻⁸ a 10⁻¹¹ en el último
punto, el cambio relativo es 0,0193 %; multiplicado por el jacobiano E²/32
de la cuadratura de cola, el cambio absoluto máximo es 1,55×10⁻⁸. Los
[datos del control](high_energy_worker_check.json) distinguen la tolerancia
del solver de la escala del observable. No se identificó una pérdida
catastrófica de precisión al evaluar los puntos grandes de esa transformación.

`tests/test_biased_harmonic_query.py` incluye ahora cuatro controles,
aprobados en 4,823 s: unión sesgada, puerto uniforme cerca de los bordes del
gap, cola firmada hasta E=12000 y cancelación real de los trabajadores tras
un fallo. El último control verifica que el error no deje ejecutándose la
cola de una campaña larga; sólo termina los procesos que pertenecen a ese
pool y conserva los resultados ya escritos.
