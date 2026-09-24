# Puente no térmico: diagnóstico de carga que decide el siguiente paso

24 de septiembre de 2026. El siguiente cálculo propuesto usa los espectros
ya obtenidos para medir una omisión concreta del cierre escalar B: si un
gradiente de ocupaciones electrónicas genera una respuesta de carga que
retroalimenta el transporte de energía. No requiere otro transiente térmico,
una cascada fotónica ni valores nuevos de parámetros materiales.

Se implementa como **respuesta lineal estática a espectro y gap congelados**.
Su salida es una dirección de respuesta, no un estado fotónico preparado
ni la adopción de una dinámica instantánea para la carga. Con él se puede
comprobar una unión entre fuerza, corriente y conversión de pares antes de
decidir la evolución temporal que falta.

## Qué se pone a prueba de B

B.13 evoluciona la ocupación simétrica electrón–hueco mediante difusión de
energía. Esa elección conserva información que una temperatura equivalente
no contiene, pero B.13 declara expresamente que no resuelve toda la cinética
de carga. B.46–B.51 justifican sus fuerzas en un espectro localmente uniforme;
la implementación espacial nueva ya permite medir dónde deja de ser cerrado
ese sector escalar.

Escribimos la distribución matricial como
$h=h_L I+h_T\tau_3$. El modo $h_L=1-2p$ representa la población simétrica;
$h_T$ mide la diferencia entre las ramas electrón y hueco. Un material con
simetría electrón–hueco no implica automáticamente $h_T=0$ en presencia
de supercorriente y gradientes de población. Las ecuaciones cinéticas
incluyen difusión longitudinal y transversal, supercorriente espectral y
un coeficiente cruzado. Véanse [Virtanen y Heikkilä, ecuaciones 7–8](https://arxiv.org/pdf/cond-mat/0405285)
y [Belzig et al., secciones 2.4.2–2.4.3](https://arxiv.org/html/cond-mat/9812297v2).

No se copia la condición de divergencia cero de un metal normal a todo
el núcleo superconductor: la conversión de carga entre excitaciones y
condensado debe permanecer en su ecuación.

## Operador lineal sobre el mismo grafo

Sean $m_i$ las masas y $c_{ij}$ las conductancias del grafo térmico,
$d_i=\Delta_i/(k_BT_c)$ y $z=\eta-iE$, con energía adimensional.
A partir de las tres componentes espectrales ya almacenadas se construye

\[
R_i=\begin{pmatrix}g_i&f_i\\\widetilde f_i&-g_i\end{pmatrix},\qquad
A_i=-\tau_3R_i^\dagger\tau_3,
\qquad
K_i=R_ih_i-h_iA_i,
\qquad
D_i=\begin{pmatrix}0&d_i\\d_i^*&0\end{pmatrix}.
\]

$D_i$ es aquí la matriz del gap, no la difusividad material. $R_i$ y $A_i$
son los propagadores retardado y avanzado; $K_i$ es su componente de
distribución. La tilde de $\widetilde f$ identifica la segunda componente
anómala y no una conjugación. Esta notación evita confundirla con la
ocupación $p$.

En la arista orientada $i\to j$, se transporta toda matriz del nodo $j$
mediante $P_{ij}=\operatorname{diag}(e^{-i\alpha_{ij}/2},e^{i\alpha_{ij}/2})$:
$\overline X_j=P_{ij}X_jP_{ij}^{\dagger}$. El flujo matricial es

\[
J^K_{ij}=\frac{c_{ij}}2
 \left(R_i\overline K_j+K_i\overline A_j
       -\overline R_jK_i-\overline K_jA_i\right).
\]

Su orientación inversa se obtiene con transporte y signo opuesto. Las dos
proyecciones reales de flujo y residuo nodal se definen por

\[
j_{L,ij}=\tfrac14\operatorname{Re}\operatorname{Tr}J^K_{ij},\qquad
j_{T,ij}=\tfrac14\operatorname{Re}\operatorname{Tr}(\tau_3J^K_{ij}),
\]

\[
Q_i=\sum_jJ^K_{ij}-\frac{m_i}{2}[D_i,K_i],\qquad
r_{L,i}=\tfrac14\operatorname{Re}\operatorname{Tr}Q_i,\qquad
r_{T,i}=\tfrac14\operatorname{Re}\operatorname{Tr}(\tau_3Q_i).
\tag{P.1}
\]

En las fórmulas exactas las proyecciones son reales. La parte imaginaria
numérica sirve para detectar inconsistencias, no para ocultar una componente
física. El término local de conversión tiene traza cero, de modo que no
crea un sumidero de energía en esta descripción estacionaria y elástica.
Sí aparece en la ecuación de carga. La **fuente del segundo miembro**
$S_T=\operatorname{Re}\operatorname{Tr}(\tau_3 m[D,K]/2)/4$ es

\[
\frac{m_i}{2}\left\{
[\operatorname{Re}(d_i\widetilde f_i)-\operatorname{Re}(d_i^*f_i)]h_{L,i}
+[\operatorname{Re}(d_i\widetilde f_i)+\operatorname{Re}(d_i^*f_i)]h_{T,i}
\right\}.
\]

El residuo P.1 contiene $-S_T$. Para gap uniforme real y fase espectral
alineada, esa contribución al residuo es $-m d\operatorname{Re}(f)h_T$.

Se conserva también el término proporcional a $h_L$. Puede ser distinto
de cero cuando la fase espectral difiere de la del gap, y equilibra una
divergencia de supercorriente espectral.

Las ecuaciones P.1 son lineales en las dos distribuciones. Se pueden ensamblar
en una matriz dispersa real aplicando las dos bases $I,\tau_3$ en cada nodo.
No requieren diferencias de grandes energías ni nuevas consultas espectrales.
En el estado normal $R=\tau_3,A=-\tau_3$ dan exactamente
$j_{L/T,ij}=c_{ij}(h_{L/T,j}-h_{L/T,i})$, lo que fija su normalización.

### Coeficientes continuos que explica la matriz

En el límite espacial continuo se recupera

\[
j_L=\mathcal D_L\nabla h_L-\mathcal T\nabla h_T+j_S h_T,\qquad
j_T=\mathcal D_T\nabla h_T+\mathcal T\nabla h_L+j_S h_L,
\]

\[
\mathcal D_{L/T}=\frac{1+|g|^2}{2}
\mp\frac{|f|^2+|\widetilde f|^2}{4},\qquad
\mathcal T=\frac{|\widetilde f|^2-|f|^2}{4}.
\]

$j_S$ es la proyección de la diferencia entre las corrientes espectrales
retardada y avanzada. La implementación de P.1 retiene todos estos términos
con los mismos enlaces; no aproxima cada coeficiente mediante un promedio
independiente. Por ejemplo, imponer $h_T=0$ todavía deja en la ecuación de
carga los gradientes de $\mathcal T\nabla h_L+j_Sh_L$ y la conversión local.

## Qué hacer con el regulador numérico eta

Los campos espectrales satisfacen

\[
\sum_jJ^{R/A}_{ij}=\frac{m_i}{2}[B_{R/A,i},R_i\text{ o }A_i],\qquad
B_R=D+(\eta-iE)\tau_3,\quad B_A=D+(-\eta-iE)\tau_3.
\]

Si se usa sin cuidado el término local $B_RK-KB_A$ en la ecuación
de distribución, aparecen tanto $[D-iE\tau_3,K]$ como
$\eta\{\tau_3,K\}$. El conmutador de energía no aporta a las
proyecciones $L/T$, pero el segundo término aporta exactamente

\[
\text{fuga artificial}_{L/T}=m_i\eta\rho_i h_{L/T,i},\qquad
\rho_i=\operatorname{Re}g_i.
\tag{P.2}
\]

Un baño físico con esa tasa necesitaría su propia distribución y autoenergía.
No existe tal baño en el modelo admitido: aquí $\eta$ sólo desplaza el
contorno de evaluación espectral. El diagnóstico P.1 **excluye explícitamente**
P.2 y registra por separado su norma sobre cada dirección $\delta h$, para
que no parezca una relajación de carga medida. Esto no elimina la sensibilidad
de los coeficientes a $\eta$; se comparan los dos valores ya calculados.

Para $h_L(E)$ uniforme en el espacio y $h_T=0$, P.1 se anula en nodos libres
por las ecuaciones espectrales. Los términos restantes son conmutadores con
$\tau_3$, cuyas dos trazas proyectadas son cero. La recuperación de esta
distribución de equilibrio es, por tanto, una identidad del grafo, salvo su
residuo espectral; no depende de atribuir a $\eta$ una tasa cinética.

## Un diagnóstico que produce una decisión

Se prescribe una dirección pequeña $\delta h_L(E,\mathbf r)$ con soporte
espacial interior. Dos perfiles bastan: una perturbación radial y otra que
rompe esa simetría, por ejemplo angular o desplazada. Son funciones de prueba
del operador, sin relación con la preparación gaussiana del fotón que sigue
abierta. El caso radial ayuda a identificar cancelaciones por simetría; en
una malla cartesiana finita no se exige que sea un cero algebraico.

Para cada espectro ya calculado se evalúa primero el residuo omitido
$r_T[\delta h_L,0]$. Después se resuelve únicamente

\[
L_{TT}\,\delta h_T=-L_{TL}\,\delta h_L,
\qquad \delta h_T=0\text{ en los nodos de distribución fijados}.
\tag{P.3}
\]

Los bordes térmicos significan que la perturbación de distribución desaparece
allí; no cambia el entorno espectral radial. El resultado compara:

- el residuo de carga que B escalar deja fuera;
- la amplitud y estructura de la respuesta $\delta h_T$;
- el cambio que esa respuesta produce en el flujo de energía y en su
  divergencia, respecto de $h_T=0$;
- la sensibilidad a $\eta$ y a las dos mallas ya disponibles;
- la fuga P.2 que se habría introducido si el regulador se confundiera con
  un baño.

Un residuo de carga no nulo demuestra que el sector cinético $h_T=0$ no
es cerrado en este ensayo de **gap, fase y potencial congelados**. No
demuestra por sí solo que falle el modelo mesoscópico completo de la
memoria: su fase KWT y su potencial normal también pueden responder.
El impacto sobre los momentos de fuerza, corriente y energía orientará
la comparación con una reducción electroquímica admisible. No se convierte un
porcentaje arbitrario nuevo en certificado de latencia. Tampoco la solución
estática P.3 demuestra que la carga se relaje instantáneamente durante un
transiente: ese régimen temporal se decidirá después.

### Direcciones lineales y soporte de Pauli

Los valores unitarios de estas funciones son derivadas de respuesta. No se
presentan como ocupaciones ya admitidas. Para convertirlas en una perturbación
finita alrededor de $h_0(E)=\tanh(E/2T_b)$, se requiere
$|h_0+\epsilon(\delta h_L\pm\delta h_T)|\le1$.
Para $\epsilon\ge0$, cada pendiente positiva $s$ limita
$\epsilon\le(1-h_0)/s$ y cada pendiente negativa limita
$\epsilon\le(1+h_0)/(-s)$. El mínimo sobre ambos signos, nodos y energías
es la cota admisible. A baja temperatura y alta energía puede ser diminuta;
si la representación redondea $h_0$ a uno, no se inventa una holgura positiva.

## Reconstrucción conjunta de fuerza y corriente

El mismo cálculo permite un control más fuerte que inspeccionar la DOS.
En las unidades de la energía térmica de grafo, los incrementos respecto
de la distribución de referencia son

\[
\delta G_{d,i}=i\,m_i\int_0^\infty\delta K_{12,i}(E)\,dE,
\qquad
\delta\overline I_{ij}=2\int_0^\infty\delta j_{T,ij}(E)\,dE.
\tag{P.4}
\]

El primero es una fuerza generalizada compleja integrada del nodo, escrita
con la misma convención de signo y escala que el gradiente térmico. No es
una fuerza por unidad de volumen. Para una distribución no térmica espacial
arbitraria todavía no se ha demostrado que provenga de una única energía
escalar: esa integrabilidad no se deduce de P.5. Al definir
$C=(f-\widetilde f^*)/(2i)$ y $S=(f+\widetilde f^*)/2$ se obtiene
$\delta G_d=-2m\int(\delta h_L C+i\delta h_T S)dE$.
Para una distribución simétrica con $\delta h_L=-2\delta p$ y $h_T=0$
recupera el peso de B.49: $4m\int C\delta p\,dE$.

Se pueden sumar estos incrementos a la fuerza y la corriente térmicas
ya calculadas. Así $\delta h=0$ devuelve exactamente esa referencia y no
se suma otra energía de vacío. Esa recuperación por sustracción de referencia
no pretende que la cuadratura retardada dispersa ya reproduzca una suma
Matsubara infinita.

Para distribuciones reales, $K_{21}=-K_{12}^*$. La conversión local de P.1
y P.4 implican la identidad discreta

\[
\boxed{
\sum_j\delta\overline I_{ij}
+\operatorname{Im}(d_i^*\delta G_{d,i})
=2\int_0^\infty\delta r_{T,i}(E)\,dE.}
\tag{P.5}
\]

Es una relación entre corriente, torque del condensado y residuo de carga.
Se conserva exactamente con una cuadratura finita común, antes de reclamar
su convergencia energética. Las reacciones de los nodos de contacto se
registran por separado. P.3 anula el segundo miembro en los nodos libres;
esto verifica que fuerza y corriente del incremento describen la misma
conversión de pares. No equivale a un balance total de energía durante el
movimiento del gap.

## Qué decisión requiere información física adicional

Este diagnóstico puede implementarse dentro del alcance ya autorizado sin
elegir una cascada, un ancho del fotón ni una tasa de relajación nueva. El
resultado orientará una comparación entre la respuesta espectral de carga
y una reducción con potencial normal y fase como la memoria. Antes de
elegir una variable dinámica adicional se compararán los momentos de fuerza,
corriente y energía relevantes, y se justificará la escala temporal de la
reducción. Un residuo espectral omitido no es por sí mismo una razón para
abandonar un cierre mesoscópico útil. Tampoco se identifica DOS nula con
conductividad transversal nula: $\mathcal D_T$ puede permanecer finita
dentro del gap. El balance de corriente total no representa todas las formas
de $h_T(E)$, pero puede ser el momento relevante de una aproximación que
se evalúe explícitamente.

Siguen pendientes el trabajo espectral al mover el gap, la deriva y transporte
de las ocupaciones y su balance con la disipación del condensado. Las
colisiones desarrolladas en B no se descartan, pero sus factores espectrales
deben usar este mismo propagador. Si la respuesta calculada exige una nueva
variable dinámica de carga o una aproximación temporal distinta, se expondrá
esa elección antes de cambiar el modelo de detector. El circuito de la
memoria y el horizonte de gatillo más margen permanecen como contratos vigentes.

## Comparación mesoscópica y cuadratura siguiente

Una comparación adicional económica restringe la respuesta de carga a la
dirección de un pequeño desplazamiento electroquímico:
$h_T(E,i)=\chi(E)v_i$, con
$\chi(E)=\partial_E\tanh[E/(2T_b)]$ en las unidades adimensionales usadas.
Se resuelve la ecuación de carga integrada en energía dentro de ese subespacio,
y se comparan sus momentos con la respuesta espectral libre. $v$ es una
coordenada electroquímica; su signo respecto del potencial electrostático de
la memoria debe fijarse por la convención de corriente. Esta proyección no
se identifica automáticamente con toda su ley óhmica o su dinámica de fase.

La malla de doce energías basta para las identidades algebraicas, pero es
gruesa para esa comparación: a $T_b=0,9$ K y con el gap de referencia actual,
$T_b/\Delta_{\rm ref}=0,05898725$. El 97,15 % de la integral de $\chi$ está
por debajo de $E/\Delta_{\rm ref}=0,25$. La regla trapezoidal sobre la malla
actual da $\int\chi\,dE=1,180225$, frente al valor exacto prácticamente uno.
No se normaliza ese error para ocultarlo.

La primera comparación ya calculada con esos doce puntos muestra, para la
dirección angular, diferencias de aproximadamente 20 % en el campo de
corriente integrado, 3,1 % en la fuerza compleja y 0,42 % en el flujo de
energía. La componente de amplitud de la fuerza cambia muy poco; la
diferencia se concentra en su respuesta de fase. Por tanto, una buena
comparación del calentamiento no decide todavía la respuesta de corriente
y fase relevante para Vout. Con la cuadratura térmica todavía gruesa, estos
datos no justifican aceptar o rechazar la reducción ni preguntar ya por una
nueva variable dinámica. La recomendación es resolver los momentos antes
de esa elección física, conservando el resultado favorable del núcleo.

Se recomienda una malla base de 31 puntos, que contiene los doce anteriores:

```text
0, .025, .05, .075, .1, .15, .2, .25, .375, .5, .625, .75,
.85, .9, .95, .975, .99, 1, 1.01, 1.025, 1.05, 1.1, 1.15,
1.25, 1.5, 1.75, 2, 2.5, 3, 4, 5
```

La malla fina añade los puntos medios de los intervalos contenidos en
$[0,0.5]$ o $[0.85,1.15]$, quedando en 50 puntos. Se concentra así el
trabajo en la ventana térmica y la estructura de coherencia, sin refinar
uniformemente la cola alta. La misma integral analítica de $\chi$ resulta
1,016162 en la base y 1,004243 en la fina. Estos números motivan la selección
de nodos; no constituyen una tolerancia universal de admisión del modelo.

La comparación relevante es entre momentos de fuerza, corriente y flujo de
energía: base frente a fina, $\eta/\Delta_{\rm ref}=0,02$ frente a 0,01,
y respuesta espectral frente a la proyección electroquímica. Una campaña
económica usa el núcleo radial 65 en ambas $\eta$ sobre la malla fina,
el radial 129 en $\eta=0,01$ sobre la fina y el asimétrico 65 en
$\eta=0,01$ sobre la base: 181 consultas espectrales antes de cualquier
reutilización con hashes verificados. Las comparaciones base/fina de los
casos finos salen de sus mismos datos anidados. Su aceptación se relacionará
con el tamaño del efecto mesoscópico comparado, no con exigir precisión
puntual arbitraria en todos los picos de la DOS.

### Interpretación predefinida de la comparación

Para cada momento integrado $M$ se registran tres cantidades en la misma
norma física: $D=\|M_{\rm potencial}-M_{\rm espectral}\|$,
$C=\|M_{h_T=0}-M_{\rm espectral}\|$ y una envolvente de variaciones
observadas $U$ construida con los cambios base/fina, de $\eta$ y de malla
disponibles. La suma de esos cambios es una referencia conservadora de la
sensibilidad observada, no una cota matemática del error verdadero.

- Si $D\le U$, la diferencia entre representaciones no queda resuelta;
  no se declara que sean exactamente equivalentes.
- Si $C\le U$, el efecto de omitir la respuesta de carga tampoco se resuelve.
  No se usa un cociente $D/C$ grande por denominador casi nulo para rechazar
  el control radial.
- Si las diferencias se resuelven, $D/C$ indica cuánto de la corrección
  ausente recupera la proyección. Si $D\ge C$ por un margen mayor que la
  variación observada, o se invierte de forma robusta el signo de un momento
  no nulo, esa proyección no mejora el ensayo y debe revisarse.
- Si la proyección reduce sistemáticamente el error del cierre escalar,
  mantiene los signos relevantes y sus momentos se estabilizan, puede
  continuar como candidata al siguiente transiente mesoscópico débil, con
  el error residual declarado. No se incorpora una variable dinámica de
  carga sólo porque el residuo espectral congelado sea distinto de cero.

El momento de energía usa $\int E j_L\,dE$, no sólo $\int j_L\,dE$.
Al comparar mallas se interpola la densidad de fuerza $\delta G_d/m$;
para corrientes se conserva la medida dual o se usa su densidad, porque el
valor integrado de una arista cambia cuando se reduce su ancho dual.
Una diferencia entre normas totales no se etiqueta como norma del error de
todo el campo. Estos criterios permiten decidir el siguiente ensayo sin
establecer una tolerancia microscópica universal ni adelantar la aceptación
del detector real.
