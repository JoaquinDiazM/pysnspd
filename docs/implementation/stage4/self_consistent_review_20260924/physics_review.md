# Decisión física después del núcleo térmico autoconsistente

Revisión del 24 de septiembre de 2026. Esta nota separa el resultado estático
ya obtenido del siguiente cambio implementable. No añade una condición de
aceptación retrospectiva ni solicita otra campaña de refinamiento térmico.

La relajación autoconsistente prueba que el funcional espacial puede producir
su propio núcleo, sin fijar su perfil interior y sin dividir por la amplitud
del condensado. Es un avance suficiente para continuar con esta representación.
El criterio medido pertenece a los nodos libres de la región de interés; no
certifica por sí mismo todos los nodos del dominio, estabilidad frente a cualquier
perturbación, nucleación de un vórtice ni dinámica del detector.

La siguiente dificultad útil es **obtener el espectro espacial en energías reales**.
Las ocupaciones electrónicas de B necesitan saber qué estados pueden ocupar.
La suma térmica actual entrega energía, fuerza y corriente de referencia, pero
no entrega directamente esos estados. Reemplazar las ocupaciones por su
temperatura equivalente eliminaría precisamente la información que B decidió
conservar. El ensayo temporal KWT a baño fijo sería legítimo como control
térmico restringido; no resolvería este obstáculo y por ello no se prioriza ahora.

## Transformación analítica que conserva el problema espectral

Se conservan el grafo, sus masas positivas $m_i$, conductancias $c_{ij}$,
enlaces $U_{ij}=e^{-i\alpha_{ij}}$, gap prescrito $d_i$ y unidades del
[contrato espacial](../spatial_energy_20260924/variational_contract.md).
El número complejo $z$ representa la frecuencia: $z>0$ en Matsubara y
$z=\eta-iE$, con $\eta>0$, en la rama retardada. Aquí energía y frecuencia
están divididas por $k_BT_c$.

El mapa térmico $g=(1+|u|^2)^{-1/2}$ no es una función holomorfa: contiene
conjugación. La continuación correcta introduce dos coordenadas complejas
independientes $a_i,b_i$:

\[
Q_i=1+a_ib_i,\qquad f_i=\frac{2a_i}{Q_i},\qquad
\widetilde f_i=\frac{2b_i}{Q_i},\qquad
g_i=\frac{1-a_ib_i}{Q_i}.
\]

Estas expresiones satisfacen exactamente $g_i^2+f_i\widetilde f_i=1$.
La tilde es el nombre de la segunda componente anómala, no una instrucción
de conjugación a igual energía. En el eje térmico se recupera $b_i=a_i^*$;
entonces $\widetilde f_i=f_i^*$, $g_i$ es real y
$a_i=u_i/(1+\sqrt{1+|u_i|^2})$ reproduce el backend existente.
El mapa requiere $Q_i\ne0$; ese denominador se controla como coordenada,
sin recortarlo ni interpretarlo como un nuevo parámetro físico.

Para cada frecuencia, la continuación algebraica de la acción ya discretizada es

\[
\begin{aligned}
\mathcal A_z={}&\sum_i m_i
  \left[2z(1-g_i)-d_i^*f_i-d_i\widetilde f_i\right]\\
&+\sum_{\{i,j\}}c_{ij}\left[
 (f_i-U_{ij}f_j)(\widetilde f_i-U_{ij}^*\widetilde f_j)
 +(g_i-g_j)^2\right].
\end{aligned}
\]

El conjugado de $U$ sólo invierte un enlace de calibre físico real; no conjuga
la incógnita espectral ni $z$. Sobre Matsubara esta acción coincide con la
parte espectral térmica implementada. Fuera de ese eje es una acción compleja:
se resuelven sus ecuaciones estacionarias, **no se minimiza su parte real**.
El contratermino $|d|^2/z$ no modifica la ecuación espectral a $d$ fijo; la
renormalización completa seguirá siendo necesaria al construir fuerzas y energía.

Esta continuación hacia la rama retardada y la separación entre espectro y
distribución son las de [Belzig et al., secciones 2.3.2–2.4.3](https://arxiv.org/html/cond-mat/9812297v2).
Las fórmulas de grafo y el residual siguientes se derivan aquí del funcional
ya adoptado; no se atribuyen literalmente a ese artículo. El uso de coordenadas
de coherencia independientes tiene como referencia adicional
[Eschrig, 2000](https://doi.org/10.1103/PhysRevB.61.9061).

## Residual y Jacobiano directamente implementables

Para cada nodo libre, con los vecinos y orientaciones correspondientes,

\[
H_i^+=m_id_i+\sum_jc_{ij}U_{ij}f_j,\quad
H_i^-=m_id_i^*+\sum_jc_{ij}U_{ij}^*\widetilde f_j,\quad
H_i^z=m_iz+\sum_jc_{ij}g_j.
\]

No hay un término propio adicional en estas sumas. La parte propia de una
arista es constante porque el vector espectral está normalizado. Las ecuaciones son

\[
R_{a,i}=H_i^-a_i^2+2H_i^za_i-H_i^+=0,\qquad
R_{b,i}=H_i^+b_i^2+2H_i^zb_i-H_i^-=0.
\]

En efecto, $\partial_{a_i}\mathcal A_z=2R_{b,i}/Q_i^2$ y
$\partial_{b_i}\mathcal A_z=2R_{a,i}/Q_i^2$. Este cálculo identifica
exactamente el problema estacionario resuelto, sin usar diferencias numéricas
de la acción como definición del modelo.

Las derivadas locales necesarias son

\[
\begin{array}{c|cc}
 &\partial_a&\partial_b\\\hline
f&2/Q^2&-2a^2/Q^2\\
\widetilde f&-2b^2/Q^2&2/Q^2\\
g&-2b/Q^2&-2a/Q^2
\end{array}.
\]

La diagonal del Jacobiano tiene
$\partial_{a_i}R_{a,i}=2H_i^-a_i+2H_i^z$ y
$\partial_{b_i}R_{b,i}=2H_i^+b_i+2H_i^z$; las otras dos derivadas
propias son cero. Para $x=a_j$ o $b_j$ en un vecino,

\[
\begin{aligned}
\partial_x R_{a,i}&=c_{ij}
 [a_i^2U_{ij}^*\partial_x\widetilde f_j
  +2a_i\partial_xg_j-U_{ij}\partial_x f_j],\\
\partial_x R_{b,i}&=c_{ij}
 [b_i^2U_{ij}\partial_x f_j
  +2b_i\partial_xg_j-U_{ij}^*\partial_x\widetilde f_j].
\end{aligned}
\]

Un Newton disperso complejo puede usar estas expresiones con reducción del
residuo y continuación de la solución. La rama normal a alta frecuencia y
la continuación desde $\operatorname{Re}z>0$ identifican la rama causal;
una raíz algebraica cualquiera no basta. La transformación de calibre es
$a_i\mapsto e^{i\beta_i}a_i$, $b_i\mapsto e^{-i\beta_i}b_i$.

## Bordes: continuar el mismo entorno del núcleo

Los nodos espectrales fijados en los ensayos anteriores usan la solución radial
con gap prescrito tipo `tanh`, radio exterior 12 y condición interna de Robin.
No son una tabla homogénea evaluada localmente. Su continuación retardada debe
resolver la misma ecuación radial, ahora compleja,

\[
\theta''+\frac{\theta'}r-
\frac{\sin\theta\cos\theta}{r^2}
-z\sin\theta+d_{\rm tanh}(r)\cos\theta=0.
\]

Se conserva la condición interior del BVP. El ángulo exterior continúa la raíz
de $d(R)\cos\theta-z\sin\theta-q(R)^2\sin\theta\cos\theta=0$.
La interpolación sobre los nodos fijados produce
$a=e^{i\phi}\tan(\theta/2)$ y $b=e^{-i\phi}\tan(\theta/2)$.
No se vuelve a conjugar el ángulo complejo. El entorno radial prescrito continúa
siendo un borde controlado, no un reservorio autoconsistente de todo el detector.

Si este BVP todavía no está disponible, las primeras pruebas pueden usar un
reservorio analítico homogéneo en una geometría declarada aparte. No pueden
presentarse como continuación del núcleo anterior con sus mismos bordes.

## Prueba inmediata y significado del siguiente pase

La implementación pequeña puede decidirse con tres resultados concretos:

1. En un material uniforme, $s=\sqrt{z^2+|d|^2}$ en la rama con
   $\operatorname{Re}s>0$ da $a=d/(z+s)$, $b=d^*/(z+s)$ y $g=z/s$.
   El estado normal debe entregar $g=1,f=\widetilde f=0$; la DOS es
   $\operatorname{Re}g$. Esto fija signos y selección de rama.
2. Para $z>0$ y los mismos bordes, el nuevo solver debe recuperar los campos
   térmicos ya obtenidos. Esta superposición entre coordenadas distintas
   prueba que la implementación conserva la acción anterior.
3. La prueba espacial retardada debe conservar la normalización compleja y
   producir una DOS causal, sin imponer una cota térmica a $|f|$. Se registra
   $\eta$ y su sensibilidad antes de interpretar espectros de baja energía.

Un $\eta$ finito aquí es un desplazamiento numérico del contorno analítico,
no un tiempo de vida medido ni un ensanchamiento material calibrado. El pase
del primer oráculo no acredita el límite $\eta\to0$ de toda la cinta.

Un conjunto pequeño de energías por debajo, alrededor y por encima del gap,
con dos valores de $\eta$ y las tres geometrías térmicas de interés, sirve para
comparar la DOS y su sensibilidad. No es todavía una cuadratura admitida de
conteo o energía. El conteo acumulado $x_i(E)=\int_0^E\rho_i(E')\,dE'$ es
monótono cuando la DOS es positiva; eso no demuestra que cada $x_i$ etiquete
un estado local adiabáticamente conservado cuando el espectro depende de sus
vecinos. Una suma de estados hasta un corte finito necesita además su cola:
no se la fuerza a coincidir con el conteo del metal mediante renormalización.

La posterior reconstrucción de la fuerza debe respetar la fase compleja. En
esta convención, la discontinuidad anómala usa la pareja
$[f^R-(\widetilde f^R)^*]/(2i)$; tomar simplemente
$\operatorname{Im}f^R$ perdería la covarianza cuando el gap tiene fase.
Este señalamiento no reemplaza la derivación de los prefactores, la resta
ultravioleta y el trabajo que necesita el cierre completo.

Antes de la prueba dinámica no térmica queda una unión concreta por derivar:
usar este mismo espectro para la fuerza del gap, la corriente y los coeficientes
de transporte. En la descripción Keldysh, $g^K=g^Rh-hg^A$ incluye la
distribución; sus componentes de energía y carga pueden acoplarse cuando hay
supercorriente y gradientes. Por ello no se supone sin comprobación que el
sector de energía de B sea cerrado en una textura arbitraria.
La referencia es [Belzig et al., ecuaciones 29–43](https://arxiv.org/html/cond-mat/9812297v2).

Tampoco se añade la energía espacial térmica como una corrección independiente
a la energía local de B: ambas contienen parte del mismo fondo condensado.
El nuevo puente debe producir un único trabajo espectral y conservar el balance
al mover el gap. Resolver únicamente la DOS no cumple todavía ese contrato.
Esta limitación identifica el próximo trabajo; no obliga a rehacer la campaña
térmica ni a aumentar sus tolerancias.

Se mantiene la malla y la proyección de continuidad de la memoria como
infraestructura reutilizable, junto con su circuito completo. El método de
avance temporal se decide después de esta unión constitutiva. Los campos
auxiliares espectrales no cambian la política del futuro fotón: mismo dispositivo,
misma precisión y un horizonte hasta el gatillo más su margen, todavía por
definir para el observable acordado.

## Revisión de la implementación preparada

Se cotejaron independientemente las expresiones de
`pysnspd/experimental/retarded_spatial_usadel.py` y el ejecutor
`sandbox/stage4_core/next_retarded_oracle.py` contra el cálculo anterior.
La acción, residual, Jacobiano y derivada de enlace coinciden con él.
El ejecutor continúa el mismo BVP radial y conserva su condición de Robin;
su Jacobiano real representa explícitamente la derivada compleja mediante
las relaciones de Cauchy–Riemann. No sustituye los contactos por una tabla
homogénea. La prueba de superposición Matsubara compara acción, campos y
derivada de corriente antes de continuar cada consulta.

El plan de 72 consultas —tres núcleos con corte 256, doce energías y dos
desplazamientos $\eta$— se considera una exploración útil del espectro.
El corte 128 no se repite: la sensibilidad del gap a ese corte ya se mide
en la campaña térmica. Seis consultas forman el piloto, con límite externo
de 240 segundos y terminación de su grupo de procesos si lo excede.
Los casos y energías comparten el presupuesto de recursos; no hay otro
pool espectral anidado. El pase de estas consultas significará que existe
un oráculo espectral ejecutable para esos campos y valores de $\eta$;
seguirá sin acreditar la cinética completa o el detector fotónico.
