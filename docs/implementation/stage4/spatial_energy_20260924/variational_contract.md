# Contrato variacional térmico para la energía espectral espacial

Fecha: 24 de septiembre de 2026. Revisión analítica independiente; propuesta de implementación experimental, sin sustitución del solver de producción.

## 1. Decisión y alcance

La siguiente implementación debe conservar los gradientes de los campos espectrales como variables auxiliares y eliminarlos por estacionariedad. Así, fuerza del condensado y corriente proceden de una misma energía espacial. No se utiliza un momento regularizado $q_\delta$ ni se añade un coeficiente constante $K_0$ a una energía uniforme.

Esta construcción es el funcional térmico de Usadel dentro de las hipótesis de A.2: singlete isotrópico, límite difusivo, acoplamiento débil, simetría electrón–hueco y DOS normal constante. «Térmico» significa una única temperatura fijada al evaluar el funcional. No demuestra un cierre fuera del equilibrio para una distribución electrónica arbitraria, no determina una movilidad temporal y no convierte automáticamente energía libre en energía interna.

La referencia continua procede de A.6, antes de su especialización uniforme A.9–A.12, y de la ecuación (20) de Virtanen, Vargunin y Silaev, *Physical Review B* **101**, 094507 (2020), [DOI](https://doi.org/10.1103/PhysRevB.101.094507). La discretización de la sección 6 es una elección consistente de esta implementación, no una identidad exacta del artículo para cualquier malla finita.

## 2. Variables, unidades y funcional continuo

Definimos

$$
z=\frac{\Delta}{k_BT_c},\quad t=\frac{T}{T_c}>0,\quad
\varepsilon_n=2\pi t(n+\tfrac12),\quad
\ell_0=\sqrt{\frac{\hbar D}{2k_BT_c}},\quad
\boldsymbol x=\frac{\boldsymbol r}{\ell_0},\quad
\boldsymbol a=\frac{2e\ell_0}{\hbar}\boldsymbol A,\quad
\mathcal D=\boldsymbol\nabla_{\boldsymbol x}-i\boldsymbol a.
$$

$N_0$ es la DOS normal por espín y $\sigma_n=2e^2N_0D$. No debe confundirse la longitud $\ell_0$ con una longitud ajustada del perfil ni con el coeficiente GL dependiente de temperatura.

Para cada frecuencia positiva se introducen un campo complejo $f_n$ y uno real $g_n$, con

$$
|f_n|^2+g_n^2=1,\qquad g_n>0.
$$

Aquí $f_n$ es el propagador anómalo de Matsubara; no es una ocupación electrónica ni una densidad de energía. Es posible escribir $f_n=\sin\vartheta_n e^{i\chi_n}$ y $g_n=\cos\vartheta_n$, pero no se debe imponer de antemano que la fase espectral $\chi_n$ coincida espacialmente con la fase de $z$. Esa coincidencia sí pertenece al ensayo radial simétrico considerado abajo.

La diferencia de energía libre superconductora–normal adimensional es

$$
\boxed{\displaystyle
\mathcal F[z,\{f_n,g_n\},\boldsymbol a]
=\int_\Omega\!\left\{
|z|^2\ln t+2\pi t\sum_{n\ge0}
\left[
\frac{|z|^2}{\varepsilon_n}+2\varepsilon_n(1-g_n)
-2\operatorname{Re}(z^*f_n)
+|\mathcal D f_n|^2+|\boldsymbol\nabla g_n|^2
\right]\right\}d^d x .}
$$

En ángulos, el término espacial es

$$
|\mathcal D f_n|^2+|\boldsymbol\nabla g_n|^2
=|\boldsymbol\nabla\vartheta_n|^2
+\sin^2\vartheta_n|\boldsymbol\nabla\chi_n-\boldsymbol a|^2.
$$

Para película de espesor constante $d_{\rm film}$ y simulación 2D, la energía física es $F=E_*\mathcal F$, con $E_*=N_0(k_BT_c)^2d_{\rm film}\ell_0^2$. En 3D, $E_*=N_0(k_BT_c)^2\ell_0^3$. Otras normalizaciones históricas requieren su conversión explícita.

La energía del metal normal omitida es $f_{\rm normal}=-\pi^2N_0(k_BT)^2/3$. No modifica las derivadas de condensado y corriente a temperatura fija, pero debe restituirse para calcular energía térmica o calor. La relación $U=F-T\partial_TF$ se evalúa a $\Delta$ física y $\boldsymbol A$ fijados, con bordes y tratamiento de cola consistentes.

## 3. Estacionariedad espectral, fuerza y corriente

Se resuelven los campos espectrales para un $z$ prescrito, sin exigir todavía autoconsistencia del gap. El teorema de la envolvente funcional permite derivar la energía reducida $\overline{\mathcal F}[z,\boldsymbol a]$ manteniendo congelados los campos espectrales estacionarios en la derivada explícita. Esto requiere que las variaciones espectrales permitidas tengan derivada nula, incluyendo los términos de borde correspondientes.

En una región donde se puedan usar ángulos, las dos ecuaciones espectrales son

$$
\begin{aligned}
0={}&\nabla^2\vartheta_n-\sin\vartheta_n\cos\vartheta_n
|\nabla\chi_n-\boldsymbol a|^2-\varepsilon_n\sin\vartheta_n
+\operatorname{Re}(ze^{-i\chi_n})\cos\vartheta_n,\\
\nabla\cdot\left[\sin^2\vartheta_n(\nabla\chi_n-\boldsymbol a)\right]
={}&\operatorname{Im}(z^*f_n).
\end{aligned}
$$

La segunda muestra por qué no debe congelarse la fase espectral a la fase local del gap para un perfil arbitrario. En variables cartesianas, la misma estacionariedad se expresa como
$\mathcal D^2f_n-(f_n/g_n)\nabla^2g_n-\varepsilon_nf_n/g_n+z=0$, junto a la normalización. El código puede usar $u_n$ para evitar la singularidad angular de $f_n=0$.

Usamos la convención real para una variable compleja:

$$
\delta\overline{\mathcal F}
=\operatorname{Re}\int G_z^*\delta z\,d^dx
-\int\overline{\boldsymbol j}\cdot\delta\boldsymbol a\,d^dx
+\text{trabajo de los bordes}.
$$

Entonces

$$
\boxed{G_z=2z\ln t+4\pi t\sum_{n\ge0}
\left(\frac z{\varepsilon_n}-f_n\right)},\qquad
\boxed{\overline{\boldsymbol j}=4\pi t\sum_{n\ge0}
\operatorname{Im}(f_n^*\mathcal D f_n)}.
$$

La derivada de Wirtinger respecto de $z^*$ es $G_z/2$. Esta distinción evita un factor dos al conectar fuerzas reales, complejas y amplitudes. La corriente física volumétrica es

$$
\boldsymbol j=\frac{2e}{\hbar}N_0(k_BT_c)^2\ell_0\,
\overline{\boldsymbol j}.
$$

En el límite uniforme recupera exactamente A.5. El gradiente del condensado ya está contenido en la solución espacial de $f_n,g_n$: añadir otro término $-2K_0\mathcal D^2z$ a esta fuerza lo contaría nuevamente.

El trabajo espectral de borde, escrito en ángulos, es

$$
4\pi t\sum_n\int_{\partial\Omega}
\left[
\partial_\nu\vartheta_n\,\delta\vartheta_n+
\sin^2\vartheta_n\,\boldsymbol\nu\cdot
(\boldsymbol\nabla\chi_n-\boldsymbol a)\,\delta\chi_n
\right]dS.
$$

Un borde de Dirichlet espectral fijado elimina estas variaciones. Recalcular su valor a partir de $z$ en cada perturbación cambia el problema variacional y añade el trabajo de ese contacto. Un espectro resuelto solo aproximadamente deja un residuo de envolvente; debe cuantificarse junto a la derivada direccional, sin identificar automáticamente la tolerancia interna del optimizador con el error físico de la fuerza.

## 4. Referencia radial reutilizable y el pequeño agujero del origen

Para un vórtice prescrito de devanado uno, $z=d(r)e^{i\phi}$, $f_n=\sin\vartheta_n(r)e^{i\phi}$, $\boldsymbol a=0$, la acción radial es

$$
\begin{aligned}
\mathcal F_{\rm rad}={}&2\pi\int_{r_0}^Rr\,dr
\left\{d^2\ln t+2\pi t\sum_n\left[
\frac{d^2}{\varepsilon_n}+2\varepsilon_n(1-\cos\vartheta_n)
-2d\sin\vartheta_n+(\vartheta_n')^2
+\frac{\sin^2\vartheta_n}{r^2}\right]\right\}
+B_{\rm inner}.
\end{aligned}
$$

Su ecuación espectral interior es la ya utilizada por la referencia:

$$
\vartheta_n''+\frac{\vartheta_n'}r
-\frac{\sin\vartheta_n\cos\vartheta_n}{r^2}
-\varepsilon_n\sin\vartheta_n+d\cos\vartheta_n=0.
$$

La fuerza de amplitud y la corriente azimutal son

$$
X_d=2d\ln t+4\pi t\sum_n
\left(\frac d{\varepsilon_n}-\sin\vartheta_n\right),\qquad
\overline j_\phi=\frac{4\pi t}{r}\sum_n\sin^2\vartheta_n.
$$

El término espectral de borde de la integral radial es

$$
8\pi^2t\sum_n[\,r\vartheta_n'\delta\vartheta_n\,]_{r_0}^R.
$$

**La condición de Robin usada en la referencia anterior necesita una acción de borde si se quiere probar la envolvente de una energía.** En el origen verdadero, regularidad significa $\vartheta_n\sim c_nr$ y $\vartheta_n(0)=0$, no $\vartheta_n'(0)=0$. En un dominio truncado, la aproximación $\vartheta_n'(r_0)=\vartheta_n(r_0)/r_0$ se obtiene variacionalmente añadiendo

$$
\boxed{B_{\rm inner}=4\pi^2t\sum_n\vartheta_n(r_0)^2.}
$$

Este es el costo espacial principal del disco omitido: se integra el perfil regular lineal $\vartheta_n(r)=\vartheta_n(r_0)r/r_0$ desde cero hasta $r_0$. Su variación cancela el flujo interior y genera precisamente Robin. Las correcciones geométricas y de potencial del pequeño disco son de orden $r_0^4$ para perfiles regulares. Robin aproxima, pero no es idéntica, a la solución regular completa: si $\vartheta=cr+br^3+\cdots$, entonces $\vartheta'-\vartheta/r=2br^2+\cdots$.

Una perturbación compacta de $d$ no garantiza $\delta\vartheta_n(r_0)=0$ después de resolver el espectro. Por tanto, no basta para justificar omitir $B_{\rm inner}$. Si se añade una cola de Matsubara también debe contemplarse la cola de este costo interior o declarar su incertidumbre. La nueva malla cartesiana que incluye el nodo del origen no tiene este agujero y no necesita este término radial artificial.

En el borde exterior, mantener fijado el ángulo de contacto hace desaparecer el trabajo. Si se vuelve a calcular el ángulo uniforme con circulación local usando $d(R)$, aparece $8\pi^2tR\vartheta_n'(R)\delta\vartheta_n(R)$. Los antiguos BVP siguen siendo referencias válidas de fuerza interior; sus resultados no acreditaron una identidad de energía reducida con ese contacto móvil.

El oráculo radial más económico conserva estos BVP y calcula su acción y derivadas con trazas espectrales fijas. No requiere buscar un vórtice autoconsistente, barreras, fotón ni un transiente. Las soluciones de frecuencias distintas siguen siendo independientes. Las comparaciones de corte y radio exterior ya disponibles orientan el costo; no sustituyen esta nueva contabilidad variacional.

## 5. Cola de Matsubara: derivar energía antes que corregir fuerza

Para campos suaves, lejos de una capa de borde incompatible con su expansión de alta frecuencia,

$$
f_n=\frac z{\varepsilon_n}
+\frac{\mathcal D^2z}{\varepsilon_n^2}
+\frac{\mathcal D^4z-|z|^2z/2}{\varepsilon_n^3}
+O(\varepsilon_n^{-4}).
$$

En la geometría radial, $\mathcal D^2z=e^{i\phi}L_1d$, con $L_1d=d''+d'/r-d/r^2$. No se sustituye por el Laplaciano de una amplitud sin fase. Definiendo

$$
S_p(N)=\sum_{n=N}^{\infty}\varepsilon_n^{-p}
=\frac{\zeta(p,N+1/2)}{(2\pi t)^p},
$$

la contribución omitida a la fuerza interior empieza por

$$
G_{z,\rm tail}=-4\pi tS_2\mathcal D^2z
+4\pi tS_3\left(-\mathcal D^4z+\frac{|z|^2z}{2}\right)+\cdots.
$$

Para variaciones compactas o con el trabajo de borde conservado, una acción asintótica que genera esos términos es

$$
\mathcal F_{\rm tail}
=2\pi t\left[
S_2\int|\mathcal D z|^2d^dx
+S_3\int\left(\frac{|z|^4}{4}-|\mathcal D^2z|^2\right)d^dx
\right]+\cdots.
$$

Antes de integrar por partes, el término de orden $\varepsilon_n^{-3}$ es
$|\mathcal D^2z|^2+2\operatorname{Re}[(\mathcal D z)^*\cdot\mathcal D(\mathcal D^2z)]+|z|^4/4$.
Por ello la forma integrada anterior lleva además
$2\operatorname{Re}\int_{\partial\Omega}(\boldsymbol\nu\cdot\mathcal D z)^*\mathcal D^2z\,dS$
dentro del coeficiente $2\pi tS_3$. No se puede descartar ese término y después usar la fórmula para bordes que varían.

El coeficiente positivo $2\pi tS_2$ es una cola calculada que tiende a cero al aumentar $N$; no es un parche $K_0$ independiente del corte. En cambio, el término negativo de cuarta derivada es solo asintótico: **no constituye un operador espacial estable para números de onda arbitrariamente grandes**. Su expansión requiere que las escalas espaciales relevantes, además de $|z|$, sean pequeñas frente a $\varepsilon_N$. No se debe instalar como una energía global sin ese control.

La primera implementación puede elegir una suma explícita finita y comparar cortes, manteniendo exactas entre sí su energía, fuerza y corriente discretas. Si añade una cola principal, debe añadir la acción $2\pi tS_2\int|\mathcal D z|^2$ y derivarla también para corriente y bordes. La cola homogénea BCS se puede evaluar con el sumando estable

$$
\frac{|z|^2}{\varepsilon}+2\varepsilon
-2\sqrt{\varepsilon^2+|z|^2}
=\frac{|z|^4}{\varepsilon(\varepsilon+\sqrt{\varepsilon^2+|z|^2})^2}.
$$

Añadir esa cola local no reconstruye los términos espaciales omitidos. Las diferencias de cortes deben darse para las magnitudes observadas, no presentarse como una cota rigurosa universal.

### Una alternativa positiva para la cola cuadrática

En el régimen lineal de amplitud, sea $L=-\mathcal D^2$ autoadjunto y positivo con condiciones homogéneas adecuadas. Para cada frecuencia,

$$
\mathcal F_n^{(2)}/(2\pi t)
=\langle z,z\rangle/\varepsilon_n
+\langle f_n,(\varepsilon_n+L)f_n\rangle
-2\operatorname{Re}\langle z,f_n\rangle.
$$

Eliminar $f_n=(\varepsilon_n+L)^{-1}z$ produce un operador positivo, y su suma de cola exacta es

$$
\mathcal F_{\rm tail}^{(2)}
=\left\langle z,
\left[\psi\!\left(N+\tfrac12+\frac L{2\pi t}\right)
-\psi(N+\tfrac12)I\right]z\right\rangle.
$$

$\psi$ es la función digamma. Esta reagrupación conserva todas las escalas espaciales del problema **cuadrático**, sin introducir una cuarta derivada negativa truncada. No es la solución no lineal completa. En una malla con masa $M$ y rigidez $K$, se aplica a $\widetilde L=M^{-1/2}KM^{-1/2}$ y $\widetilde z=M^{1/2}z$.

Con contactos espectrales inhomogéneos se requiere un levantamiento de sus trazas: en los nodos interiores,
$(\varepsilon M_I+K_{II})f_I=M_Iz_I-K_{IB}f_B$.
Una traza $f_B$ dependiente de frecuencia impide aplicar sin más la fórmula homogénea de digamma. Además, su corriente debe derivarse del mismo operador dependiente de enlaces; sustituirla por una corriente GL independiente rompería el contrato. Por sencillez, esta optimización de cola puede esperar al control finito explícito.

## 6. Acción sobre la malla Delaunay–Voronoi

Sean $V_i>0$ las áreas duales adimensionales y $w_{ij}$ las conductancias dual/primal. Cada enlace no orientado se cuenta una sola vez en la energía. Para la orientación elegida $i\to j$, definimos

$$
\alpha_{ij}=\int_i^j\boldsymbol a\cdot d\boldsymbol l,
\qquad R_{ij}=R_z(-\alpha_{ij}).
$$

$R_z$ rota las dos componentes anómalas y deja la tercera invariante. Se usa la parametrización cartesiana regular

$$
u_{ni}=f_{ni}/g_{ni},\qquad
\boldsymbol n_{ni}=\frac{(\operatorname{Re}u_{ni},\operatorname{Im}u_{ni},1)}
{\sqrt{1+|u_{ni}|^2}}=(\operatorname{Re}f_{ni},\operatorname{Im}f_{ni},g_{ni}).
$$

La acción finita sin cola es

$$
\begin{aligned}
\mathcal F_h={}&\sum_iV_i|z_i|^2\ln t\\
&+2\pi t\sum_{n=0}^{N-1}\left\{
\sum_iV_i\left[\frac{|z_i|^2}{\varepsilon_n}
+2\varepsilon_n(1-g_{ni})-2\operatorname{Re}(z_i^*f_{ni})\right]
+\sum_{\{i,j\}}w_{ij}|\boldsymbol n_{ni}-R_{ij}\boldsymbol n_{nj}|^2
\right\}.
\end{aligned}
$$

No aparece un factor $1/2$ adicional en esta suma de enlaces contados una sola vez, usando las conductancias convencionales que aproximan $\int|\nabla v|^2$. Esa convención debe comprobarse con la geometría real. La energía presupone pesos compatibles con el volumen finito ortogonal; no justifica recortar silenciosamente pesos negativos ni masas de borde. Celdas duales de borde y condiciones de contacto requieren su tratamiento geométrico explícito.

Esta acción de diferencias de vectores unitarios converge al término espacial continuo al refinar una malla adecuada. **No equivale exactamente en una malla finita** a discretizar primero las ecuaciones angulares ni a cualquier operador anterior llamado «Usadel». La independencia del parámetro $\delta$ proviene de la variable cartesiana, no de una equivalencia a resolución finita.

### Derivadas analíticas

Sean $s_i=\sqrt{1+|u_i|^2}$, $P=(\boldsymbol e_x,\boldsymbol e_y)$ y

$$
J_i=\frac{\partial\boldsymbol n_i}{\partial(\operatorname{Re}u_i,\operatorname{Im}u_i)}
=\frac{(I-\boldsymbol n_i\boldsymbol n_i^T)P}{s_i}.
$$

Para cada modo, con $\boldsymbol b_i=(\operatorname{Re}z_i,\operatorname{Im}z_i,\varepsilon_n)$,

$$
\boldsymbol G_{n_i}=4\pi t\left[
\sum_jw_{ij}(\boldsymbol n_i-R_{ij}\boldsymbol n_j)-V_i\boldsymbol b_i
\right],\qquad
\boldsymbol G_{u_i}=J_i^T\boldsymbol G_{n_i}.
$$

Para la orientación inversa se usa $R_{ji}=R_{ij}^{T}$. La condición espectral es $\boldsymbol G_{u_i}=0$ en nodos libres. Cuando $g_i\to0$, $u_i$ puede ser grande y su Jacobiano pequeño; por ello también se debe inspeccionar el residuo tangente físico $(I-\boldsymbol n_i\boldsymbol n_i^T)\boldsymbol G_{n_i}$, no solamente el gradiente en coordenadas $u_i$.

La derivada del gap, en la misma convención compleja real de la sección 3, es

$$
G_{z_i}=2V_i z_i\ln t+4\pi tV_i\sum_{n=0}^{N-1}
\left(\frac{z_i}{\varepsilon_n}-f_{ni}\right).
$$

Es una fuerza integrada del nodo. La fuerza por volumen dual es $G_{z_i}/V_i$; comparar normas sin distinguir ambas magnitudes genera una dependencia espuria de la cantidad de nodos.

La corriente integrada orientada del enlace es

$$
\boxed{\overline I_{ij}=-\frac{\partial\mathcal F_h}{\partial\alpha_{ij}}
=4\pi t\,w_{ij}\sum_n
\operatorname{Im}\left(f_{ni}^*e^{-i\alpha_{ij}}f_{nj}\right).}
$$

Es positiva de $i$ a $j$ para una fase que aumenta en esa dirección y $\boldsymbol a=0$. Para película 2D, $I_{ij}^{\rm phys}=(2e/\hbar)E_*\overline I_{ij}$. Una cola añadida a la acción aporta también sus derivadas a $G_z$ y a $\overline I$.

### Calibre y reacciones de los contactos

La transformación $z_i\mapsto e^{i\lambda_i}z_i$, $u_{ni}\mapsto e^{i\lambda_i}u_{ni}$ y $\alpha_{ij}\mapsto\alpha_{ij}+\lambda_j-\lambda_i$ deja la acción invariante. En un ensayo de calibre se transforman también las trazas espectrales prescritas.

Si la incidencia $B$ tiene $+1$ en la cola y $-1$ en la cabeza, la identidad local, incluso antes de resolver el espectro, es

$$
\operatorname{Im}(z_i^*G_{z_i})
+\sum_n\operatorname{Im}(u_{ni}^*G_{u_{ni}})
+(B\overline I)_i=0.
$$

La segunda contribución desaparece en los nodos espectrales libres estacionarios. En contactos fijados representa una reacción y no debe descartarse al reclamar conservación de corriente. La condición de gap estacionario tampoco está impuesta en un perfil prescrito; su torque puede equilibrar una divergencia de corriente. Esto separa una identidad de calibre de la afirmación física más fuerte de continuidad en un estado totalmente estacionario.

## 7. Pruebas que deciden el siguiente paso

1. **Normal y homogéneo.** $z=u=0$ tiene energía de condensación y corriente nulas. Una solución uniforme sin campo reproduce A.3–A.5, con idéntica convención de corte y de cola.
2. **Límite lineal independiente.** Para amplitud pequeña, el problema espectral satisface $(\varepsilon M+K)f=Mz$ con sus contactos. Una solución analítica de Fourier apropiada o una solución independiente del operador lineal prueba el límite del sistema no lineal y sus factores geométricos.
3. **Derivadas antes de eliminar.** Diferencias direccionales de la acción para $z$, $u$ y enlaces prueban las derivadas explícitas. No precisan resolver cientos de transientes.
4. **Envolvente después de eliminar.** Volver a resolver el espectro para $z\pm h\eta$ o enlaces perturbados, con las mismas trazas espectrales físicas, y comparar la derivada de la energía reducida. Si se usa el oráculo radial con agujero, incluir $B_{\rm inner}$.
5. **Calibre con y sin estacionariedad.** Invariancia de energía, covarianza de fuerza e identidad local completa, conservando el torque espectral y las reacciones de borde donde corresponda.
6. **Núcleo con gap nulo.** Atravesar $z=0$ en coordenadas cartesianas debe dar campos y derivadas finitos sin calcular $\arg z$, dividir por $|z|$ ni escoger un $\delta$.
7. **Corte y dominio.** Comparar energía, fuerza y corriente con cortes distintos y con la referencia radial en su región interior. Si hay cola, comprobar sus tres derivadas y separar su aproximación de los errores del optimizador y de malla.
8. **Convergencia espacial útil.** Después de las identidades anteriores, comparar la acción de malla contra el oráculo continuo para el mismo perfil y contactos. La diferencia de representaciones a malla finita debe medirse, no declararse inexistente.

No se establece una tolerancia universal nueva. Cada aceptación debe relacionar error de estacionariedad, derivada, corte y malla con el cambio físico observado y con el costo. El primer resultado útil es una energía espacial común que produzca fuerzas y corrientes compatibles y que se acerque a la referencia independiente del núcleo. Todavía no valida barreras de activación, dinámica con fotón, un vórtice estable ni las tasas disipativas del dispositivo real.

## 8. Próxima acción concreta

Implementar primero la acción finita cartesiana con campos auxiliares, trazas de contacto fijadas y derivadas analíticas. Reutilizar el BVP radial como dato independiente de inicialización y comparación interior, sin atribuirle retrospectivamente una validación de energía que no realizó. Resolver las identidades y el corte en estos ensayos estáticos acotados; solo después decidir la aproximación de cola o el método de aceleración que haga falta. El transporte cinético no térmico y la dinámica temporal mantienen contratos separados hasta derivar su conexión con esta energía.
