---
title: "E. Cuaderno de aprendizaje del modelo"
subtitle: "Conceptos, problemas resueltos y evaluación · Documento 5 de 5"
date: "8 de septiembre de 2026 · Revisión 0.3"
lang: es
---

\setlength{\parskip}{3.5pt plus 1pt}

# E.0. Cómo trabajar con este cuaderno

Este cuaderno conecta herramientas ya familiares de mecánica lagrangiana, modos normales y dispositivos electrónicos con las ideas que utilizan A–D. El objetivo es poder reconstruir un argumento, reconocer sus hipótesis y detectar una interpretación incorrecta. No se evalúa memorizar nombres ni reproducir frases del texto.

Cada tema contiene una introducción con figura, un problema resuelto, una actividad evaluada y una conexión concreta con A–D. Las figuras son construcciones pedagógicas reproducibles: no representan datos de NbN salvo que se indicara expresamente. Se pueden consultar las explicaciones al responder. Conviene escribir el razonamiento, las unidades y cualquier incertidumbre en el Markdown, debajo del identificador correspondiente.

| ID estable | Tema | Estado | Evidencia para esta revisión |
|:--|:--|:--|:--|
| E01 | Campos, partículas y Modelo Estándar | **ABIERTO** | Tema nuevo; respuestas pendientes |
| E02 | Nambu y espín: dos índices distintos | **ABIERTO** | Tema nuevo; respuestas pendientes |
| E03 | Variaciones, funcionales y envolvente | **ABIERTO** | Tema nuevo; respuestas pendientes |
| E04 | Electrones y huecos: del diodo al superconductor | **ABIERTO** | Tema nuevo; respuestas pendientes |
| E05 | De modos normales a DOS y DFPT | **ABIERTO** | Tema nuevo; respuestas pendientes |

**Criterio de cierre.** Cada tema vale 10 puntos y requiere al menos 8, además de superar la condición esencial indicada. Un error aritmético pequeño con método correcto recibe crédito parcial; una confusión central mantiene el tema abierto aunque otras respuestas sean buenas. “CERRADO / APROBADO” se asignará únicamente después de revisar respuestas. Si se incorpora contenido nuevo o una respuesta posterior revela una dificultad fundamental, el tema vuelve a abierto, conservando el historial. Abierto significa pendiente de aprendizaje o evaluación, no reprobado.

En la siguiente iteración se conservarán las respuestas, se añadirá la pauta completa de los ejercicios entregados y se registrará qué evidencia permite cerrar o qué ejercicio breve conviene repetir. Esta versión contiene soluciones de **ejemplos distintos**; no incluye por anticipado la solución de las actividades evaluadas.

Orden sugerido: E01 para ubicar el lenguaje; E04 si se prefiere comenzar desde semiconductores; E02 para entender las matrices de A; E03 para reconstruir la variación; E05 para interpretar las entradas materiales de B. Cada tema puede estudiarse en una sesión separada.

# E.1. Campos, partículas y Modelo Estándar

**ID E01 · Estado: ABIERTO · Objetivo:** distinguir una partícula elemental, una partícula compuesta y una excitación colectiva, e interpretar masa, carga y espín sin tratarlos como una descripción completa de un estado físico.

## E.1.1. Introducción: de un modo a una excitación

Una cadena de masas permite escribir desplazamientos como suma de modos normales. En un campo, la variable se asigna a cada punto del espacio y también puede descomponerse en modos. Al cuantizar, la energía de un modo se intercambia en cuantos. Esta conexión ayuda a entender una idea de teoría cuántica de campos: una partícula es una excitación del campo correspondiente. Es un puente conceptual, no una deducción del electrón a partir de resortes [T1].

![Figura E.1. Un perfil espacial y dos reglas de ocupación de un estado cuántico. Para un estado bosónico pueden existir muchas ocupaciones; para un estado fermiónico completamente especificado, incluidos sus índices internos, sólo 0 o 1. La onda dibujada sirve de analogía de un modo; no es la trayectoria de una bolita ni una imagen literal del campo electrónico.](figuras/E_01_campos_ocupacion.png){width=100%}

La masa relaciona energía y momento de una partícula libre mediante

$$
E^2=p^2c^2+m^2c^4.
\tag{E.1}
$$

La carga eléctrica determina cómo participa en la interacción electromagnética. El espín describe cómo transforma su estado bajo rotaciones: es momento angular intrínseco, no una esfera girando sobre sí misma. Para un electrón, $s=1/2$ y una medición de $S_z$ produce $\pm\hbar/2$; esos resultados no cambian el valor de $s$. El momento angular orbital es otro grado de libertad [T2].

Masa y espín clasifican cómo se transforman excitaciones libres bajo simetrías del espacio-tiempo; las cargas describen su respuesta a simetrías e interacciones internas. Esta organización responde al tipo de ecuaciones que podemos construir, no a la imagen de bolitas con propiedades pegadas.

**Masa, carga y espín no bastan para describir todo.** También hacen falta el estado de movimiento, la ocupación, correlaciones y otros números cuánticos. Por ejemplo, los quarks tienen carga de color y participan en la interacción fuerte; los leptones no tienen esa carga. El Modelo Estándar organiza los campos elementales conocidos y sus interacciones electromagnética, débil y fuerte [P].

| Familia elemental | Dos ejemplos | Carga eléctrica en unidades de $e>0$ | Espín | Qué distingue a la familia |
|:--|:--|:--|:--|:--|
| Quarks, fermiones | up $u$, down $d$ | $+2/3$, $-1/3$ | $1/2$ | Tienen color; forman hadrones compuestos |
| Leptones, fermiones | electrón $e^-$, neutrino electrónico $\nu_e$ | $-1$, $0$ | $1/2$ | Sin color; el electrón también participa electromagnéticamente |
| Bosones de calibre | fotón $\gamma$, $W^+$ | $0$, $+1$ | $1$ | Asociados a interacciones electromagnética y débil |
| Bosón escalar | Higgs $H$ | $0$ | $0$ | Excitación del campo de Higgs |

Quarks y leptones son subfamilias de fermiones; “bosón” es una clasificación por espín y estadística. También existen bosones compuestos, por lo que la palabra no significa necesariamente elemental. Hay tres generaciones de fermiones; la tabla muestra sólo ejemplos, no el inventario completo. El protón y el neutrón son compuestos. Un fonón es una excitación colectiva de un sólido, no una nueva partícula elemental de esta tabla [P, CERN].

El mecanismo de Higgs permite masas de fermiones cargados y de $W/Z$ dentro de la teoría; no significa que toda la masa de un protón sea la suma de masas de sus quarks. Gran parte corresponde a la dinámica de quarks y gluones. Tampoco se identifica masa con carga eléctrica: pueden existir partículas masivas neutras. La gravedad no está incluida en el Modelo Estándar y la materia oscura no está identificada por él; las masas de neutrinos requieren extender su versión mínima. Éste es un marco extraordinariamente exitoso, no una teoría final de todo [H, CERN].

## E.1.2. Problema resuelto: qué dice y qué no dice una composición

**Problema.** A partir del contenido de quarks de valencia $uud$ del protón y $udd$ del neutrón, y las cargas de la tabla, obtener la carga de ambos y de un átomo de hidrógeno. Decidir si también puede obtenerse su masa y espín sumando las etiquetas de los constituyentes.

**Paso 1: carga.** Es aditiva:

$$
Q_p/e=2(2/3)-1/3=1,\qquad
Q_n/e=2/3-2(1/3)=0,\qquad
Q_{\rm H}=Q_p-e=0.
\tag{E.2}
$$

**Paso 2: masa.** La neutralidad del neutrón no implica masa nula. La energía interna y las interacciones contribuyen a la masa del sistema. La operación “sumar cargas” no se traslada sin cambios a la masa.

**Paso 3: espín.** Tres números $1/2$ no autorizan a asignar $3/2$ al protón. Los momentos angulares se acoplan y la estructura interna importa. Las etiquetas de los constituyentes restringen los estados posibles, pero no seleccionan por sí solas el estado compuesto observado.

**Lectura para el detector.** Se puede estudiar un fonón como cuanto de un modo colectivo sin resolver de nuevo los quarks de cada núcleo. La elección de variables efectivas aprovecha una separación de escalas; no contradice la descripción más fundamental.

## E.1.3. Actividad evaluada

**E01.1 — Clasificación razonada, 3 puntos.** Clasificar electrón, fotón, protón y fonón como elemental, compuesto o colectivo. Para dos de ellos, indicar qué información adicional a masa, carga y espín se necesita para describir un estado concreto. Justificar por qué “neutro” no equivale a “sin interacciones”.

**Respuesta E01.1:** 
    - El electron, al ser un fermion, se lo considera una particula elemental segun el modelo estandar, bajo la misma logica el foton, al ser un boson, es elemental. Un proton en cambio es una particula compuesta, en particular de 3 particulas elementales, dos up y un down. Finalmente el fonon se lo considera una particula (Aqui acepto que estoy estirando el significado de la palabra particula, pero a falta de una mejor palabra en el español creo que esa es la mas general) colectiva, dado que su existencia esta estrictamente ligada a la interaccion de multiples iones (Particulas compuestas) de la red cristalina de un solido (He leido de la existencia de fonones o similares para describir lo que pasa en superfluidos de astros donde estamos fuera de la categoria de solido, pero estoy restringiendo mi conociemiento a lo concerniente a SNSPD). 
    - Para caracterizar completamente un proton, por ejemplo, debemos entender el origen subyacente de la entidades que unen o provocan la interaccion de las particulas elementales que lo conforman, es decir los gluones. Para el fonon, concepto mucho mas alejado de las particulas elementales, se podrian dar muchos ejemplos de cosas que faltan para describir su estado completamente desde el modelo estandar, pero si o si hace falta hablar de su energia, $\Omega$.
    - Se le asigna la propiedad de neutro a una particula exclusivamente cuando la carga es nula, es decir no se puede percivir como una perturbacion del campo electromagnetico, mientras que sin interacciones es algo que no se puede percibir como la perturbacion de cualquier campo. 


**E01.2 — Balance y límites, 4 puntos.** Un núcleo idealizado contiene dos protones y dos neutrones, y hay tres electrones alrededor. Obtener la carga total. Explicar por qué ese cálculo no permite obtener automáticamente la masa ni el espín total. No se pide determinar si ese ion sería estable.

**Respuesta E01.2:** 
    - Cada proton tienen 2 up y 1 down, como la carga total de una particula compuesta cumple superposicion podemos decir que cada proton de la particula fundamental suma 1, es decir hasta ahora llevamos un +2 de carga. Los neutrones al tener 2 down y 1 up no aportan nada asi sean 100, seguimos en +2. Cada electron trae un -1 desde el modelo estandar, por lo que el nuevo y ultimo total de carga es -1. Supongo que un nucleo idealizado es termina siendo un ion de helio y su estabilidad me importa poco tambien.
    - Solo la carga cumple con esa suerte de superpocicion de cantidad sobre las particulas elementales que componen a las de mas alto nivel. Calcular la masa y el spin requiere conceptos mas avanzados (En los cuales no quiero inmiscuirme todavia, hasta que te pida "Campos, partículas y Modelo Estándar II").

**E01.3 — Analizar una figura, 3 puntos.** En la figura E.1, un compañero ocupa dos veces el mismo estado fermiónico porque “pueden ser dos electrones de espín distinto”. Precisar cuándo la afirmación es correcta y cuándo contradice la definición de estado utilizada. Explicar por qué un fonón no exige añadir una nueva especie al Modelo Estándar.

**Respuesta E01.3:** 
    - La afirmacion esta correcta bajo el precepto o descripcion del estado de una particula compuesta o de mayor nivel, por ejemplo en un atomo tiene multiples orbitales que pueden ser ocupados, algunos en el mismo nivel pueden ser ocupados por dos electrones, fermiones, mientras tengan diferente spin, por ello la descripcion quimica de la materia, a diferencia de la descripcion del modelo estandar, requiere mas cantidades para representar completamente el estado de un atomo o molecula, como el numero cuantico principal, secundario, magnetico, etc.
    - Ahora, si hablamos desde la definicion de estado del modelo estandar la frase esta errada, aqui se tiene como principio fundamental que un estado fermionico puede ser ocupado o no por un solo un fermion. Me parece que lo suelen llamar principio de exclusion o principio de Pauli, pero no estoy seguro. En suma, un estado fermionico incluye de por si al spin, decir que tenemos dos electrones en el mismo estado pero con diferente spin contradice, o deja abirta la ambiguedad de, la definicion de estado en el modelo estandar que de por si YA contempla al spin.
    - Un fonon puede construirse a partir de iones en una red cristalina y su dinamica, es decir a partir de conceptos que ya se pueden definir del modelo estandar, o incluso la descripcion quimica de la materia en estado solido. La tarea no es expandir el modelo mas intrinsecto y microscopico de la fisica, es usar lo que ya tiene y sus sub-productos para describir nuevos conceptos como los fonones.

**Condición esencial:** distinguir estado de especie, y excitación colectiva de partícula elemental. La tabla y la aritmética deben acompañarse de esa distinción.

## E.1.4. Conexión con A–D y fuentes

En **A.1–A.2**, espín es un índice que interviene en trazas y normalizaciones. En **B.1–B.3**, fermiones y bosones explican factores de ocupación diferentes. En **C.0**, el campo de orden es una descripción colectiva, no el campo de una nueva partícula fundamental. **D** utiliza esta separación para fijar el alcance del modelo efectivo.

Fuentes: [CERN, Modelo Estándar](https://home.cern/science/physics/standard-model/) [CERN]; [Particle Data Group, clasificación](https://pdg.lbl.gov/chris/museum_version/standard_model.html) [P]; [D. Tong, introducción a campos](https://www.damtp.cam.ac.uk/user/tong/qft/qfthtml/S0.html) [T1]; [D. Tong, partículas y espín](https://www.damtp.cam.ac.uk/user/tong/pp/pp2.pdf) [T2]; [CERN, cinco años del Higgs y masa compuesta](https://home.cern/sites/default/files/file/scientists/CCJulAug17_HIGSAT5.pdf) [H]. Las analogías y el problema resuelto son construcciones de este cuaderno.

\newpage

# E.2. Nambu y espín: dos índices distintos

**ID E02 · Estado: ABIERTO · Objetivo:** construir la matriz Nambu por espín de un singlete sencillo, interpretar su traza y explicar por qué duplicar componentes no duplica los electrones físicos.

## E.2.1. Introducción: una tabla con dos índices

Una matriz puede organizar simultáneamente dos clasificaciones. En circuitos, un vector puede llevar un índice de nodo y otro de componente; cuatro entradas no significan cuatro circuitos completos. Aquí el índice de espín tiene dos componentes y Nambu organiza operadores de destrucción y creación. La combinación es un producto tensorial de dos espacios de dimensión dos.

En una base adaptada a inversión temporal, compatible con el singlete utilizado en A, se escribe

$$
\Psi_{\mathbf k}=\begin{pmatrix}
c_{\mathbf k\uparrow}\\c_{\mathbf k\downarrow}\\
c^\dagger_{-\mathbf k\downarrow}\\-c^\dagger_{-\mathbf k\uparrow}
\end{pmatrix}.
\tag{E.3}
$$

Los operadores $c$ destruyen electrones y $c^\dagger$ los crean. El sector inferior es $i\sigma_y(c^\dagger_{-\mathbf k\uparrow},c^\dagger_{-\mathbf k\downarrow})^T=(c^\dagger_{-\mathbf k\downarrow},-c^\dagger_{-\mathbf k\uparrow})^T$. Esta elección convierte la estructura singlete en una identidad de espín en la base utilizada. Las dos entradas inferiores usan el par relacionado por inversión temporal, con un signo convencional. No son dos nuevas especies de antipartículas. Otra base puede cambiar los signos y colocar $i\sigma_y$ en el término de emparejamiento; los resultados físicos no cambian al transformar consistentemente toda la matriz [N, A].

![Figura E.2. Las cuatro entradas de una base Nambu por espín. Las filas distinguen destrucción y sector conjugado; las columnas, los dos índices de espín de la base elegida. Las matrices tau actúan entre filas y sigma entre columnas. El bloque conjugado invierte orden y signo de forma explícita; no se añaden cuatro electrones independientes.](figuras/E_02_nambu_spin.png){width=92%}

Para una banda independiente del espín, sin acoplamiento espín–órbita ni campo magnético de espín, sin superflujo ($q=0$), con $\xi_{-\mathbf k}=\xi_{\mathbf k}$ y $\Delta$ real,

$$
H_{\rm BdG}=(\xi_{\mathbf k}\tau_3+\Delta\tau_1)\otimes\sigma_0,
\qquad \xi_{\mathbf k}=\varepsilon_{\mathbf k}-\mu.
\tag{E.4}
$$

$\tau_i$ y $\sigma_i$ son matrices de Pauli en espacios diferentes. Por ejemplo, $\tau_3\otimes\sigma_0$ distingue los dos sectores Nambu; $\tau_0\otimes\sigma_3$ distingue los dos índices de espín. La representación BdG completa tiene redundancia electrón–hueco: energías positivas y negativas están relacionadas. El prefactor $1/2$ al escribir un Hamiltoniano cuadrático completo, sumando todo el espacio de momentos, en la base duplicada evita contar la misma física dos veces; los términos constantes deben conservarse cuando se calcula energía [N].

## E.2.2. Problema resuelto: diagonalizar y contar

**Problema.** En unidades arbitrarias de energía, usar $\xi=3$, $\Delta=4$. Encontrar la energía positiva y los pesos Nambu; después calcular $\operatorname{tr}_{N,s}(\tau_3\hat g)$ para $\hat g=(c\tau_3+s\tau_1)\otimes\sigma_0$.

El bloque de un par es

$$
\begin{pmatrix}3&4\\4&-3\end{pmatrix}
\begin{pmatrix}u\\v\end{pmatrix}
=E\begin{pmatrix}u\\v\end{pmatrix},\qquad
E^2=3^2+4^2=25.
\tag{E.5}
$$

Para $E=5$, la primera fila exige $v=u/2$. La normalización $|u|^2+|v|^2=1$ da $|u|^2=4/5$, $|v|^2=1/5$. Se trata de una mezcla coherente de amplitudes; no de una moneda que escoge dos partículas clásicas independientes.

Dentro de una traza sobre ambos espacios, la abreviatura $\tau_i$ significa $\tau_i\otimes\sigma_0$. Para la traza se usa $\operatorname{tr}_N\tau_i\tau_j=2\delta_{ij}$ y $\operatorname{tr}_s\sigma_0=2$:

$$
\operatorname{tr}_{N,s}(\tau_3\hat g)
=\big(c\operatorname{tr}_N\tau_3^2+s\operatorname{tr}_N\tau_3\tau_1\big)
\operatorname{tr}_s\sigma_0=4c.
\tag{E.6}
$$

El cuatro procede de la traza en ambos espacios. No puede añadirse de nuevo como “cuatro partículas” al integrar una DOS que ya incorpore convenciones de espín y excitaciones. Ésta es la operación concreta usada al reducir el funcional de A.

## E.2.3. Actividad evaluada

**E02.1 — Cálculo, 4 puntos.** Repetir la diagonalización para $\xi=5$ y $\Delta=12$ en las mismas unidades. Dar ambas energías y los pesos normalizados de la rama positiva. Mostrar al menos una verificación con la matriz original.

**Respuesta E02.1:** _Escribir aquí._

**E02.2 — Índices, 3 puntos.** Calcular las matrices diagonales $\tau_3\otimes\sigma_0$ y $\tau_0\otimes\sigma_3$ en la base E.3. Indicar qué parejas de entradas diferencia cada una. Explicar el significado del signo de la cuarta componente.

**Respuesta E02.2:** _Escribir aquí._

**E02.3 — Auditoría de un argumento, 3 puntos.** Un cálculo usa $N_0$ por espín, toma una traza Nambu/espín completa y después multiplica todo por cuatro “porque hay cuatro partículas en el espinor”. Identificar qué comprobaciones de normalización exigirías y por qué esa justificación es incorrecta. No se pide adivinar el prefactor final sin conocer la definición de la integral.

**Respuesta E02.3:** _Escribir aquí._

**Condición esencial:** reconocer la redundancia Nambu y mantener separados los espacios de espín y Nambu. No interpretar los huecos de esta base como positrones libres.

## E.2.4. Conexión con A–D y fuentes

**A.2.2, A.6b–A.6d:** base, trazas y gradiente covariante. **B.2:** el espectro anómalo lleva información de coherencia que la DOS por sí sola no contiene. **C.0:** la fase y amplitud del condensado corresponden a una descripción colectiva derivada; no son componentes del espín. **D:** la normalización debe conservarse al trasladar resultados entre capas.

Fuentes: [curso de TU Delft, representación BdG y redundancia](https://topocondmat.org/test/w1_topointro/0d.html) [N]; [Virtanen et al., expresión quasiclásica](https://arxiv.org/pdf/1909.00992v1), ecuación (20) [A]. La base de E.3 es la adaptación explícita usada aquí; no se copian sin transformar signos de bases diferentes. La diagonalización y el cálculo de trazas anteriores son ejercicios originales de esta revisión.

\newpage

# E.3. Variaciones, funcionales y teorema de la envolvente

**ID E03 · Estado: ABIERTO · Objetivo:** calcular una primera variación con condiciones de borde y distinguir variar un campo de eliminar una variable estacionaria.

## E.3.1. Introducción: una coordenada en cada punto

Una función $F(y)$ recibe números. Un funcional $\mathcal F[y]$ recibe una función completa y devuelve un número. En mecánica lagrangiana ya aparece un ejemplo: la acción recibe una trayectoria. La energía de una cuerda deformada recibe su perfil espacial. La derivada funcional responde cuánto cambia el resultado si se modifica el perfil cerca de cada posición.

En una malla con coordenadas $y_1,\ldots,y_N$, la energía es una función de $N$ variables. El continuo organiza sus derivadas como una densidad:

$$
\delta\mathcal F[y;\eta]
=\left.\frac{d}{d\epsilon}\mathcal F[y+\epsilon\eta]\right|_{\epsilon=0}
=\int\frac{\delta\mathcal F}{\delta y(x)}\eta(x)\,dx,
\tag{E.7}
$$

con los términos de borde tratados según las variaciones permitidas. $\eta$ es una dirección de prueba; $\epsilon$ controla su tamaño. $\delta\mathcal F/\delta y$ no es la energía ni la derivada de una curva de equilibrio, sino el coeficiente de cada perturbación local.

![Figura E.3. Izquierda: un perfil y una perturbación local admisible. Derecha: la variable interna z se acomoda en el valle z=2y de una energía de dos variables. Seguir el valle ilustra eliminar una variable estacionaria; la pendiente y la curvatura deben distinguirse.](figuras/E_03_variacion_envolvente.png){width=87%}

## E.3.2. Problema resuelto: perfil estático y variable eliminada

**Parte I: funcional.** En unidades adimensionales y con extremos fijos $y(0)=y(\pi)=0$, considerar

$$
\mathcal F[y]=\int_0^\pi\left[\frac{(y')^2}{2}+\frac{y^2}{2}-2y\sin x\right]dx.
\tag{E.8}
$$

Insertar $y+\epsilon\eta$ y conservar términos lineales produce

$$
\begin{aligned}
\delta\mathcal F&=\int_0^\pi[y'\eta'+y\eta-2\sin x\,\eta]dx\\
&=[y'\eta]_0^\pi+\int_0^\pi[-y''+y-2\sin x]\eta\,dx.
\end{aligned}\tag{E.9}
$$

Como los extremos de $y$ están fijados, $\eta(0)=\eta(\pi)=0$: desaparece el término de borde. La arbitrariedad de $\eta$ da $-y''+y=2\sin x$. La solución $y=\sin x$ satisface ecuación y fronteras. Además, $\delta^2\mathcal F=\int[(\eta')^2+\eta^2]dx>0$ para toda perturbación no nula: aquí el estacionario es un mínimo. La primera variación nula por sí sola no habría demostrado esa estabilidad.

Si los extremos fuesen libres, no se puede borrar el término de borde por costumbre. Para este funcional, la condición natural sería $y'=0$ en esos extremos. Es la misma lógica que distingue imponer un desplazamiento e imponer una fuerza en mecánica.

**Parte II: envolvente.** Considérese ahora la función ordinaria

$$
F(y,z)=\frac{y^2}{2}+\frac{(z-2y)^2}{2}.
\tag{E.10}
$$

La condición interna $F_z=z-2y=0$ da $z^*(y)=2y$. La función reducida es $\bar F(y)=F(y,z^*(y))=y^2/2$. Por regla de la cadena,

$$
\frac{d\bar F}{dy}=F_y+F_z\frac{dz^*}{dy}
=y-2(z^*-2y)+0=y.
\tag{E.11}
$$

El término de la derivada de $z^*$ desaparece porque multiplica $F_z=0$, no porque $z^*$ sea constante. Es la forma de la envolvente que utiliza A.2.3: los ángulos espectrales ya satisfacen su ecuación, mientras la amplitud todavía puede variar.

**Una precaución de segundo orden.** No se debe congelar otra vez la variable al calcular una segunda derivada. En este ejemplo $F_{yy}=5$ a $z$ fijo, pero $\bar F''=1$. Para una variable interna no degenerada,

$$
\bar F''=F_{yy}-\frac{F_{yz}^2}{F_{zz}}=5-\frac{(-2)^2}{1}=1.
\tag{E.12}
$$

La relajación interna reduce aquí la rigidez. Esto importa para capacidades y estabilidad: una identidad de primeras derivadas no autoriza a omitir la respuesta espectral en todas las derivadas posteriores. En un extremo degenerado o un cambio de rama se deben revisar las hipótesis de diferenciabilidad.

## E.3.3. Actividad evaluada

**E03.1 — Variación y fronteras, 4 puntos.** Para

$$
\mathcal G[y]=\int_0^\pi\left[(y')^2+\frac{y^2}{2}-3y\sin x\right]dx,
\qquad y(0)=y(\pi)=0,
\tag{E.13}
$$

obtener la primera variación, conservar el término de borde antes de usar las condiciones y encontrar una solución proporcional a $\sin x$. Explicar si el estacionario es estable mediante la segunda variación.

**Respuesta E03.1:** _Escribir aquí._

**E03.2 — Error de eliminación, 3 puntos.** Para $H(y,z)=y^2+(z-3y)^2/2$, encontrar $z^*(y)$ y comparar $dH(y,z^*)/dy$ con la derivada a $z$ fijo evaluada en esa rama. Después repetir la comparación cuando alguien utiliza por error $z=3y+0.1$. Indicar cuál hipótesis del teorema dejó de cumplirse.

**Respuesta E03.2:** _Escribir aquí._

**E03.3 — Lectura física, 3 puntos.** Explicar por qué se pueden eliminar los ángulos estacionarios de A sin imponer $G=0$ durante todo un evento. Indicar qué cantidad adicional se necesita para convertir la pendiente de energía en una velocidad de evolución. Usar la figura E.3 para justificar la diferencia entre pendiente y rigidez.

**Respuesta E03.3:** _Escribir aquí._

**Condición esencial:** no eliminar términos de borde sin justificarlo y no confundir la estacionariedad espectral con el equilibrio instantáneo de la amplitud.

## E.3.4. Conexión con A–D y fuentes

**A.2.3, A.8a–A.8b:** aplican la variación y la envolvente a [Virtanen et al., ecuación (20)](https://arxiv.org/pdf/1909.00992v1). **B.5 y B.11:** fijar espectro u ocupaciones produce derivadas distintas. **C.1–C.4:** fuerzas, movilidad y trabajo completan la dinámica. **D:** las identidades prueban coherencia interna, no validación experimental. Los ejemplos E.7–E.13 se derivan aquí.

\newpage

# E.4. Electrones y huecos: del diodo al superconductor

**ID E04 · Estado: ABIERTO · Objetivo:** calcular carga y energía respecto de una referencia ocupada y distinguir un hueco de semiconductor, una antipartícula y una componente Nambu.

## E.4.1. Introducción: llevar la cuenta de lo que falta

En un semiconductor, describir por separado todos los electrones de una banda casi llena resulta poco práctico. Se elige la banda llena como referencia y se sigue la escasa población de vacantes. Si un electrón ocupa la vacante vecina, la vacante se desplaza en sentido contrario. Es una nueva contabilidad del mismo sistema electrónico [M1].

En un diodo se comparan poblaciones de electrones de conducción y huecos de valencia. En equilibrio, difusión y arrastre por el campo de la unión se compensan; eso no significa que cada portador esté inmóvil. Bajo polarización cambian las poblaciones y aparecen corrientes netas. El campo y las bandas del diodo hacen familiar la elección de una referencia, pero el gap de un semiconductor no procede del emparejamiento de Cooper [M2].

![Figura E.4. Izquierda: un electrón ocupa una vacante vecina y el hueco efectivo se desplaza en sentido contrario. Derecha: comparación entre un electrón añadido por encima de mu y uno retirado por debajo; ambas operaciones producen excitaciones de energía positiva. Los niveles son un dibujo de contabilidad, no un perfil de bandas de un diodo real.](figuras/E_04_electron_hueco.png){width=100%}

Aquí “energía de excitación” se mide respecto del equilibrio a potencial químico fijado: es el cambio de $U-\mu N$. Si cambia el número electrónico, $\delta U=\mu\,\delta N+E_{\rm excitaciones}$; la energía interna absoluta no se confunde con esa referencia.

Si el electrón retirado tenía energía $\varepsilon<\mu$ y carga $-e$, la excitación relativa al mar lleno cambia la carga en $+e$ y cuesta $\mu-\varepsilon>0$. La energía negativa relativa $\xi=\varepsilon-\mu$ del estado electrónico no es una energía negativa de la excitación. En un metal normal, cerca de la superficie de Fermi, se puede llevar esta cuenta separando electrones añadidos y vacantes.

En superconductividad, la excitación mezcla ambos canales. Para un bloque BCS uniforme, $E=\sqrt{\xi^2+|\Delta|^2}$ y

$$
|u|^2=\frac12\left(1+\frac{\xi}{E}\right),\qquad
|v|^2=\frac12\left(1-\frac{\xi}{E}\right).
\tag{E.14}
$$

En la contabilidad de carga electrónica de ese modo, $Q_{\rm qp}=-e(|u|^2-|v|^2)=-e\xi/E$. Es una expectativa de carga relativa de la excitación, no una nueva especie con carga fundamental fraccionaria. El condensado y la respuesta de pantalla completan la conservación eléctrica del sistema. En $\xi=0$, los dos pesos son iguales: una excitación puede contener energía aun cuando esa expectativa de carga se anule [N].

El hueco de un sólido no es un positrón: requiere una referencia electrónica ocupada y propiedades de banda. La componente conjugada de Nambu tampoco añade una banda de valencia física. Es una organización útil para escribir el emparejamiento. Simetría electrón–hueco y ausencia de desequilibrio de carga son hipótesis de un sector cinético; no obligan a que la población sea térmica ni a que no exista una corriente superconductora.

## E.4.2. Problema resuelto: energía positiva con carga neta nula

**Problema.** Con referencia al estado ocupado a $T=0$, añadir un electrón a $\mu+3$ meV y retirar otro de $\mu-2$ meV. Obtener número neto de electrones, carga y energía de excitación. Después compararlo con una cuasipartícula BCS en $\xi=0$.

La primera operación aporta $\delta N=+1$, $\delta Q=-e$ y 3 meV. La segunda aporta $\delta N=-1$, $\delta Q=+e$ y 2 meV. Por tanto,

$$
\delta N_{\rm total}=0,\qquad \delta Q_{\rm total}=0,
\qquad E_{\rm excitaciones}=5\,\mathrm{meV}.
\tag{E.15}
$$

Un calorímetro puede registrar la energía aunque no se haya inyectado carga neta. Dos poblaciones de energías diferentes pueden mantener esa neutralidad. La misma distinción es importante al absorber un fotón en un sólido.

Para una cuasipartícula BCS con $\xi=0$, E.14 da $|u|^2=|v|^2=1/2$ y $E=|\Delta|$. Aquí hay una superposición coherente dentro del modo; no es el mismo estado que dos excitaciones normales clásicamente separadas. La coincidencia de carga neta no identifica los dos estados.

## E.4.3. Actividad evaluada

**E04.1 — Contabilidad, 4 puntos.** Añadir un electrón a $\mu+4$ meV y retirar dos electrones, uno a $\mu-1$ meV y otro a $\mu-2$ meV. Obtener $\delta N$, $\delta Q$ y energía total de excitación. Explicar por qué no se suman los valores de $\xi$ con sus signos para obtener esa energía.

**Respuesta E04.1:** 
    - $\delta N = -1$, $\delta Q = +e$ y $E_{\rm excitaciones}=7$. 
    - Porque hay que tener en cuenta el piso de energia para las QPs en un superconductor, el gap, no se pueden crear aquellas con energia menor al gap en condiciones normales/ideales. 

**E04.2 — Analizar el dibujo, 3 puntos.** Si la vacante del panel izquierdo se mueve a la derecha, indicar el sentido del salto electrónico y el sentido de corriente convencional asociado a ese salto. Explicar qué parte de esta analogía sirve para la unión p–n y qué parte no basta para describir una cuasipartícula superconductora.

**Respuesta E04.2:** 
    - Para que la vacante se mueva en ese sentido un electron adyacente debe haber tomado el antiguo lugar, izq, y dejar libre el lugar al que llega la vacante o hueco, der, es decir se movio a la izquierda. La corriente convencional va en sentido del movimiento de carga positiva (A pesar de que las cargas que se mueven son en realidad los electrones y estos tienen carga negativa), en este caso el hueco, por lo que la corriente va a la derecha.
    - En esta analogia podemos describir el proceso elemental para la creacion de la zona de depresion/polarizacion en un diodo el intercambio de electrones con huecos, a un nivel de mayor escala de este proceso se puede recuperar el comportamiento del diodo o una juntura p-n. Eso respecto a semiconductores, respecto a superconductores la cosa es mas compleja dada la existencia del gap y exitaciones que tiene que tienen que interactuar con el condensado. Por poner un ejemplo es facil añadir un electron con un nivel energetico de 1meV en un semiconductor, pero si se trata de un superconductor como el NbN (Tc aproxx 10 Kelvin) ese electron no alcanza a calificar como QPs porque su energia es demasiado baja y se va a incorporar al condensado de manera espontanea.

**E04.3 — Separar conceptos, 3 puntos.** Evaluar estas afirmaciones con una corrección breve: “un hueco es un positrón dentro del metal”; “energía de excitaciones positiva exige carga neta distinta de cero”; “simetría electrón–hueco permite reemplazar cualquier distribución por Fermi–Dirac”.

**Respuesta E04.3:** 
    - A esa frase le sobra el termino positron, lo que realmente es un hueco es la ausencia de carga que, dentro del mar de electrones en el cristal o semiconductor, se percibe como una zona local relativamente positiva en cuanto a carga electrica. 
    - Esta frase se cae al exigir carga carga neta distinta de cero, por contradiccion, si añado 1 electron con energia n meV y saco uno con m meV, y la n > m, la carga neta es nula a pesar de que la energia de las excitaciones es positiva.
    - La frase trata de unir dos conceptos a travez de una implicancia falsa. La simetria electron-hueco habla de como el desplazamiento de la carga en posiciones fijas se puede interpretar, mientras que la distribucion Fermi-Dirac nos dice que una nube o conjunto de electrones distribuye de cierta manera estadistica dada una temperatura y potencial quimico.

**Condición esencial:** definir la referencia antes de asignar signos de energía y carga, y distinguir neutralidad, coherencia y termalidad.

## E.4.4. Conexión con A–D y fuentes

**A.2:** la contribución normal cuenta electrones y vacantes; explica factores que reaparecen en las energías. **B.1 y B.11:** el sector simétrico todavía conserva distribuciones no térmicas. **C.0–C.4:** corriente y conservación de energía son ecuaciones relacionadas, pero distintas. **D:** una sola temperatura equivalente no reemplaza esos balances ni todas las poblaciones.

Fuentes: [MIT 6.012, portadores en semiconductores](https://live.ocw.mit.edu/courses/6-012-microelectronic-devices-and-circuits-fall-2009/a37c73e7fcff4bcb5c5f1177e12d0891_MIT6_012F09_lec01.pdf) [M1]; [MIT 6.012, transporte y unión p–n](https://web.mit.edu/6.012/www/) [M2]; [TU Delft, representación BdG](https://topocondmat.org/test/w1_topointro/0d.html) [N]. E.14 sigue también de diagonalizar la matriz de E.4; la contabilidad numérica es un ejemplo propio.

\newpage

# E.5. De modos normales a DOS y DFPT

**ID E05 · Estado: ABIERTO · Objetivo:** construir una DOS a partir de modos, distinguirla de una ocupación y de un espectro ponderado de interacción, y explicar qué calcula DFPT para alimentar B.

## E.5.1. Introducción: cuántos modos hay y cuánto interactúan

Para $N$ masas con un desplazamiento escalar cada una, la teoría de oscilaciones pequeñas entrega $N$ modos normales. Escribir $N$ modos no implica que todos estén excitados. Los autovectores describen las formas; los autovalores, las frecuencias; las ocupaciones indican cuánta energía hay en cada modo. Son tres piezas distintas.

Una DOS reúne las frecuencias en una función que cuenta estados. Para una cadena finita,

$$
F_\omega(\omega)=\sum_{\nu=1}^{N}\delta(\omega-\omega_\nu),
\qquad \int_0^\infty F_\omega(\omega)\,d\omega=N.
\tag{E.16}
$$

La delta de Dirac concentra un área de un modo. Su altura puntual no es una población infinita. En una figura se puede reemplazar por un pico estrecho de área conservada; ese ensanchamiento gráfico no representa automáticamente una vida media física.

![Figura E.5. Cadena de dos masas entre paredes y sus dos modos: en fase y en oposición. Las frecuencias alimentan dos líneas de DOS con área unitaria. Los mismos modos pueden tener pesos de acoplamiento diferentes; contar modos y medir cuánto interactúan responde a preguntas distintas.](figuras/E_05_modos_dos.png){width=100%}

En un cristal tridimensional con $s$ átomos por celda hay $3s$ ramas por vector de onda $\mathbf q$. Una convención por celda y por energía es

$$
F_\Omega(\Omega)=\frac{1}{N_q}\sum_{\mathbf q\nu}
\delta(\Omega-\hbar\omega_{\mathbf q\nu}),\qquad
\int F_\Omega\,d\Omega=3s.
\tag{E.17}
$$

La normalización por átomo dividiría este conteo por $s$. Si $N_i$ cuenta celdas por volumen, debe multiplicar una DOS por celda; no una DOS ya normalizada por átomo. Además,

$$
F_\Omega(\Omega)=\frac1\hbar F_\omega(\Omega/\hbar).
\tag{E.18}
$$

Cambiar el eje sin transformar la altura cambia el número de modos. Si se usa frecuencia ordinaria $\nu$ en THz, la energía es $h\nu$, no $\hbar\nu$. La DOS electrónica normalizada $\rho$ de A es otra cantidad: debe multiplicarse por su escala electrónica, no por $N_i$.

**Dónde entra DFPT.** Las frecuencias provienen de la curvatura de la energía electrónica más iónica respecto de desplazamientos atómicos. DFPT calcula la respuesta electrónica lineal autoconsistente a una perturbación; con ella se obtienen constantes de fuerza y acoplamientos electrón–fonón. Así reemplaza un resorte supuesto por una respuesta del material calculada dentro de DFT [DF].

## E.5.2. Problema resuelto: de resortes a entradas materiales

**Parte I: dos masas.** Dos masas iguales $m$ están unidas por tres resortes iguales $K$, con paredes fijas. Para desplazamientos $u_1,u_2$,

$$
V=\frac K2[u_1^2+(u_2-u_1)^2+u_2^2],\qquad
m\ddot{\mathbf u}=-K\begin{pmatrix}2&-1\\-1&2\end{pmatrix}\mathbf u.
\tag{E.19}
$$

Los autovectores normalizados son $(1,1)/\sqrt2$ y $(1,-1)/\sqrt2$. Los autovalores de la matriz adimensional son 1 y 3: $\omega_1=\sqrt{K/m}$ y $\omega_2=\sqrt{3K/m}$. No aparecen cuatro modos por tener dos masas y dos resortes exteriores; cuentan los grados de libertad independientes.

**Parte II: población.** Elegir $\hbar\sqrt{K/m}=1$ en una unidad de energía. Si las ocupaciones son $n_1=2$ y $n_2=0$, la energía de excitaciones es 2 unidades. La energía de punto cero se omite de esta referencia y no cambia ese resultado. La DOS contiene dos modos antes y después de excitar: lo que cambió fue $n_\nu$.

**Parte III: interacción.** Como ejemplo espectral, sea

$$
\alpha^2F(\Omega)=W_1\delta(\Omega-\Omega_1)
+W_2\delta(\Omega-\Omega_2),\qquad
\lambda=2\left(\frac{W_1}{\Omega_1}+\frac{W_2}{\Omega_2}\right).
\tag{E.20}
$$

Con eje de energía y $\alpha^2F$ adimensional, $W_i$ tiene unidades de energía. Si $\Omega_1=1$, $\Omega_2=\sqrt3$, $W_1=0.10$ y $W_2=0.20$ en esa unidad, $\lambda=0.43094$. Los pesos no son ocupaciones ni áreas de la DOS sin ponderar. Si se duplican ambos $W_i$, se duplica $\lambda$ sin cambiar las frecuencias ni el número de modos. Este ejemplo ilustra la integral de A.32, no un ajuste de NbN.

**Parte IV: lo que sustituye al resorte en DFPT.** Cerca de una estructura de equilibrio,

$$
E_{\rm BO}(\{u\})=E_0+\frac12\sum_{IJ}\Phi_{IJ}u_Iu_J+O(u^3),
\qquad \Phi_{IJ}=\frac{\partial^2E_{\rm BO}}{\partial u_I\partial u_J}.
\tag{E.21}
$$

$I,J$ incluyen átomo, celda y dirección. La matriz dinámica transforma estas curvaturas mediante las masas: $D_{IJ}=\Phi_{IJ}/\sqrt{M_IM_J}$, con la transformada espacial apropiada para cada $\mathbf q$. Sus autovalores son $\omega^2$. No confundir esta matriz $D_{IJ}$ con el coeficiente escalar de difusión $D$ de pySNSPD.

En un ejemplo electrónico no degenerado, la ecuación de respuesta tiene la estructura

$$
(H_{\rm KS}-\varepsilon_n)\,\delta\psi_n
=-(\delta V_{\rm KS}-\delta\varepsilon_n)\psi_n.
\tag{E.22}
$$

La densidad cambia con $\delta\psi_n$, y $\delta V_{\rm KS}$ cambia con esa densidad: hay que resolver la respuesta autoconsistentemente y fijar la libertad de fase del estado. En metales y subespacios degenerados se necesitan ocupaciones y proyectores adecuados; E.22 enseña la estructura, no es una receta completa de implementación [DF].

Finalmente, la perturbación del potencial proyectada sobre un modo conecta estados electrónicos. De forma esquemática,

$$
g_{mn\nu}(\mathbf k,\mathbf q)
\sim\left\langle\psi_{m,\mathbf k+\mathbf q}\left|
\delta_{\mathbf q\nu}V_{\rm KS}\right|\psi_{n,\mathbf k}\right\rangle
\times\text{amplitud cuántica del modo}.
\tag{E.23}
$$

En E.23, $\delta_{\mathbf q\nu}V_{\rm KS}$ designa la derivada del potencial respecto de la coordenada normal del modo, todavía sin multiplicar por su amplitud cuántica. Su normalización debe especificar las masas; si la amplitud ya estuviera incluida en la perturbación, no se multiplicaría otra vez.

El espectro $\alpha^2F$ combina pesos $|g|^2$, modos fonónicos y estados electrónicos cerca del nivel de Fermi. La documentación de PHonon define sus prefactores y normalización; E.23 sólo separa las piezas físicas [QE]. Los factores de ocupación y coherencia superconductora de B se aplican después: DFPT no calcula por sí sola una cascada no térmica, una burbuja inicial ni la gTDGL del detector.

![Figura E.6. Cadena de dependencias: estructura y respuesta electrónica producen curvaturas y elementos de interacción; de ellos se obtienen modos, DOS y espectro ponderado. B añade ocupaciones y factores de coherencia para formar tasas y balances. Las flechas muestran dependencias, no una equivalencia entre todos los objetos.](figuras/E_06_dfpt_a_cinetica.png){width=100%}

## E.5.3. Actividad evaluada

**E05.1 — Modos y DOS, 4 puntos.** Para tres masas iguales con cuatro resortes iguales y extremos fijos, formar la matriz de rigidez. Sin necesidad de hallar fórmulas exactas para las tres frecuencias, justificar cuántos modos hay y cuánto vale la integral de la DOS total. Explicar cómo cambiaría esa integral si se normalizara por oscilador y, por separado, por unidad de masa física de toda la cadena. Indicar las unidades en ambos casos. Comprobar la dimensión de la matriz dinámica.

**Respuesta E05.1:** _Escribir aquí._

**E05.2 — Interpretar resultados, 3 puntos.** Dos cálculos tienen idénticas frecuencias y DOS fonónica, pero el segundo duplica $|g|^2$ para todos los procesos. Decidir qué cambia entre número de modos, $\alpha^2F$, $\lambda$ y ocupación inicial de fonones. Explicar qué información falta para predecir una potencia electrón–fonón.

**Respuesta E05.2:** _Escribir aquí._

**E05.3 — Unidades y frontera entre métodos, 3 puntos.** Una DOS está tabulada por celda y por THz. Describir cómo convertir tanto el eje como los valores a una DOS por volumen y por joule, identificando la densidad de celdas necesaria. Explicar en qué paso del esquema E.6 entran los factores de coherencia de B.2 y las ocupaciones de B.3.

**Respuesta E05.3:** _Escribir aquí._

**Condición esencial:** no confundir DOS, ocupación y peso de interacción; conservar el número de modos al cambiar eje y normalización.

## E.5.4. Conexión con A–D y fuentes

**A.1 y A.5:** convenciones y significado de $\lambda$; la microscopía de equilibrio y la cinética deben declarar sus aproximaciones. **B.2–B.3:** las entradas $F$ y $\alpha^2F$ no sustituyen los espectros superconductores ni sus factores de ocupación. **C:** estas tasas afectan las fuerzas y el intercambio, pero no fijan automáticamente la movilidad de orden. **D:** la auditoría material exige saber de qué estructura, normalización y aproximación proviene cada tabla.

Fuentes: [Baroni, de Gironcoli, Dal Corso y Giannozzi, DFPT](https://arxiv.org/abs/cond-mat/0012092), revisión publicada en Reviews of Modern Physics **73**, 515 (2001) [DF]; [Quantum ESPRESSO, coeficientes electrón–fonón](https://www.quantum-espresso.org/Doc/ph_user_guide/node19.html) [QE]; [PHonon, capacidades y estructura](https://www.quantum-espresso.org/Doc/ph_user_guide/node5.html); [EPW, interpolación de acoplamientos](https://arxiv.org/abs/1005.4418). La cadena de masas, las DOS discretas, sus pesos y las actividades son ejemplos originales; no se ejecutó DFT/DFPT del material para confeccionarlos.

# E.6. Registro de evaluación para la siguiente iteración

| Tema | Respuestas recibidas | Puntaje | Condición esencial | Próximo estado |
|:--|:--|:--|:--|:--|
| E01 | Pendientes | Sin evaluar | Sin evaluar | Abierto |
| E02 | Pendientes | Sin evaluar | Sin evaluar | Abierto |
| E03 | Pendientes | Sin evaluar | Sin evaluar | Abierto |
| E04 | Pendientes | Sin evaluar | Sin evaluar | Abierto |
| E05 | Pendientes | Sin evaluar | Sin evaluar | Abierto |

**Comentarios de aprendizaje para la revisión 0.4:** _Escribir aquí qué explicación funcionó, qué paso quedó oscuro o qué analogía convendría cambiar. Estos comentarios no sustituyen los ejercicios y no afectan negativamente la evaluación._

**Pauta de las actividades evaluadas:** pendiente de la próxima revisión con respuestas. Se incorporará al mismo tema y conservará los identificadores E01–E05 para poder seguir los avances sin perder versiones anteriores.
