# Qué reutilizar de la memoria y cómo representar la energía espacial

Revisión documental y algebraica del 24 de septiembre de 2026. No ejecuta
cálculos ni modifica producción. Su alcance inmediato es el bloque **térmico
espacial** de referencia, con el mismo material difusivo, singlete isotrópico,
acoplamiento débil y ausencia de campo propio ya declarados. No lo convierte
en un transiente no térmico completo.

La simplificación útil consiste en cambiar las coordenadas espectrales y
eliminar variables estacionarias. No consiste en volver a cerrar la física
del núcleo con una división por la amplitud del condensado.

## 1. Lo que realmente hace el solver de la memoria

Se cotejaron los anexos C y D de
`tmp/pdfs/modelo_v0_2/sources/memoria_02.pdf`, con su texto asociado, y los
módulos de producción siguientes:

| Pieza | Implementación y alcance |
|---|---|
| Malla y operadores | `pysnspd/mesh/pytdgl_like.py` y `pysnspd/gtdgl/tdgl_operators.py`: Delaunay–Voronoi, áreas de control, gradiente covariante en aristas y divergencia compatible. Memoria D.3. |
| Fuerza histórica | `pysnspd/gtdgl/allmaras.py:compute_allmaras_forcing_dimensionless`: laplaciano complejo con `xi_mod²`, reacción interpolada GL y corrección imaginaria de divergencias Usadel–GL. Memoria C.3 y D.4–D.6. |
| Corriente histórica | `pysnspd/gtdgl/usadel_current.py`: tabla uniforme de rigidez en temperatura, amplitud y gradiente de fase, multiplicada por el flujo regular de pares de la arista. Memoria D.5. |
| Paso local | `pysnspd/solver/core.py:solve_for_psi_squared`: fuerza evaluada en el tiempo anterior y nueva amplitud implícita mediante una ecuación cuadrática. Memoria D.2–D.3. |
| Potencial y circuito | Proyección de continuidad con matriz dispersa y circuito de tres estados. La topología y parámetros circuitales vigentes siguen siendo los de la memoria. |

La ley de corriente de producción es

\[
j_{ij}=K_{ij}^{\rm tabla}
\frac{\operatorname{Im}(\Delta_i^*U_{ij}\Delta_j)}{\ell_{ij}}.
\]

No divide ese flujo por `|Delta|²`: desaparece cuando un extremo tiene
amplitud cero. Sin embargo, la corrección de fase sí contiene el cociente de
la diferencia de divergencias por `|psi|²`. `PhaseDriveContinuationSolver`
prolonga armónicamente el campo no resuelto usando el mismo grafo. La memoria
describe esta operación como un procedimiento numérico controlado; no aporta
por esa vía un nuevo funcional microscópico de núcleo.

El apéndice B.3 de Allmaras explica también por qué usar partes real e
imaginaria del condensado evita la singularidad de sus coordenadas polares.
Su texto distingue esa representación de las aproximaciones físicas del
cierre. Cambiar coordenadas no demuestra que una fuerza y una corriente
arbitrarias sean derivadas de la misma energía.

**Reutilización segura:** geometría, pesos, enlaces de calibre, ensamblaje
disperso, proyección de corriente normal, trazabilidad de unidades, gestión
de pasos aceptados y circuito. La actualización local KWT puede estudiarse
como integrador de una fuerza nueva una vez fijadas sus unidades y
coeficientes; el paso algebraico no la valida físicamente. No se suma la
corrección Usadel–GL histórica a una fuerza que ya deriva de la energía
espacial completa: se duplicaría parte del acoplamiento. Tampoco se heredan
automáticamente los contactos normales de producción si el ensayo vigente
requiere reservorios superconductores.

## 2. Una representación espectral sin dividir por Delta

El punto de partida es el funcional difusivo de Virtanen, Vargunin y Silaev,
ecuación 20, conservando `Delta` independiente antes de imponer
autoconsistencia. La reducción uniforme ya está desarrollada en A.6–A.12.
La formulación siguiente mantiene sus gradientes espaciales y explicita una
transformación algebraica para implementarlos. No es una nueva aproximación
microscópica introducida para ajustar las curvas.

Para cada frecuencia positiva `epsilon_n=pi*kB*T*(2n+1)` se representa el
propagador mediante un campo anómalo complejo `f_n` y una componente normal
real `g_n`, con

\[
|f_n|^2+g_n^2=1,\qquad g_n>0,
\qquad D_A=\nabla-i\frac{2e}{\hbar}\mathbf A.
\]

La `f_n` de esta sección es un **propagador espectral**, no la función de
ocupación electrónica. La energía libre relativa al metal normal, para
temperatura uniforme impuesta, es

\[
\begin{aligned}
\mathcal F/N_0={}&\int dV\,|\Delta|^2\ln(T/T_c)\\
&+2\pi k_BT\sum_{n\ge0}\int dV\left[
\frac{|\Delta|^2}{\epsilon_n}+2\epsilon_n(1-g_n)
-2\operatorname{Re}(\Delta^*f_n)
+\frac{\hbar D}{2}\left(|D_Af_n|^2+|\nabla g_n|^2\right)
\right].
\end{aligned}
\]

`N0` es por espín, como en A, y `sigma_n=2*e²*N0*D`. La resta normal y la
renormalización por `Tc` se realizan de manera común antes de retirar el
corte espectral. No se añade a esta expresión otro término `K0*|grad Delta|²`:
la respuesta espacial ya está en los propagadores.

En una región donde `f=sin(Theta)*exp(i*chi)` y `g=cos(Theta)`, la identidad

\[
|D_Af|^2+|\nabla g|^2
=|\nabla\Theta|^2+\sin^2\Theta
\left|\nabla\chi-\frac{2e}{\hbar}\mathbf A\right|^2
\]

recupera la forma angular habitual. En un núcleo basta trabajar con el
campo complejo; no hace falta construir `arg Delta`, ni imponer un mínimo a
`|Delta|`, ni definir `q_delta`.

### Dos coordenadas equivalentes en la rama térmica

La elección directa `g=sqrt(1-|f|²)` es regular alrededor de `f=0`. Su
gradiente tiene la forma

\[
|D_A f|^2+
\frac{[\operatorname{Re}(f^*D_Af)]^2}{1-|f|^2}.
\]

El denominador se acerca a cero cuando `g` se acerca a cero, no cuando
desaparece el condensado. Por tanto no debe describirse esta coordenada
como uniformemente bien condicionada en todo el dominio.

Otra elección exacta es el campo complejo `u=f/g`:

\[
f=\frac{u}{\sqrt{1+|u|^2}},\qquad
g=\frac{1}{\sqrt{1+|u|^2}}.
\]

Todo `u` finito satisface la normalización y `g>0`. Su métrica espacial es

\[
\frac{|D_Au|^2}{1+|u|^2}
-\frac{[\operatorname{Re}(u^*D_Au)]^2}{(1+|u|^2)^2}.
\]

Por Cauchy–Schwarz es al menos
`|D_Au|²/(1+|u|²)²`: el término de gradiente conserva signo. En el estado
homogéneo de corriente nula `u=Delta/epsilon_n`; con `T>0` el estado
`Delta=0` es una coordenada ordinaria. Las frecuencias pequeñas pueden
producir `|u|` grande, así que siguen siendo necesarios escalado y control
del residuo. No se promete una condición numérica uniforme al tomar `T→0`.

## 3. La fase espectral no se bloquea artificialmente a la del gap

En 2D general, `chi_n=arg(f_n)` debe poder depender de la frecuencia y ser
distinta de `arg(Delta)`. La variación de su fase exige

\[
\frac{\hbar D}{2}\nabla\cdot
\operatorname{Im}(f_n^*D_Af_n)
=\operatorname{Im}(\Delta^*f_n).
\]

Fijar `chi_n=arg Delta` anularía el segundo miembro y exigiría que cada
corriente espectral tuviera divergencia nula por separado para esa fase
común. Eso no se cumple en una textura arbitraria prescrita. En el caso
uniforme y en el vórtice radial simétrico hay simplificaciones particulares;
no se transfieren como identidad a toda la cinta 2D. Las ecuaciones
angulares con `Delta` complejo de Harwin y colaboradores, ecuaciones 1a–1b,
dan una comprobación independiente de esta separación de fases.

Al estacionar los propagadores **manteniendo Delta fijo**, el teorema de la
envolvente permite evaluar la fuerza compleja sin derivar la solución del
solver espectral en cada dirección:

\[
h_\Delta=2N_0\left[\Delta\ln(T/T_c)
+2\pi k_BT\sum_{n\ge0}
\left(\frac{\Delta}{\epsilon_n}-f_n\right)\right],
\qquad
\delta\mathcal F=\int dV\,\operatorname{Re}(h_\Delta^*\delta\Delta).
\]

No se impone `h_Delta=0` durante una consulta de fuerza fuera del equilibrio
del condensado. La corriente es

\[
\mathbf j_s=\frac{2\pi\sigma_n k_BT}{e}
\sum_{n\ge0}\operatorname{Im}(f_n^*D_Af_n).
\]

Con estas convenciones, la invariancia de calibre implica

\[
\nabla\cdot\mathbf j_s
=-\frac{2e}{\hbar}\operatorname{Im}(\Delta^*h_\Delta).
\]

La corriente superconductora sola no necesita ser solenoidal si el
condensado está fuera de autoconsistencia. Su acoplamiento a corriente
normal y potencial debe conservar la continuidad total.

## 4. Discretización compatible con Delaunay–Voronoi

La siguiente construcción es una propuesta de implementación derivada
aquí. Reutiliza el grafo de la memoria con volúmenes positivos `V_i` y pesos
`w_ij=A_dual,ij/ell_ij`; `A_dual` incluye el espesor de la película. Cada
arista no orientada se cuenta una vez. Para
`alpha_ij=(2e/hbar)*integral_i^j A·dl` y `U_ij=exp(-i*alpha_ij)`, el
término espacial discreto es

\[
\mathcal F_{h,\nabla}
=2\pi N_0k_BT\frac{\hbar D}{2}
\sum_{n\ge0}\sum_{\{i,j\}}w_{ij}
\left[|U_{ij}f_{n,j}-f_{n,i}|^2+(g_{n,j}-g_{n,i})^2\right].
\]

El potencial local anterior se suma con `V_i`. Si se usan incógnitas `u`,
primero se construyen `f(u)` y `g(u)` en los nodos y después se evalúa esta
energía. Aplicar una supuesta regla de la cadena exacta entre diferencias
discretas de `u` y `f` generaría otro esquema; las derivadas numéricas deben
salir de la energía realmente evaluada.

La transformación `f_i→exp(i*beta_i)f_i`,
`Delta_i→exp(i*beta_i)Delta_i`,
`alpha_ij→alpha_ij+beta_j-beta_i` deja la energía invariante exactamente.
Derivando cada enlace se obtiene la corriente orientada de `i` a `j`:

\[
I_{ij}=-\frac{2e}{\hbar}\frac{\partial\mathcal F_h}{\partial\alpha_{ij}}
=\frac{2\pi\sigma_n k_BT}{e}w_{ij}
\sum_{n\ge0}\operatorname{Im}(f_{n,i}^*U_{ij}f_{n,j}).
\]

En un campo uniforme con gradiente de fase pequeño recupera A.5 por
`Im(conj(f_i)*U*f_j)≈|f|²*q*ell`. La divergencia debe usar esa misma corriente,
con los mismos pesos y orientaciones. La identidad de calibre discrimina
errores de signo, de factor dos y de espesor mejor que comparar dos
fórmulas interpoladas independientemente.

Una malla que no tenga los pesos y volúmenes requeridos debe identificarse
como tal; no se reparan pesos negativos mediante recortes. En los bordes
aislantes se impone flujo espectral normal nulo. En un reservorio se fija
el propagador de su estado declarado; un contacto normal usa `f=0,g=1`,
pero no sustituye un reservorio superconductor del contrato vigente. Una
interfaz resistiva necesita su término de borde correspondiente: no se
añade por analogía con una resistencia circuital macroscópica.

## 5. Ruta mínima y reducciones que conservan el problema

1. Mantener `Delta` prescrito y temperatura uniforme en una geometría
   pequeña; resolver únicamente las incógnitas complejas espectrales por
   frecuencia. El problema admite paralelismo entre frecuencias sin
   aproximarlas como independientes espacialmente.
2. Obtener energía, fuerza y corriente de la misma expresión discreta y
   contrastar el límite uniforme, C.29–C.31, un cambio puro de calibre y el
   núcleo radial ya disponible. No se necesita un transiente para ese
   contraste.
3. Conservar unas pocas frecuencias bajas explícitas y tratar las altas
   con una expansión controlada solo donde se compruebe su error. Reutilizar
   patrones dispersos, precondicionadores e inicializaciones cercanas no
   cambia la física. El número de modos y su cola deben quedar registrados.
4. Resolver primero el espectro a `Delta` fijo es una eliminación de
   variables estacionarias; imponer simultáneamente gap autoconsistente
   sería otro problema. Para un Hessiano del funcional reducido, la
   dependencia espectral aparece mediante un complemento de Schur, no
   mediante la Hessiana con propagadores congelados.
5. Probar una fuerza de esta referencia en un integrador KWT constituye
   otro control, con la movilidad efectiva explícita. No convierte por sí
   solo KWT en una dinámica microscópica ni identifica cómo se reparte el
   calor producido.

Esta ruta permite aprovechar el núcleo geométrico de la memoria y mejorar
la representación física que se ensaya. No requiere simplificar el circuito,
cambiar el dispositivo ni reducir de antemano la dimensión espacial. La
política de parar después del disparo de Vout afecta al horizonte futuro,
no a estas ecuaciones.

## 6. Límite preciso para acoplar poblaciones no térmicas

El campo `g=sqrt(1-|f|²)` es una coordenada de la rama de **Matsubara térmica**.
No puede continuarse ingenuamente sustituyendo la frecuencia por una
energía retardada: el módulo y la conjugación compleja no son funciones
holomorfas. En el problema retardado la componente normal suele ser
compleja; la normalización involucra dos componentes anómalas relacionadas
por las simetrías pertinentes, no una cota universal `|f_R|<1`.

Fuera del equilibrio hay que resolver la rama retardada causal y las
distribuciones apropiadas, con su transporte y condiciones de borde. La
parametrización Keldysh `g^K=g^R h-h g^A` separa espectro y distribución;
no se deduce sustituyendo una temperatura equivalente en la suma Matsubara.
La revisión de Belzig y colaboradores da el marco y sus ecuaciones cinéticas.

B.11 adoptó un espectro adiabático **localmente uniforme**, con variables
de conteo y fuerzas propias. Una solución espacial del espectro puede
depender de campos vecinos y de los contactos. Se debe justificar de nuevo
la correspondencia entre conteos, energía, fuerzas y corrientes, así como
el transporte a igual energía y el trabajo espectral. Añadir la corrección
térmica espacial a la energía de B no garantiza esa correspondencia y puede
contar dos veces contribuciones del vacío o del gradiente.

Por tanto, este bloque exacto dentro de sus hipótesis es una referencia
útil inmediata y un posible componente posterior. Su promoción al modelo
no térmico exige una derivación específica; no requiere afirmar que se ha
resuelto ya toda la dinámica microscópica ni desechar la cinética desarrollada.

## 7. Siguiente cálculo físico: núcleo térmico autoconsistente estático

Una vez evaluada la campaña de campos prescritos, el paso útil es dejar que
el propio funcional determine la amplitud y la fase interiores del
condensado. Esto produce un núcleo y una corriente que responden al mismo
problema térmico espacial. Todavía no introduce tiempo físico, movilidad
KWT, fotón, poblaciones cinéticas ni una barrera de nucleación.

### La actualización del gap es un mínimo cuadrático exacto

En las unidades de `implementation_options.md`, `d=Delta/(kB*Tc)`,
`t=T/Tc` y `epsilon_n=2*pi*t*(n+1/2)`. Con los propagadores y las trazas
de borde fijos, la dependencia del funcional de corte finito en cada nodo
interior es

\[
F_N(d;f,g)=\sum_i m_i\left[C_N|d_i|^2
-4\pi t\operatorname{Re}\!\left(d_i^*\sum_{n<N}f_{ni}\right)\right]
+\text{términos independientes de }d,
\quad
C_N=\ln t+\sum_{n=0}^{N-1}\frac1{n+1/2}.
\]

Para `C_N>0`, su mínimo único en ese bloque es

\[
d_i^*=\frac{2\pi t}{C_N}\sum_{n<N}f_{ni},\qquad
F_N(d;f,g)-F_N(d^*;f,g)
=\sum_{i\;\mathrm{interior}}m_iC_N|d_i-d_i^*|^2.
\]

Los cortes 128 y 256 a `t=0,9/8,65` tienen `C_N>0`; la implementación debe
comprobar esta condición antes de emplear la actualización. La suma conserva
la parte compleja completa: no proyecta la fase a un vórtice radial ni fija
un cero interior. El gradiente integrado del gap es precisamente
`G_i=2*m_i*C_N*(d_i-d_i*)`.

Alternar una resolución espectral con `d` fijo y esta actualización del
gap busca las ecuaciones estacionarias del **mismo** `F_N`. No es un paso de
tiempo ni una nueva aproximación de cierre. El bloque del gap no necesita
mezcla para asegurar descenso: la identidad anterior lo demuestra. Para
el bloque espectral, el descenso depende de que el solver entregue su rama
estacionaria mediante pasos que reduzcan la energía, como exige su búsqueda
actual. Esa propiedad no demuestra un mínimo global ni una convergencia
rápida de la alternancia. Pueden existir ramas o estados metastables.

No se agrega una cola únicamente a la fuerza. Una corrección espectral debe
provenir de una energía común; puede cambiar el problema cuadrático del gap
y dejar de permitir esta actualización local exacta.

### Campaña concreta, conservando las geometrías existentes

Se proponen cuatro problemas: `65×65/N=128`, `65×65/N=256`,
`129×129/N=256` y un perfil inicialmente asimétrico `65×65/N=256`,
todos sobre el cuadrado existente de
`12*ell0` de lado, con `T=0,9 K`, `Tc=8,65 K` y `D=0,5 cm²/s`.
Se mantienen las trazas complejas de `d` y las trazas espectrales radiales
del ensayo prescrito anterior. Por tanto, sólo se libera el condensado
interior. Es un control local con borde sujetado, no el dispositivo completo
de 80 nm ni un nuevo reservorio autoconsistente.

El contorno del condensado mantiene una vuelta de fase. Esa condición
permite estudiar un estado con vorticidad impuesta; no modela su entrada
desde un borde ni determina una barrera de activación. Las trazas BVP usadas
en el borde ya no son una solución de referencia para el interior cuando
se relaja `d`. La comparación física será entre los núcleos obtenidos, no
contra el perfil `tanh` inicial.

El perfil asimétrico ya disponible es una inicialización útil para no
ocultar las libertades de fase 2D mediante una simetría inicial perfecta.
No se conserva esa deformación como restricción durante la relajación.
Las soluciones gruesas y de menor corte pueden inicializar las siguientes;
eso no cambia sus ecuaciones ni sus bordes. Cada corte requiere su propia
relajación: truncar a 128 modos un estado obtenido con 256 no lo convierte
en una solución autoconsistente del corte menor.

Un único pool puede resolver todas las frecuencias pendientes, compartido
entre casos. Dentro de cada caso se termina la suma de una iteración antes
de modificar `d`. No se actualiza el gap con una mezcla de frecuencias
calculadas sobre condensados distintos. El coste de una iteración completa,
la memoria y la contracción observada fijarán el presupuesto y la ETA;
no se promete de antemano un número pequeño de iteraciones.

### Cierre operativo y resultados que permiten decidir

Como objetivo inicial de la relajación se propone

\[
r_d=\frac{\left[\sum_{i\in\mathcal C}
m_i|d_i-d_i^*|^2\right]^{1/2}}
{d_{\mathrm{ref}}\left[\sum_{i\in\mathcal C}m_i\right]^{1/2}}
\le10^{-3},
\]

donde $\mathcal C$ contiene únicamente los nodos interiores libres con
$r\le4\ell_0$, y `d_ref` es la amplitud homogénea de referencia empleada al construir
los mismos bordes. Es un objetivo numérico inicial, no una precisión física
adquirida, una condición universal ni una declaración de estacionariedad
en todo el dominio. Se registra además el máximo nodal global y
la evolución de los observables para detectar un defecto localizado que
la norma promediada pudiera ocultar. Sólo se exige mayor precisión si el
residuo de relajación impide interpretar las diferencias entre cortes o
entre las dos mallas. La tolerancia espectral ya usada se conserva como
punto de partida; su residuo y cualquier fallo se informan.

El residuo de gap se calcula **después** de resolver los espectros para el
`d` que se va a informar. Si ya satisface el objetivo, se guarda ese par
coherente `d,f,g`; no se presenta `d_new` junto con los espectros de `d_old`
como una solución final. Un presupuesto agotado produce un estado
incompleto con checkpoint; no equivale a convergencia ni provoca un
reinicio automático.

Los productos centrales del cálculo serán:

- Mapa de amplitud y fase, localización del núcleo y vueltas de fase,
  mostrando si la solución se desplazó o desarrolló más estructura.
- Radios a media amplitud respecto de `d_ref`, medidos desde el núcleo
  en distintas direcciones, y curva de recuperación de amplitud. No se
  fuerza una media radial cuando oculta una asimetría apreciable.
- Corriente circulante y sus perfiles, con longitud en nm y corriente
  normalizada o en SI sólo cuando se declare toda la escala material.
- Energía `F_N`, disminución respecto del estado inicial del mismo caso,
  residuo espectral, residuo de gap e identidad de corriente interior.
  La energía no se etiqueta como barrera ni como energía de nucleación.
- Cambios de esos observables al variar el corte y la malla existentes,
  separados del error de la relajación. Sin una cola consistente, el
  resultado sigue siendo una familia de sumas finitas explícitas.

Esta campaña entrega un objeto físico adicional —el núcleo térmico
determinado por el funcional— y permite decidir si la representación
espacial merece acoplarse a una dinámica posterior. No completa por sí
sola los requisitos cinéticos ni el transiente del dispositivo.

## Fuentes cotejadas

- Memoria local, anexos C y D: `tmp/pdfs/modelo_v0_2/sources/memoria_02.pdf`,
  pp. PDF 152–165 para la relación condensado/corriente y su discretización.
- [Allmaras, tesis 2020](https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf),
  apéndice B.3, pp. impresas 241–244: representación compleja y reorganización
  del término temporal. La copia textual local fue cotejada.
- [Virtanen, Vargunin y Silaev, preprint v1](https://arxiv.org/pdf/1909.00992v1),
  ecuación 20, p. 3; [publicación 2020](https://doi.org/10.1103/PhysRevB.101.094507).
  La especialización y las transformaciones de coordenadas anteriores se
  exponen como desarrollo algebraico de esta nota.
- [Harwin y colaboradores, 2018](https://arxiv.org/pdf/1805.06215), ecuaciones
  1a–1b y 2: amplitud, fase espectral y gap complejo. Se usan como contraste
  de la estructura; no se importan su geometría TES ni su modelo de lectura.
- [Belzig y colaboradores, 1999](https://arxiv.org/pdf/cond-mat/9812297),
  sección de sistemas fuera del equilibrio: propagadores retardados/avanzados,
  componente Keldysh y distribuciones. No se importan nuevas hipótesis
  cinéticas a producción en esta entrega.
