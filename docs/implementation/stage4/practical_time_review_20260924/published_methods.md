# Ruta numérica práctica después del transiente térmico

Revisión del 24 de septiembre de 2026. Fuentes primarias consultadas en esa
fecha; los archivos numéricos históricos no se modifican. Esta nota selecciona
el siguiente trabajo y documenta las conversiones necesarias. No promueve el
modelo experimental a producción.

## Decisión

La integración no lineal llegó a 1 ps. La interrupción posterior corresponde
al certificado comparativo: 95 de sus 96 comprobaciones pasan y falla una
señal cruzada casi extinguida, el torque de fase inducido por la sonda de
amplitud. No es un fallo de Newton, de Krylov ni de continuidad de la solución.
El diagnóstico cuantitativo y los checkpoints recuperados de esta entrega
conservan el criterio original y su fallo; no se cambia retrospectivamente a
un pase.

No se solicita otra campaña térmica hasta 1 ps ni un nuevo integrador para
perseguir esa cola. La evidencia ya justifica trasladar la acción y sus flujos
a la geometría Delaunay–Voronoi que se utilizará en el dispositivo. Se mantiene
el transiente térmico completado como referencia de desarrollo con precisión
limitada de esa señal secundaria. El siguiente control debe resolver una
pregunta nueva: si los mismos operadores conservativos funcionan en la malla
de producción y con sus medidas físicas correctamente convertidas.

## Qué se reutiliza de implementaciones publicadas

La [documentación oficial de pyTDGL](https://py-tdgl.readthedocs.io/en/latest/background.html)
describe volúmenes finitos Delaunay–Voronoi, enlaces covariantes y un paso de
Euler con resolución algebraica local de la amplitud. Su fuerza espacial se
evalúa en el estado anterior: no es Heun ni una discretización temporal de
segundo orden. El
[código oficial](https://github.com/loganbvh/py-tdgl/blob/main/tdgl/solver/solver.py)
y el [artículo del autor](https://arxiv.org/abs/2302.03812) permiten comprobar
esa distinción. La publicación del método respalda su reutilización; no
acredita automáticamente los cierres adicionales de este modelo a 0,9 K.

Ya existen en este repositorio las piezas necesarias:

| Pieza | Implementación existente | Uso inmediato |
|:--|:--|:--|
| Triangulación y dual | `pysnspd/mesh/pytdgl_like.py` y `mesh/finite_volume/` | Generar la geometría con el backend que se prevé usar. |
| Flujos y derivadas | `pysnspd/gtdgl/tdgl_operators.py` | Reutilizar áreas, longitudes duales y orientación de aristas. |
| Acción espacial Usadel | `pysnspd/experimental/thermal_spatial_usadel.py` | Cambiar únicamente la construcción de `ThermalGraph`; mantener acción, raíces, fuerza y corriente. |
| Paso local KWT | `pysnspd/solver/core.py:solve_for_psi_squared` | Su argumento `forcing_dimensionless` ya sustituye toda la fuerza GL. La conversión de abajo permite estudiar su reutilización sin sumar nuevamente el gradiente antiguo. |
| Circuito de tres estados | `pysnspd/circuit/readout.py` y adenda `circuito_memoria_20260922.md` | Conservar topología, parámetros y tres variables cuando se conecten los puertos. |

El oráculo térmico actual ya usa un grafo conservativo de volúmenes de control
rectangulares; no requiere reescribir Usadel como otro método espacial. El
primer banco GLL de la etapa 4 y este grafo son antecedentes distintos.

## Mapa de la malla: qué cambia y qué permanece

Sean `A_i` el área física del volumen de control, `l_ij` la longitud de una
arista Delaunay y `s_ij` la longitud de su cara Voronoi. Con
`ell0 = sqrt(hbar D / (2 kB Tc))`, el grafo de la acción recibe

\[
m_i=A_i/\ell_0^2,\qquad c_{ij}=s_{ij}/l_{ij},\qquad
\mathbf X_i=\mathbf r_i/\ell_0.
\]

`m_i` integra el área adimensional y `c_ij` multiplica el flujo entre dos
nodos. No se debe introducir otra área al evaluar la fuerza integrada.
La misma arista aporta signos opuestos a las dos celdas vecinas: la cancelación
de flujos interiores es algebraica. El contorno físico se integra mediante
sus caras y sus condiciones de borde; los contactos fijos del ensayo de núcleo
no se convierten por ello en contactos del dispositivo.

La dual permite adaptar la resolución a bordes y zonas de interés, evita una
divergencia construida con pesos incompatibles y comparte la geometría con
producción. No elimina el error constitutivo de un cierre local, no determina
la preparación fotónica ni resuelve el reparto de trabajo espectral. Tampoco
elimina por sí sola la rigidez temporal: al reducir la separación de nodos
aparecen escalas difusivas más cortas. Un campo inicial suave puede reducir su
excitación, pero los núcleos y el frente del hotbelt deben seguir resolviéndose.

## Conversión exacta para reutilizar el paso KWT de la memoria

Este apartado deriva una compatibilidad de ecuaciones continuas a temperatura
impuesta. No afirma que Euler y ETD2 den el mismo resultado con un paso finito.
Define cada escala porque pasar el RHS entre ambos códigos sin convertirlo
cambiaría la movilidad o aplicaría esta dos veces.

El oráculo usa `d = Delta/(kB Tc)`, el tiempo `tau = t/tD` con
`tD = hbar/(2 kB Tc)` y el potencial `v = 2 e phi/(kB Tc)`. Su gradiente complejo
integrado `G_i` se define por
`delta F = Re sum_i conj(G_i) delta d_i`. Para `theta = T/Tc`, definimos

\[
C_T=\frac{\pi}{2}\sqrt{\frac{1+\theta}{2}},\qquad
\kappa=4\left(\frac{k_B T_c\tau_\psi}{\hbar}\right)^2,
\qquad S_i=\sqrt{1+\kappa|d_i|^2}.
\]

`tau_psi` es el tiempo físico de relajación KWT conservado por el modelo;
`ThermalKWTNormal` lo obtiene sumando las tasas electrón–electrón y
electrón–fonón heredadas. Su inversión local equivale a

\[
\frac{C_T}{S_i}\left[
\left(\partial_\tau+\frac{i v_i}{2}\right)d_i
+\frac{\kappa}{2}d_i\partial_\tau|d_i|^2\right]
=-\frac{G_i}{m_i}.
\tag{P.1}
\]

La producción usa `psi = Delta/Delta0`, `t_prime = t/tau0`, donde
`tau0 = pi hbar/(8 kB Tc)`. Sea `d_star = Delta0/(kB Tc)`, una escala fija de
normalización del material, distinta del gap local y del valor térmico del
ensayo. Sustituir `d = d_star psi` y aplicar la regla de la cadena en P.1 da
exactamente el formato del paso local con

\[
u=1,\qquad \gamma^2=\kappa d_\star^2,
\qquad \mu_i=\frac{\tau_0}{2t_D}v_i,
\qquad
\texttt{forcing}_i=-\frac{\tau_0}{t_D}
\frac{G_i}{m_i C_T d_\star}.
\tag{P.2}
\]

En particular, `gamma = 2 Delta0 tau_psi/hbar`, `tau0/tD = pi/4` y
`mu = 2 e phi tau0/hbar`, coherentes con las escalas ya presentes en
`material.py` y `mesh/device.py`. La cantidad de P.2 es la fuerza **antes**
de invertir la movilidad. No es `velocity`, ni lleva el término de potencial
que ya aplica el enlace temporal. Sustituye por completo el bracket antiguo;
no se suman la reacción GL, un gradiente `K0` o la corrección Usadel–GL histórica.

La corriente también debe proceder de la misma acción. La ecuación de Poisson
de producción admite un ensamblador reutilizable, pero su corriente nativa
`Im(conj(psi) grad psi)` no representa automáticamente la nueva corriente
espectral. En el control actual la caída física es
`phi_i - phi_j = (kB Tc/(2e)) (v_i-v_j)`; para conectar el circuito debe usarse
el voltaje pasivo de puerto y la corriente **total**, siguiendo la adenda de la
memoria. No se agrega al circuito la inductancia de la región ya resuelta.

Estas conversiones permiten una comprobación ligera y concreta antes de
promover el paso heredado: comparar el límite de paso pequeño de la actualización
local con `ThermalKWTNormal.response` en un mismo estado, con potencial y
movilidad no nulos. Una sola comparación de este tipo descubre factores de dos,
signos y dobles movilidades; no exige otro transiente completo de referencia.

## Trabajo que ya no es necesario para el siguiente avance

No se abre ahora una competición de integradores. Los algoritmos publicados
[phipm](https://arxiv.org/abs/0907.4631) y
[KIOPS](https://arxiv.org/abs/1804.05126), con
[código de sus autores](https://gitlab.com/stephane.gaudreault/kiops), ofrecen
adaptación de las acciones exponenciales si su coste se vuelve un problema
medido. El fallo observado no fue de esa clase. Asimismo,
[BDF de SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.BDF.html)
es una opción establecida, pero su interfaz requiere una matriz jacobiana o
aproximarla. Eliminar Usadel y Poisson produce una respuesta espacial global;
inventar una jacobiana reducida con solo vecinos sería incorrecto. Introducir
otro backend ahora añadiría trabajo sin resolver el motivo del certificado.

Tampoco se necesitan más extremos artificiales para admitir el traslado de
malla. Son suficientes las identidades discretas de unidades, derivadas y
flujos; después, un estado suave sobre la geometría relevante y una comparación
de resolución de los observables que se usarán. La energía libre térmica medida
es una referencia útil, pero no sustituye el balance de energía interna cuando
se evolucionen poblaciones. El acoplamiento no térmico permanece como trabajo
físico explícito de la etapa 4, no como otra tolerancia a apretar.

## Contrato de las próximas figuras

Cada figura debe declarar en su título o pie: magnitud física, fórmula o
procedimiento de extracción, unidad, normalización y referencia sustraída.
Una norma espacial debe identificar su peso de área; una diferencia relativa,
su denominador. Debe distinguirse la sonda de amplitud de la sonda angular y
el torque que cada una induce. La cola cercana a cero se presenta también en
escala absoluta: un porcentaje grande de una cola minúscula no se confunde
con un cambio grande del condensado. Los mapas compartirán escala cromática
cuando se comparen tiempos. Ninguna curva térmica del núcleo se rotulará como
pulso, latencia o `V_out` del dispositivo.
