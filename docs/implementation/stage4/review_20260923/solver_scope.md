# Alcance del núcleo numérico y siguiente contraste físico

Revisión del 23 de septiembre de 2026. Esta nota distingue el solver de
producción de los controles experimentales de etapa 4A. No cambia código ni
elige el backend definitivo para los transientes del nuevo modelo.

## Qué conserva producción

La generación de malla y los operadores de producción conservan la estructura
adaptada de pyTDGL: triangulación Delaunay, áreas de control Voronoi, gradientes
en aristas y divergencias ponderadas por sus longitudes duales. Las variables
de enlace mantienen las derivadas covariantes. Es una adaptación dentro de
pySNSPD, con cierres propios, no la ejecución intacta de toda la biblioteca
pyTDGL. Véanse `pysnspd/mesh/pytdgl_like.py`,
`pysnspd/gtdgl/tdgl_operators.py` y `pysnspd/solver/core.py`.

El condensado conserva la actualización algebraica local de pyTDGL,
`solve_for_psi_squared`, y el control `adaptive_euler_step`. La documentación
oficial llama al método *implicit Euler*. Su descripción precisa aquí es
**Euler de primer orden con tratamiento implícito local de la amplitud y
fuerza espacial explícita**. Se resuelve una ecuación algebraica para
`|psi_new|^2`; no hay dos evaluaciones de fuerza que constituyan Heun o punto
medio. En el límite `gamma=mu=0`, el código se reduce exactamente a
`psi_new = psi + dt*forcing(psi)/u`, Euler explícito. La adaptación del paso
no aumenta el orden del método. El controlador heredado de pyTDGL tiene
además ajustes propios, como el límite de crecimiento del siguiente paso.

El circuito sí tiene un paso RK2 de punto medio en
`pysnspd/circuit/readout.py:step_circuit_rk2`. Su voltaje de entrada se mantiene
fijo durante ese paso. El controlador de transientes encadena el bloque
mesoscópico y el circuito; ese RK2 aislado no acredita orden dos del problema
acoplado completo.

Fuentes primarias: [método espacial y temporal de pyTDGL](https://py-tdgl.readthedocs.io/en/latest/background.html),
[implementación oficial del solver](https://github.com/loganbvh/py-tdgl/blob/main/tdgl/solver/solver.py).
El orden afirmado arriba se comprueba también directamente en la actualización
local de este repositorio; no se infiere del nombre de una función.

## Qué usaron las campañas 4A

Los controles espaciales usan `RectangularSpatialFunctional`: elementos
nodales Legendre–Gauss–Lobatto (GLL), con cuadratura tensorial en un rectángulo
2D y valores independientes en sus cuatro bordes. Esta representación
experimental **es diferente de Delaunay–Voronoi**. La energía, la fuerza y la
corriente se obtienen del mismo funcional discreto.

Las campañas ejecutadas prescriben campos y poblaciones. Calculan energía,
derivadas, estabilidad y disipación instantánea; no integran un transiente.
La velocidad KWT evaluada en un estado no equivale a haber seguido una
trayectoria. El inicio de 4A (`da1d9cb`) no modificó los módulos de producción
`solver`, `mesh` ni `gtdgl`. Tampoco adoptó un nuevo integrador temporal.

Conviene conservar ambas afirmaciones explícitas en los informes: la
producción mantiene su núcleo heredado y el banco experimental ya tiene una
discretización espacial diferente. La futura promoción requiere comparar
observables y costes; no queda decidida por estos controles estáticos.

## Interpretación de los siguientes controles

Para un perfil prescrito, la comparación entre su derivada analítica y la
derivada de la malla permite identificar una excitación espuria del gradiente
antes de atribuir un signo negativo de D.36 al modelo continuo. Deben
conservarse ambos resultados: el obtenido por la malla y el de la referencia
analítica. No se reemplaza silenciosamente uno por el otro ni se recorta el
signo. Refinar el campo de amplitud reducida es pertinente si modifica ese
diagnóstico; repetir únicamente el balance energético no resuelve el error
espacial.

Fijar los grados de libertad del borde permite separar la disipación interior
de la respuesta del borde libre. La reacción que cancela exactamente la
fuerza en un borde inmóvil es la reacción de esa restricción de Dirichlet.
Como su velocidad es cero, su trabajo instantáneo es cero. Esto no constituye
una referencia independiente para la fuerza ni convierte el borde en un
reservorio superconductor con inyección de corriente. La reacción, el trabajo
y la geometría restringida deben declararse al presentar el resultado.

## Un contraste físico independiente y barato

Se realizó como contraste nuevo **la rigidez de una modulación
infinitesimal de amplitud**, a temperatura impuesta y corriente de fondo nula,
con la referencia espacial Usadel linealizada de C.29–C.31. No necesita un
vórtice, una cascada fotónica ni adoptar una nueva malla de producción.

Se usa el estado uniforme térmico autoconsistente para `Tc=8.65 K` y
`T=0.9 K`. Con `t=T/Tc`, `d=|Delta|/(kB*Tc)`, `epsilon_n=pi*t*(2n+1)`,
`E_n=sqrt(epsilon_n^2+d^2)` y `kbar=k*ell0`, se compara

\[
H_U(k)-H_U(0)=4\pi t\sum_{n\ge0}
\frac{(\epsilon_n/E_n)^2\,\bar k^2}{E_n(E_n+\bar k^2)},
\qquad
H_{K_0}(k)=H_U(0)+\frac{\pi}{2}\bar k^2.
\]

Aquí `H` es la variación de la fuerza respecto de la amplitud de la pequeña
modulación. La primera expresión permite que el espectro responda
espacialmente; la segunda es el gradiente local del modelo candidato.
El mismo estado y la misma normalización deben emplearse en ambas.

La serie corta `kbar=0.15, 0.30, 0.50, 1.00` muestra
el crecimiento del sesgo. El resultado principal es el cociente de rigideces.
Si se muestran tiempos, debe aplicarse la **misma** movilidad KWT a las dos
respuestas y llamarlos tiempos lineales condicionales. No son latencias del
detector ni una validación de las constantes de relajación. La suma puede
comprobarse duplicando una vez el número de frecuencias; no requiere una
campaña de transientes. Las longitudes usadas deben permanecer dentro del
régimen difusivo de la referencia; la falta de un libre recorrido confirmado
impide convertir esta serie por sí sola en un margen experimental.

Este contraste utiliza el funcional térmico a temperatura impuesta. No debe
confundirse con una segunda variación a poblaciones congeladas: cambiar de
conjunto de restricciones cambia la susceptibilidad. En particular, una
perturbación puramente real a fase constante tiene `q_delta=0`; no diagnostica
la singularidad de fase del núcleo ni valida el valor de `delta`.

El contraste existe como antecedente documental en
`docs/modelo_v0_4/C_condensado_corriente_y_senal_v0_4.md`, C.8. Sus cifras
anteriores no cuentan como una nueva ejecución de etapa 4. La nueva repetición
independiente queda en `linear_usadel_reference.md` y su JSON: el exceso de
rigidez local crece de 0.763 % a 30.881 % en esa serie. Este sesgo físico
puede persistir aunque el refinamiento elimine una inestabilidad aparente.

## Lo que sigue pendiente

D.4.4 pide contrastar núcleo y disipación frente a una referencia física, no
solo identidades del candidato. Ni positividad de D.36 ni convergencia de una
malla valida barreras, eventos de fase, movilidad KWT o reparto de calor. El
ensayo lineal anterior acota una parte del error de gradiente; los perfiles de
núcleo, la cinética absoluta y los primeros eventos requieren trabajo
adicional. Los transientes exigirán además condiciones de borde realmente
usadas, el circuito de la memoria, transporte y balance acoplado, y
convergencia temporal medida. La etapa 4 completa no puede declararse cerrada
solo con las campañas estáticas de 4A.
