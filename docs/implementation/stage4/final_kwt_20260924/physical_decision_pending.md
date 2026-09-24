# Consulta física en investigación: cuándo puede ajustarse instantáneamente la distribución de carga

La etapa 4 permanece abierta. El control térmico sobre la malla dual y el
control longitudinal débil permiten avanzar sin repetirlos. Ante la consulta
sobre **conservar la evolución del desequilibrio electrón–hueco o aproximar
su ajuste como instantáneo**, se solicitó investigar primero si la física
del dispositivo justifica esa separación de tiempos. No se ha elegido
automáticamente ninguna de las dos rutas ni cambiado el modelo ejecutable.

La [evaluación cuasiclásica y de escalas](quasiclassical_assessment.md)
distingue cuatro aproximaciones: cuasiclásica, difusiva, espectro adiabático
y eliminación temporal del modo de carga. Recomienda conservar las dos
primeras y las poblaciones dinámicas, aprovechando el espectro adiabático
como aproximación separada. No exige un solver de dos tiempos completo.
Tampoco resuelve automáticamente el trabajo del espectro móvil ni los
términos de calibre pendientes. La elección concierne a la física; no es
una preferencia entre integradores.

## Qué representa la variable consultada

El modo de energía hL describe la ocupación simétrica de las excitaciones.
El modo de carga hT describe la diferencia de ocupación entre las ramas de
tipo electrón y de tipo hueco. No es por sí solo el voltaje, la carga de un
capacitor ni otra corriente del circuito. Una película eléctricamente neutra
puede tener un desequilibrio de ramas, y su relación con el potencial y la
fase debe imponerse de forma consistente.

Esta distinción aparece en [Golub, ecuaciones 9–10](https://www.jetp.ras.ru/cgi-bin/dn/e_044_01_0178.pdf),
donde el momento de distribución y el potencial invariante de calibre
intervienen conjuntamente en la densidad de carga. Su aproximación cerca de
Tc no se importa como calibración para NbN a 0,9 K. La expansión temporal
de los productos cuasiclásicos requiere conservar sus términos de fase y
potencial; véase [Larkin y Ovchinnikov](https://www.jetp.ras.ru/cgi-bin/dn/e_041_05_0960.pdf).

## Las dos rutas de la consulta

| Ruta | Significado físico | Trabajo concreto después de elegirla |
|---|---|---|
| Conservar la evolución de hT | Las distintas energías pueden tardar en ajustar su desequilibrio; se retiene esa memoria durante el transiente. | Unir la ecuación cinética temporal con los términos del gap y potencial móviles, la neutralidad y las condiciones de los contactos. Reutilizar los operadores espaciales existentes, sin confundir hT con capacitancia. |
| Ajuste instantáneo espectral provisional | En cada instante se calcula la distribución de carga completa, suponiendo que se adapta antes de que cambien apreciablemente las demás variables. | Resolver la reducción algebraica por energía y declarar su banda de validez. Conservar la forma energética completa y comprobar la separación de tiempos en el régimen que realmente se utilice. |

**Ambas rutas conservan el oráculo espectral espacial.** Elegir la segunda no
autoriza volver al catálogo local, eliminar todo hT ni aproximarlo mediante
una temperatura o un único potencial. Tampoco autoriza cambiar el circuito
de la memoria, la geometría o los tiempos materiales.

La [revisión de momentos](../moment_review_20260924/physics_decision.md)
distinguió dos cuestiones que no deben confundirse. Comprimir la forma
energética de hT a una sola dirección térmica produjo diferencias resueltas
en corriente y torque del ensayo de vórtice congelado. En cambio, la
eliminación algebraica **con todas las energías** reprodujo bien el diagnóstico
dinámico lento. Éste sólo respalda la banda lenta comprobada, no la respuesta
de pocos picosegundos de Korzh. El ensayo de vórtice congelado tampoco descarta
por sí solo la fase y el potencial móviles del modelo de la memoria.

## Lo que esa respuesta no resuelve automáticamente

1. **El trabajo de un espectro que cambia.** B.46–B.51 deriva fuerzas a partir
   de un catálogo local y de sus ocupaciones a conteo fijo. El nuevo espectro
   espacial depende del entorno. La [revisión del puente no térmico,
   P.4](../self_consistent_review_20260924/bridge_review.md) declara que la
   fuerza Keldysh no térmica espacial no está acreditada como derivada de una
   energía escalar para cualquier distribución. La elección temporal de hT
   no completa ese cambio de representación.
2. **La conexión entre fase, potencial y neutralidad.** El diagnóstico
   `frozen_charge_response.py` se construyó con gap y espectro congelados.
   Añadir su masa DOS al sistema móvil no proporciona los términos de
   calibre que faltan. La identidad de conversión de carga ayuda a verificar
   la conexión; no la determina por sí sola.
3. **El calor de orden cuadrático.** B.41 ya especifica una deposición
   espectral efectiva de la potencia KWT y de Joule. Se conserva esa elección,
   pero su potencia debe proceder de un balance común con el trabajo
   reversible. Añadir a una gráfica el calor perdido por una perturbación
   no demuestra que se haya transferido a las poblaciones dinámicas.
4. **Los datos materiales y la transferencia fotónica.** La normalización de
   tasas, el tiempo cinético material y la preparación de la cascada siguen
   sus decisiones de 3.5. No se deducen de una tolerancia más flexible.

## Qué aportan directamente las fuentes de implementación

[Vodolazov 2017, ecuación 1](https://arxiv.org/pdf/1611.06060) contiene el
trabajo radial local asociado al cambio de amplitud. Su cascada inicial usa
gap fijo y sus modelos espaciales posteriores utilizan aproximaciones
térmicas. No proporciona por simple copia la unión espacial no térmica
general que aquí falta.

El [solver.m de los autores de 2025](https://github.com/qnngroup/proj-KE-solver/blob/main/solver.m)
lee espectros de archivos, mantiene el gap y evoluciona las ocupaciones con
Euler. Después calcula la corrección espectral de fuerza. Es una fuente útil
de cinética y convenciones; no contiene el ensamblaje de ambos modos de
distribución, condensado móvil y circuito de este proyecto.

La reutilización práctica consiste en conservar los operadores y códigos
probados, y completar sólo las uniones que corresponden a la ruta elegida.
No corresponde anunciar una implementación final mediante la suma de los
controles independientes, ni solicitar otra campaña larga antes de que ese
ensayo tenga ecuaciones y observables definidos.
