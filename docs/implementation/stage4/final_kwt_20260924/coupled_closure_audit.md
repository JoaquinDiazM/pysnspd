# Unión no térmica: qué está especificado y qué elección física falta

Esta revisión evita convertir otra corrida en sustituto de una ecuación que
todavía no está definida. No propone otro barrido estático ni un integrador.

## Lo que ya está determinado

- La geometría dual, el paso KWT heredado y el circuito de tres estados de la
  memoria tienen implementaciones concretas.
- B.41 ya seleccionó el reparto espectral efectivo del calor. No hace falta
  volver a preguntar cómo distribuir QDelta: debe reutilizarse esa elección
  con sus unidades y sus poblaciones, una sola vez.
- En la referencia real uniforme, el trabajo radial que entra en la ecuación
  longitudinal queda fijado por el puente recíproco. El control dinámico ya
  realizado lo verifica y distingue disponibilidad de energía interna.
- La fuerza y la corriente térmicas de la acción Usadel espacial son comunes
  y están implementadas. El cambio de malla no exige cambiar esa acción.

## La unión general no es sólo conectar estas funciones

B.46–B.51 define la energía adiabática fuera del equilibrio mediante niveles
de un catálogo **local** y ocupaciones a conteo fijo. El nuevo oráculo espacial
resuelve espectros influidos por los vecinos. La revisión previa del
[puente no térmico, P.4](../self_consistent_review_20260924/bridge_review.md)
advierte que su fuerza Keldysh general no está acreditada como derivada de una
energía escalar para una distribución espacial arbitraria. Por tanto no se
puede afirmar que B.46 haya pasado sin cambios al nuevo oráculo.

Existe además una diferencia entre responder a una perturbación con gap fijo
y evolucionar el condensado. El operador de carga congelado satisface una
identidad de la forma

```
div(delta I_s) + Im(conj(d) delta G) = 2 integral(r_T dE).
```

En el diagnóstico estático se impone rT=0. Si se lo sustituye por una masa de
almacenamiento y una derivada temporal de hT, su segundo miembro ya no es
cero. Agregar ese estado a la ecuación de potencial óhmico de la memoria no
proporciona automáticamente la neutralidad ni los términos de trabajo
asociados al potencial y a la fase móviles. `frozen_charge_response.py`
declara expresamente que no contiene esa conexión.

La distinción tiene antecedente publicado: [Golub, ecuaciones 9–10](https://www.jetp.ras.ru/cgi-bin/dn/e_044_01_0178.pdf)
separa el momento de la distribución de carga del potencial invariante de
calibre y aplica neutralidad. Sus aproximaciones próximas a Tc no son una
calibración automática para la película a 0,9 K. La cinética completa usa
productos temporales cuya expansión adiabática exige retener coherentemente
los términos de potencial y gap; véase [Larkin y Ovchinnikov](https://www.jetp.ras.ru/cgi-bin/dn/e_041_05_0960.pdf).

## Qué se puede reutilizar de las fuentes propuestas

[Vodolazov 2017, ecuación 1](https://arxiv.org/pdf/1611.06060) aporta el trabajo
radial local proporcional a la velocidad de la amplitud y a la derivada de
la ocupación respecto de energía. La cascada inicial se estudia con el gap
fijo; para la dinámica espacial posterior se presentan aproximaciones
térmicas. No es una implementación del cierre espacial no térmico general
que aquí se intentaría añadir.

El [solver.m del repositorio de 2025](https://github.com/qnngroup/proj-KE-solver/blob/main/solver.m)
carga la DOS y R2 de archivos, fija el gap, evoluciona ocupaciones con Euler
y calcula después la corrección de fuerza Phi. Es reutilizable para cinética,
datos y convenciones. No contiene una trayectoria acoplada de un gap espacial
móvil con ambos modos de distribución y el circuito de este proyecto.

Estas observaciones describen el alcance visible de las fuentes consultadas;
no sostienen que no exista otra formulación publicada utilizable.

## Dos rutas físicas distintas

| Ruta | Qué se conserva y qué se aproxima | Consecuencia |
|---|---|---|
| Ajuste instantáneo espectral de la distribución de carga | Se conserva el espectro espacial y se resuelve algebraicamente la forma energética completa de hT. No se la comprime a un solo potencial ni se vuelve al catálogo local. | Menos estados temporales. Su justificación requiere una separación de tiempos que sólo se ha comprobado en el diagnóstico lento. No resuelve por sí sola el trabajo de un espectro móvil. |
| Evolución espacial de la distribución de carga | Se conserva el espectro espacial y la distribución de carga guarda su evolución temporal, junto con la de energía. | Requiere completar los términos dinámicos de calibre, neutralidad y trabajo del espectro antes de incorporar esa capacidad al detector. No exige por definición equivalencia microscópica exacta, pero sí declarar una expansión temporal y un cierre consistentes. |

El contraste previo del 20,55 % de corriente y 77,53 % de torque se obtuvo
al congelar un vórtice y comprimir la forma energética de la distribución de
carga a una sola dirección térmica. **No rechaza la eliminación algebraica
que conserva todas las energías.** Tampoco acredita esa ruta en una hotbelt, un vórtice o
un transiente ultrarrápido. Una geometría suave elimina algunas singularidades
geométricas, pero no los términos de carga y trabajo que dependen del estado.

La elección es física, no una preferencia entre Euler y otro método ni una
tolerancia más o menos exigente. Las fuentes y documentos vigentes no
determinan inequívocamente qué aproximación autoriza la sustitución del
catálogo local por el oráculo espacial en todas las ecuaciones no térmicas.

## Por qué no basta otra prueba débil a corriente cero

Alrededor del equilibrio uniforme sin corriente, amplitud/energía y
fase/carga se desacoplan a primer orden por simetría. Superponer sus controles
produce otro control lineal legítimo, pero no verifica los términos cruzados
que aparecen a segundo orden o con polarización. El calor disipativo también
es de segundo orden. Añadir su integral al registro no equivale a evolucionar
la energía interna de las poblaciones.

Por ello no se propone otra cantidad de consultas espectrales como condición
de cierre. Para un ensayo general útil debe estar seleccionada la ruta
constitutiva y especificada su unión; el siguiente cálculo entonces comprobará
esa unión en un sistema suave y polarizado. No hará falta volver a certificar
el núcleo térmico ni el control longitudinal ya aprobados.
