# Paso heredado y aceptación práctica, revisión 2

24 de septiembre de 2026. Esta revisión aplica la instrucción de adoptar
tolerancias útiles para el dispositivo y cerrar discrepancias pequeñas en vez
de acumularlas como trabajo futuro. Es una decisión de ingeniería del proyecto,
no una tolerancia que las publicaciones atribuyan al material.

## Decisión sobre el ensayo ya terminado

**El contraste temporal térmico se acepta con una precisión práctica del 2 %
respecto a la escala inicial de cada señal.** Las 96 comparaciones cumplen,
incluidas las 24 que comparan secuencias temporales distintas. Se cierra la
discrepancia del torque de fase inducido por la sonda de amplitud; no queda
una campaña pendiente para esa cola ni se exige medir su tiempo de relajación.

El error máximo es 8,3046528×10⁻⁸ en la norma espacial del torque adimensional,
equivalente al 1,50001697 % de su propia respuesta inicial 5,53637256×10⁻⁶.
El límite nuevo es 1,10727451×10⁻⁷, exactamente el 2 % de esa referencia.
La cola restante es mucho menor que el error y, por tanto, no proporciona
un tiempo de relajación preciso. Esa información describe la resolución de
la curva; ya no bloquea la aceptación del ensayo.

El [certificado original](../practical_time_review_20260924/decision.json)
y sus datos conservan su resultado fallido con el umbral anterior. La nueva
[decisión legible por máquina](acceptance_policy.json) identifica la fuente,
su huella, el umbral revisado y el resultado. No se modifica el cálculo pasado
ni se atribuye a aquel plan un criterio que no tenía.

## Qué clase de paso temporal se utiliza

Se utiliza la implementación real `TDGLSolver.solve_for_psi_squared` de la
memoria, mediante `inherited_kwt_step`, conservando la nueva fuerza y corriente
espectrales. No se suma nuevamente la fuerza GL ni se convierte la población
electrónica en una temperatura para poder llamar al paso.

La [documentación oficial de pyTDGL](https://py-tdgl.readthedocs.io/en/latest/background.html#implicit-euler-method)
denomina el procedimiento Euler implícito. Sus ecuaciones 13–21 evalúan la fuerza
espacial en el nivel temporal anterior, resuelven una ecuación cuadrática para
el nuevo módulo al cuadrado y recuperan el campo complejo. El
[código oficial](https://github.com/loganbvh/py-tdgl/blob/main/tdgl/solver/solver.py)
es la referencia de implementación.

**El grado de la ecuación algebraica es dos; el orden temporal es uno.**
Esto no descalifica el método: permite tratar el término KWT local sin una
iteración no lineal global. Una comprobación sencilla se obtiene tomando
γ=0 y potencial cero: el código reduce exactamente a
ψⁿ⁺¹=ψⁿ+Δt F(ψⁿ)/u, es decir, Euler hacia delante para la fuerza restante.
Dos operaciones para obtener módulo y campo complejo no forman un corrector
de segundo orden. El enlace de fase exponencial tampoco eleva el orden de
los demás términos congelados durante el paso.

La [prueba local guardada](../practical_time_review_20260924/kwt_bridge/receipt.json)
ya mide cocientes de error 2,00003 y 1,99893 al dividir Δt por dos. Esta
observación es consistente con primer orden y verifica la conversión de
unidades. La selección de Δt de una trayectoria debe resolver sus observables;
el discriminante positivo de la cuadrática asegura una raíz admisible,
pero no es por sí solo un estimador de precisión temporal.

## Criterios del siguiente ensayo sin fotón

Se declara antes de ejecutar la geometría, preparación suave, horizonte,
observables primarios y escalas. No se repite el barrido de colas anteriores.

| Comprobación | Criterio práctico |
|---|---|
| Campos y señales principales | Diferencia de dos pasos temporales ≤2 % de la escala inicial declarada del mismo observable. Para una señal inicial nula se declara antes del cálculo una escala característica no nula del mismo observable. No se divide por el valor de una cola que tiende a cero. |
| Balance térmico integrado | Valor absoluto de `F(t)-F(0)+integral(P_loss-P_in)dt` ≤2 % de la energía libre inicialmente disponible de la preparación. El término de entrada incorpora trabajo de borde cuando corresponda. Se conserva el signo de cada potencia. |
| Diagnósticos secundarios | Se presentan en unidades absolutas y respecto a la escala inicial del canal. Una cola residual por debajo del presupuesto absoluto no causa un rechazo general. |
| Raíces y campos | Todos los campos finitos, ramas espectrales continuas y condiciones de borde satisfechas. Un discriminante inadmisible permite reducir el paso del intento; no se oculta ni se transforma en raíz arbitraria. |
| Disipación y conservación | Se informa el residuo instantáneo como diagnóstico. No se exige un nuevo umbral relativo de 10⁻⁹ cuando numerador y denominador son cancelaciones pequeñas. Un cambio de signo físico o crecimiento de energía libre sin entrada que exceda el 2 % del presupuesto, y persista al reducir el paso, obliga a corregir. |

La energía libre disponible es una diferencia con una referencia térmica
compatible con los mismos bordes, no el enorme valor absoluto de una energía
con cero arbitrario. Si esa diferencia no está acreditada o es nula, el plan
debe declarar otro presupuesto energético de la excitación antes de ejecutar.
No se fabrica un denominador a posteriori. Tampoco se cuenta dos veces como
calor una pérdida que ya aparece en el intercambio espectral.

La tolerancia de Newton no necesita cambiar para corregir el rechazo pasado:
todas las raíces ya cumplieron 10⁻⁷. El criterio nuevo modifica la aceptación
de observables; conserva de momento esa tolerancia interna que funcionó.
Las identidades algebraicas de orientación y unidades son controles baratos
de implementación. Una tolerancia práctica no admite inversiones de signo,
valores no finitos o una población fuera de su soporte físico de manera
significativa; tampoco justifica recortes ocultos.

## Qué concluye cada ejecución

El próximo transiente KWT sobre la malla dual puede cerrar la transferencia
espacial y temporal del **núcleo térmico**. Su balance usa energía libre a
temperatura impuesta. No se etiqueta como balance de energía interna total
del detector ni como un transiente no térmico completo.

El estado actual del código distingue tres piezas:

- `FrozenKineticOperator` calcula transporte y fuerza para distribuciones
  sobre un espectro espacial congelado; no es un integrador de poblaciones.
- `LongitudinalReciprocalBridge` admite una respuesta lineal de amplitud y
  población alrededor de gap real y fase fija. Su disponibilidad cuadrática
  es energía libre de la perturbación; no incorpora colisiones, modo de carga,
  electrostática y circuito en una trayectoria no lineal.
- `electrical_ports` ya contiene potencial, potencias de puerto y las tres
  variables del circuito de la memoria, pero requiere corrientes y conductancias
  físicas suministradas por el modelo espacial.

Por tanto, para cerrar el contrato completo de etapa 4 hace falta un ensayo
suave **acoplado**: poblaciones, trabajo del espectro al cambiar el gap,
potencial y circuito, con balance integrado identificado. No basta sustituir
una etiqueta de etapa o relajar un umbral de error. El ensamblaje debe
reutilizar los bloques existentes; esta observación no pide otro barrido de
integradores, casos singulares ni una recuperación de nanosegundos.

La etapa 5 añade el fotón y la comparación del detector con el experimento.
La transferencia inicial específica y las tasas materiales no admitidas no se
rellenan con valores inventados para cerrar etapa 4. Si se usa un control de
tasas explícitamente prescritas para verificar el ensamblaje, se identifica
como control del modelo, sin presentarlo como predicción calibrada de NbN.
