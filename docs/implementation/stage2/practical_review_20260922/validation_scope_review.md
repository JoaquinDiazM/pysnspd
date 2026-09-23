# Revisión del alcance de la validación

**22 de septiembre de 2026. Revisión y preparación; no implementación espacial nueva.**

El fallo recuperado no corresponde a una trayectoria interrumpida por poblaciones
inválidas o por pérdida del balance de energía. Corresponde al presupuesto
suplementario de una comparación entre dos pasos temporales. Se conserva su
dictamen original `FAIL_PAIR_DIAGNOSTIC`; esta revisión no cambia tolerancias,
no convierte resultados fallidos en aprobados y no declara cerrada la etapa 2.

La finalidad del trabajo es mejorar la capacidad de modelamiento de SNSPD.
Conviene distinguir tres preguntas: si el programa implementa las ecuaciones
elegidas, con qué error las aproxima y qué precisión física ofrecen esas
ecuaciones para el dispositivo. Una identidad algebraica puede necesitar una
comprobación cercana al redondeo, aunque todavía exista una incertidumbre física
grande en el núcleo, la movilidad o los datos materiales. Cumplir la primera
pregunta no resuelve las otras dos; exigir todas las comprobaciones antes de
cualquier trabajo independiente tampoco las resuelve más deprisa.

## Qué se observó

El [resultado original del par fonónico de 2049 nodos](raw/two_guarded_ph2049_pair_320_640.json)
compara 320 y 640 pasos en 321 tiempos comunes. Su diferencia máxima es
`4.565195781e-5`, frente al presupuesto suplementario `2.5e-5`. El máximo aparece
en la población fonónica, en `t = 0.76875`, con diferencia absoluta
`2.651832074e-6` y norma de referencia `0.058088025`.

La norma utilizada suma las diferencias absolutas de las ocupaciones después
de multiplicarlas por las capacidades fonónicas; mide una diferencia en la
población representada, no una comparación entre índices sin peso físico. La
referencia no es casi nula. Los controles de balance integrado, balance
instantáneo y poblaciones aprobaron en ambas trayectorias.

El resultado equivale a una diferencia relativa de aproximadamente **0.004565 %**.
Supera en un factor de aproximadamente 1.826 el presupuesto auxiliar; equivale
al 4.565 % del presupuesto de comparación de mallas `1e-3`. Esta proporción pone
el fallo en contexto, pero no demuestra el error respecto de la solución
continua: dos pasos temporales pueden compartir un sesgo. Tampoco acredita la
convergencia de las mallas que aún no se compararon.

El presupuesto `2.5e-5` se introdujo como una cuarta parte del criterio temporal
`1e-4`, para ayudar a separar error temporal y error de malla. Es una elección
de validación numérica registrada, no una ley física ni una precisión
experimental exigida por los datos de un SNSPD.

Los [criterios congelados de etapa 2](../acceptance_criteria.json) identifican
`agreed_by: ["root", "release_audit"]`. Esos nombres corresponden a agentes que
acordaron el protocolo; **no documentan una aprobación explícita del usuario de
cada cifra**. El usuario sí pidió verificaciones independientes y anticipar el
comportamiento de los transitorios. Esta distinción permite revisar el alcance
del protocolo abiertamente, conservando los veredictos obtenidos bajo sus
reglas originales.

## Qué debe controlar cada comprobación

| Comprobación | Pregunta que responde | Papel en un primer ensayo espacial estático |
|:--|:--|:--|
| Unidades, procedencia, multiplicidades y geometría de eventos; misma transferencia en los balances | ¿Se implementa el modelo declarado, sin factores o energía duplicados? | Conservar los controles pertinentes al funcional; las identidades cinéticas siguen siendo regresiones de etapa 2. |
| Equilibrio, positividad y límites de ocupación | ¿Los operadores respetan sus estados físicos y sus cancelaciones previstas? | Las poblaciones de entrada deben ser físicas. La dinámica de los operadores que no se ejecutan no es una dependencia de un ensayo estático. |
| Energía, fuerzas y corriente obtenidas de un potencial común; soporte del catálogo | ¿Las cantidades evaluadas son consistentes dentro del dominio utilizado? | Necesarios para los estados y campos del ensayo. |
| Derivadas espaciales, interfaces y condición D.36 | ¿La discretización espacial representa la energía elegida y tiene el dominio de estabilidad declarado? | Son controles propios del ensayo espacial; deben comprobarse directamente. |
| Refinamiento temporal, balance durante la evolución y actividad suficiente | ¿La dinámica acoplada está resuelta y realmente ejercita sus términos? | No son dependencias de un cálculo sin evolución temporal ni colisiones. |
| Refinamiento de poblaciones electrónicas/fonónicas y pares temporales auxiliares | ¿La dinámica cambia al aumentar su resolución? | Necesarios para admitir posteriormente esa dinámica; no bloquean por sí solos una prueba estática con ocupaciones congeladas. |
| Calibración NbN, núcleo, disipación, umbrales y latencias | ¿El modelo reproduce el dispositivo en el régimen de interés? | Posteriores; un resultado estático favorable no las acredita. |

Las tolerancias pequeñas de las identidades algebraicas no significan que el
dispositivo se conozca con doce cifras. Sirven, por ejemplo, para detectar un
factor de dos o una transferencia aplicada sólo a uno de los dos sectores.
Son verificaciones normalmente baratas y no conviene abandonarlas. Los
presupuestos de error de una trayectoria cumplen una función distinta y deben
relacionarse con los observables que esa trayectoria pretende resolver.

## Posible contrato separado para una etapa 3A estática

El [contrato de entrada preparado para etapa 3](../../stage3/entry_contract.json)
describe como primer trabajo un funcional espacial con ocupaciones congeladas,
sin integración cinética, fotón ni circuito temporal. Ese trabajo depende de
la energía y sus derivadas, del soporte del catálogo y de D.36. No depende de
que hayan terminado las comparaciones dinámicas de todas las mallas fonónicas.
El bloqueo actual de toda implementación de etapa 3 es, por tanto, más amplio
que las dependencias científicas de ese primer trabajo.

Una apertura **exploratoria** requeriría un contrato separado, explícito y
registrado antes de obtener sus resultados:

1. Limitar el trabajo a una tira pequeña y sintética: energía discreta D.8,
   fuerzas cartesianas, corriente, estados uniformes y gradientes débiles, con
   ocupaciones térmicas y no térmicas congeladas.
2. Identificar la representación espectral y verificar de forma independiente
   energía, fuerzas, corriente y soporte en los estados efectivamente usados.
   Comprobar D.36 antes de admitir un estado para futura evolución. Los estados
   inestables pueden conservarse como controles negativos identificados.
3. Mantener `stage2.numerical_admission = false` mientras no se cumpla su
   certificado estricto. Conservar el fallo del par y todos los demás
   resultados históricos, con sus métodos y alcances originales.
4. Registrar los resultados nuevos como implementación exploratoria y, cuando
   corresponda, validación estática del dominio probado. Excluir pulsos,
   latencias, umbrales de detección, integración acoplada del circuito y
   promoción a producción.

**Este documento sólo prepara esa posibilidad. No modifica el dictamen global
de etapa 2 ni el contrato de entrada de etapa 3, y no inicia su implementación.**

## Política proporcional para la continuación

Primero debe declararse qué diferencia interesa resolver: recuperación del
condensado, transferencia de energía, forma del pulso, umbral o latencia. A
partir de esa pregunta se establece un presupuesto numérico prospectivo, con
sus escalas, referencias y alcance. Cerca de un umbral, la estabilidad de la
clasificación puede ser más relevante que un porcentaje sobre un observable
que se anula. Este documento no elige una precisión final nueva.

Un fallo estructural —unidades incompatibles, trabajo contado dos veces,
cantidades no finitas, poblaciones inválidas o salida del dominio admitido—
debe detener el cálculo afectado. Un presupuesto de precisión incumplido debe
mantener su dictamen y su incertidumbre visibles; no implica automáticamente
que todo trabajo independiente de ese observable deba detenerse. Cualquier
plan futuro que distinga esos casos debe declararlo antes de ejecutarse.

Cambiar un informe o un orquestador tampoco obliga por sí solo a repetir la
física. Los datos originales pueden reutilizarse y volver a evaluarse cuando
sus fuentes físicas, parámetros, método, condiciones iniciales y archivos
sean compatibles y verificables. Un cambio real del mapa de integración o de
la representación numérica exige sus comprobaciones pertinentes; no permite
transferir aprobaciones de forma automática.

La continuación propuesta conserva la validación y coloca cada control en el
trabajo que depende de él: avanzar en el funcional espacial estático mientras
se mantiene abierto el cierre numérico de las celdas, y exigir la validación
dinámica y del dispositivo antes de atribuir fiabilidad a transitorios o
promover el modelo a producción.
