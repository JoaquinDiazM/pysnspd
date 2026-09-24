# Horizonte de los futuros transientes fotónicos

Estado: política de observación acordada, pendiente de implementación en el
integrador experimental. No cambia producción en esta entrega. Se aplica a
las futuras comparaciones del nanohilo de 80 nm de Korzh; los controles
estáticos sin fotón de etapa 4 no necesitan un disparo eléctrico.

## Decisión de alcance

La prioridad es estudiar el retardo relativo entre 775 y 1550 nm. Se seguirá
la evolución hasta un disparo confirmado de **Vout**, más un intervalo
posterior registrado. No se exigirá simular la recuperación completa durante
nanosegundos para obtener ese observable. La recuperación larga queda
pospuesta hasta establecer la comparación de retardos pertinente.

Se cambia **únicamente el horizonte de integración y observación**. Hasta el
instante de parada permanecen las mismas ecuaciones físicas, el circuito de
la memoria con sus tres estados, la geometría, condiciones de borde,
parámetros y acoplamientos. No se vuelve cuasiestático el circuito, no se
eliminan inductancias o capacitancias, no se congela su evolución y no se
reinicia artificialmente la cola del pulso. Se conserva el estado completo
final para poder continuar una corrida posteriormente si hace falta.

## Definición del evento y de la parada

`Vout` es el voltaje de salida de la carga del circuito, distinto del voltaje
mesoscópico dentro de la zona simulada. Se fija una referencia previa al
evento y una polaridad común. El disparo se refiere al primer cruce ascendente
de un umbral registrado por la señal de salida respecto de esa referencia.

Sean `t_i,t_(i+1)` dos tiempos aceptados que encierran el cruce y `v_i,v_(i+1)`
las excursiones de voltaje con esa polaridad. El tiempo interpolado es

\[
t_\times=t_i+
\frac{V_{\rm umbral}-v_i}{v_{i+1}-v_i}(t_{i+1}-t_i).
\]

El cruce deberá confirmarse con la duración e histéresis que se registren
para el criterio operacional. Su tiempo de detección sigue siendo
`t_cross`, no el instante posterior en que se terminó de confirmarlo. La
confirmación evita contar una oscilación numérica como un disparo; no añade
un retardo físico a la comparación de colores. Solo se usan estados aceptados
y muestras suficientemente densas para resolver el cruce y su confirmación.
Si ya se empieza por encima del umbral o falta una referencia previa válida,
no se inventa retrospectivamente un cruce.

Una vez confirmado, el destino de parada será

\[
t_{\rm fin}=\max(t_\times+\Delta t_{\rm extra},t_{\rm confirmación}).
\]

Se registran el instante interpolado, el instante de confirmación y el instante
real de parada por separado. El último paso debe ajustarse al destino o
declarar cualquier exceso por resolución del integrador. No se necesita
esperar al máximo del pulso ni a su recuperación completa para esta regla.
Las medidas que sí dependan del máximo o de toda la cola no se declararán
completas con una trayectoria truncada.

También habrá un techo temporal finito registrado para corridas sin disparo.
Si termina una trayectoria numéricamente válida sin cruce confirmado, el
resultado queda **censurado por la ventana temporal**: no se asigna retardo
cero, infinito ni igual al techo como si hubiera ocurrido una detección. Si
hay un fallo numérico, se informa el fallo por separado; no se interpreta
como ausencia física de detección. Si el cruce ya es válido pero el techo
impide completar el margen posterior, pueden conservarse el cruce y su
latencia, indicando que el margen quedó incompleto.

## Comparación entre colores

775 y 1550 nm compartirán señal, referencia operacional, umbral, polaridad,
confirmación, histéresis, margen posterior y criterio del techo temporal.
No se ajustará un umbral distinto para hacer coincidir cada color. Se
registrará, con signo explícito,

\[
\Delta t_{\rm modelo}=
(t_{\times,1550}-t_{0,1550})-(t_{\times,775}-t_{0,775}).
\]

El reloj `t0` sigue siendo el de transferencia al modelo. La posible
diferencia entre retardos ópticos omitidos sigue explícita; acortar el
horizonte no la determina. Si uno de los dos cruces está censurado, la
diferencia no se presenta como un valor medido completo. La incertidumbre
temporal numérica se estima sobre el cruce relevante, no imponiendo una
tolerancia universal ajena al observable.

## Situación del código y decisiones pendientes

Producción ya contiene un cruce interpolado con confirmación en
`pysnspd/analysis/timing.py:_confirmed_crossing_time`. Sin embargo, su
`early_stop_mode="latency"` **no implementa esta nueva regla**: espera un
máximo confirmado y un margen desde ese máximo. Véanse
`timing.py:_detection_result` y `solver/transient.py:CoupledTransientConfig`.
Por tanto no basta con cambiar el nombre del modo de parada. La nueva
semántica deberá implementarse y comprobarse expresamente en la futura ruta
experimental, sin alterar retrospectivamente las corridas archivadas.

Los valores por defecto existentes, entre ellos el umbral de 100 microvoltios,
no se adoptan automáticamente como valores físicos para Korzh. En este
contrato quedan abiertos el umbral, la polaridad, el margen posterior, la
confirmación, la histéresis y el techo temporal de ausencia de disparo. Antes
de preparar la primera campaña fotónica se intentará fijarlos con las fuentes
y la cadena de lectura que se represente. Si esa información no basta, se
consultará al usuario explicando:

- El umbral: excursión de voltaje de salida que se contará como detección.
- El margen: tiempo adicional observado después del cruce, para comprobar su
  persistencia y capturar la respuesta inmediata de interés.
- El techo: tiempo máximo observado si el cruce todavía no ocurre.
- La confirmación e histéresis: persistencia mínima y tolerancia a un descenso
  breve alrededor del umbral.

No hace falta resolver esas elecciones para avanzar ahora en los controles
de núcleo y disipación sin fotón. Esta política tampoco elimina los
pendientes de preparación fotónica, tasas materiales o admisión dinámica.

## Integración documental propuesta

La portada vigente debe enlazar esta política como adenda de horizonte
temporal. `MODELO_VIGENTE.md` debe indicar que el sistema físico y el circuito
siguen iguales; `SECUENCIA_VIGENTE.md` debe distinguir la comparación temprana
de retardos y la recuperación larga posterior. Esta entrega conserva esos
archivos históricos y sus hashes; no los reescribe desde esta nota.
