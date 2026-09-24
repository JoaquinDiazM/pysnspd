# Modelo previsto y ruta de implementación vigente

Actualizado el 24 de septiembre de 2026. El tag `v1.0.0` y el solver de producción
conservan la implementación de la memoria. Lo siguiente identifica el modelo
previsto y sus pruebas experimentales posteriores.

La base es el [sistema continuo 0.4](../modelo_v0_4/D_sintesis_plan_y_verificaciones_v0_4.md),
con esta sustitución normativa solicitada por el usuario:

- **Circuito:** [red de polarización y lectura de la memoria](../modelo_v0_4/actualizaciones/circuito_memoria_20260922.md).
  Sus tres variables son las corrientes de polarización y del detector y el
  voltaje del capacitor. Sustituye el circuito ideal simplificado de C/D,
  sus balances y su inicialización. La adenda contiene las ecuaciones completas,
  las convenciones de signo y la partición de inductancia que evita doble conteo.
- **Implementación:** el [cierre de desarrollo de etapa 2](stage2/closure_20260922/closure_decision.json),
  autorizado por el usuario, acepta la evidencia disponible para continuar con
  la secuencia espacial, sus bordes y el circuito. No cambia retrospectivamente
  las tolerancias ni convierte fallos en pases.

La etapa 2 queda cerrada como desarrollo. El
[informe final](stage2/closure_20260922/Informe_cierre_etapa_2.md)
([PDF](../../output/pdf/implementation/Informe_cierre_etapa_2.pdf)) reúne las
13 trayectorias del último lote, los pases temporales de una/dos celdas y las
limitaciones medidas. Su [certificado numérico estricto](stage2/stage2_admission.json)
sigue incompleto: se conserva el fallo auxiliar y no se atribuye convergencia
a las mallas dinámicas pendientes.

La etapa 3 se cerró como desarrollo con pendientes dinámicos explícitos, y la
investigación 3.5 se cerró dentro de su alcance autorizado. La
[etapa 4 sin fotón](stage4/README.md) está en ejecución: los controles de campo
prescrito permiten separar sensibilidad de malla y de cierre constitutivo.
El cierre local regularizado falló el contraste radial y se conserva como
antecedente rechazado para esa dinámica. La nueva representación espacial
recupera dicho contraste y produce sus propios núcleos térmicos. La
[revisión autoconsistente](stage4/self_consistent_review_20260924/README.md)
acepta los cuatro casos y obtiene su espectro a energías reales y su respuesta
cinética lineal de carga y energía.
Siguen pendientes las capacidades dinámicas antes de interpretar el detector.
La transferencia fotónica, su ancho y la normalización de tasas NbN no se
admiten mediante este avance estático.

## Ventana inicial de observación del mismo dispositivo

La prioridad acordada es el retardo relativo de 775 frente a 1550 nm en la cinta
de 80 nm de Korzh (2020). Los futuros transientes se seguirán hasta un cruce
confirmado del gatillo en Vout más un margen declarado; no se exigirá primero
la recuperación completa durante nanosegundos. La
[política de horizonte](stage4/followup_20260923/horizon_policy.md) define el
observable, el tratamiento de corridas sin cruce y las decisiones aún abiertas.

**Sólo cambia el horizonte de integración/observación.** Ecuaciones, material,
depósito, geometría, malla, bordes, circuito completo y política de precisión
permanecen iguales entre una corrida corta y su extensión. No se eliminan
variables lentas ni se acortan constantes de tiempo. Se conserva el estado
final para continuar. Umbral y margen se especificarán para el ensayo; no se
adoptan automáticamente los valores predeterminados de producción.

El modo heredado `latency` espera un máximo confirmado y todavía no implementa
la regla gatillo más margen. Esta adenda registra la política futura; no cambia
el solver de producción. El reloj sigue comenzando en la transferencia al modelo,
con el retardo óptico previo desconocido por separado.

Antes de una predicción dinámica se fija qué diferencia del dispositivo se
quiere resolver y se registra su presupuesto numérico. Las identidades de
energía, unidades, dominios y signos siguen siendo controles necesarios; una
precisión numérica alta no acredita por sí sola el modelo físico del detector.

Los C/D y PDF históricos se mantienen íntegros. Cuando se solicite su siguiente
edición integrada, deberá incorporar esta adenda; no debe recuperarse el circuito
simplificado como contrato vigente por leer un PDF anterior.
La [revisión previa del alcance](stage2/practical_review_20260922/validation_scope_review.md)
se conserva como antecedente de esta decisión; su preparación limitada a 3A
no sustituye la secuencia completa vigente.

## Referencia térmica espacial regular

La [entrega espacial](stage4/spatial_energy_20260924/README.md) implementa
un oráculo térmico con campos espectrales complejos, sin completado δ ni
gradiente K0 añadido. Recupera el contraste radial y deriva fuerza y corriente
de la misma energía. Los cuatro casos de la campaña autoconsistente cumplen el criterio declarado.
No sustituye el cierre de poblaciones no térmicas: la rama retardada y el
operador cinético congelado ya están implementados, pero falta el balance
dinámico de trabajo espectral, transporte y disipación. Producción sigue intacta.


La [decisión actual](stage4/self_consistent_review_20260924/acceptance_decision.json)
cierra únicamente el problema estático térmico al nivel medido, sin otra
refinación. Su continuación retardada usa coordenadas complejas independientes
y el mismo entorno radial de borde. Es una transformación de la acción
espectral, no una sustitución de distribuciones por temperatura. El
[puente físico pendiente](stage4/self_consistent_review_20260924/physics_review.md)
ya obtiene incrementos de fuerza y corriente con una identidad común de
conversión de carga. La [revisión de momentos y frecuencia](stage4/moment_review_20260924/README.md)
resuelve la discrepancia de la proyección a un potencial por nodo: 20,55 % en
corriente y 77,53 % en torque de fase, frente a variaciones numéricas de 2,87 %
y 3,64 %. No se admite esa compresión energética para este ensayo. La eliminación
algebraica espectral completa conserva el detalle energético y aproxima bien la
respuesta lenta medida; no acredita todavía el régimen ultrarrápido del detector.
Tampoco se rechaza por este control congelado toda la fase y el potencial de la
memoria. El operador temporal térmico ya se contrastó con 2304 raíces espectrales;
la [trayectoria afín completada](stage4/time_review_20260924/README.md) mantiene
la deriva de referencia y pasa el contraste temporal. El posterior contraste
no lineal también pasa tras corregir cancelación de redondeo en Newton. Sigue
la integración térmica no lineal ETD2 con la misma fuerza, KWT y potencial normal. Falta unir trabajo y transporte no térmicos antes de evaluar
un fotón; el circuito completo y el horizonte acordado siguen vigentes.
