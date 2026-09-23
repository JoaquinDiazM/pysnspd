# Modelo previsto y ruta de implementación vigente

Actualizado el 22 de septiembre de 2026. El tag `v1.0.0` y el solver de producción
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

La [etapa 3](stage3/README.md) queda **preparada, no iniciada**, con su
[contrato de entrada](stage3/entry_contract.json). La secuencia comienza por
energía espacial y estabilidad con ocupaciones congeladas; continúa con bordes
y reservorios, conservación de carga, circuito de tres variables y después
dinámica débil y deposición localizada sintética. Cada ensayo registrará sus
controles propios antes de ejecutarse. El cierre de desarrollo no constituye
una promoción a producción ni una validación de transitorios del detector.

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
