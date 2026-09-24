# Etapa 4: respuesta de carga resuelta y preparación dinámica

La [revisión actual](moment_review_20260924/README.md) reúne la campaña
del usuario y el diagnóstico temporal ligero. La etapa 4 continúa: ya no hace
falta otra campaña estática para decidir la compresión del modo de carga,
pero todavía falta completar la dinámica acoplada y su balance de energía.

La campaña de momentos terminó en 198,12 s: 181 espectros, 362 respuestas a dos
sondas y siete proyecciones. Reducir el defecto de cuadratura de 18,02 % a
0,424 % no elimina la diferencia de corriente del potencial único, 20,55 %.
La diferencia del torque de fase es 77,53 %; ambas superan la variación
numérica observada. Son normas de respuestas de prueba, no errores de latencia.

La respuesta armónica reutilizó esos espectros: 24 casos en 21,19 s. A la
frecuencia lenta ν=0,001, conservar toda la dependencia energética y eliminar
algebraicamente el modo de carga cambia la corriente sólo 0,0258 % y el torque
0,0422 %. Esto respalda esa reducción en el diagnóstico lento, sin añadir un
estado temporal por defecto. Las frecuencias altas prueban el cierre aproximado;
no acreditan la respuesta AC subpicosegundo del detector.

El operador para la evolución térmica débil con gap móvil ya pasó en 32,94 s,
con la misma acción espacial, KWT y potencial normal. Las respuestas tangentes
reutilizan 256 factores espectrales; sus derivadas se contrastaron con 2304
raíces independientes. Se conserva el residuo real del núcleo y su corrección
temporal del 1,1 %. Sigue la trayectoria afín de todos los nodos hasta 1 ps,
con referencia y perturbaciones débiles separadas. Se distingue energía libre
térmica del trabajo espectral no térmico pendiente. Véase la
[decisión física](moment_review_20260924/physics_decision.md).

El [cuaderno vigente](../../GEMINGA_COMMANDS.md) identifica la siguiente prueba,
salidas, recursos y comandos ligeros. Los cálculos previstos de más de cinco
minutos son manuales y se entregan también en el chat. El presupuesto compartido
máximo es 28 de 32 hilos, reservando dos núcleos físicos completos en Geminga.

## Evidencia anterior conservada

- [Núcleos autoconsistentes y primeros espectros](self_consistent_review_20260924/README.md):
  cuatro núcleos aceptados, sin otra relajación térmica requerida.
- [Energía espacial](spatial_energy_20260924/README.md): recupera la fuerza radial
  al mismo corte y deriva fuerza y corriente de la misma energía.
- [Diagnóstico del cierre local](followup_20260923/README.md): explica por qué
  no bastaba refinar la malla del cierre anterior.
- [Controles iniciales](review_20260923/Informe_avance_etapa_4A.md): conservan
  valores, conclusiones y límites originales.

Producción y v1.0.0 siguen intactos. Se conserva el circuito completo de la
memoria. Para el futuro fotón, sólo cambia la ventana hasta Vout más margen.
Las incertidumbres de transferencia fotónica y tasas materiales permanecen
explícitas; la etapa 5 no ha comenzado.
