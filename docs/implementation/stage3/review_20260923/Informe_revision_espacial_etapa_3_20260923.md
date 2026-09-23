# Etapa 3: revisión del lote y siguiente ejecución

Resultados del ensayo espacial | 23 de septiembre de 2026

## El lote terminó; falta precisión espacial

Los <b>18 casos terminaron en 15.1 minutos</b>. Pasaron las derivadas de energía, los controles de fase, los signos locales y los contrastes espectrales muestreados. El dictamen final fue precisión espacial insuficiente: no hubo un fallo de integración temporal ni una ejecución interrumpida.

| Respuesta en la malla de 32 celdas | Estimación de error |
| --- | --- |
| Amplitud térmica: energía / fuerza / corriente | 1,86 % / 2,00 % / 69,78 % |
| Fase térmica: energía / fuerza / corriente | 5,89 % / 8,66 % / 6,09 % |
| Fase no térmica: energía / fuerza / corriente | 5,87 % / 8,76 % / 6,07 % |
| Control de gradiente aislado | 0,321 %: aprobado |

![Izquierda: estimaciones condicionales al refinar 8/16/32 celdas; la corriente nula del vacío no lleva porcentaje. Derecha: corriente uniforme frente al límite continuo, con población congelada. El punto verde es el nuevo piloto; aún no es una curva de convergencia nodal.](figures/spatial_review.png)

La causa principal es geométrica: promediar dos números complejos con distinta fase acorta su amplitud, aunque ambos nodos tengan el mismo módulo. En el estado prescrito esa reducción introduce corriente numérica adicional. Una expansión analítica explica el 98,23 % del sesgo uniforme de vacío a 32 celdas.

Los 11 fallos espaciales se conservan. Los 10 contrastes electrónicos 630/1260 mostraron cambios de matriz de hasta 4,20e-9; no aportan evidencia de que aumentar globalmente el espectro sea el siguiente paso útil.

## La siguiente corrida cambia la discretización

Se implementó una energía <b>nodal de cuarto orden</b>: consulta la amplitud en cada nodo y usa operadores covariantes compatibles para gradiente y corriente. Fuerza y corriente siguen siendo derivadas de la misma energía. El modelo continuo, el catálogo, la regularización y el objetivo del 1 % permanecen iguales.

| Verificación nueva | Resultado |
| --- | --- |
| Suite del esquema nodal y controles históricos | 102 pruebas aprobadas en 7.26 s |
| Piloto nodal de ocho celdas | Aprobado en 33.2 s |
| Corriente uniforme del piloto frente al continuo | Diferencia 0.535 % |
| Convergencia espacial del esquema nodal | Pendiente del lote del usuario |

La energía de gradiente también controla oscilaciones alternantes entre nodos. Usar solamente una derivada central de orden alto dejaría ese modo sin penalización; el operador positivo complementario evita ese defecto.

El nuevo lote mantiene los seis perfiles y las mallas de 8/16/32 celdas. Las diferencias independientes de energía se ejecutan para cada caso en 16 celdas, además del piloto de ocho; no se repiten las 24 consultas derivativas en cada malla. Los controles de fase, referencias uniformes, signos y respuestas espaciales se mantienen en las tres mallas.

<b>Próximo paso:</b> ejecutar el bloque activo de /home/jdiaz/GEMINGA_COMMANDS.md. La terminal muestra barras, tiempo transcurrido y ETA; progress.jsonl guarda el avance. Se reutiliza únicamente el piloto nodal con fuentes y estados verificados.

Este lote decide el cierre del bloque estático 3A. Después quedan bordes y reservorios, potencial eléctrico, las tres ecuaciones del circuito de la memoria y la dinámica débil con depósito sintético. La etapa 3 completa sigue abierta; todavía no se ha calculado una señal de detección.
