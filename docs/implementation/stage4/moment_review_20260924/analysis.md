# Revisión de momentos: resultado y decisión

La campaña completó las 181 consultas espectrales, sus 181 respuestas cinéticas y siete proyecciones en **198,12 s**. No falló un solver ni hace falta repetir la relajación térmica. El resultado distingue una limitación de la proyección estática de potencial de un error de resolución.

Se comprobaron 369 archivos completos, las fuentes y el plan ejecutados. Los siete mapas integrados se descargaron y los momentos se recalcularon localmente. No se ejecutaron nuevas soluciones espectrales.

## Qué significan D, C y U

`D` mide la diferencia entre la proyección de potencial y la respuesta completa; `C`, la corrección que se pierde al omitir el modo de carga. `U_D` suma los cambios observados al refinar energía, reducir eta y refinar la malla en ambas respuestas comparadas; `U_C` hace lo mismo para la omisión. Se conserva el criterio previo: una diferencia se resuelve si `D>U_D`; la proyección mejora la omisión si `D+U_D<C−U_C`. U no es una cota rigurosa de error.

## La sonda angular discrimina el cierre

Referencia con los tres controles: núcleo radial 65², eta/Delta=0,01, malla fina de 50 energías. Las cifras siguientes se normalizan por la norma completa de cada observable.

| Observable | D / norma completa | U_D / norma completa | D / C | Interpretación |
|---|---:|---:|---:|---|
| Corriente de carga | 20.54644 % | 2.87139 % | 0.40875 | Diferencia resuelta; mejora la omisión |
| Flujo ponderado por energía | 0.43817 % | 0.19597 % | 0.98229 | Diferencia resuelta; mejora incierta |
| Fuerza compleja | 3.22860 % | 1.53305 % | 0.88083 | Diferencia resuelta; mejora incierta |
| Fuerza de amplitud | 0.00153 % | 1.52018 % | 3.19027 | Diferencia no resuelta |
| Torque de fase | 77.52784 % | 3.63645 % | 0.77528 | Diferencia resuelta; mejora la omisión |

El 3,23 % de diferencia de fuerza compleja oculta dos situaciones muy distintas: la fuerza de amplitud apenas cambia, mientras el torque de fase difiere un 77,53 % respecto de su propia norma. El torque proyectado tiene una norma 1,717 veces la completa. Esto no implica un error de latencia del 77,53 % ni invalida automáticamente una dinámica con fase y potencial acoplados.

En la sonda radial todas las diferencias quedan sin resolver frente a U. Sus pequeñas componentes de carga/fase no sirven por sí solas para aceptar el cierre, ni sus cocientes pueden trasladarse a la sonda angular. El caso asimétrico conserva datos útiles, pero sólo tiene 31 energías: no dispone del mismo presupuesto independiente de controles.

## Qué cambió al resolver mejor la energía

La integral conocida de la susceptibilidad pasa de un exceso de 18,02 % con 12 puntos a 1,616 % con 31 y 0,4243 % con 50. No se renormalizó esa función. La discrepancia angular de corriente permanece alrededor del 20,55 % después de mejorar la cuadratura. No procede atribuirla enteramente a los 12 puntos iniciales.

| Control del resultado completo, sonda angular | Corriente | Fuerza compleja | Torque de fase | Flujo energético |
|---|---:|---:|---:|---:|
| 31 a 50 energías | 0.30514 % | 0.04434 % | 0.42264 % | 0.00758 % |
| eta 0,02 a 0,01 | 0.98608 % | 0.65909 % | 0.44289 % | 0.03043 % |
| 65² a 129² | 0.14804 % | 0.05976 % | 0.24939 % | 0.06198 % |

Se comparan densidades de fuerza sobre coordenadas comunes y corrientes por ancho de cara dual; nunca corrientes brutas de enlaces de diferente tamaño. El flujo energético lleva el peso E en la integral. Todas las cuadraturas de fuerza, corriente y residuo son las mismas.

## Decisión para continuar

Este ensayo de resolución queda terminado. La proyección de un solo perfil térmico de carga por potencial no puede presentarse como sustitución exacta de la respuesta cinética congelada: deja diferencias resueltas en corriente y fase. Sí mejora omitir por completo la respuesta de carga.

La siguiente decisión requiere estudiar la unión dinámica entre fase del condensado, potencial y poblaciones, incluyendo su balance energético. No se añade todavía un nuevo estado dinámico de carga, no se rechaza automáticamente el esquema completo de la memoria y no se altera producción. La etapa 4 sigue abierta; el transiente de fotón y la latencia del dispositivo siguen pendientes.

Los archivos raw/ conservan los resúmenes, recursos, siete mapas integrados y el recibo de integridad remoto. analysis.json contiene magnitudes absolutas, cocientes, componentes de U y verificaciones.
