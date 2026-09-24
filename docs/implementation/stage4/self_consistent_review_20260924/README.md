# Núcleo autoconsistente, espectro y respuesta de carga

Los cuatro casos de la continuación cumplen el criterio anunciado. Se acepta
la referencia térmica autoconsistente. También se calcularon 72 espectros a
energías reales, 144 respuestas cinéticas y seis comparaciones con un potencial
electroquímico por nodo. El siguiente paso comprueba sus integrales energéticas.
No se solicita otra relajación térmica ni se endurece retrospectivamente su
tolerancia. La etapa 4 sigue abierta; todavía no corresponde un transiente
con fotón ni declarar implementado el modelo completo.

El [informe ilustrado](../../../../output/pdf/implementation/Informe_etapa_4_nucleo_autoconsistente.pdf)
resume la decisión y los resultados. El [análisis reproducible](analysis.md)
identifica medidas, normas e integridad; la [revisión física](physics_review.md)
deduce la continuación del mismo problema espectral y sus límites.

## Resultado aceptado

La corrida del usuario duró 484,52 s y resolvió 12.672 problemas espectrales.
Al incluir el piloto, los casos necesitaron 15, 16, 16 y 17 barridos; esos índices
son pasos de minimización, no tiempo del detector. Los residuos RMS del núcleo
quedaron entre 0,084 % y 0,098 %, dentro del 0,1 % declarado. El máximo nodal
global se conserva aparte y alcanza 0,195 %.

Las energías disminuyen en los 64 registros completos. Los 23 checkpoints
completos se comprobaron por SHA-256; los cuatro estados finales se verificaron
con su propio espectro. El campo `next_d` se mantiene separado: es una propuesta
de actualización y no una solución espectral aceptada.

| Comparación de estados relajados | Gap, norma L2 | Densidad de corriente, norma L2 |
|---|---:|---:|
| Malla 65² a 129², 256 frecuencias | 0,0384 % | 0,0798 % |
| Corte 128 a 256, malla 65² | 0,0258 % | 0,6151 % |

Son sensibilidades observadas, no incertidumbres experimentales. Se comparan
coordenadas comunes, áreas duales y densidades de corriente. La corriente
integrada de un enlace cambia con el ancho de su cara; compararla directamente
entre mallas daría una diferencia geométrica espuria.

El vórtice asimétrico conserva una vuelta de fase y se desplaza unos 0,46 nm.
El cero queda entre nodos: un mínimo nodal distinto de cero no significa que
el vórtice desapareció. La ubicación se calcula con interpolación bilineal
del campo complejo, con el alcance espacial de esa interpolación.

## Espectro y respuesta de carga implementados

El módulo experimental `retarded_spatial_usadel.py` resuelve el mismo problema
espectral con dos campos complejos independientes, manteniendo el exterior
radial de los bordes. Pasaron las 72 consultas en 88,48 s. La normalización,
la rama causal y la correspondencia con la acción térmica se comprobaron.
La diferencia espacial máxima de DOS es 0,454 %; cambiar el desplazamiento
numérico η de 0,04 a 0,02 modifica la DOS hasta 9,92 %. η no es una tasa física.
Véase el [análisis espectral](retarded_oracle/analysis.json).

`frozen_kinetic_usadel.py` construye el transporte de energía y carga desde
esos mismos campos. Las 144 direcciones lineales, calculadas en 5,13 s,
comprueban conjuntamente fuerza, corriente y conversión de pares. La identidad
integrada tiene residuo máximo 1,14×10⁻¹⁶. La perturbación angular cambia el
flujo espectral de energía hasta 18,10 % frente a congelar la respuesta de
carga; la radial, hasta 1,52 %. Son valores por energía y para gap fijo,
distintos de los momentos integrados que siguen. No constituyen ocupaciones
fotónicas finitas ni prueban que falle el cierre completo de fase y potencial
de la memoria. [Evidencia cinética](kinetic_response/analysis.json).

La [proyección electroquímica](potential_projection/summary.json) restringe
esa respuesta a una dirección térmica conocida multiplicada por un potencial
por nodo. Tardó 2,20 s y conserva el balance de carga integrado. En las sondas
angulares difiere de la respuesta completa en 19,92–20,57 % de corriente,
3,11–3,17 % de fuerza y 0,415–0,430 % de flujo de energía ponderado por energía.
Son diferencias de campos en las normas declaradas, no sólo diferencias de
sus magnitudes globales. Los contactos mantienen la misma perturbación nula.

La cuadratura actual de doce energías sobrestima un peso térmico de integral
conocida en 18,02 %. No se corrige renormalizándolo ni se atribuye toda la
diferencia de corriente al modelo. La [campaña siguiente](resolution_campaign/README.md)
usa 31 y 50 energías anidadas, concentradas en las escalas que faltan. Compara
momentos de fuerza, corriente y energía, junto con la sensibilidad a η y a la
malla espacial, para distinguir una reducción física útil de un error de
cuadratura. No repite la relajación del núcleo ni impone una tolerancia diminuta
a cada punto de DOS.

La [revisión analítica](bridge_review.md) define los operadores, las unidades
y los límites de la comparación. Antes de cerrar la etapa 4 falta el trabajo
espectral al mover el gap y su unión dinámica con poblaciones, potencial,
fase y disipación. Una identidad de carga estática no sustituye ese balance
de energía. Tampoco se admite automáticamente una variable dinámica nueva.

La malla, el método temporal y el circuito de la memoria permanecen como
infraestructura de producción. El trabajo nuevo está en módulos experimentales.
La futura ventana hasta Vout más margen continúa siendo sólo una elección del
horizonte de observación del mismo sistema.

## Datos conservados

`raw/` contiene identidad, plan ejecutado, resumen, historial y cuatro mapas
compactos. `raw/extraction_receipt.json` acredita la comprobación remota de campos
completos; el postproceso tardó 1,31 s y no volvió a resolver espectros.
Las matrices espectrales completas permanecen en
`/home/jdiaz/scratch/pysnspd_stage4_self_consistent_20260924` y los checkpoints
del piloto anterior en su ruta original. El comando y los recursos de la
siguiente prueba se publican en [la libreta](../../../GEMINGA_COMMANDS.md)
y se copian explícitamente en el chat.

La carpeta `previous_delivery_exact/` conserva los documentos de estado de la
entrega anterior. Las fuentes ejecutadas de la campaña térmica no se modifican.
