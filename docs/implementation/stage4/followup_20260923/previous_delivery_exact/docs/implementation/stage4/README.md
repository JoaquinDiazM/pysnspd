# Etapa 4: resultados 4A y continuación enfocada

Los 40 controles locales y seis estados 2D terminaron. El
[informe de resultados](review_20260923/Informe_avance_etapa_4A.md)
([PDF](../../../output/pdf/implementation/Informe_avance_etapa_4A.pdf)) separa
estabilidad constitutiva, errores de borde y resolución del núcleo.

El perfil suave cambia sólo 1,50 % en disipación interior entre mallas. La
duplicación del calor total original proviene de mover libremente sus bordes
prescritos. Se reconstruyó la respuesta con cargas conjugadas que los mantienen
fijos, sin nuevas consultas espectrales. Este control no representa los bordes
ni los contactos del dispositivo.

El núcleo suprimido todavía cambia 28,9 % en disipación interior. Su único signo
espacial negativo no se reproduce con la derivada analítica del perfil; requiere
resolver el gradiente antes de atribuirlo a la física. Se conservan por separado
los cuatro controles locales de gradiente fuerte que sí tienen D.36 negativo.

La continuación tiene [sólo dos estados nuevos](review_20260923/next_campaign_plan.json):
δ = 0,05 en 8×4 elementos y δ = 0,1 en 16×8, ambos de grado 4. Se compararán
derivadas, fuerza volumétrica, calor interior y corriente sumada por sección.
El comando manual, con barras y ETA, está en
[GEMINGA_COMMANDS.md](../../GEMINGA_COMMANDS.md). Estimación: 40-65 min.

La [referencia Usadel linealizada](review_20260923/linear_usadel_reference.md)
actualizada a Korzh cuantifica un sesgo constitutivo de rigidez de 3,03 %
a 98,4 nm y 30,88 % a 29,5 nm. Es una referencia térmica independiente, no un
transiente a ocupación fija ni una validación de núcleo físico.

**La etapa 4 sigue abierta.** Faltan resolución del núcleo, su referencia física
y evolución débil con las dependencias dinámicas correspondientes de etapa 3.
No hay fotón, hotbelt, latencia ni promoción a producción. El circuito previsto
sigue siendo el de la memoria.

[Procedencia numérica](review_20260923/solver_scope.md): producción conserva
Delaunay-Voronoi y Euler adaptativo de primer orden; el banco 4A usa GLL 2D
experimental y no integra el tiempo. No se ha escogido un reemplazo definitivo.

La [preparación original](start_20260923/README.md) y sus estados de pendiente
se conservan como historial. La revisión preserva resultados, fuentes y entrega
anteriores; se verifica con `sandbox/stage4_core/verify_review.py`.
