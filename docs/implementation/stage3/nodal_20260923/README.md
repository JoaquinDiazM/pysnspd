# Etapa 3A: corrección de la discretización y siguiente lote

**Estado: piloto nodal aprobado; campaña espacial nodal pendiente del usuario.**

El lote anterior completó 18 casos en 15,15 minutos. Pasaron las derivadas de
energía y los signos muestreados; 11 respuestas espaciales incumplieron el
objetivo del 1 %. El promedio cartesiano de una hélice reduce su amplitud y
produce un sesgo de corriente. La expansión analítica explica el 98,23 % de
ese sesgo en el uniforme de vacío de 32 celdas.

La [auditoría](../review_20260923/audit_notes.md) conserva resultados, fallos,
fuentes y recibos. No hay indicios en los contrastes 630/1260 muestreados de
que ampliar globalmente el espectro sea el siguiente paso útil. No se repite
la etapa 2.

## Qué cambia en la implementación

`spatial_nodal.py` consulta amplitudes y poblaciones en los nodos. La energía
usa una derivada covariante de cuarto orden y un Laplaciano positivo compatible;
la fuerza y la corriente proceden de esa misma energía. La corriente de una
ruta de dos enlaces se reparte en las caras físicas correspondientes.
El término positivo evita el modo alternante nulo de una derivada central sola.

Se mantienen el modelo continuo, el catálogo R2, la regularización y los seis
perfiles físicos. La energía y las fuerzas tienen orden formal cuatro; la
corriente conjugada de caras se contrasta conservadoramente con orden esperado
dos, además del orden medido. El objetivo del 1 % no cambia. La
[prerregistración](registration.json) y sus aclaraciones previas al cálculo
quedan preservadas.

## Resultados nuevos y cálculo pendiente

- 102 pruebas aprobaron en 7,26 s en Geminga.
- Piloto de ocho celdas aprobado en 33,2 s.
- Diferencias absolutas de fuerza y corriente frente a variaciones de energía:
  1,35e-11 y 1,14e-12.
- El sesgo de corriente uniforme del piloto es 0,535 % respecto del límite
  continuo. Esto no certifica aún los perfiles no uniformes ni su convergencia.

El [informe con figuras](../review_20260923/Informe_revision_espacial_etapa_3_20260923.md)
([PDF](https://github.com/JoaquinDiazM/pysnspd/blob/c80c0f8612d8fd5a2b390938ca5ed4b9d58d3caa/output/pdf/implementation/Informe_revision_espacial_etapa_3_20260923.pdf))
separa el lote anterior, el piloto nuevo y lo pendiente.

El bloque activo de `/home/jdiaz/GEMINGA_COMMANDS.md` ejecuta los seis casos
en 8/16/32 celdas, reutilizando únicamente el piloto nodal verificado.
Duración estimada: unos siete minutos, orientativamente 6-12; un proceso CPU,
menos de 1 GB RAM estimado. El agente no lo lanzó. Hay barras de caso y lote,
tiempo transcurrido, ETA y registro `progress.jsonl`.

Las diferencias independientes de energía se ejecutan por caso en N=16 y en
el piloto N=8; las otras mallas declaran que ese control no se repitió. Se
mantienen en todas las mallas las referencias uniformes, fase, signos y
respuestas espaciales. El contraste espectral usa puntos representativos N=32.

## Qué falta para culminar la etapa 3

Este lote decide la precisión del bloque estático 3A. Después siguen 3B
(bordes y reservorios), 3C (potencial y conservación de carga), 3D (circuito de
tres estados de la memoria) y 3E (dinámica débil con depósito sintético).
No hay todavía señal del detector ni promoción a producción.

La entrega anterior completa está archivada, incluida su libreta. El nuevo
`verify_nodal_delivery.py` verifica ese archivo histórico, el cierre de etapa 2,
las fuentes y el piloto actual. `check_nodal_handoff.py` comprueba el reuso y
la carpeta de salida sin calcular física. Ninguno modifica resultados.
