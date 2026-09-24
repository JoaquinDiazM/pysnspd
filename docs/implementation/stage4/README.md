# Etapa 4: núcleo, referencia física y ejecución paralela

La [continuación actual](followup_20260923/README.md) recoge las dos corridas
finalizadas y la comprobación del nuevo ejecutor paralelo. **La etapa 4 sigue
abierta; la referencia física ya calculada exige diagnosticar el cierre del núcleo.**

El signo negativo del núcleo grueso con δ=0,05 desapareció al resolver mejor
el gradiente. Para δ=0,10, pasar de 561 a 2.145 nodos cambia la fuerza interior
un 0,572 % y el calor interior del par KWT Korzh un −2,551 %. En cambio, variar
δ de 0,10 a 0,05 a igual malla cambia ese calor un +21,24 %. No se propone seguir
refinando indefinidamente este perfil; la sensibilidad del cierre necesita una
referencia física.

El piloto serial/paralelo conserva exactamente los resultados y los 33 arrays
de cada caso: 19,5554 frente a 8,71745 s con cuatro trabajadores. El ejecutor
comparte las consultas espectrales entre casos y admite hasta 27 trabajadores
más el coordinador en Geminga, reservando dos núcleos físicos y respetando
los límites efectivos de CPU y memoria. Conserva barras, ETA y procedencia.

La [referencia radial Usadel](followup_20260923/radial_reference_README.md)
resolvió 512 problemas espectrales y seis evaluaciones del candidato en
4,852 s con el pool compartido. La fuerza del candidato δ=0,10 difiere
158,09 % en norma L² radial frente a la referencia de 256 frecuencias;
la variación por corte espectral es 1,87 %. El núcleo físico no queda
admitido: el paso útil es examinar el cierre antes de lanzar una dinámica
costosa o repetir refinamientos. No se ajusta δ ni se cambia automáticamente
el modelo. El transiente débil sin fotón y sus balances siguen pendientes.
Las condiciones de borde fijas de los controles anteriores no se presentan
como contactos del dispositivo.

La [política de horizonte](followup_20260923/horizon_policy.md) para futuros
disparos fotónicos está autorizada, todavía sin implementar. No altera las
ecuaciones ni el circuito de la memoria. No hay predicción admitida de hotbelt,
latencia o jitter, ni promoción a producción.

El banco GLL 2D sigue siendo experimental; producción conserva
Delaunay–Voronoi y su integración existente. Se preservan la
[procedencia numérica](review_20260923/solver_scope.md), el
[informe anterior](review_20260923/Informe_avance_etapa_4A.md) y la
[preparación inicial](start_20260923/README.md).
