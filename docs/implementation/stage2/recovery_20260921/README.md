# Recuperación numérica de la etapa 2

**Estado: preparado para relanzamiento autorizado. La etapa 2 sigue abierta.**

El lote anterior completó las comparaciones temporales RK4 de una y dos celdas
en la malla candidata de 630 estados electrónicos y 1025 fonónicos. La
[auditoría de archivos y resultados](completed_batch_audit.json) verificó los
errores registrados, las condiciones iniciales, las fuentes y la reutilización
de trayectorias. Sus referencias de 320 pasos se conservan completas en
`audit_inputs/`; el plan nuevo usa esas copias archivadas.

## Por qué se cambia el integrador

La trayectoria con 2520 estados electrónicos falló durante una etapa interna
de RK4: una ocupación de la cola electrónica salió del dominio físico. El
diagnóstico está en [first_step_diagnosis.json](first_step_diagnosis.json).
El resultado fallido y su estado inicial permanecen disponibles.

Se probó después reducir el paso por bisección. El ensayo agotó 24 evaluaciones
del lado derecho sin aceptar ningún avance. El último rechazo registrado ocurrió
con paso 9,765625 × 10⁻⁵ y una ocupación de −1,57 × 10⁻²⁷ en el último estado
electrónico de la segunda celda. Esta cantidad es muy pequeña, pero el método
necesita poblaciones admisibles para evaluar sus productos de ocupación. No se
convirtió el dato negativo en cero. El
[registro de bisección](positive_first_interval_probe.json) documenta el límite
del ensayo; no demuestra que toda reducción posible de paso vaya a fallar.

## Ensayo registrado

La nueva propuesta combina SSPRK3, un método explícito de tres etapas, con un
factor de disponibilidad θ para cada evento. Si un evento consume electrones,
huecos o fonones, el factor tiene en cuenta todas las cantidades que puede
consumir. El mismo θ multiplica sus aportes a ambos sectores y a su registro
de energía. Así se mantiene la identidad de intercambio por evento.

Esto cambia el procedimiento de integración. Los núcleos físicos originales,
el catálogo R2 y los umbrales de aceptación permanecen intactos. No hay recorte
posterior de poblaciones ni corrección artificial del balance energético.
La conservación por evento tampoco garantiza por sí sola la precisión temporal
de la trayectoria completa: debe medirse la convergencia y registrarse cuándo
actúa el limitador.

Las tres pruebas cortas contra DOP853 pasaron: el error máximo del caso fino
fue 1,5798 × 10⁻⁶, con reducciones de 3,74 y 4,13 al refinar el paso. La
prueba del primer intervalo fino también pasó: mantuvo poblaciones físicas,
con error del balance acumulado de 6,76 × 10⁻¹⁰. El defecto de flujo ponderado
por energía que introdujo el limitador en ese intervalo fue 8,92 × 10⁻²³.
Estos resultados cubren un intervalo corto; no anticipan la aceptación del
lote completo. La regresión independiente del repositorio pasó sus 509 pruebas.

El [plan de relanzamiento](limited_acceptance_plan.json) ordena el trabajo:

1. Comparar tres resoluciones SSPRK3 para una y dos celdas con las referencias
   RK4 de 320 pasos ya archivadas, verificando que comparten el problema físico.
2. Ensayar la malla electrónica fina y su paso reducido. Este control pareado
   se identifica como diagnóstico suplementario; no sustituye tres niveles
   temporales ni emite por sí solo una admisión.
3. Comparar las tres mallas electrónicas y las tres fonónicas manteniendo el
   mismo integrador y paso. La propia malla candidata debe pasar su límite.
4. Conservar las trayectorias de equilibrio y fronteras de población y revisar
   los campos y el soporte en al menos diez tiempos físicos por celda.

El lote para ante el primer fallo, sin reintentos automáticos ni borrado de
salidas. Su finalización tampoco cierra automáticamente la etapa: quedan la
revisión de resultados y el dictamen independiente. No hay un informe final
nuevo ni una validación del transiente completo del detector en esta entrega.

## Ejecución y alcance de la autorización

El usuario autorizó expresamente al agente a relanzar este cálculo largo en
`screen code_000`, desacoplar la sesión y responder sin esperar ni sondear su
estado. La estimación preparada es de 3–6 horas, incierta, con un hilo y una
reserva de 12 GiB. El comando y las rutas de salida deben quedar en
`/home/jdiaz/GEMINGA_COMMANDS.md`; la siguiente revisión leerá esas salidas.
Esta excepción corresponde a la recuperación actual, no a otros cálculos
largos del proyecto.

El ejecutor
[`run_limited_acceptance_batch.py`](../../../../sandbox/stage2_cells/recovery_20260921/run_limited_acceptance_batch.py)
funciona en primer plano dentro de la sesión. Por defecto sólo imprime el
plan (`--dry-run`); ejecutar exige el plan revisado y los hashes de fuentes y
referencias congelados. Sus pruebas de control sintéticas están registradas
en [limited_batch_offline_tests.json](limited_batch_offline_tests.json).

El alcance sigue siendo numérico y sintético: una o dos celdas. La calibración
absoluta de NbN, los gradientes espaciales del condensado, el circuito y los
transientes completos de SNSPD requieren sus propias validaciones.
