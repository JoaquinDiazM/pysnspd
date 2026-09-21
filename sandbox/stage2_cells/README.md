# Diagnósticos de la etapa 2

El estado y la evidencia están en
[docs/implementation/stage2](../../docs/implementation/stage2/README.md).
La etapa está pendiente de la referencia temporal; no iniciar una secuencia
espacial ni asumir que una ejecución terminada implica aceptación.

`run_coupled.py` registra parámetros, poblaciones iniciales y hashes antes de
evaluar. El estado `COMPLETED_NOT_YET_ADJUDICATED` requiere evaluación posterior;
`INCOMPLETE` nunca pasa un criterio. El valor predeterminado fonónico es 1025.

`run_reference_progress.py` observa la misma ecuación y escribe progreso cada
100 evaluaciones. Está preparado para la ejecución manual pendiente indicada
en `/home/jdiaz/GEMINGA_COMMANDS.md`. No fue ejecutado como cálculo largo por el
agente. `assess_time_refinement.py` y `assess_grid_refinement.py` contrastan
archivos identificados por hash y todos los checkpoints comunes.

Los scripts de figuras y ensamblaje del informe quedan preparados. El informe
final se construirá después de completar las puertas pendientes; no se entrega
un PDF preliminar como si fuera el cierre de la etapa.

El directorio `history/` conserva fuentes de variantes anteriores. Para ejecutar
un runner histórico cuya raíz se calcula mediante `parents[2]`, copiarlo a un
nombre nuevo directamente en `sandbox/stage2_cells/`, sin reemplazar el actual.
Los resultados y sus configuraciones indican qué versión se utilizó.
