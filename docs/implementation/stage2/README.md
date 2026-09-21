# Etapa 2: celdas con cinética acoplada

**Estado: pendiente de validación temporal; no cerrada ni promovida a producción.**

La implementación experimental incorpora eventos electrón-fonón conservativos,
transporte a energía fija entre espectros distintos, movilidad KWT, BGK, calor y
escape fonónico. La malla candidata tiene 630 estados electrónicos y 1025
fonónicos por celda. El catálogo R2 y el tag v1.0.0 permanecen intactos.

Resultados estáticos y ensayos realizados:

- Transporte: error máximo de respuesta de 0,00395 % con 630 estados; disminuye
  al refinar a 1260/2520. La versión de 180 estados fallaba el límite de 0,1 %.
- Colisiones: comparación electrónica independiente por debajo de 0,1 % en los
  casos registrados. El control de distribución fonónica exige 1025 nodos;
  513 no pasó todos los controles. Ver `EVENTS_STATUS.md` para los alcances
  distintos de los cortes infrarrojos y las referencias.
- Energía y fuerzas: la malla se verifica con referencias independientes,
  incluidas la respuesta a Gamma y una banda estrecha próxima al borde.
- Los ensayos preliminares de una y dos celdas conservan energía y mantienen
  ocupaciones físicas. Usan 513 nodos fonónicos y **no certifican la candidata
  final de 1025**. Se conservan íntegros en `trajectories/` y `pilots/`.
- Una referencia adaptativa DOP853 quedó incompleta bajo la ejecución limitada
  a 240 segundos. No se repitió ni se fraccionó para eludir el límite.

## Continuación necesaria

Ejecutar el comando manual con progreso de [command_addendum.md](command_addendum.md),
también añadido a `/home/jdiaz/GEMINGA_COMMANDS.md`, y devolver la salida con
`reinicia`. El registro exacto de la ejecución incompleta está en
[reference_handoff.json](reference_handoff.json).

Faltan el dictamen temporal frente a la referencia, los refinamientos dinámicos
en la configuración final, los controles de soporte sobre esas trayectorias y
el informe final. El dictamen completo de puertas está en
[stage2_admission.json](stage2_admission.json). No avanzar a la etapa espacial
mientras siga pendiente.

## Evidencia y reproducción

- `acceptance_criteria.json`: umbrales fijados antes de los resultados.
- `review/complementary_grid_amendment.json`: cambio explícito de la malla,
  manteniendo todos los umbrales y la evidencia fallida.
- `NUMERICAL_DECISIONS.md`: errores encontrados y decisiones adoptadas.
- `event_summary.json`, `transport_results.json`, `closure_results.json`:
  evidencia cinética, transporte y cierres con hashes de fuentes.
- `review/`: referencias independientes, controles negativos y validación de
  los evaluadores de convergencia.
- `regression_result.json`: pruebas del repositorio, independientes del cálculo
  temporal que falta completar.
- `NEXT_STAGE.md`: contrato preparado para la etapa espacial; su ejecución
  depende del cierre de esta etapa.

Los scripts en `sandbox/stage2_cells/` son diagnósticos optativos. Los operadores
están en `pysnspd.experimental`; el solver de producción no los importa. Las
entradas Debye son sintéticas y no acreditan tasas o latencias absolutas de NbN.
