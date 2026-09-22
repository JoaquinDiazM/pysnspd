# Etapa 2: celdas con cinética acoplada

**Estado: preparado para relanzamiento autorizado; etapa no cerrada ni promovida a producción.**

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
  a 240 segundos. La ejecución manual posterior terminó en 1164,55 s. Tres
  ensayos cortos RK4 pasaron contra ella; el fino tardó 24,94 s y alcanzó un
  error máximo de 1,57 × 10⁻⁶. La salida de screen y las trayectorias se
  conservan en [la continuación](resume_20260921/README.md).
- En esa continuación también pasaron la comparación continua al corte 0,005
  y los controles de un piloto completo de dos celdas con 1025 nodos fonónicos.
  El lote posterior aprobó la convergencia temporal RK4 de una y dos celdas
  en la malla candidata. Los archivos, las condiciones iniciales y la
  reutilización se verificaron en
  [completed_batch_audit.json](recovery_20260921/completed_batch_audit.json).
- El refinamiento electrónico de 2520 estados falló una condición de población
  física en una etapa interna de RK4. Un diagnóstico con bisección agotó su
  presupuesto sin avanzar. Ambos resultados se conservan; no se recortaron
  poblaciones ni se reinterpretó el fallo como una aprobación.

## Continuación necesaria

La [recuperación preparada](recovery_20260921/README.md) registra un ensayo
SSPRK3 con un factor común para los flujos de cada evento. Primero debe superar
la comparación temporal frente a las referencias RK320 archivadas; después se
evalúan las mallas, las fronteras de población y los campos visitados. El lote
se detiene ante cualquier fallo y conserva la evidencia.

Para esta recuperación, el usuario autorizó expresamente lanzar el cálculo
largo en `screen code_000`, desacoplar la sesión y responder sin esperar ni
sondear el proceso. Los comandos se registran en `/home/jdiaz/GEMINGA_COMMANDS.md`.
La regla habitual de ejecución manual sigue vigente para otros cálculos largos.

Faltan el dictamen temporal del nuevo integrador, los refinamientos dinámicos,
los controles sobre esas trayectorias y el informe final. El dictamen de puertas
conservado como punto de partida está en
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
