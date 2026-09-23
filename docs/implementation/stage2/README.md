# Etapa 2: celdas con cinética acoplada

**Estado al 22 de septiembre: desarrollo cerrado; certificado numérico estricto
de mallas incompleto. Etapa 3 preparada, no iniciada.**

La [decisión de cierre](closure_20260922/closure_decision.json), autorizada por
el usuario, acepta la evidencia disponible para continuar la implementación.
El [informe final](closure_20260922/Informe_cierre_etapa_2.md)
([PDF](../../../output/pdf/implementation/Informe_cierre_etapa_2.pdf)) reúne los
resultados, los plots y los límites que se trasladan a la siguiente etapa.

El último lote protegido completó 21 tareas, incluidas 13 trayectorias físicas.
Los controles temporales propios de una/dos celdas pasan en la malla candidata.
El máximo defecto energético escalado es 5,43766 × 10⁻⁸ y el máximo residuo
instantáneo registrado es 1,63498 × 10⁻¹⁴. El par de 2049 nodos fonónicos tiene
una diferencia de 0,004565 %, frente al presupuesto auxiliar de 0,0025 %;
conserva su FAIL original. Los
[datos auditados](practical_review_20260922/saved_results_audit.json) y el
[dictamen numérico](stage2_admission.json) no se reinterpretan como convergencia
completa de mallas ni como admisión de producción.

La implementación experimental incorpora eventos electrón-fonón conservativos,
transporte a energía fija entre espectros distintos, movilidad KWT, BGK, calor y
escape fonónico. La malla candidata tiene 630 estados electrónicos y 1025
fonónicos por celda. El catálogo R2 y el tag v1.0.0 permanecen intactos.

## Resultados estáticos y antecedentes

Los siguientes ensayos documentan la evolución de la implementación. Sus
informes y README conservan los estados y pendientes de cada fecha; la decisión
de cierre enlazada arriba define el estado de desarrollo actual.

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

## Siguiente etapa preparada

**No se recomienda otro lote largo ahora.** La [etapa 3](../stage3/README.md)
queda preparada para implementar energía espacial y estabilidad, bordes y
reservorios, conservación de carga y el circuito de tres variables de la
[memoria](../../modelo_v0_4/actualizaciones/circuito_memoria_20260922.md), antes
de pasar a dinámica débil y deposición localizada sintética. La secuencia no
se ha iniciado. Cada ensayo debe registrar sus entradas, observables y controles
propios; antes de una predicción dinámica se fija la precisión que necesita.

El certificado completo de mallas sigue pendiente. Los controles estructurales
de energía, unidades, signos y dominio se mantienen; el fallo auxiliar no se
usa como bloqueo general de trabajos con otras dependencias. La
[revisión anterior](practical_review_20260922/validation_scope_review.md) explica
esa distinción y permanece histórica. La [libreta vigente](../../GEMINGA_COMMANDS.md)
ofrece un resumen verificable del cierre y conserva los comandos anteriores.

## Evidencia y reproducción

- `acceptance_criteria.json`: umbrales fijados antes de los resultados.
- `review/complementary_grid_amendment.json`: cambio explícito de la malla,
  manteniendo todos los umbrales y la evidencia fallida.
- `NUMERICAL_DECISIONS.md`: errores encontrados y decisiones adoptadas.
- `event_summary.json`, `transport_results.json`, `closure_results.json`:
  evidencia cinética, transporte y cierres con hashes de fuentes.
- `review/`: referencias independientes, controles negativos y validación de
  los evaluadores de convergencia.
- `regression_result.json`: regresiones del checkpoint correspondiente;
  no sustituyen los resultados y alcances del informe final.
- `closure_20260922/closure_decision.json`: decisión vigente de cierre de
  desarrollo y límites del certificado numérico.
- `../stage3/entry_contract.json`: contrato de la siguiente secuencia preparada;
  `NEXT_STAGE.md` conserva la referencia de transición de etapa 2.

Los scripts en `sandbox/stage2_cells/` son diagnósticos optativos. Los operadores
están en `pysnspd.experimental`; el solver de producción no los importa. Las
entradas Debye son sintéticas y no acreditan tasas o latencias absolutas de NbN.
