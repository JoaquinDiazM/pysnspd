# Etapa 3: empalme espacial y acoplamiento instantáneo

La campaña abierta ejecutada por el usuario terminó sus seis casos en 271,28 s.
La auditoría `open_audit.json` permite continuar: los errores de corriente y
gradiente terminal disminuyen con la malla. Todas las soluciones partieron y
terminaron en la hélice uniforme; no se atribuye a ese lote una relajación dinámica.

El nuevo empalme identifica la traza transversal uniforme del condensado entre
un rectángulo 2D y dos continuaciones 1D. Conserva las cuadraturas espectrales de
cada sector. Las poblaciones no se identifican por índice: la continuidad cinética
a igual energía y el reparto temporal de fonones siguen pendientes.

El piloto de 45 grados de libertad y 55 cuadraturas pasó en 78,34 s. Evalúa fuerza,
corriente, potencial, KWT, depósito de calor y el circuito de tres estados de la
memoria en un mismo instante. Las diferencias finitas se comprueban sólo en el
perfil perturbado, manteniendo sus poblaciones fijas. El balance CM.9 cierra con
un residuo absoluto máximo de 3,02e-24 W. No integra una trayectoria.

## Decisión sobre los extremos

Los extremos libres originan el 99,999956 % del calor del condensado en la hélice.
No se solicita una campaña larga para converger ese calor: el perfil con corriente
no satisface la condición natural de carga nula. Los datos se conservan como
control de las identidades. El siguiente diagnóstico usa una carga de fase
prescrita por la corriente del reservorio y registra su trabajo explícitamente.
La amplitud terminal fija es una condición de este diagnóstico instantáneo;
no sustituye la rama dependiente de la corriente durante un transiente.

El diagnóstico con carga pasó: el calor normalizado de la hélice baja de 1,01146
a 6,19216e-7 y el del perfil perturbado queda en 8,01523e-5. El residuo CM.9
registrado es como máximo 4,85e-27 W; la recomputación independiente da 5,66e-27 W.
La corriente normal en la arista interior adyacente es -0,041935 % de Iref:
se informa sin confundirla con el flujo externo. Se reprodujo en Geminga en
0,046 s sin consultas espectrales nuevas. Cinco pruebas y seis subpruebas
adicionales del módulo de reservorio pasaron allí en 0,09 s.

## Decisión numérica

Se acepta `SeededCountCatalog`: una tabla proporciona la estimación inicial y
el mismo solver causal corrige la energía, antes de evaluar las derivadas
implícitas. Nueve puntos de contraste dan una aceleración mediana de 1,904 veces
por kernel y diferencias relativas de energía de hasta 1,35e-14. No se infiere
una aceleración del RHS o de la trayectoria completa.

Dos aproximaciones de energía anteriores se conservan rechazadas: el polinomio
directo cruzaba niveles; reconstruir incrementos positivos preservaba el orden
pero no las fuerzas. Nunca se utilizan esas energías ni fuerzas en el piloto.
El archivo del polinomio rechazado se utiliza exclusivamente como semilla interna
de Newton. La caja de campos es explícita y no se extrapola.

## Evidencia y reproducción

- `raw/stage3b_open_full_20260923/`: lote abierto original, completo.
- `raw/patch_check/` y `raw/monotone_check/`: alternativas rechazadas y fuentes exactas.
- `seeded_check/`: comparación del acelerador causal y tabla utilizada como semilla.
- `raw/mixed_pilot/`: dos estados mixtos, fuerzas, espectros, calentamiento y balances.
- `mixed_registration.json`: criterios y alcance fijados antes del piloto.
- `coupled_tests.log`: 102 pruebas y 17 subpruebas aprobadas en 0,59 s en Geminga.
- `reservoir_postprocess_v2/`: respuesta con carga y balance del trabajo de borde.
- `reservoir_reproduction_geminga/`: reproducción algebraica en el entorno de ejecución.
- `reservoir_tests.log`: cinco pruebas y seis subpruebas adicionales en Geminga.
- `mixed_pilot_review.json` y `reservoir_postprocess_review.json`: auditorías desde los datos.
- `reservoir_postprocess/`: intento detenido al comparar masas entre plataformas
  bit por bit; conserva sus fuentes. La recuperación usa las masas archivadas
  y sólo admite redondeo al contrastarlas con la geometría reconstruida.
- `Informe_empalme_y_acoplamiento_etapa_3_20260923.md`: informe con figuras.
- `delivery_manifest.json`: inventario reproducible de esta entrega.

El comando vigente se encuentra en `/home/jdiaz/GEMINGA_COMMANDS.md`, reflejado en
`docs/GEMINGA_COMMANDS.md`. La libreta anterior se conserva íntegra en
`GEMINGA_COMMANDS_before_coupled.md`. No hay instrucciones de gestión de sesiones.

La etapa 3 sigue en desarrollo: faltan el ensayo espacial con reservorios,
el transporte cinético del empalme y la evolución débil acoplada. El siguiente
contraste temporal debe concentrarse en la señal inducida y el balance integrado,
sin repetir identidades algebraicas ya verificadas. El modelo experimental no
se activa en producción y el tag `v1.0.0` permanece intacto.
