# Auditoría independiente del lote mixto completo

Dictamen: **VERIFIED_SIX_SNAPSHOT_CAMPAIGN_SCOPE_LIMITED**. El lote terminó en 355.651 s; contiene seis estados instantáneos, dos perfiles sobre tres mallas. La revisión usa sólo archivos guardados y no ejecuta física nueva.

## Integridad e identidades

Se cotejaron las fuentes del corredor (6), del cálculo espacial (18) y del postproceso (12), cinco artefactos de la semilla y 38 entradas del postproceso. Los casos coinciden con sus resúmenes; el inventario contiene 61 archivos. Incidencias: 0.

Campo, poblaciones y masas permanecen idénticos entre la instantánea original y su lectura con reservorio. Se recompusieron momentos energéticos, carga de fase, trabajo de borde, continuidad, reacción radial, ecuaciones circuitales y CM9. Los criterios originales se mantienen. El control de diferencias finitas se ejecutó sólo en m0 perturbada, como estaba registrado; no se atribuye a los otros cinco casos.

## Comparación de las mallas

| Perfil/malla | Energía [bar] | QDelta [bar] | Velocidad máxima [bar] | Vdev [V] | max abs(Jn)/Iref adyacente |
|---|---:|---:|---:|---:|---:|
| m0_helix | -42.2815831992 | 6.19215888e-07 | 0.0011960404 | 9.04496123e-07 | 0.000419349496 |
| m0_smooth_perturbation | -42.2815663219 | 8.01522847e-05 | 0.0140151015 | 9.15340354e-07 | 0.000419349496 |
| m1_helix | -42.2815786379 | 1.00783946e-09 | 3.44242738e-05 | 3.04579618e-09 | 6.03483039e-06 |
| m1_smooth_perturbation | -42.2815635833 | 7.787294e-05 | 0.0125388243 | 3.6945357e-09 | 6.03483039e-06 |
| m2_helix | -42.2815786225 | 3.4088131e-12 | 1.03180245e-06 | 1.11272743e-11 | 9.04413102e-08 |
| m2_smooth_perturbation | -42.2815635517 | 7.79071892e-05 | 0.0124186646 | 6.60197464e-10 | 9.04413102e-08 |

Las etiquetas [bar] indican magnitudes normalizadas del modelo; no representan la unidad de presión. La hélice cargada reduce fuertemente su respuesta residual al refinar. La perturbación retiene una respuesta interior finita. Para evitar que el error de la hélice o una gran energía de fondo oculten esa respuesta, se compara además perturbación menos hélice:

| Exceso inducido | m0 | m1 | m2 | Cambio m1→m2 / valor m2 |
|---|---:|---:|---:|---:|
| energy_bar | 1.687728056e-05 | 1.50545884e-05 | 1.507088756e-05 | 0.10815% |
| condensate_heat_loaded_bar | 7.953306883e-05 | 7.787193218e-05 | 7.790718576e-05 | 0.0452507% |
| Vdev_V | 1.084423128e-08 | 6.487395159e-10 | 6.490701901e-10 | 0.0509458% |

El menor margen D.36 observado es 1.57079604; el residuo máximo de CM9 recompuesta es 2.044e-25 W. El trabajo de reservorio queda incluido con su signo. El JSON conserva los valores de trabajo, calor terminal, incertidumbre de D.36 y las diferencias por malla.

Estos cambios son sensibilidad observada. El registro no fijó un umbral de precisión multimalla y aquí no se inventa uno. Una disminución del defecto de corriente en la cara interior próxima al extremo tampoco acredita por sí sola la traza normal externa de D.27. Las energías de fondo se informan, pero no se usan para ocultar el tamaño relativo del exceso inducido.

## Cierre de desarrollo y pendientes

La evidencia respalda el cierre de desarrollo acotado que el usuario autorizó y el inicio de 3.5 como investigación de rangos y márgenes. Ese cierre conserva los resultados y sus límites. No equivale a completar una validación temporal, toda D.27 o la interfaz cinética; los campos stage3_closed=false de los datos originales permanecen intactos.

Antes de atribuir resultados a transientes completos siguen pendientes el radio de reservorio dependiente de Is y su derivada, el intercambio cinético y energético en la unión 2D–1D, la condición normal externa apropiada y una trayectoria débil acoplada con balance integrado y comparación temporal. El depósito localizado necesita su propio ensayo. Estas obligaciones se registran para las etapas dependientes; no se ejecutan barridos ni se bloquea la investigación 3.5 por no disponer aún de ellos.

El contrato de entrada histórico describía la preparación inicial y contiene flags de etapa no iniciada. No debe citarse como estado actual ignorando la evidencia posterior. La decisión actual debe referenciar esta auditoría y el cierre de desarrollo separado.

Reproducción: `python sandbox/stage3_spatial/closure_20260923/audit_full_reservoir.py`. Sólo aritmética e integridad sobre archivos guardados; rechaza sobrescribir esta revisión.
