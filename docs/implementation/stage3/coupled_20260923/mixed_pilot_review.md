# Revisión independiente del piloto mixto

El piloto conserva su resultado de identidades instantáneas. Los 18 hashes de fuentes, los 5 artefactos de la semilla y los diez artefactos referenciados por los dos resultados coinciden. El registro guardado coincide con el vigente. Las matrices son finitas, y el campo y las poblaciones iniciales permanecen idénticos. No se ejecutaron espectros ni dinámica nuevos.

La conservación de corriente, Noether, la potencia campo–circuito, el depósito de calor y las dos derivadas direccionales verificadas son compatibles con los criterios originales. La fuerza y corriente de la perturbación tienen errores FD de 6.76e-12 y 5.43e-11. La hélice no recibió un nuevo control FD: así lo declara el registro. El margen mínimo D.36 permanece positivo. Los contrastes causales directos se limitan a tres cuadraturas por caso.

## El calor inicial procede de los extremos

| Magnitud guardada | Hélice | Perturbación suave |
|---|---:|---:|
| QDelta integrado | 1.01145942556 | 1.01153895863 |
| QDelta en los dos terminales | 1.01145898406 | 1.01145898406 |
| Fracción terminal | 99.999956351% | 99.992093773% |
| QDelta interior | 4.41495867e-07 | 7.99745647e-05 |
| Máxima velocidad material interior | 0.000360268155 | 0.0140151015 |

Estas tasas y velocidades usan las unidades normalizadas del registro. La velocidad máxima de ambos perfiles, 2.85332871, se encuentra en los extremos. El aporte terminal es exactamente el mismo en ambos estados guardados. La respuesta interior sí cambia con la perturbación; queda oculta en el calor global por el relajamiento terminal.

El dato inicial impone una hélice con q=0.1, mientras los extremos tienen carga conjugada nula y todas sus velocidades materiales quedan libres. Por ello no satisface el flujo natural terminal nulo. A grado polinómico fijo, el esfuerzo terminal integrado es de orden J y la masa temporal GLL de orden h: se espera velocidad inicial de orden 1/h y calor terminal integrado de orden 1/h. Esta es una inferencia analítica, no una ley medida con varias mallas. No implica que una trayectoria a tiempo positivo o su energía disipada integrada diverjan: puede aparecer una capa de relajación inicial.

## Decisión y siguiente ensayo

Conservar el PASS algebraico del piloto y sus datos. No recomiendo gastar el lote de seis instantáneas libres como comprobación de precisión del dispositivo: su calor inicial estaría dominado por esa incompatibilidad terminal. Primero debe fijarse y registrarse la carga de reservorio, distinguir la corriente superconductora prescrita de la corriente total del circuito y contabilizar el trabajo de las cargas y reacciones. Una reacción radial impuesta no debe ocultarse ni ajustarse a posteriori para forzar una identidad.

Después, seis instantáneas con ese contrato pueden aportar sensibilidad espacial del estado cargado y de la perturbación inducida. El caso libre sigue siendo un control histórico válido. Las instantáneas, incluso conservativas, no acreditan un transiente ni la etapa 3 completa.

En un transiente débil se deben guardar corrientes y voltajes circuitales, amplitud y q/Gamma, velocidad material, calor interior/terminal, energía electrónica y temperatura equivalente, junto al balance integrado dominio+circuito con trabajo de reservorio. Los controles de población, dominio del catálogo, estabilidad y comparación temporal deben seguir las mismas variables efectivamente integradas; no se fijan aquí tolerancias nuevas.

## Límites de esta revisión

Las energías individuales del stencil FD y todos los espectros instantáneos no están archivados. Se verificaron las estimaciones FD guardadas, sus errores y presupuesto, y la aritmética registrada del momento electrónico; no se recalcularon esas evaluaciones físicas. La coincidencia térmica inicial es una identidad de preparación e inversión, no una demostración de relajación. La conservación del depósito Joule no certifica su distribución local. No se alteraron resultados, criterios ni registros previos.

El archivo `mixed_pilot_review.json` conserva el inventario SHA256, las métricas recalculadas y las limitaciones. Reproducción: `python sandbox/stage3_spatial/coupled_20260923/audit_mixed_pilot.py` (rechaza sobrescribir la revisión existente).
