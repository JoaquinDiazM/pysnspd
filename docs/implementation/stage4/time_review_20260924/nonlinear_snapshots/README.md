# Contraste constitutivo no lineal de la trayectoria térmica

La campaña corregida terminó en **22,74 s**, con 2048 raíces espectrales nuevas. Los contrastes registrados entre la ley no lineal y su tangente cumplen el margen del 1 % de la mayor escala entre la respuesta inicial y la respuesta restante. Este resultado admite usar la tangente como referencia débil y avanzar a la integración térmica no lineal; no valida aún el trabajo no térmico ni un evento de detección.

Se evaluaron los estados guardados a 0, 0,01 y 1 ps para la base, la perturbación de amplitud y la perturbación angular. La base inicial se reutilizó: quedaron ocho estados nuevos por cada una de las 256 frecuencias de Matsubara. En cada tiempo se restó la **base no lineal del mismo tiempo** antes de comparar respuestas. Se conservaron la malla, los contactos espectrales, la acción, la movilidad KWT, el potencial normal, todas las constantes y la tolerancia espectral `1e-7`.

| Diferencia a 1 ps, respecto de la respuesta inicial | Amplitud | Fase angular |
|---|---:|---:|
| Velocidad del gap, RHS | 0,4056 % | 0,002785 % |
| Corriente | 0,07091 % | 0,002204 % |
| Torque de fase | 0,6171 % | 0,002925 % |

La cola de torque de amplitud llega a una diferencia del 124,6 % respecto de lo que queda a 1 ps, pero sólo al 0,6171 % del impulso inicial. Por eso se publican ambas normalizaciones y no se impone precisión relativa a una cola casi nula. La identidad instantánea de potencia de la ley no lineal deja un residuo absoluto máximo inferior a `1.8e-12`; representa energía libre térmica, no un balance completo de energía interna del detector.

## Fallo conservado y corrección numérica

La primera campaña se detuvo en la búsqueda de descenso de Newton. Para la frecuencia `n=9`, estado angular a 1 ps, el residuo era `1.007e-7`, todavía por encima de la tolerancia `1e-7`. El paso Newton realmente reducía la acción en aproximadamente `3.153e-15`, pero restar dos energías próximas a `-62.296` producía por redondeo una diferencia positiva `1.421e-14`. El criterio de Armijo rechazaba entonces un descenso correcto.

El nuevo backend `thermal_stable_newton.py` evalúa la **misma diferencia de acción mediante una identidad algebraica estable**, antes de sumar o restar el fondo grande. Conserva el residuo, el Jacobiano, el criterio de Armijo, la tolerancia y los contactos. No se modificó el backend histórico. El caso que reproducía el fallo está guardado en [roundoff_fixture.npz](roundoff_fixture.npz), junto con su [diagnóstico](n09_diagnosis.json). Las cuatro pruebas de diferencias exactas y Newton estable pasaron; véanse [el recibo](test_receipt.json) y [el registro](tests.log).

| Campaña | Plan y código | Evidencia original en Geminga |
|---|---|---|
| Intento detenido, conservado | [plan.json](plan.json), `thermal_nonlinear_snapshots.py` | `/home/jdiaz/scratch/stage4_nonlinear_snapshots_20260924`: fallo y 45 modos terminados |
| Campaña completa con diferencia estable | [stable_plan.json](stable_plan.json), `thermal_nonlinear_stable_snapshots.py` | `/home/jdiaz/scratch/stage4_nonlinear_stable_20260924`: 256 modos, 2048 raíces nuevas |

## Evidencia y alcance

El [resumen](summary.json) contiene los resultados por estado, las dos normalizaciones y los residuos de cada raíz. Los [campos compactos](nonlinear_comparison_fields.npz) conservan fuerzas, corrientes, energía, velocidad y potencial; la [identidad](identity.json) registra insumos, fuentes y recursos. El [certificado](verification_receipt.json) verificó los 256 archivos espectrales completos (268 275 142 bytes), los 265 insumos y las diez fuentes, y acredita también el intento anterior. Los campos espectrales completos permanecen en Geminga con sus hashes.

Estas son evaluaciones constitutivas sobre estados de una trayectoria ya existente. No son una nueva integración temporal. El siguiente control integra la ley térmica no lineal con ETD2, conserva todos los nodos y compara dos discretizaciones temporales; el Jacobiano guardado se usa sólo como parte del integrador exponencial. El trabajo espectral fuera del equilibrio, el fotón y el circuito con puertos reales siguen fuera de este contraste.
