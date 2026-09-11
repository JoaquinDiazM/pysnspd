# Revisión independiente — etapa 1, iteración 2

**Dictamen final: `PASS_SAMPLED_ELECTRONIC_DOMAIN`; material NbN: `REJECTED`.** El certificado externo es `../catalog_admission.json`. Sus 49 puertas enlazan el NPZ exacto, los criterios y la evidencia. Las referencias se calcularon antes de seleccionar el catálogo final; las tolerancias permanecieron sin cambios.

El catálogo final supera 234 comparaciones físicas, 150 identidades de derivación, 476 consultas de soporte que incluyen los 326 fallos históricos y 6.855 consultas densas nuevas. No queda ningún rechazo en esas pruebas. El control de la cúbica en todo cada intervalo de la coordenada de depairing, a 65 amplitudes muestreadas, revisa 9.260 celdas: su menor incremento es `6.99e-10`, positivo. Cuatro cotas suficientes de Bernstein son negativas; los mínimos reales de esas cúbicas son positivos, por lo que no son fallos del espectro.

Los máximos errores finales de energía, fuerza de amplitud y respuesta Γ son `4.97e-7`, `7.60e-6` y `5.03e-4` escalados. Los errores relativos máximos de las fuerzas son 0,0265 % y 0,0513 %. En el caso crítico Γ=0, la respuesta es `0.678800868` frente a `0.678584013`: error 0,03196 %. El siguiente paso que estos resultados permiten es el ensayo sintético de una o dos celdas.

El contrato congelado es `../acceptance_criteria.json`, SHA-256 `ef8164fe9f7d7d27f4ecbba1a318fc8642ea3dcc9723dd2e1a9d9c1447dc6bfb`. La energía debe tener error escalado ≤`1e-5`; las dos fuerzas, error escalado y relativo para referencias no nulas ≤`1e-3`. Se evalúa el error total contra el límite causal, además de separar interpolación, cuadratura, corte y regulador. El dictamen se aplica al dominio y poblaciones ensayados; no es una cota uniforme para cualquier distribución.

## Referencias que no importan el catálogo

- `causal_reference.py/json`: 150 casos fijos, cinco poblaciones y 30 pares de campos. Integra una parametrización real de la solución causal con η=0 y un jacobiano analítico de la cuenta de estados. No usa la ecuación cuártica compleja ni un regulador. Diferencia máxima al refinar corte y cuadratura: `1.12e-14`; tiempo local: 1.41 s.
- `gamma_zero_reference.py/json`: límite BCS exacto y dependencia con η. En Γ=0 se incluye la contribución integrada del borde espectral a la respuesta Γ. El valor puntual del núcleo en η=0 no sustituye esa contribución.
- `thermal_reference.py/json`: nueve casos FD. Una integral causal real y sumas de Matsubara con extrapolación ultravioleta independiente coinciden dentro de `3.16e-13`. La variación de la extrapolación es ≤`5.21e-12`; tiempo local: 0.68 s. Se compara energía libre tras restar entropía y fuerzas con las ocupaciones fijas durante la variación.
- `tail_reference.py/json`: 318 integrales adaptativas de las colas hasta infinito, para cortes 12 y 20. La mayor contribución restante es `3.74e-11`; cada cola satisface el presupuesto de error de referencia del contrato.

Los JSON contienen los puntos, resultados, errores estimados y hashes del código y del contrato. Los tiempos son mediciones locales; no representan tiempos de un transiente.

## Controles del código del candidato

`assess_candidate.py` recibe el NPZ exacto, registra su hash y evalúa 234 casos físicos: 150 fijos, 60 adicionales fuera de nodos, nueve FD y 15 normales. Conserva todos los fallos. Comprueba además soporte, Pauli, factores SI, paridad de corriente, identidades de derivación y Maxwell, y guardado/carga sin resolver espectros nuevos. El control de derivadas para Γ diminuto usa paso complejo de una representación polinómica local equivalente; no resta energías casi idénticas con pasos de `1e-14`. Este control verifica la diferenciación de la representación, mientras la comparación causal verifica su precisión física.

`spectral_crosscheck.py/json` vuelve a evaluar la ecuación de Usadel sin elevar al cuadrado, la normalización, la positividad y la inversión de cuenta para 90 combinaciones de campos y η, con puntos próximos al borde espectral. Los fallos del comparador escalar se registran separadamente de las puertas del constructor por lotes: ausencia de comparación nunca se presenta como acuerdo. Los resultados deben regenerarse si cambia el módulo.

La consulta escalar actual delega en el mismo algoritmo por lotes; su igualdad comprueba la API, no independencia física. La revisión causal real y Matsubara conserva esa independencia. El control espectral actualizado supera las 90 combinaciones, incluidas η=`3e-9`, con residuo de Usadel ≤`4.16e-16` e inversión de cuenta ≤`7.97e-13`; los resultados anteriores que detectaron regresiones quedan archivados.

## Pilotos rechazados y regresiones conservadas

`direct_17_assessment.json` encontró 248 rechazos de soporte en 2.352 consultas densas: 140 energías no positivas y 108 pérdidas de orden entre estados. Al alinear el eje interno con Γ/|Δ|, `ratio_17_assessment.json` dejó 40 pérdidas de orden en 4.820 consultas, todas dentro del sector con brecha. La alineación eliminó los valores negativos muestreados, pero por sí sola no cerró la admisión.

`ratio_bulk_17_assessment.json` ya supera los 234 valores físicos contra el límite causal, incluidos FD y el caso crítico de Γ=0. Sus máximos de error escalado son `4.97e-7`, `7.60e-6` y `5.03e-4` en energía, fuerza de amplitud y respuesta Γ, respectivamente. Sin embargo, todavía rechaza 60 de 6.519 consultas densas y 43 puntos históricos. Su estado sigue siendo **FAIL**. Una energía integrada precisa no justifica admitir un espectro que pierde el orden de sus estados.

`support_regression_points.json` conserva cada punto donde se observó un rechazo. Esos puntos se repiten después de cambiar la malla o la representación, además de las nuevas muestras en cada celda. No se desplazan ni eliminan para obtener aprobación. El refinamiento localizado posterior resolvió todos los rechazos observados; `final_assessment.json` contiene la revisión del archivo final.

## Alcance de una posible aprobación

La aprobación de las puertas electrónicas permite preparar el ensayo sintético de una o dos celdas previsto en D.4.2. Los datos fonónicos de NbN permanecen rechazados y el catálogo continúa separado del solver de producción. Las dificultades independientes del núcleo positivo y del símbolo principal identificadas en v0.4 siguen abiertas.

La admisión se expresa en un certificado externo que referencia los hashes del NPZ y del contrato. El NPZ conserva su estado de construcción pendiente de ese certificado, incluidas unidades, soporte, procedencia y método. No se reescribe el archivo después de medirlo sólo para cambiar su estado y romper la trazabilidad.
