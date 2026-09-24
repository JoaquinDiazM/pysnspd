# Resultado físico de la referencia radial térmica

La diferencia del núcleo **no se resuelve refinando la malla ni ajustando el reloj de KWT**. La referencia espacial independiente terminó en 4,85 s y muestra una discrepancia grande del cierre local regularizado, incluso cuando se comparan fuerzas a la misma temperatura y sobre el mismo perfil. No se cambió todavía el modelo físico ni producción.

## Qué se comparó y qué se comprobó

Se prescribe d(r)=d_eq tanh(r/ℓ₀), fase azimutal y potencial vector nulo. Ambos lados usan energía libre térmica a 0,9 K. La referencia resuelve el espectro espacial de Usadel; el candidato conserva C.9 y su regularización finita. Esto no es un vórtice autoconsistente, una barrera de nucleación ni la cinta completa de Korzh.

Se verificaron los 521 archivos registrados, el plan y las fuentes. Los dos cortes de Matsubara cambian la fuerza de referencia 1,868 % en norma L² sin estimación de cola y 0,0800 % al incluir por separado el término asintótico principal. Aumentar el radio exterior de 8ℓ₀ a 12ℓ₀ cambia la fuerza interior aproximadamente 1,06×10⁻⁶ en términos relativos. Son sensibilidades observadas, no cotas rigurosas ni intervalos de confianza física.

También se comprobó de forma independiente que la fuerza del candidato es la derivada de su propia energía. Se perturba la amplitud con funciones de soporte compacto y se compara la variación numérica de energía con la integral de X·δd. Los tres regularizadores y una perturbación más extensa coinciden hasta 6,47×10⁻¹⁰ relativo. Por tanto, este resultado no apunta a un signo mal programado o a una derivada omitida: señala una limitación de la energía elegida.

## Magnitud y sentido del desacuerdo

Las normas usan el peso radial 2πr en 0,001≤r/ℓ₀≤4.

| δ/Δ₀ | Diferencia L² frente a suma Usadel N=256 | Frente a referencia con estimación principal de cola |
|---|---:|---:|
| 0,05 | 232,95 % | 228,79 % |
| 0,10 | 158,09 % | 155,20 % |
| 0,20 | 75,37 % | 73,46 % |

La observación más importante no es únicamente una diferencia de magnitud. Para δ=0,10, a r=0,1161ℓ₀≈0,546 nm, el candidato da X=−11,879 y la referencia espacial da X=+0,20265, o +0,20907 con la estimación principal de cola. X es la derivada de energía; con movilidad radial positiva la evolución va en el sentido −X. Así, sobre ese mismo perfil, los cierres indican sentidos opuestos de relajación local. Multiplicar la movilidad por otro factor positivo no puede corregirlo.

El libre recorrido medio no está acreditado y los máximos subnanométricos no deben interpretarse como fuerzas microscópicas medidas. El diagnóstico tampoco depende exclusivamente de esos máximos: para δ=0,10, excluyendo todo r<0,5ℓ₀≈2,35 nm todavía queda una diferencia L² del 36,35 %; excluyendo r<ℓ₀≈4,70 nm queda 27,21 %. Estas dos cifras emplean la referencia con estimación principal de cola.

## Qué parte de la energía produce el pico

Como control algebraico se separó una energía local BCS con gradiente cartesiano K₀ ordinario. Su fuerza es X_ordinaria=X_BCS(d,0)−2K₀L₁d. No se propone adoptar ese control como modelo.

En el máximo negativo de δ=0,10:

| Contribución | Fuerza normalizada |
|---|---:|
| Referencia espacial, con estimación de cola | +0,2091 |
| BCS local más gradiente K₀ ordinario | +0,0530 |
| Corrección al introducir la respuesta uniforme con q_δ | −11,9320 |
| Candidato completo | −11,8790 |

El pico nace de incorporar el cierre uniforme de corriente mediante q_δ dentro de un núcleo, donde la respuesta espectral es espacial. A la vez, quitar completamente esa incorporación tampoco resuelve el problema: el control BCS+K₀ tiene 83,68 % de discrepancia L² respecto de Usadel. La suma de normas de estas contribuciones no es el error total, porque hay términos cruzados; el JSON conserva esa distinción.

## Acción recomendada

Se debe mantener el candidato actual como referencia de comparación y abrir la revisión de la energía de núcleo antes de promover transientes físicos con ceros del condensado. No corresponde ajustar δ para aproximar este único perfil, retocar KWT para ocultar la fuerza ni retornar automáticamente a BCS+K₀.

El siguiente cambio de implementación concreto es incorporar una **respuesta espacial térmica de referencia** en la interfaz experimental de energía, fuerza y corriente, conservando el mismo régimen estadístico en las comparaciones. La solución radial ya calculada ofrece un primer caso independiente, económico y reproducible. Cualquier nueva energía efectiva de gradiente debe contrastarse simultáneamente con ese caso finito y con la respuesta lineal espacial C.29–C.31; sus fuerzas y corrientes deben derivarse de esa misma energía. Modificar sólo la fuerza volvería a perder la coherencia que motivó esta actualización.

No se ha elegido ni implementado aquí una nueva formulación no térmica de Usadel. La compatibilidad de una alternativa con las poblaciones cinéticas de B y su coste para el dispositivo sigue siendo una decisión posterior del modelo. El resultado actual permite tomar esa decisión con evidencia: la etapa 4 continúa abierta por fidelidad física del cierre, no por exigir más cifras de precisión a los cálculos existentes.

Procedencia reproducible: `radial_reference_analysis.json`, datos originales bajo `raw/stage4_radial_reference_20260923` y `sandbox/stage4_core/analyze_radial_reference.py`. El control de variación de energía tardó 7,45 s y no resolvió nuevos problemas espectrales espaciales.
