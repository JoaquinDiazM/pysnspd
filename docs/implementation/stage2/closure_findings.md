# Cierres locales: resultado y alcance

Los cierres KWT y la entrada Debye sintética pasan sus controles algebraicos y
de cuadratura. Las 61 regresiones del módulo terminaron en 0,12 s; el diagnóstico
de Geminga tardó 0,477 s. La evidencia y los hashes están en
[closure_results.json](closure_results.json).

La velocidad y el calor calculados en unidades normalizadas coinciden con una
evaluación independiente en SI: los errores relativos escalados máximos son
3,69·10⁻¹⁵ y 7,20·10⁻¹⁵. La movilidad es positiva, el límite de amplitud cero es
finito y el calor coincide con el trabajo perdido por el condensado. Se conservó
el factor 1/2 entre la fuerza respecto de la amplitud real y la fuerza compleja.
El tiempo de relajación cinética es un parámetro distinto y no interviene en esta
movilidad.

La entrada Debye es una prueba sintética declarada: ΩD/Δ₀ = 4, λ = 0,1 y
nátomo = 10 N₀Δ₀. No es una calibración de NbN. Su normalización de modos y su
acoplamiento integrado tienen errores menores que 4,74·10⁻¹⁶; con 129 nodos, la
energía y la capacidad térmica coinciden con la cuadratura adaptativa independiente
a 1,05·10⁻¹⁵ relativo en las temperaturas ensayadas.

El corte infrarrojo cambia explícitamente la entrada, sin renormalizarla. Su
efecto térmico depende de la temperatura:

| Ωmín/Δ₀ | Energía Debye omitida a kBT/Δ₀ = 0,03 | A kBT/Δ₀ = 0,12 |
|---:|---:|---:|
| 0,020 | 1,1743 % | 0,022312 % |
| 0,010 | 0,16740 % | 0,0028787 % |
| 0,005 | 0,022312 % | 0,00036554 % |

Para un dominio que incluya la temperatura 0,03, el corte 0,005 satisface el
presupuesto de 0,1 % de energía térmica en esta prueba. Este control no sustituye
la convergencia de los procesos electrón–fonón ni autoriza un corte para cualquier
distribución no térmica.

El lector de la forma fonónica de NbN verifica la trazabilidad de la tabla
derivada y exige factores independientes para convertir el eje y la DOS.
Estos factores se declaran sintéticos: el lector no resuelve la procedencia de
las unidades, la base átomo/celda ni la densidad absoluta pendiente. Por ello
su estado sigue siendo condicional y no admite tasas SI de NbN.

La comprobación KWT valida la implementación de la ley heredada. No demuestra
su precisión fenomenológica lejos de Tc ni predice tiempos de respuesta de un
material real.
