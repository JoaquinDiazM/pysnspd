# Diagnóstico de la monotonía del catálogo local

Los dos interpolantes archivados quedan **sin admisión**. Se evaluaron únicamente sus coeficientes guardados: no se consultó la fuente física, no se cambiaron coeficientes y no se ejecutó dinámica.

El primer fallo del recorrido es el grado4, probe2, a=0.9904, Gamma=0.00818 (desplazamiento +2e-4 del stencil). La pareja de índices 19→20 (base cero) tiene E=0.9303258878771008 y 0.9303243512854876: incremento **-1.5365916e-06**. Los conteos correspondientes son 5.9540013e-08 y 7.2082764e-08.

El grado6 también falla: sus centros de probe rechazados son [2, 3]. En particular ya pierde monotonía en el centro del probe2, antes de su stencil. Omitir el grado4 no arreglaría el candidato.

| Grado | Probes + stencils no monótonos | Muestras de caja no monótonas | Peor incremento muestreado |
|---:|---:|---:|---:|
| 4 | 22/63 | 2693/6561 | -4.7851689e-05 |
| 6 | 21/63 | 2827/6561 | -2.0236947e-05 |

La malla de diagnóstico81×81 sólo mide fallos muestreados; no certifica toda la caja. Las energías permanecen positivas en este muestreo. La pérdida de orden se concentra cerca de la transición de borde a conteo muy pequeño. Los saltos negativos exceden por muchos órdenes de magnitud el espaciado de float64: no son redondeo que deba recortarse.

La reconstrucción de los nodos de entrenamiento conserva energías ordenadas. El problema aparece entre nodos: ajustar cada energía con un polinomio independiente no preserva que las curvas vecinas estén ordenadas. Aumentar el grado reduce algunos errores, pero no aporta esa garantía.

Una representación adecuada para la siguiente iteración es interpolar los logaritmos de los incrementos positivos: `d0=E0`, `dj=Ej−E(j−1)`, `Lj=log(dj)`, y reconstruir `Ej=sum(k≤j) exp(Lk)`. Las derivadas se obtienen de esa misma energía: `∂tEj=sum exp(Lk)∂tLk`, con t=a o Gamma. Esto conserva positividad y orden por construcción sin clipping, ordenamiento posterior ni cambio de poblaciones.

La propuesta aún debe registrarse y contrastarse con la fuente física usando el mismo presupuesto de precisión. Sus derivadas pueden seguir necesitando una caja menor o una representación por tramos cerca del borde. Mantener la monotonía no basta para admitir momentos, curvaturas ni dinámica. También deben rechazarse desbordamientos o incrementos que pierdan resolución numérica al reconstruir.

Trazabilidad completa, índices, pesos, todos los puntos fallidos, hashes y fórmulas: `patch_monotonicity_diagnosis.json`. Reproducción sin sobrescribir: `python sandbox/stage3_spatial/coupled_20260923/diagnose_patch_monotonicity.py --verify-only`.
