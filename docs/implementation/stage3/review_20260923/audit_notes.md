# Auditoría de resultados guardados de 3A

Se verificaron los 18 resultados y sus archivos de estado, recibos y siete fuentes. La campaña completó los casos, pero no cumplió la precisión espacial registrada. Esta auditoría no evaluó espectros ni fuerzas nuevos, y no integró trayectorias.

Integridad: PASS. Archivos inventariados: 76. Enlaces comprobados: 336. Duración registrada del lote: 15.149 min.

## Lo comprobado

- Los 18 casos aprobaron las identidades discretas. Diferencias máximas de fuerza y corriente: 1.375e-11 y 1.125e-11; invariancia/covariancia de fase: 1.556e-15.
- El menor valor propio es 1.57079633; la incertidumbre de diferenciación máxima es 2.774e-05.
- Los 10 puntos representativos 630/1260 conservan signo positivo. Desplazamiento máximo de la matriz: 4.202e-09; diferencia máxima escalada de momentos: 5.639e-09. Esto describe los estados y resoluciones muestreados, no una cota global del espectro continuo.
- El control negativo D.36 permanece rechazado: valor propio -0.30661962. No se interpreta su rechazo esperado como fallo del código.
- El error del gradiente analítico baja de 5.0359% a 1.2785% y 0.3209%. Pasa el control de 1% en 32 celdas; la corriente nula de este caso pasa sólo su control absoluto.

## Precisión espacial que falta

Los siguientes porcentajes son estimaciones condicionales de extrapolación, no errores frente a una solución exacta. El presupuesto original sigue siendo 1%.

| Caso | Energía de exceso | Fuerza inducida | Corriente inducida |
| --- | ---: | ---: | ---: |
| weak_amplitude_thermal | 1.8592% | 1.9968% | 69.7780% |
| weak_phase_thermal | 5.8938% | 8.6638% | 6.0911% |
| weak_phase_nonthermal | 5.8707% | 8.7560% | 6.0679% |
| weak_amplitude_vacuum | 1.0097% | 1.0101% | Cero analítico; control absoluto |

Son 11 observables fuera del objetivo. Las diferencias decrecen al refinar, pero eso no basta para afirmar precisión del 1%. El caso de amplitud en vacío queda apenas por encima; su fallo se conserva.

## Por qué una identidad exacta puede coexistir con ese error

Para una hélice nodal de amplitud a y flujo q, el promedio cartesiano consulta `a_link = a cos(qh/2)` y `|d_link| = 2a sin(qh/2)/h`. La amplitud del punto medio depende así del flujo y del tamaño de celda. La fuerza local no es cero en el estado prescrito a=0.9: no es un equilibrio autoconsistente. Por ello, derivar correctamente esa energía discreta incluye un sesgo de corriente de orden h².

Con `m=a²/(a²+delta²)`, `c=m/4-1/6`, `b=-1/12-m²(2c-1/4)` y `Gamma=(mq)²/r`, la expansión es:

`e_h=e+h²[-a*u_a*q²/8+2*u_Gamma*Gamma*q²*c+kappa*a²*q⁴*b]+O(h⁴)`.

Su derivada a ocupaciones y amplitud fijas da `J_h=J+h²*C_J+O(h⁴)`, con:

`C_J=-a*q*(u_a+u_aGamma*Gamma)/4+4*c*Gamma*q*(u_GammaGamma*Gamma+2*u_Gamma)+4*kappa*a²*q³*b`.

Se usa aquí sólo la expresión analítica ya documentada del vacío con gap. Los tres aportes y las medidas originales quedan separados en el JSON de auditoría.

En el uniforme de vacío a 32 celdas, la corriente discreta difiere del límite continuo en 5.089%. Su desplazamiento medido es 0.00809645; el término h² predice 0.00795307, es decir 98.23% del desplazamiento. Este acuerdo identifica un sesgo de reconstrucción espacial, aunque no demuestra que todos los errores no uniformes provengan exclusivamente de ese término.

El error en respuestas inducidas puede ser mucho mayor que en la corriente total, porque se resta una referencia uniforme y la señal restante es pequeña. No corresponde ocultarlo normalizando con la energía total del fondo ni reinterpretarlo como fallo de la corriente física.

## Qué conservar y qué volver a comprobar

- Conservar los 18 resultados y sus fallos espaciales; no repetir identidades sin cambios sólo para recrear evidencia.
- Esta evidencia no exige repetir trayectorias de etapa2 ni aumentar globalmente la resolución espectral.
- Una reconstrucción espacial nueva necesita comprobaciones variacionales y de respuesta espacial propias con los mismos perfiles físicos, además de D.36 en sus estados muestreados.
- Si sólo se añaden mallas espaciales más finas con el mismo operador, reutilizar esta evidencia bajo sus fuentes y entradas exactas y registrar previamente la ampliación.
- No atribuir precisión al siguiente ensayo espacial usando solamente conservación de energía o error relativo sobre una gran energía de fondo.

La evidencia sostiene el desarrollo del funcional, sus signos y la estabilidad local ensayada. Falta resolver o delimitar la precisión de la reconstrucción espacial antes de atribuir exactitud a la respuesta del siguiente ensayo. Bordes, reservorios, potencial y circuito conservan sus verificaciones propias; la etapa 3 completa sigue abierta.

## Reproducción y fuentes

`python sandbox/stage3_spatial/review_20260923/audit_saved.py --raw <directorio_del_lote> --verify-only`

El programa lee exclusivamente JSON/NPZ, verifica hashes y reproduce comparaciones algebraicas. Para generar una copia de la auditoría, indicar un directorio nuevo mediante `--output-root`. El inventario por archivo y los hashes de las fuentes están en `saved_results_audit.json`.

- Modelo continuo: `docs/modelo_v0_4/D_sintesis_plan_y_verificaciones_v0_4.md`, D.8–D.10 y D.36.
- Vacío analítico: `pysnspd/experimental/energy_catalog.py`, `vacuum_state`.
- Reconstrucción discreta: `pysnspd/experimental/spatial_functional.py`.
- Registro previo: `docs/implementation/stage3/iteration_20260923/registration.json`.
- Resultados originales: manifiesto, 18 casos y `summary.json` del directorio del lote.
