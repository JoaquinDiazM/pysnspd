# Diagnóstico del costo de la referencia temporal

La referencia DOP853 completada por ejecución manual conserva la energía, pero su costo no justifica trasladar directamente ese procedimiento a todos los ensayos. El diagnóstico estático favorece como explicación la baja suavidad del operador discretizado al moverse el condensado y el control muy estricto de cada ocupación. No demuestra ausencia de rigidez en todos los modos ni a lo largo de toda la trayectoria.

No se cambió ninguna ecuación, kernel o tolerancia. Se evaluaron 17 estados próximos al estado inicial guardado de la referencia: 630 ocupaciones electrónicas y 1025 fonónicas. No se integró el tiempo. La ejecución tardó **2.65 s**, con **982409 eventos** y aproximadamente **0.144 s por RHS**. Los hashes de fuentes, criterios y condiciones iniciales están en [rhs_smoothness.json](rhs_smoothness.json).

## Qué muestran los números

Las diferencias centrales modifican únicamente la amplitud, manteniendo fijas las poblaciones. La norma de poblaciones es la L1 ponderada por las capacidades del catálogo, como exige el criterio de aceptación.

| Incremento de amplitud h | Segunda diferencia del RHS / norma del RHS | Segunda derivada estimada en norma L1 | Cambio relativo de la pendiente respecto al h anterior |
|---:|---:|---:|---:|
| 1e-3 | 5.87e-4 | 473 | — |
| 1e-4 | 2.52e-5 | 2031 | 15.6% |
| 1e-5 | 6.40e-7 | 5157 | 4.65% |
| 1e-6 | 9.89e-9 | 7970 | 1.35% |
| 1e-7 | 2.60e-11 | 2098 | 0.209% |
| 1e-8 | 7.96e-13 | 6419 | 0.00423% |

La curvatura estimada cambia mucho al reducir la escala. Las mayores segundas diferencias para h entre 1e-3 y 1e-6 aparecen en modos fonónicos altos. A h de 1e-8, la mayor segunda diferencia de una ocupación individual es 1.19e-9 en un nodo electrónico de muy pequeña capacidad; su contribución ponderada es mínima. Esto ilustra por qué una norma sobre cada componente puede gastar trabajo en detalles de peso físico muy pequeño.

Al mantener fija la amplitud y perturbar las poblaciones en la dirección de su evolución inicial, la escala local obtenida es **1.5114 por unidad de tiempo**, estable entre incrementos de 1e-5 y 1e-6. La diagonal fonónica más negativa del Jacobiano de eventos es **−0.4900801**, corroborada por diferencias finitas; la norma ponderada de esa columna, dividida por su capacidad, es **2.0237**. Estas medidas no muestran una escala extremadamente rápida en las direcciones examinadas. No son un cálculo del espectro completo del Jacobiano.

Las diferencias del RHS anteriores no son errores de una trayectoria. No certifican por sí mismas la tolerancia temporal de 1e-4; esa comparación debe hacerse con estados integrados y la referencia disponible.

## Origen en el código

- `coupled_cells.py:132–135` construye el campo electrónico y toda la lista de eventos en cada evaluación. Reutilizar una lista congelada durante el movimiento del condensado cambiaría las energías de reacción y no es una optimización admisible sin una nueva formulación y comprobación.
- `kinetic_events.py:380–414` cambia la partición de integración con las energías instantáneas. Resuelve los tramos electrónicos y paneles uniformes de energía de anchura 0.125, pero no incluye todos los nodos fonónicos como límites de integración.
- `kinetic_events.py:491–493, 519–526` asigna cada evento mediante funciones baricéntricas de dos nodos. Al cruzar una energía nodal, cambia su pendiente. El cambio es continuo dentro del soporte, pero no tiene la suavidad que un método de orden alto aprovecha en un problema analítico.
- `kinetic_events.py:528–546` utiliza potencias geométricas coherentes con esas mismas asignaciones. Conservan las identidades de equilibrio, pero también heredan los cambios de pendiente. En ocupaciones exactamente nulas pueden dejar de ser Lipschitz; esto no explica directamente el ensayo iniciado con poblaciones interiores positivas.
- `kinetic_events.py:548–562` acumula muchas contribuciones y divide por capacidades que pueden ser tan pequeñas como 2.02e-9 para electrones y 1.34e-10 para fonones. Eso aumenta la relevancia de los errores de redondeo si se observa solamente la derivada de una ocupación individual.
- `run_coupled.py:170–172` solicita DOP853 con tolerancia relativa 1e-9 y absoluta 1e-13 iguales para todas las componentes. Su `max_step` es un límite superior, no una garantía de pocas evaluaciones.
- Hay costos menores repetidos: `coupled_cells.py:137` invierte la energía para obtener la temperatura y `cell_validation.py:115` repite esa inversión dentro del calentamiento. Resolverlos no elimina los casi un millón de eventos por evaluación.

## Secuencia numérica recomendada

1. **Usar la referencia ya terminada como control independiente del intervalo corto.** Comparar tres pasos de RK4 con las mismas fuentes, entradas y mallas; revisar poblaciones ponderadas, energías, amplitud, transferencias y todos los puntos temporales comunes. No aceptar solamente el residuo total de energía. Mantener los umbrales congelados.
2. **Decidir a partir de la convergencia medida.** Si RK4 satisface los criterios con costo razonable, conservarlo en esta etapa y registrar el costo de DOP853. Si el error no converge, no reemplazar esa evidencia por una tolerancia más laxa.
3. **Si hace falta mejorar la regularidad**, evaluar una variante de cuadratura que incluya los nodos fonónicos reales entre los límites de integración en energía. El integrando y las asignaciones se integrarían por sus tramos, en vez de permitir que una muestra virtual cruce una esquina de la base fonónica durante un paso temporal. Mantener coeficientes positivos, las mismas actividades y la misma estequiometría permite conservar las identidades discretas. Cambia la cuadratura, por lo que requiere repetir los controles independientes y la convergencia espacial; también puede aumentar el número de eventos. Primero medir un único estado y su costo.
4. **Optimizar sin modificar el operador:** reutilizar constantes y reglas de Gauss, evitar la segunda inversión térmica, reducir copias y fusionar el cálculo de tasas con las sumas por nodo. Verificar equivalencia de RHS, equilibrio y balance antes de volver a integrar. No interpolar o congelar arbitrariamente las tasas entre amplitudes.
5. **No adoptar un integrador implícito solo por el número de llamadas.** Un Jacobiano numérico denso de unas 1660 variables exige miles de RHS y no elimina esquinas del operador. Un método implícito queda justificado si una medición posterior identifica modos rígidos; un control temporal basado en las normas físicas necesita validación contra la referencia, límites de positividad y los mismos umbrales de aceptación.

La decisión inmediata es medir convergencia temporal contra la referencia existente. La mejora posterior de cuadratura y costo se apoya en un defecto numérico identificable, sin reinterpretar el modelo físico o el dispositivo.
