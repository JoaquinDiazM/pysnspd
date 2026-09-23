# Canales aislados y orden suave del escape

Resultado: **los canales de reacción y el escape tienen evidencia aprobada en
su alcance; el transporte aislado revela un defecto numérico reproducible y
permanece rechazado**. No se modificaron los kernels, el integrador congelado ni
los criterios de aceptación. Todos los experimentos se prerregistraron y se
ejecutaron con un hilo y timeout total de 240 s por bloque autorizado.

Este dictamen corresponde al SSP anterior, cuyo fallo se conserva. El prototipo
posterior y sus verificaciones acotadas se documentan por separado en
`isolated_guarded_findings.md`; todavía requieren validación global propia.

Se usó el catálogo seleccionado de 630 estados electrónicos, 1025 nodos
fonónicos, IR=0.005 y el mapa SSP original. Los campos se mantuvieron fijos
explícitamente. Para separar emisión, absorción, recombinación y creación, el
adaptador selecciona la actividad forward o reverse original y la familia de
eventos correspondiente, conservando coeficientes, geometría y estequiometría.
No corta el signo de una tasa neta. Estos suboperadores dirigidos no representan
el equilibrio del operador reversible completo.

Las cuatro reacciones usan amplitud 0.5. La banda 1.1<x<1.5 comienza llena,
con p=1, y su exterior vacío, p=0; creación comienza con todos los electrones
vacíos. Emisión y recombinación comienzan con n=0; absorción y creación, con
n=0.05. Cada trayectoria tiene duración 0.02 y 21 tiempos guardados. BGK,
movimiento del condensado, calor y escape están deshabilitados explícitamente.
No se afirma haber ensayado un estado electrónico globalmente lleno dentro de
la inversión de temperatura positiva.

| Caso | Resultado | Error de energía | Error de conteo conservado |
|---|---|---:|---:|
| Emisión | PASS de los estados guardados, recuperados sin reintegración | 1.380e-15 | 1.249e-15 |
| Absorción | PASS | 1.381e-15 | 1.249e-15 |
| Recombinación | PASS | 1.380e-15 | 1.249e-15 |
| Creación | PASS | 1.480e-15 | 1.480e-15 |
| Transporte | Rechazo de soporte reproducido | No se admite la trayectoria | No se admite la trayectoria |

Para scattering se conserva el número electrónico; para recombinación/creación
se comprueba Q+2N_ph. Todos los estados guardados de las reacciones permanecen
en sus dominios físicos y presentan evolución en la dirección esperada. La
absorción activa el limitador en 133694 evaluaciones; su defecto integrado
ponderado por energía es 1.8250e-6. Esto demuestra soporte y balance, **no
precisión temporal** de ese caso. Recombinación y creación no activaron el
limitador en estos ensayos.

El primer intento completó emisión, pero falló al serializar un `np.bool_` en el
JSON. Se conservan intactos `isolated_emission.npz`, el JSON parcial, el log y
`isolated_results.json` con estado de fallo. Un segundo script convierte
exclusivamente escalares NumPy a sus equivalentes Python y escribe en otro
directorio. Recupera emisión desde NPZ sin integración ni RHS nuevo. Los
contadores del limitador y los valores numéricos instantáneos que no llegaron a
escribirse se declaran perdidos; no se reconstruyen ficticiamente. Los otros
tres casos aprobados fueron ejecutados una sola vez.

## Defecto localizado en transporte

El caso de transporte usa amplitudes 0.5 y 1, una banda llena en la celda
izquierda, electrones vacíos en la derecha y fonones vacíos en ambas. Todas las
reacciones están deshabilitadas. El primer ensayo rechazó una etapa. El
diagnóstico posterior, autorizado y prerregistrado, añadió exclusivamente un
observador que guarda la entrada y el estado rechazado; no modificó operaciones,
flujos, resultados ni excepciones del mapa.

El rechazo ocurre en el paso macro 17, etapa SSP 2, llamada FE 50, a tiempo de
etapa 0.017 con h=0.001. Sólo se detecta un electrón negativo: celda 0, índice
297, x=1.6757735026918965. No hay fonones negativos.

Sea q=4.9406564584124654e-324, el menor float64 positivo:

| Magnitud | Valor observado |
|---|---:|
| Ocupación de entrada | 70q = 3.46e-322 |
| Capacidad | 0.040000000000000036 |
| Producto exacto capacidad por ocupación | aproximadamente 2.8q |
| Inventario almacenado en float64 | 3q = 1.5e-323 |
| Cociente de disponibilidad tras margen 0.9 | 2.158998389522019e-10 |
| Incremento calculado de ocupación | -75q = -3.7e-322 |
| Ocupación rechazada | -5q = -2.5e-323 |

El margen del presupuesto se pierde al redondear productos y actualizaciones
en el régimen subnormal. El inventario redondeado sobreestima el producto exacto;
la actualización `h*count_change/capacity` agota más que la ocupación disponible.
La suma decimal exacta de los operandos float y el estado rechazado se conservan
en el diagnóstico. Los avisos de división por cero/desbordamiento se observaron
en otros cocientes; no se usan por sí solos como explicación de este nodo,
cuyo cociente registrado es finito.

Este resultado identifica aritmética de presupuestos y actualización, no una
escala física rápida que se haya medido. Su peso energético diminuto no autoriza
aceptar la población negativa. La corrección debe preservar los flujos comunes
y el soporte **antes** de aceptar el paso; no se aplicó clipping, reparación
posterior de energía, otro método ni tolerancia más laxa. El caso queda abierto
para una versión numérica corregida y su correspondiente verificación.

## Escape: orden tres medido

El ensayo independiente de escape, que no se había alcanzado antes del fallo
de transporte, se ejecutó con el mismo mapa SSP, campo fijo, tau=0.7 y duración
0.2. Se compara contra la solución exacta
`n_bath+(n_initial-n_bath)*exp(-t/tau)`; no hay reacciones, transporte ni BGK.

| Pasos | Error L1 ponderado por capacidades | Reducción |
|---:|---:|---:|
| 10 | 9.149185718e-9 | — |
| 20 | 1.130643268e-9 | 8.092 |
| 40 | 1.405248725e-10 | 8.046 |

Las reducciones superan 6, el mínimo preregistrado correspondiente a la razón
esperada 8 con margen del 25%. Las tres trayectorias mantienen poblaciones
físicas; el error energético máximo es 1.943e-16. El orden suave del escape queda
comprobado. No concede precisión al transporte rechazado ni cierra la etapa.

## Evidencia

- `isolated_registration.json`: primer prerregistro y fuentes congeladas.
- `isolated_recovery/isolated_registration.json`: recuperación autorizada de
  serialización, hashes del intento inicial y de sus archivos parciales.
- `isolated_recovery/isolated_emission_recovered.json`: invariantes recuperadas
  sin reintegrar, con estadísticas perdidas declaradas.
- `isolated_recovery/isolated_absorption.json`, `isolated_recombination.json` e
  `isolated_creation.json`: trayectorias, identidades y actividad del limitador.
- `isolated_transport_diagnostic/isolated_transport_diagnosis.json` y
  `isolated_transport_rejected_state.npz`: etapa, inventarios, factores, estados
  y magnitud del rechazo reproducido.
- `isolated_transport_diagnostic/isolated_escape_order.json`: errores, razones
  y hashes del control de orden, junto con sus tres JSON/NPZ.

Los tiempos medidos de los bloques numéricos fueron 3.53 s hasta el error de
serialización, 7.32 s para la recuperación y los canales pendientes hasta el
rechazo, y aproximadamente 0.21 s para el diagnóstico observacional más escape,
sin incluir importación y transporte de archivos. No hubo cálculo largo.
