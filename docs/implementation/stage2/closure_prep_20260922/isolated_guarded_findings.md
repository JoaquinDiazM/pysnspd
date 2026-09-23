# Protección numérica de inventarios subnormales

El prototipo nuevo supera las comprobaciones acotadas. **No constituye todavía
admisión temporal o de malla del modelo completo.** Los kernels físicos, el SSP
anterior y sus resultados permanecen intactos. La nueva fuente es
`limited_ssp_guarded.py`, SHA256
`9b77b91ff7b4e68786ea07a7e915f119d330f1c957ca22d2dd7b07028b00fad4`.

La rama ordinaria llama exactamente al mapa antiguo. El integrador reutiliza su
mismo bytecode con un diccionario global privado que sustituye únicamente el
despachador de pasos; no modifica el módulo original. Cuando un inventario
positivo de electrones, huecos o fonones cae por debajo de
`1024*float64.tiny`, la rama protegida calcula inventarios, demandas, factores
comunes, acumulación y suma de poblaciones en `longdouble`. La misma extensión
de evento actualiza todos sus participantes y sus libros de energía.

Se verifica el soporte antes y después de convertir el resultado a float64.
No hay clipping, reparación energética ni asignación forzada de ceros. El
redondeo ordinario de cantidades positivas por debajo del mínimo representable
puede producir cero al convertir el dtype; el prototipo cuenta esos casos.
La rama usa `np.add.at`, porque `np.bincount` no acepta las sumas extendidas.

En Geminga, `longdouble` dispone de almacenamiento de 128 bits, mantisa de 63
bits explícitos y exponente mínimo -16382. En el entorno Windows observado es
equivalente a float64 y el módulo rechaza explícitamente su ejecución. No se
afirma portabilidad de esta implementación a ese entorno.

## Resultados medidos

La QA principal duró **6.655 s**, con 14 comprobaciones aprobadas y un único
límite de 240 s. El caso del error original se recupera sin aceptar un negativo:

| Comprobación | Resultado |
|---|---|
| Celda seleccionada 630/1025, condensado móvil, duración 0.01, cinco pasos | Todos los bytes de tiempos y estados coinciden con el SSP anterior; 15/15 etapas delegadas |
| Etapa que antes produjo -2.5e-323 | Produce +3.5e-323; rama extendida activa; sin clipping |
| Transporte fijo completo, duración 0.02 y 20 pasos | PASS; 27 etapas extendidas y 33 ordinarias |
| Balance energético de ese transporte | Error máximo 2.017e-15 |
| Conteo electrónico de ese transporte | Error máximo 1.249e-15 |
| Reacciones y transporte competidores en un sistema de cuatro estados | Soporte físico y energía conservada, limitador común activo |
| Escape, 10/20/40 pasos | Estados bit a bit iguales a los tres archivos anteriores; reducciones 8.092 y 8.046 |

El transporte completo registra 77 conversiones de cantidades positivas a cero
por falta de representación en float64. El defecto energético integrado de
los flujos limitados durante sus etapas extendidas es 2.3744e-6. Este último
valor mide actividad del limitador, **no un residuo de conservación**, y no
certifica por sí mismo precisión temporal. La convergencia debe medirse con el
prototipo cuando se quiera admitir esa trayectoria como resultado dinámico.

Una revisión independiente solicitó cubrir también tasas negativas, que consumen
fonones. Se añadieron dos microtrayectorias de absorción y creación, cada una
con 20 pasos y 21 tiempos. El fonón de energía 0.7 comienza en 1e-320 y participa
realmente en las tasas y en su consumo. Ambas activan 25 etapas extendidas de
60, mantienen poblaciones físicas y conservan energía y la estequiometría
correspondiente:

| Caso | Cambio del conteo fonónico | Cambio del conteo electrónico | Error energético | Error del conteo conservado |
|---|---:|---:|---:|---:|
| Absorción | -0.0010033771 | -3.886e-16 | 2.776e-16 | 3.886e-16 |
| Creación | -0.0339093791 | +0.0678187582 | 1.943e-16 | 1.665e-16 |

El suplemento duró **0.059 s**. Su alcance es un sistema espectral normal de
cuatro estados y seis fonones con la misma ley proyectada y el mismo mapa.
Una primera preparación fue rechazada antes de integrar porque un único nodo
electrónico ocupado, rodeado de ceros, anula las actividades geométricas de los
intervalos. Ese intento de 0.019 s y su prerregistro se conservan. El fixture
corregido tiene poblaciones electrónicas interiores para absorción; creación
conserva su condición inicial de electrones vacíos. No se cambió el prototipo,
ningún kernel ni ningún umbral para obtener estos resultados.

## Evidencia y siguiente límite

- `isolated_guarded/isolated_results.json`: QA principal, capacidades de
  precisión, activación y defectos del limitador.
- `isolated_guarded_negative/isolated_results.json`: precheck rechazado antes
  de integración; conservado como historial.
- `isolated_guarded_negative_v2/isolated_results.json`: suplemento de absorción
  y creación con consumo de fonones bajo la rama extendida.
- `isolated_guarded_inventory.json`: inventario externo de **24 archivos**,
  incluidos `isolated_normal_equivalence.npz` e
  `isolated_stage_correction.npz`, sin reescribir los resultados congelados.
  Registra Python 3.10.20, NumPy 2.2.6 y SciPy 1.15.3, observados en la misma
  instalación de Geminga después de la QA. Los 24 hashes se verificaron también
  tras copiar los archivos al repositorio local.

La igualdad bit a bit se midió para el intervalo normal corto y las tres
trayectorias de escape indicadas; no se extrapola a las trayectorias completas
ONE/TWO. El módulo protegido requiere su propia validación temporal y de mallas,
salvo para cualquier trayectoria concreta cuya equivalencia se pruebe de manera
completa y reproducible. No se ejecutó ningún cálculo largo ni se modificó el
plan manual de mallas desde este subtask.
