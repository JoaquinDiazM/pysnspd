# Casos especiales SSP: equilibrio y poblaciones límite

Los tres casos registrados pasaron sus controles propios, con 630 estados
electrónicos y 1025 fonónicos. Las fuentes físicas, SSP y tolerancias permanecen
iguales. Esta evidencia tiene alcance de casos cortos; no cierra la etapa 2.

| Caso | Intervalo / pasos | Defecto energético máximo | Resultado |
|:--|:--|:--|:--|
| Dos celdas en equilibrio, Gamma 0 y 0,1 | 0,1 / 10 | 0 | PASS |
| Una celda desde vacío fonónico | 0,2 / 40 | 1,91e-9 | PASS |
| Una celda con ocupaciones electrónicas iniciales 0 y 1 | 0,2 / 40 | 2,82e-10 | PASS |

El equilibrio se verifica en toda la trayectoria, no sólo mediante energía:
las amplitudes permanecieron constantes; la mayor variación de las poblaciones,
en norma L1 ponderada, fue 2,93e-18 frente al presupuesto registrado 1,02e-12.
La fuerza inicial máxima fue 7,37e-17. El vacío fonónico produjo fonones y las
poblaciones electrónicas límite se relajaron hacia el interior, conservando
Pauli y positividad. El limitador estuvo inactivo en estos tres escenarios.

La referencia independiente en los campos reales del equilibrio comparó
momentos sobre los 630 estados en 11 tiempos por celda. Para Gamma 0,1 se usó la
parametrización causal real ideal, que permite ver conjuntamente la diferencia
de representación y del regulador: error escalado máximo 5,35e-9. Una segunda
comprobación con 60 dígitos verificó el núcleo a regulador finito en nueve
estados de ese campo, con error máximo 4,37e-16. La identidad derivativa tuvo
defecto 3,61e-11. El muestreo de este segundo control se declara explícitamente;
por sí solo no certifica todos los momentos.

Los campos de ambas trayectorias de frontera pasan la referencia BCS y su cota
de absorción hacia energías superiores no resueltas es 8,89e-10 relativa. El
equilibrio con Gamma 0,1 necesita su control de soporte separado; no se le
atribuye el resultado de Gamma cero ni se supone una población externa arbitraria.

## Evidencia y reproducción

- [special_plan.json](special_plan.json): registro anterior a los resultados.
- [special_summary.json](special_summary.json) y [CSV](special_summary.csv):
  resumen con enlaces y hashes de los resultados.
- [special_raw/special_assessment.json](special_raw/special_assessment.json):
  controles de las tres trayectorias y registro original completo.
- [special_postprocessed/special_postprocess.json](special_postprocessed/special_postprocess.json):
  referencias de campos y soporte de las trayectorias de frontera.

El bloque inicial tardó 104,32 s y se detuvo después de completar las tres
trayectorias, al importar una dependencia de la referencia. Se conservó ese
fallo. `mpmath` ya estaba disponible en el directorio de dependencias de revisión;
se terminó exclusivamente el postprocesamiento con esa ruta, en 8,56 s. No se
repitió ninguna integración. Ambas ejecuciones tuvieron límite de 240 s y un hilo.

Los scripts nuevos están en `sandbox/stage2_cells/closure_prep_20260922/`:
`special_cases.py` registra y ejecuta el bloque; `special_postprocess.py` sólo
lee sus trayectorias y calcula referencias. La primera versión de alcance,
archivada como `special_*_initial_unexecuted`, nunca se ejecutó: el control de
soporte de Gamma finito se asignó a otra revisión antes de correr el bloque.

Para reproducir sólo el postprocesamiento en una carpeta **nueva**, se utiliza
`PYTHONPATH=/home/jdiaz/pysnspd/tmp/stage1_r2_review/deps`, el intérprete `snspd`,
`special_postprocess.py`, el plan registrado y
`--input-root tmp/stage2_special_20260922`. No hace falta relanzar las tres
trayectorias. Estos casos tampoco acreditan la convergencia de malla ni la
precisión temporal de un régimen con limitador activo.
