# Estado vigente de la etapa 3

La [iteración de bordes y circuito](ports_20260923/README.md) contiene el estado
actual, el informe y el siguiente cálculo. La campaña nodal terminó sus 18 casos:
3A queda cerrado para desarrollo en los ensayos registrados. Diez estimaciones
relativas pasan; una diferencia de energía conserva su dictamen sin certificado
relativo. No se reclasificaron fallos históricos ni se relajaron tolerancias.

Ya están implementados y probados el funcional longitudinal abierto, la rama
superconductora del reservorio, su intercambio electrónico, el potencial y el
circuito de tres estados de la memoria. Pasaron 168 pruebas y 23 subpruebas. El
piloto físico abierto pasó en 20,9 s. Las seis comparaciones abiertas de longitud
y malla están preparadas para ejecución manual en Geminga, con progreso y ETA.

La etapa 3 sigue abierta: faltan esa comparación, el empalme 2D-1D y el ensayo
débil acoplado en el tiempo. El control circuital usa una resistencia prescrita,
no una señal calculada por el material. Producción y v1.0.0 no cambian.

La [iteración nodal](nodal_20260923/README.md) y la
[primera entrega 3A](iteration_20260923/README.md) son históricas. Sus manifiestos
se verifican con las instantáneas exactas de las libretas y entradas anteriores;
el verificador vigente es `sandbox/stage3_spatial/ports_20260923/verify_delivery.py`.

`README.md` y `entry_contract.json` de esta carpeta conservan la propuesta
del 22 de septiembre, incluida en la entrega inmutable del cierre de etapa 2.
La secuencia física propuesta sigue vigente; su estado «no iniciada» es histórico.
