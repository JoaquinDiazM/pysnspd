# Etapa 4: evolución térmica y siguiente transiente no lineal

La [revisión actual](time_review_20260924/README.md) acepta la trayectoria afín
térmica ejecutada por el usuario: 1 ps físico en 71,10 s de cómputo. El máximo
desplazamiento es 0,222 % del gap; la comparación temporal registrada pasa.
Se resta la deriva de referencia al evaluar las perturbaciones.

El contraste no lineal de los campos guardados también pasa: 2048 raíces en
22,74 s. Un fallo de Newton se resolvió evaluando de forma estable el mismo
incremento energético, sin alterar tolerancia, ecuación ni contactos.
El [informe](../../../output/pdf/implementation/Informe_etapa_4_evolucion_termica.pdf)
muestra evolución, mapas, disipación y diferencias no lineales con su escala.

Sigue la [trayectoria térmica no lineal](time_review_20260924/nonlinear_time/README.md):
ETD2 reutiliza el operador rígido, mientras cada fuerza usa espectros actuales.
Se retienen todos los nodos y las derivadas KWT/potencial normal. Las funciones
phi conservan algebraicamente el término constante sin un estado auxiliar.

También se preparó el [puente longitudinal](time_review_20260924/longitudinal_coupling/README.md)
para unir amplitud y poblaciones. Sus pruebas comprueban reciprocidad y balance
de disponibilidad en fase constante y corriente nula. El trabajo espectral
general y el balance de energía interna siguen pendientes. La etapa 4 no se cierra.

El [cuaderno vigente](../../GEMINGA_COMMANDS.md) identifica el único cálculo
largo solicitado, salidas y recursos. Se limita el presupuesto compartido al
90 %; en Geminga quedan libres dos núcleos físicos completos. No se pide
repetir los cálculos completados ni crear sesiones de terminal.

## Evidencia anterior conservada

- [Momentos de carga, respuesta armónica y operador térmico](moment_review_20260924/README.md).
- [Núcleos autoconsistentes y espectros](self_consistent_review_20260924/README.md).
- [Energía espacial](spatial_energy_20260924/README.md).
- [Diagnóstico del cierre local](followup_20260923/README.md).
- [Controles iniciales](review_20260923/Informe_avance_etapa_4A.md).

Producción y v1.0.0 siguen intactos. Se conserva el circuito completo de la
memoria. Para el futuro fotón, sólo cambiará la ventana hasta Vout más margen.
La transferencia fotónica y las tasas materiales permanecen abiertas; la
etapa 5 no ha comenzado.
