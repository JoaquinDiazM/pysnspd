# Etapa 4: preparada la unión dinámica polarizada

La [entrega vigente](coupled_closure_20260924/README.md) contiene el nuevo
comando: dos referencias estacionarias con corriente y cinco respuestas
acopladas de espectro, poblaciones, fase, potencial y circuito completo.
Se ejecutaron los pilotos ligeros; el lote largo queda para ejecución manual.
Mantiene el margen práctico del 2 % y el presupuesto global de 28 hilos.
La etapa todavía no se declara cerrada.

## Resultados anteriores conservados

El [ensayo actual](final_kwt_20260924/README.md) completó **1 ps en 153,81 s**
con el paso KWT real de la memoria, la malla dual de 1712 nodos y 256 frecuencias.
Las dos sondas y sus dos resoluciones temporales cumplen el presupuesto del
2 %. El máximo defecto integrado es 1,147 % para la sonda de fase principal;
se reduce a 0,577 % al dividir el paso. No hubo crecimiento de energía.

La reutilización del Jacobiano espectral uniforme evita reconstrucciones:
893440 predicciones cumplieron el residuo no lineal completo, sin necesitar
Newton adicional en esta campaña. No se sustituyó la fuerza ni se congeló
el espectro. Euler sigue siendo de primer orden temporal aunque su amplitud
se obtenga de una ecuación cuadrática. El cálculo mantuvo dos núcleos físicos
libres y un hilo numérico por proceso.

El [control longitudinal](final_kwt_20260924/longitudinal/README.md) también
terminó: intercambio débil no térmico y disponibilidad integrada en un sector
invariante exacto sobre la misma malla. No es una evolución no lineal general
de poblaciones, calor, carga y circuito.

La [revisión de tolerancias](final_kwt_20260924/acceptance_policy.json) acepta
el ensayo anterior al 2 % sin recalcularlo. Su desacuerdo máximo del 1,500017 %
queda cerrado bajo ese criterio; se conserva íntegro el certificado original.

Antes de preparar el lote actual, a petición del usuario se investigó la
[separación de escalas](final_kwt_20260924/quasiclassical_assessment.md). Usadel
difusivo sigue siendo el marco efectivo apropiado; no implica que las
poblaciones ni el desequilibrio electrón–hueco sean instantáneos. La ruta
recomendada conserva su evolución y trata el espectro instantáneo como una
aproximación adiabática separada. Falta completar la unión de trabajo
espectral, fase, potencial y calor antes de cerrar toda la etapa 4.

El [informe de resultados](../../../output/pdf/implementation/Informe_etapa_4_Euler_y_acoplamiento.pdf)
identifica cada campo, norma, referencia y normalización. El
[cuaderno](../../GEMINGA_COMMANDS.md) conserva los comandos ya completados y
las reproducciones opcionales. Producción, v1.0.0 y el circuito de tres estados
permanecen intactos; etapa 5 y la predicción fotónica todavía no se inician.

## Evidencia anterior conservada

- [Diagnóstico de la interrupción y adaptadores](practical_time_review_20260924/README.md).
- [Trayectoria afín, corrección de Newton y ETD2](time_review_20260924/README.md).
- [Momentos de carga y respuesta armónica](moment_review_20260924/README.md).
- [Núcleos autoconsistentes](self_consistent_review_20260924/README.md).
- [Energía espacial](spatial_energy_20260924/README.md).
- [Diagnóstico del cierre local](followup_20260923/README.md).
