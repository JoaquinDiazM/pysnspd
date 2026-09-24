# Etapa 4: evolución térmica y corrección del Newton espectral

La trayectoria térmica afín solicitada terminó correctamente en **71,10 s**.
Los 37 checkpoints se recuperaron y verificaron. El desplazamiento máximo fue
**0,222 %** del gap, dentro del vecindario débil declarado del 2 %. Restando
la deriva de referencia, a 1 ps permanece el **73,59 %** de la perturbación
de amplitud y el **1,751 %** de la perturbación angular (normas del campo completo).
La comparación temporal registrada pasa; no se pide repetirla.

El contraste posterior evalúa las ecuaciones no lineales en los campos guardados
a 0, 0,01 y 1 ps. Un primer lote falló porque la resta de energías casi iguales
rechazaba un descenso verdadero de Newton. El nuevo backend evalúa la misma
diferencia mediante una identidad algebraica estable. Conserva residual,
Jacobiano, acción, condiciones de borde y tolerancia espectral 10⁻⁷.
Las fuentes y el fallo anteriores se preservan.

Con esa corrección, **2048 raíces espectrales terminaron en 22,74 s** y pasan
todos los márgenes del 1 % relativos a la respuesta inicial/remanente. A 1 ps
la discrepancia de velocidad de amplitud es **0,406 % del impulso inicial**;
el residuo instantáneo del balance no lineal de potencia no supera 1,8×10⁻¹²
en unidades internas. Son evaluaciones constitutivas en campos guardados,
no una trayectoria no lineal integrada ni un balance de energía interna.

- [Informe ilustrado](https://github.com/JoaquinDiazM/pysnspd/blob/c80c0f8612d8fd5a2b390938ca5ed4b9d58d3caa/output/pdf/implementation/Informe_etapa_4_evolucion_termica.pdf).
- [Análisis reproducido de la trayectoria](analysis.md) y [datos](analysis.json).
- [Newton estable y contraste no lineal](nonlinear_snapshots/README.md).
- [Siguiente trayectoria térmica no lineal](nonlinear_time/README.md).
- [Prototipo de intercambio longitudinal con poblaciones](longitudinal_coupling/README.md).
- [Decisión de alcance](decision.json).

La cola del torque casi extinguida no tiene precisión relativa demostrada.
La energía cuadrática de la trayectoria afín y su defecto de Taylor se distinguen
de la energía no lineal evaluada después. La pequeña deriva de la coordenada
auxiliar constante del Arnoldi anterior queda registrada; el nuevo paso usa
funciones phi sobre las coordenadas físicas y conserva la fuerza constante
algebraicamente. No se retocan los datos anteriores.

La etapa 4 permanece abierta: siguen el transiente térmico no lineal, el trabajo
espectral no térmico con poblaciones y los controles acoplados pertinentes.
El prototipo longitudinal demuestra una identidad de disponibilidad cuadrática
en el sector de fase constante y corriente nula; no acredita el sistema general
fuera del equilibrio ni una simulación del detector. La etapa 5 no comienza.
Producción, v1.0.0 y el circuito de tres estados de la memoria se conservan.
