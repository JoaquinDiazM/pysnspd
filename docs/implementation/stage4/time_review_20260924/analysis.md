# Etapa 4: trayectoria térmica débil sobre todos los nodos

La corrida completó 1 ps en 71.10 s, con 246 aplicaciones por lotes del operador. La revisión verifica 37 checkpoints (16 observaciones y 21 pasos aceptados), las fuentes ejecutadas y 257 insumos del operador. No ejecuta nuevas ecuaciones físicas.

## Resultado que puede admitirse

Pasan las comparaciones de refinamiento registradas para este sistema afín térmico. El mayor desplazamiento total es 0.221713 % del gap de referencia, frente al corte declarado de 2 %. Ese corte es un dominio de trabajo lineal, no un rango medido del material. La etapa 4 sigue abierta; aún no se acredita el transiente no lineal.

| Perturbación, tras restar la base | Norma restante a 0,03 ps | Norma restante a 1 ps |
|---|---:|---:|
| Amplitud | 99.0481 % | 73.5926 % |
| Fase angular | 48.4260 % | 1.7511 % |

Son normas espaciales respecto a la misma perturbación inicial, no velocidades de respuesta del detector. La fase angular cruza la mitad entre 0,01 y 0,03 ps; el muestreo no permite dar un tiempo exacto de cruce.

## Refinamiento y colas

| Sonda | Observable | Máxima diferencia / señal inicial | Máxima diferencia / señal restante |
|---|---|---:|---:|
| amplitude | displacement | 0.00193278 % | 0.00262632 % |
| amplitude | current | 0.00140385 % | 0.00145268 % |
| amplitude | force_density | 0.00632199 % | 0.00873127 % |
| amplitude | phase_torque_density | 0.0940221 % | 15.8809 % |
| angular_phase | displacement | 0.0383158 % | 1.03842 % |
| angular_phase | current | 0.0425061 % | 0.793649 % |
| angular_phase | force_density | 0.113258 % | 17.5003 % |
| angular_phase | phase_torque_density | 0.115472 % | 160.316 % |

La aceptación usa la tolerancia previamente registrada: 1e-8 más 0,5 % de la mayor norma entre inicial y restante. Por eso no autoriza afirmar igual precisión relativa en las colas. El torque angular tardío puede diferir más que su valor restante; no se le asigna un tiempo de relajación preciso.

## Deriva, energía y límites

La base residual también se mueve: su desplazamiento final es 0.015577733 en norma de área, 5.199 veces la perturbación inicial de amplitud. Se resta esa trayectoria para identificar cada respuesta; un mapa total tardío de fase está dominado por esta deriva.

La energía cuadrática disminuye en todas las observaciones y su tasa calculada es negativa. La disipación aproximada KWT más normal es positiva. El mayor defecto instantáneo tasa más disipación es 0.0769197 % de la disipación; el mayor defecto constitutivo de velocidad es 0.231574 % de la velocidad total. Son defectos de Taylor evaluados sobre los estados guardados, no una verificación de la energía no lineal exacta. Las integrales trapezoidales de ocho observaciones se conservan como diagnóstico de muestreo y no certifican un balance temporal.

La coordenada auxiliar constante de Arnoldi deriva hasta 2.20994e-05 respecto de uno (0.00220994 %): es una pequeña variación numérica del forzamiento afín. Se registra sin imponer un umbral retrospectivo ni repetir la campaña. Puede preservarse exactamente en una futura versión del propagador.

El siguiente contraste útil es evaluar la respuesta térmica no lineal exacta en los estados ya guardados a 0, 0,03 y 1 ps. El punto de 0,03 ps conserva aproximadamente la mitad de la perturbación de fase; el final concentra la deriva y la truncación. No requiere repetir la trayectoria ni incorporar un fotón.

## Procedencia

`raw/extraction_receipt.json` certifica los originales conservados en Geminga. `raw/trajectory_compact.npz` conserva todos los campos observados, retirando únicamente la repetición de la geometría. `analysis.json` conserva las normas, verificaciones reproducidas, sensibilidades y hashes.
