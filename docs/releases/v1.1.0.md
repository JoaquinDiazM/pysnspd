# pySNSPD 1.1.0: cierre de la etapa 5

12 de octubre de 2026. Se completan la evolución conjunta del modelo neutral efectivo BCS/KWT, su balance independiente y los controles de paso, malla, longitud y espectro. El estado y los supuestos están en [MODELO_VIGENTE](../implementation/MODELO_VIGENTE.md); la [admisión final](../implementation/stage5/effective_validation/final_20261012/admission.json) identifica datos y decisiones.

A 775 nm se observa una hotbelt transversal a 15,3 ps y el cruce eléctrico de 0,1 mV a 31,9413 ps, seguido de 5 ps. A 1550 nm no hay cruce en 200 ps. Los cuatro refinamientos de 775 nm pasan en 0–36,6 ps, con confirmaciones cubiertas hasta 33 ps. Todas las colas externas conservan auditoría individual, sin certificar convergencia de cuatro ejes. El margen espectral de señal del 5 % se adoptó después de los datos; se conserva el fallo previo al 2 %.

La implementación usa una energía común BCS+gradiente+fonones+circuito, colisiones y transporte conservativos, neutralidad eléctrica, tres estados circuitales y punto medio implícito. La API es pysnspd.effective; las pipelines 05 físicas y visuales están disponibles en [comandos presentes](../GEMINGA_COMMANDS_PRESENT.md). La referencia exploratoria utiliza 54 estados; 80 es un control opcional.

Los nominales de Allmaras (2020), τee(Tc)=5 ps y τep(Tc)=24,7 ps, son calibradores libres. Material acústico, depósito postcascada y rigidez radial/KWT son hipótesis explícitas. La admisión no acredita calibración experimental universal, jitter estadístico, latencia óptica absoluta ni reset total.

El [informe físico final](../../output/pdf/implementation/Informe_final_etapa5_modelo_efectivo_20261011.pdf) y los documentos [Principal](../../Matematica/Principal.pdf), [Auxiliar](../../Matematica/Auxiliar.pdf), [Precomputo](../../Matematica/Precomputo.pdf) están compilados y revisados. Los originales numéricos, logs, builds y renders permanecen en scratch; Git conserva fuentes, certificados pequeños y PDF finales.

La memoria permanece congelada en v1.0.0. El candidato de primer orden Moyal rechazado, los resultados previos y las políticas anteriores permanecen históricos.
