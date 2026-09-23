# Revisión independiente del piloto abierto y del circuito

La evidencia guardada es consistente y sus fuentes coinciden con los hashes congelados: **PASS_REVIEWED_SAVED_EVIDENCE**. Se verificaron nueve fuentes del piloto abierto y tres del diagnóstico eléctrico, los artefactos eléctricos y la identidad del resultado abierto con su resumen. La revisión sólo realizó lectura, hashes y aritmética de datos guardados.

El piloto abierto de 360 nm y 17 nodos cumple los criterios registrados: residuo de estacionariedad **1.93553e-5**, error relativo de corriente **0.000603483%**, error de flujo regularizado convertido a q terminal **0.0713006%**, y diferencia de amplitud **1.11e-16**. D.36 conserva margen positivo. Duró **20.86 s**.

**El optimizador hizo cero iteraciones y una evaluación.** La hélice inicial, construida con la rama uniforme del reservorio, ya satisfacía las tolerancias. Esto comprueba consistencia estática de ese estado; no demuestra relajación desde una perturbación ni un dominio de convergencia del solucionador.

La referencia del reservorio da amplitud 0.995218, corriente **8.63351 µA**, inductancia resuelta **0.166621 nH** y exterior **9.833379 nH**, cuya suma es 10 nH. El diagnóstico del circuito usó por separado **7 nH** como entrada exterior ilustrativa. Estos dos resultados no constituyen una simulación acoplada ni deben mezclarse como una sola identificación del dispositivo.

El grafo eléctrico respeta corriente total, gauge y descomposición de potencia. La integración del circuito frente a la referencia exponencial matricial tiene error máximo de estado escalado **4.29526e-11**, frente al presupuesto 1e-8. La diferencia máxima de Vout es **6.44518e-14 V** y el residuo instantáneo guardado de potencia es **1.69407e-21 W**. La reconstrucción algebraica independiente del balance a partir de los estados también cumple el presupuesto registrado. No se calculó aquí otra trayectoria.

El máximo muestreado de **1.425 mV** corresponde a una resistencia impuesta que salta de 0 a 1000 Ω. **No es un pulso de detección ni una latencia SNSPD.** El balance comprobado es instantáneo; no sustituye el balance de energía integrado del futuro sistema acoplado.

El log registra **168 pruebas y 23 subpruebas aceptadas en 9.87 s**. Se preservó su hash; el log compacto no enumera la invocación ni cada prueba. El productor abierto guardó hashes de fuentes y registro, pero no recibos individuales de sus salidas: esta revisión incorpora sus hashes observados sin atribuirles un recibo previo inexistente.

El alcance queda separado:

- **3A:** cierre de desarrollo en los ensayos periódicos registrados, conservando una respuesta espacial `ROUND_OFF_LIMITED` sin pase relativo.
- **3B:** piloto uniforme abierto e identidades algebraicas de interfaces. Sigue pendiente el lote de seis casos de longitud y resolución, preparado para ejecución del usuario; tampoco se han demostrado fronteras transparentes ni relajación no uniforme.
- **3D:** circuito ante entrada resistiva prescrita y álgebra del puerto de potencial. Falta el acoplamiento con la evolución espacial del condensado y las poblaciones.

No hay fundamento en estos resultados para cerrar la etapa 3 global, afirmar simulación del detector completo ni promover el modelo a producción. Las limitaciones son de alcance; no invalidan los pases registrados.

Los hashes, comprobaciones y métricas quedan en `validation_review.json`. Ningún registro, manifiesto, resultado histórico o umbral fue modificado.
