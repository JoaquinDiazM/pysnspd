# Etapa 4: momentos resueltos y respuesta lenta de carga

La ejecución del usuario terminó correctamente en **198,12 s**: 181 consultas
espectrales, 362 respuestas cinéticas a dos sondas y siete comparaciones de
momentos. La nueva integración permite decidir qué representación de carga
conviene seguir estudiando, sin otra relajación del núcleo.

| Resultado | Valor | Decisión |
|---|---:|---|
| Defecto de la integral del peso térmico: 12 → 31 → 50 energías | 18,02 → 1,62 → 0,424 % | La malla fina resuelve la limitación anterior para esta comparación |
| Diferencia del potencial único, corriente | 20,55 %; variación numérica 2,87 % | Discrepancia resuelta; no admitir esa compresión como respuesta cuantitativa |
| Diferencia del potencial único, torque de fase | 77,53 %; variación numérica 3,64 % | La amplitud casi coincidente ocultaba un cambio importante de fase |
| Respuesta espectral instantánea frente al diagnóstico armónico, ν=0,001 | 0,0258 % corriente; 0,0422 % torque | Continuar estudiando la eliminación algebraica que conserva todas las energías |

Los porcentajes describen normas de respuestas infinitesimales, no errores
de latencia. La variación numérica suma cambios observados de cuadratura,
contorno y malla; no es una cota rigurosa. No se descarta la dinámica completa
de fase y potencial de la memoria a partir de este ensayo con gap fijo.

El [informe ilustrado](../../../../output/pdf/implementation/Informe_etapa_4_respuesta_de_carga.pdf)
presenta los resultados. La [decisión física](physics_decision.md) desarrolla
la reducción de Schur, el alcance del almacenamiento adiabático y el siguiente
control térmico con condensado móvil. La [secuencia vigente](../../SECUENCIA_VIGENTE.md)
conserva los pendientes del dispositivo.

## Evidencia y trazabilidad

- [Análisis de momentos](analysis.md) y [datos procesados](analysis.json):
  diferencias, sensibilidades y procedencia. `raw/` contiene metadatos y siete
  mapas integrados; 369 mapas originales fueron verificados en Geminga.
- [Diagnóstico de frecuencia](charge_frequency/README.md): 24 trabajos en
  21,19 s, sin nuevas resoluciones espectrales. Incluye fuentes, plan, 61 tests,
  balances y regresión de frecuencia cero. Los 24 mapas completos permanecen
  en scratch con hashes; el repositorio conserva un compacto reproducible.
- [Revisión independiente](frequency_review.md): reproduce normas del compacto,
  verifica signos y distingue la banda lenta de las frecuencias exploratorias.
- `previous_delivery_exact/` conserva los cuatro documentos vigentes de la
  entrega anterior; sus resultados y fuentes numéricas no se modificaron.

En ν=0,001 la escala angular 1/ω es 442 ps; el período es 2π/ω.
Los puntos ν≥0,01 son exploratorios para este cierre con espectros fijos.
Una buena respuesta lenta no acredita automáticamente la respuesta AC de pocos
picosegundos. Tampoco identifica el desequilibrio electrónico con el potencial
electrostático del circuito.

## Próximo avance

El [operador térmico de gap móvil](thermal_weak/README.md) pasó sus contrastes
en 32,94 s: 256 frecuencias, 2304 raíces espectrales y cinco pruebas nuevas.
La diferencia máxima del Jacobiano temporal frente a variaciones finitas de
escala 0,01 es 0,00485 %. Diferenciar implícitamente Usadel permite reutilizar
sus factorizaciones para diferentes perturbaciones sobre todos los nodos.
Se conserva la fuerza residual del núcleo y su derivada de movilidad: esa
corrección aporta aproximadamente 1,1 % del RHS. No se redefine la referencia
como equilibrio exacto ni se exige otra minimización térmica.

Sigue la [trayectoria afín térmica](thermal_time/README.md) hasta 1 ps: una evolución base y perturbaciones
débiles independientes de amplitud y fase, con KWT y potencial normal. Su
espacio de aproximación temporal representa todos los nodos; no impone unos
pocos perfiles como grados de libertad físicos. El residuo temporal, su
contraste refinado y el límite de validez débil se medirán explícitamente.

El [cuaderno de comandos](../../../GEMINGA_COMMANDS.md) describe la preparación,
su coste y las salidas. Una preparación del operador no equivale todavía a un
transiente del dispositivo. La etapa 4 sigue abierta por la evolución acoplada
y por el trabajo espectral no térmico. Etapa 5, producción y v1.0.0 permanecen
sin cambios; se conserva el circuito completo y el futuro horizonte Vout más
margen.
