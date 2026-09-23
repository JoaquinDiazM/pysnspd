# Auditoría de los seis ensayos abiertos

**PASS_SAVED_OPEN_AUDIT.** Las nueve fuentes y el registro coinciden con los hashes congelados. Los 21 archivos contienen exactamente los seis casos previstos, con sus registros de iteración y reservorio. Los resultados individuales y el resumen son idénticos. Se verificaron las métricas mediante aritmética de los campos guardados, sin ejecutar física nueva.

La campaña registrada duró **271.28 s** (4.521 min). Los tres dominios tienen dos resoluciones: elementos de 90 y 45 nm.

| Longitud | Nodos | Residuo estacionario | Error de corriente | Error de q terminal |
|---:|---:|---:|---:|---:|
| 360 nm | 17 | 1.93553e-05 | 0.000603483% | 0.07130059% |
| 360 nm | 33 | 6.31939e-07 | 9.04413e-06% | 0.004868235% |
| 720 nm | 33 | 1.93553e-05 | 0.000603483% | 0.07130059% |
| 720 nm | 65 | 6.31939e-07 | 9.044132e-06% | 0.004868235% |
| 1080 nm | 49 | 1.93553e-05 | 0.000603483% | 0.07130059% |
| 1080 nm | 97 | 6.31939e-07 | 9.044133e-06% | 0.004868235% |

Los errores de corriente y q terminal cumplen el presupuesto registrado de 1%; el residuo estacionario queda bajo 1e-4 y la amplitud bajo 0.001. Al duplicar los elementos, los errores disminuyen en las tres longitudes. La corriente de referencia es **8.63351 µA** y la amplitud **0.995218**.

**Los seis optimizadores hicieron cero iteraciones y una sola evaluación.** La hélice inicial ya satisfacía los criterios. La campaña comprueba consistencia estática y refinamiento de esa hélice; no demuestra relajación desde una perturbación. La diferencia de amplitud cercana al redondeo no representa recuperación dinámica.

El cambio de q terminal pasa aproximadamente de **0.0713006% a 0.00486823%** y el de corriente de **0.000603483% a 0.00000904413%**. Las dos resoluciones describen reducción de error frente a esta referencia: no justifican por sí solas un orden de convergencia universal ni un certificado de tres mallas.

A igual tamaño de elemento, la corriente media y las densidades de energía apenas cambian al alargar el dominio uniforme. Este resultado es coherente con una hélice uniforme; no mide reflexiones ni permite inferir transparencia ante una excitación localizada.

| Longitud | Inductancia resuelta | Inductancia exterior |
|---:|---:|---:|
| 360 nm | 0.166621 nH | 9.833379 nH |
| 720 nm | 0.333243 nH | 9.666757 nH |
| 1080 nm | 0.499864 nH | 9.500136 nH |

La partición diferencial suma 10 nH y se identifica en la misma rama antes de cualquier evolución. Se verificaron su escala lineal con la longitud y sus unidades. Las márgenes registradas del reservorio y de D.36 son positivas.

El alcance estático uniforme abierto queda aceptado para continuar el desarrollo acoplado. **La validación general de bordes y la etapa 3 completa siguen abiertas.** Faltan perturbaciones no uniformes, intercambio de poblaciones, balances del sistema acoplado y evolución espacial con circuito. No hay certificado de interfaz 2D–1D ni predicción de detección.

Límites de esta auditoría: D.36 y la estacionariedad se conservaron como métricas y registros, sin las matrices ni gradientes completos; se comprobaron sus márgenes y coherencia, pero no se recalcularon físicamente. El productor guardó hashes de fuentes y registro, no recibos individuales de salidas. El inventario SHA256 de esta revisión identifica los archivos observados sin atribuirles recibos previos.

Reproducción: `python sandbox/stage3_spatial/coupled_20260923/audit_open.py --verify-only`. Sin esa opción, el auditor crea `open_audit.json` y `open_audit.md` únicamente si todavía no existen. Las fuentes, resultados anteriores y tolerancias no se modificaron.
