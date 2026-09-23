# Pruebas restantes tras la convergencia temporal de una celda

Estado: **diseño para prerregistro; no es un plan ejecutable ni un certificado**.
Fecha: 2026-09-22. No se ejecutaron RHS ni nuevas trayectorias para elaborar esta propuesta.
Se conserva `acceptance_criteria.json`, SHA256
`48a56a76b64b58b81b176325a6a34cc2e6c690bfde3b4eb739a8cd2723c702cf`.

La convergencia SSP de una celda ya permite continuar. El lote histórico de 23
tareas no debe relanzarse: repite resultados disponibles, sus pasos 40/80 fallaron
o fueron descartados con justificación, y no contiene todos los diagnósticos
necesarios para cerrar literalmente los criterios de fronteras y equilibrio.
El siguiente lote debe contener todas las fases pendientes, con evaluación entre
fases y parada al primer fallo, sin modificar umbrales ni reintentar a escondidas.

## Evidencia que se conserva

Los cuatro resultados nuevos tienen el mismo operador, catálogo, malla y estado
inicial. El evaluador temporal devuelve PASS. La siguiente tabla transcribe los
archivos `raw/one_ssp_*.json` y `raw/one_ssp_time_assessment.json`; la auditoría
independiente de procedencia pertenece al certificado final.

| Pasos, duración 2 | Error máximo respecto a SSP1280 | Residuo máximo de energía | Tiempo observado |
|---:|---:|---:|---:|
| 160 | 5.415703e-5 | 5.437656e-8 | 98.43 s |
| 320 | 1.569984e-5 | 6.741160e-9 | 195.84 s |
| 640 | 3.741810e-6 | 8.391935e-10 | 390.22 s |
| 1280 | Referencia separada | 1.045992e-10 | 776.06 s |

Las reducciones medidas son 3.45 y 4.20. El paso de 160 satisface su propio
presupuesto temporal de 1e-4 y el balance de 1e-7. Esto no demuestra que 160 pasos
sean suficientes para dos celdas, otras mallas o estados iniciales. El presupuesto
suplementario de 2.5e-5 del plan de recuperación es más estricto; ONE160 no lo
satisface y ONE320 sí. Ese presupuesto no debe confundirse con el criterio
temporal principal, ni suprimirse silenciosamente donde se lo haya prerregistrado.

Se retienen, dentro de sus escenarios declarados, los controles estáticos de
energía y fuerzas, KWT/Debye, equilibrio de eventos, caras físicas, entropía,
transporte continuo, referencia continua al corte IR=0.005, refinamiento
electrónico separado, interpolación fonónica, omisiones infrarrojas y superiores,
y controles negativos. Las fuentes están congeladas. El ensayo SSP de orden
tres para BGK y la equivalencia de su mapa no limitado con el RHS están en
`recovery_20260921/test_limited_ssp.py` y en su recibo de 41 pruebas aprobadas.
No se sustituyen esos resultados por el simple balance de una trayectoria.

## Inventario exacto de lo que falta

| Puerta | Evidencia existente | Validación todavía necesaria |
|---|---|---|
| Tiempo de dos celdas con SSP | RK4 candidato aprobado; SSP sólo un intervalo en malla fina | SSP160/320/640 frente a SSP1280, duración 2, parámetros TWO congelados; todos los observables y poblaciones |
| Malla electrónica dinámica | 630 y 1260 estados con RK4; 2520 falló positividad; prueba SSP corta 2520 | Tres trayectorias completas SSP con 630/1260/2520, fonones 1025; precisión del candidato 630 y descenso de errores |
| Malla fonónica dinámica | Interpolación estática aprobada | Tres trayectorias SSP con 1025/2049/4097, electrones 630; precisión del candidato 1025 y descenso de errores |
| Contaminación temporal de mallas | Sólo tiempo en candidato; SSP fino activa limitador | Control de paso propio de las mallas modificadas, especialmente 2520; no inferir su precisión del candidato |
| Equilibrio completo | RHS inicial y equilibrio estático | Trayectoria completa estacionaria, y evaluación explícita del cambio de amplitudes, energías y poblaciones |
| Fronteras del sistema acoplado | Pilotos antiguos en otra malla; caras iniciales del RHS | Casos 630/1025 de vacío fonónico y banda electrónica exactamente 0/1, con todas las invariantes |
| Fronteras de canales aislados | Test SSP pequeño con reacción y transporte simultáneos; no aislamiento por canal | Trayectorias de emisión, absorción, recombinación, creación y transporte aislados, desde las caras aplicables |
| Campos y soporte reales | Piloto TWO40, 22 muestras, PASS | Al menos diez tiempos físicos por celda, extremos incluidos, en las trayectorias finalmente admitidas |
| Campo del equilibrio con Gamma=0.1 | Referencia causal estática existe | Usarla explícitamente en los puntos visitados; el evaluador de campos a Gamma=0 rechaza este caso |
| Actividad SSP y balance de baño/fuente | Actividad RK4 aprobada | Recalcular en la trayectoria SSP aceptada; amplitud >=0.02 y ambas transferencias >=1% de la escala de excitación registrada |
| Orden suave del método actual | BGK SSP de orden tres; escape sólo medido antes con RK4 | Pequeño control SSP de escape frente a la exponencial exacta; no requiere red de eventos grande |
| Cierre reproducible | Regresión anterior y hashes | Auditoría final de archivos, mapas de criterios y límites de alcance; suite acotada si cambian auxiliares |

La prueba aislada no debe forzar un estado electrónico totalmente lleno dentro
de una inversión BGK de temperatura positiva. Se fijan los campos y se desactivan
los canales ajenos al ensayo, usando la misma ley de eventos y el mismo mapa de
extensiones. El estado lleno sí es admisible para el operador de colisiones
aislado. Un esquema docente de cuatro nodos aporta una prueba unitaria; no
reemplaza la evidencia con el catálogo declarado cuando ésta se invoque para
cerrar el modelo de celda.

## Secuencia propuesta para un solo lote completo

1. **Admisión previa sin integración.** Verificar hashes y resultados ONE ya
   guardados; reutilizarlos sin ejecutar. Preparar y probar todos los evaluadores
   que todavía faltan antes de pedir el lote largo. Los controles nuevos de
   software deben rechazar fuentes, métodos, mallas, estados iniciales o duración
   diferentes. Una pareja temporal no puede figurar como admisión de tres niveles.

2. **Controles cortos que pueden descubrir fallos temprano.** Equilibrio TWO:
   duración 0.1, 10 pasos, escape 15 y calentamiento 0. Fronteras ONE: duración
   0.2, 40 pasos, escape infinito y calentamiento 0, escenarios
   `phonon_vacuum` y `sparse_electrons`. Son los parámetros ya registrados en el
   plan anterior. Añadir los cinco canales aislados anteriores, con al menos
   diez tiempos guardados y paso declarado. Comprobar poblaciones, energía,
   balance de número aplicable y evolución interior real; no sólo evaluar el
   estado inicial. Evaluar el equilibrio contra el estado estacionario, no sólo
   contra el balance total. Si el paso de fronteras no cumple las invariantes,
   el lote se detiene y conserva la evidencia.

3. **Tiempo TWO.** Ejecutar 160, 320, 640 y 1280 pasos hasta tiempo 2. Mantener
   exactamente `case=two`, `escape=15`, `heating=0.01`, electrones 630, fonones
   1025, IR=0.005 y el resto de argumentos del caso TWO histórico. Exigir PASS
   del evaluador SSP frente a 1280 y las invariantes de cada trayectoria. Con los
   resultados medidos, seleccionar el menor N del conjunto que cumpla el
   presupuesto temporal previamente definido para la comparación de mallas.
   No elegirlo ahora usando sólo el resultado ONE. El resultado N ya existe:
   no se repite para iniciar cada refinamiento.

4. **Fonones antes del cálculo electrónico caro.** Ejecutar 2049 y 4097 nodos a
   N y 2N, electrones 630; cuatro trayectorias nuevas. Cada pareja controla la
   contaminación temporal en su propia malla. Usar el presupuesto suplementario
   ya registrado de 2.5e-5, manteniendo todas las invariantes. Comparar después
   1025/2049/4097 a N sobre los mismos tiempos físicos, con el candidato como
   primera fila. Esta fase precede a 2520 para descubrir un posible fallo
   fonónico antes de gastar horas en la otra familia.

5. **Electrones.** Ejecutar 1260 y 2520 estados a N y 2N, fonones 1025; cuatro
   trayectorias nuevas. Aplicar los mismos controles pareados y luego el
   contraste 630/1260/2520 a N. Registrar el número de eventos limitados, el
   defecto integrado de flujo y su cambio con el paso. Un factor mínimo pequeño
   no demuestra error macroscópico, y una masa pequeña no permite aceptar una
   población negativa. La convergencia medida es la decisión.

6. **Posprocesamiento y auditoría final.** Comprobar campos, fuerzas, derivadas,
   soporte, actividad y entrada/escape en las trayectorias seleccionadas; conservar
   el análisis del equilibrio de Gamma finito separado. Publicar el mapa completo
   de criterios con enlaces y hashes. Sólo emitir admisión si todas las puertas
   obligatorias medidas pasan. El alcance sigue siendo celdas sintéticas; no
   material NbN absoluto, circuito, problema espacial ni transiente completo.

Esta matriz base añade **15 trayectorias acopladas**: 3 especiales, 4 temporales,
4 fonónicas y 4 electrónicas; además de los ensayos aislados cortos y evaluadores.
No añade tres trayectorias infrarrojas largas: el criterio de omisión IR ya tiene
una referencia analítica independiente. Tampoco repite las cuatro ONE aprobadas.

El evaluador de mallas actual exige igualdad de todos los parámetros temporales.
Por eso la ruta anterior usa un N común. Si se decide aprovechar distintos N por
malla, antes del lote se necesita un **nuevo contrato explícito**: tiempos comunes
sin interpolación temporal, una pareja propia por malla, presupuesto de
contaminación medido por observable y separación de ese error del de malla.
Quitar simplemente la comprobación de `steps` del evaluador sería incorrecto.
La alternativa no debe incorporarse después de observar un fallo sin un nuevo
registro. Las parejas suplementarias no sustituyen la secuencia temporal TWO.

## Coste y recursos: estimaciones, no mediciones futuras

El nuevo lote ONE tardó 1462 s de pared según el manifiesto; la suma de tiempos
de sus cuatro trayectorias es 1460.54 s. La tasa asintótica es aproximadamente
0.606 s por paso de una celda. Los archivos RK4 completos dan 174.66 s para
TWO160 y 769.77 s para TWO1260/160. El primer paso SSP de TWO2520 costó 32.11 s
hasta el callback y 43.88 s incluyendo diagnósticos: sólo es una medida inicial,
no una tasa garantizada para toda la trayectoria.

| Bloque | Estimación con N=160 | Base de la estimación |
|---|---:|---|
| TWO160/320/640/1280 | 45–65 min | Dos celdas aproximadamente duplican el coste ONE, más transporte |
| Parejas fonónicas 2049 y 4097 | 15–30 min | Mismo número de estados electrónicos; no crecimiento cuadrático de la red por nodos fonónicos |
| Pareja electrónica 1260, N/2N | 35–55 min | Relación de costes RK4 medidos y sobrecoste SSP de inventarios |
| Pareja electrónica 2520, N/2N | 3–4.5 h | Extrapolación prudente del paso SSP fino medido; domina el lote |
| Especiales, aislados y auditorías | 5–20 min | Estimación; algunos diagnósticos deben acotarse por timeout |

Total orientativo si N=160: **5–7.5 horas, un hilo**. Si las mallas requieren
N=320, su bloque aproximadamente se duplica; no prometer el mismo coste. El caso
fino puede variar al evolucionar el condensado y reconstruir la red. No hay motivo
para repetir una referencia DOP853 larga.

Reservar **12 GiB** para el lote con 2520 estados, coherente con el plan anterior.
No se ha medido un pico RSS nuevo. El algoritmo conserva una sola gran red de
reacciones viva cada vez, pero construye vectores temporales de índices, tasas,
coeficientes y demanda. Una reserva de 2 GiB sólo se había propuesto para las
trayectorias de 630 estados. Ejecutar secuencialmente con OMP, OpenBLAS y MKL a
un hilo; no lanzar refinamientos en paralelo.

## Trabajo autónomo ligero y preparación previa al lote

- Verificación de hashes, manifiesto, formas, poblaciones, balances almacenados,
  tiempos comunes, actividad y limitador: lectura y álgebra de JSON/NPZ; segundos.
- Campos a Gamma=0: el control de 22 muestras anterior tardó 0.074 s. Puede
  hacerse ahora sobre ONE y luego sobre cada trayectoria seleccionada.
- Soporte superior: el control TWO630 de 22 muestras tardó 2.48 s. Reconstruye
  eventos y calcula actividades, aunque no integra; no confundirlo con lectura
  pura. En las mallas grandes usar timeout240 y registrar cualquier incompleto.
- El evaluador causal estático independiente tardó 1.81 s; adaptar su lectura de
  estados reales para el equilibrio, con comprobación de fuentes. La existencia
  del archivo estático no certifica automáticamente esos puntos.
- Orden suave SSP de escape y guardas del nuevo orquestador: problemas pequeños,
  sin transientes grandes. Los tests ya aprobados se reutilizan por hash.
- Preparar previamente los controles aislados y de estacionariedad, el mapa de
  puertas, los recibos y la reanudación por **identidad completa** de resultados.
  Conservar archivos fallidos; no sobrescribir. Si el lote se detiene, los
  resultados aprobados siguen disponibles y no deben calcularse otra vez.

No se propone ejecutar por el agente ninguno de los bloques largos, ni dividirlos
para eludir el límite. El comando único del lote completo se escribirá en
`/home/jdiaz/GEMINGA_COMMANDS.md` cuando exista el plan congelado y probado.
