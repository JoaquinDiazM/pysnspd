# Respuesta uniforme: comparación independiente de la campaña

Se completaron 720 consultas, con 0 fallos. Son respuestas armónicas lineales a fondo uniforme sin corriente, sin fotón ni transiente de calor.

La referencia de cuadratura usa orden 24 por panel, transición E*=64 kBTc y 24 nodos en cada cola infinita. La transformación x=E*/|E| integra el mismo integrando hasta infinito; no ajusta contraterminos. Se compara orden 16 frente a 24 y transición 32 frente a 64. La variación de eta se informa por separado: eta es un regulador espectral, no una tasa física de relajación.

| Comparación numérica | Máximo | Supera 2 % |
|:--|--:|--:|
| Admitancia compleja, norma relativa | 1.14303e-10 % | 0/60 |
| Vout complejo, norma relativa | 2.94154e-15 % | 0/60 |
| Matriz de respuesta candidata, norma Frobenius relativa | 1.01355e-05 % | 0/300 |
| Núcleo, norma Frobenius relativa | 0.0019954 % | 0/300 |
| Diferencia de disipación / gamma KWT | 9.25419e-05 % | 0/300 |
| Diferencia Re Y / conductancia normal | 7.30981e-11 % | 0/60 |

La comparación relativa de Re Y o de gamma resuelta consigo mismas puede amplificar diferencias de una señal casi nula. El JSON conserva esas comparaciones, además de las escalas normal y KWT anteriores; no se ocultan signos ni se recortan valores.

| Omega = hbar omega/(kBTc) | gamma resuelta / gamma KWT, mínimo–máximo entre cinco modos; eta/Delta=0,002 |
|--:|--:|
| 0.02 | 1.142487e-06 – 9.580831e-06 |
| 0.05 | 3.573531e-07 – 2.000564e-06 |
| 0.1 | 2.426661e-07 – 7.522453e-07 |
| 0.2 | 2.162209e-07 – 3.575278e-07 |
| 0.5 | 2.210012e-07 – 2.304179e-07 |
| 1 | 2.877377e-07 – 3.304202e-07 |
| 2 | 3.191293e-05 – 4.194038e-05 |
| 3 | 0.000153279 – 0.0002027868 |
| 4 | 0.02418561 – 0.03270608 |
| 6 | 0.03399713 – 0.03769117 |

| Sensibilidad eta=0,005/0,01 frente a 0,002, cuadratura fina | Máximo |
|:--|--:|
| Admitancia compleja, norma relativa | 1.70675 % |
| Vout complejo, norma relativa | 6.71297e-05 % |
| Matriz de respuesta candidata, norma Frobenius relativa | 1.67164 % |
| Núcleo, norma Frobenius relativa | 2.16274 % |
| Diferencia de disipación / gamma KWT | 0.0812608 % |
| Diferencia Re Y / conductancia normal | 1.15659 % |

El aumento de disipación cerca de Omega≈2Delta/(kBTc)=3,52785 coincide con la apertura del canal de rotura de pares de la referencia BCS. Por debajo, una parte de las pequeñas colas disipativas sigue dependiendo de eta; no se utiliza para identificar una vida media material.

El máximo residuo cinético es 2.809e-15; el de la ecuación de neutralidad resuelta es 3.603e-13. El balance medio del circuito tiene defecto máximo 2.48e-33 W.

## Diagnóstico y corrección de la cola

En la campaña inicial, el residuo bruto de calibre con corte 64 no tendía a cero al refinar sólo los paneles: queda el término de extremo. Para el puerto uniforme, su valor esperado es `1 − integral(C−Omega/2,C+Omega/2)[tanh(E/2T) rho(E)] dE / Omega`. Con este corte ronda −0,000380, aproximadamente −0,038 %. Es distinto del residuo algebraico de la neutralidad acoplada. Una corrección de cola debe salir de esa integral, nunca de ajustar una constante a la respuesta. La matriz de respuesta de fase/potencial amplificaba esta omisión hasta un 54,9977 % al comparar los cortes 32 y 64, aunque el núcleo cambiaba menos del 0,07 %. La columna radial permanecía estable. Ese diagnóstico original se conserva en `finite_cutoff_review/uniform_analysis.md` y en `uniform_results/`.

La campaña actual integra ambas colas con la transformación recíproca. El piloto independiente comparó 12/24 nodos por cola y transiciones 32/64: el cambio máximo de la respuesta fue 1,44e−9 relativo. El defecto grande desaparece al recuperar la contribución faltante; no se relajó la tolerancia para aceptarlo. Las cifras de las tablas corresponden a las 720 consultas nuevas con colas, no al piloto ni a resultados retocados.

La resta térmica estática y el incremento retardado comparten la misma referencia. El circuito conserva los tres estados y sus constantes de la memoria. La pequeña identidad de potencia circuital no demuestra todavía el balance total de energía de la película.

## Alcance de la comparación polarizada pendiente

`biased_coupled_response.py` proyecta las respuestas completas del grafo alrededor de una rama con corriente, resuelve amplitud, fase y potencial juntos y devuelve residuos también fuera del subespacio modal. Eso permite medir los términos cruzados ausentes por simetría en el ensayo uniforme. Sin embargo, calcula respuesta de primer orden y potencia de puerto de segundo orden; no evoluciona poblaciones/fonones con B.41 ni reconstruye el calor de segundo orden. El balance circuital no cierra ese trabajo interno por sí solo. Su pase autorizaría el acoplamiento débil dentro del alcance registrado, no un transiente no lineal con fotón.

Las derivadas y constantes de los datos brutos se conservaron. El JSON registra hashes y cada comparación; no se ejecutó de nuevo el modelo para producir este análisis.
