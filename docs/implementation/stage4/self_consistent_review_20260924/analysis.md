# Núcleo térmico autoconsistente: revisión de resultados

Los cuatro casos alcanzaron el criterio anunciado. Se acepta este control estático y no se solicita otra corrida estática más estricta. La etapa 4 continúa abierta para el enlace espectral retardado, cinético y disipativo.

La continuación terminó en 484.52 s y resolvió 12672 consultas espectrales. Se comprobaron los SHA-256 de 23 checkpoints completos y la coherencia algebraica de los cuatro estados finales. El postproceso no volvió a resolver espectros.

| Caso | Barridos totales | Residuo RMS del núcleo | Máximo libre | Energía final |
|---|---:|---:|---:|---:|
| radial_65_N128 | 15 | 0.08442% | 0.12633% | -175.19740542 |
| radial_65_N256 | 16 | 0.09742% | 0.14669% | -174.96655194 |
| radial_129_N256 | 16 | 0.09760% | 0.14661% | -174.94387083 |
| asymmetric_65_N256 | 17 | 0.08952% | 0.19533% | -174.96381306 |

Todos los bloques redujeron la misma energía y conservaron el arrollamiento de fase +1. Se verificó el estado final con su propio espectro: el campo next_d es una propuesta distinta y no se mezcló con el espectro anterior para fingir convergencia. Los máximos libres se muestran aparte porque el criterio era RMS en r≤4ℓ₀.

## Sensibilidad útil para continuar

| Cambio | Diferencia del gap (L2) | Diferencia de densidad de corriente (L2) | Cambio del radio de medio gap |
|---|---:|---:|---:|
| cutoff_N128_to_N256 | 0.0258% | 0.6151% | 0.1599% |
| mesh_65_to_129 | 0.0384% | 0.0798% | 0.0144% |

Las comparaciones usan estados relajados separadamente, áreas duales y coordenadas comunes. La corriente de cada enlace se divide por el ancho de su cara dual antes de reconstruir una densidad nodal: comparar directamente corrientes integradas de mallas distintas produciría un cambio geométrico espurio. Las diferencias son sensibilidades observadas, no cotas rigurosas de error.

El vórtice asimétrico se desplaza a aproximadamente (-0.097364, 0.000000)ℓ₀. Su arrollamiento permanece +1. El mínimo nodal de |Δ| no debe interpretarse como desaparición del vórtice: el cero queda entre nodos; su ubicación se obtuvo del interpolante bilineal complejo.

## Decisión

El problema estático que motivó esta corrida queda resuelto al nivel declarado. No se vuelve a ajustar el cierre local ni se exige convergencia arbitrariamente más dura. Se continúa con la correspondencia entre este oráculo térmico y la descripción espectral retardada/cinética que necesita el transiente. La suma de Matsubara representa equilibrio; no autoriza sustituir las poblaciones no térmicas por una temperatura ni extender esta parametrización compleja por sustitución ingenua de frecuencia.

No se ha validado aún el transiente de un fotón, la transferencia de energía inicial, la latencia de Korzh ni el circuito dinámico. Los índices de barrido son iteraciones de minimización, nunca tiempos físicos.
