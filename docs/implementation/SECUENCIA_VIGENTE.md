# Secuencia vigente de implementación e investigación

Actualización del 24 de septiembre de 2026. Las etapas 1-3 se cerraron como
desarrollo con los límites registrados; 3.5 se cerró como investigación. La
etapa 4 continúa sin fotón. Los contratos y resultados previos se preservan.

| Etapa | Estado | Resultado o propósito |
|---|---|---|
| 1. Datos y catálogo | Desarrollo cerrado | Admisión explícita de fuentes, unidades y catálogo dentro del alcance registrado. |
| 2. Cinética de una y dos celdas | Desarrollo cerrado con límites | Acoplamientos conservativos y trayectorias registradas; certificado dinámico de malla incompleto. |
| 3. Infraestructura espacial, bordes y circuito | Desarrollo cerrado con pendientes explícitos | Energía y corriente comunes, empalme de campo 2D–1D, potencial, cargas de reservorio y circuito de la memoria; evidencia estática e instantánea. |
| **3.5. Dominio de confianza físico y numérico** | **Investigación cerrada con límites explícitos** | Establecer qué parámetros, estados, geometrías y ventanas de observación tienen fundamento físico, y qué márgenes numéricos necesitamos para estudiarlos. |
| 4. Núcleo y disipación | Euler dual y control longitudinal completados; escalas investigadas, unión dinámica pendiente | Cerrar el acoplamiento dinámico de trabajo espectral, poblaciones, carga, calor y puertos. |
| 5. Transientes completos y comparación experimental | No iniciada | Comparar el dispositivo y sus observables con el experimento, incluyendo incertidumbre y cadena de lectura. |

La [decisión de cierre de etapa 3](stage3/closure_20260923/closure_decision.json)
no declara completado todo el contrato histórico D.4.3. La [etapa 3.5](stage3_5/README.md)
se cerró como investigación y conserva las incertidumbres de cascada y tasas materiales.

La [referencia térmica espacial](stage4/spatial_energy_20260924/README.md)
recupera el contraste del núcleo al mismo corte espectral. La [continuación
autoconsistente](stage4/self_consistent_review_20260924/README.md) terminó con los
cuatro casos aceptados; no constituye todavía un cierre fuera del equilibrio
ni altera los requisitos dinámicos siguientes.

## Pendientes que no desaparecen con el cierre

Antes de interpretar resultados dinámicos de las etapas 4–5 deben completarse:

- La evolución espacial débil acoplada al potencial, las poblaciones y el circuito,
  con contraste temporal sobre los observables relevantes y balance integrado.
- La rama de reservorio dependiente de la corriente, su trabajo y las condiciones
  externas completas; el radio fijo del diagnóstico no resuelve este problema.
- El transporte cinético a igual energía en el empalme y una representación única
  de cada población fonónica compartida.
- La justificación dinámica de cualquier reducción 1D, la protección del dominio
  de interés frente a los extremos y la permanencia en el dominio espacial admitido.

Estos pendientes son requisitos de entrada al ensayo que los necesite. No se
convierten en una afirmación de validación mediante el cambio de número de etapa.

## Regla física para 3.5

El experimento objetivo es Korzh y colaboradores (2020), sobre resolución
temporal sub-3 ps; la tesis de Allmaras (2020) es una fuente complementaria.
Se seleccionará una configuración experimental identificable, evitando mezclar
geometrías, materiales o cadenas de lectura de distintas muestras.

La simulación espacial será 2D y, únicamente cuando se justifique, 1D. El espesor
sigue siendo un parámetro físico de la película, no una tercera coordenada de
la malla. El intervalo ilustrativo 1,5–6 para longitud 2D/ancho no es un rango
adoptado. Se investigará una ventana que permita observar la hotbelt y controlar
la influencia de los extremos sin gastar cómputo innecesario.

3.5 estudia procedencia, restricciones, márgenes y límites de interpretación de
los parámetros. No es un barrido de simulaciones para comprobar si todo pasa,
ni un ajuste de parámetros destinado únicamente a reproducir una curva. Una
resolución temporal experimental sub-3 ps tampoco fija por sí sola el paso de
tiempo del solver o la latencia determinista del modelo.

## Horizonte inicial y ejecución paralela

La [política de observación](stage4/followup_20260923/horizon_policy.md) prioriza
el gatillo Vout más un margen para la latencia relativa de la cinta de 80 nm.
Sólo se acorta la ventana: el sistema físico completo es el mismo, incluido el
circuito de la memoria. La recuperación larga se pospone. Un cruce no observado
antes del techo temporal queda censurado, sin latencia inventada.

Los casos y consultas independientes se distribuyen con un presupuesto común
máximo del 90 % de los recursos detectados. Geminga tiene actualmente 16 núcleos
físicos y 32 hilos; el límite es 28 hilos, reservando dos núcleos completos.
Se limita la anidación BLAS/OpenMP y la memoria disponible. El paralelismo no
exime de entregar manualmente cálculos previstos de más de cinco minutos.

La [revisión actual de etapa 4](stage4/final_kwt_20260924/README.md) contiene los
resultados, el contraste físico siguiente y el estado real de implementación.
