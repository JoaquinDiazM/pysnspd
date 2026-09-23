# Secuencia vigente de implementación e investigación

Actualización del 23 de septiembre de 2026. Esta entrada añade la etapa 3.5 y
registra el alcance del cierre de desarrollo de la etapa 3 autorizado por el
usuario. Los contratos y resultados anteriores se conservan con sus dictámenes.

| Etapa | Estado | Resultado o propósito |
|---|---|---|
| 1. Datos y catálogo | Desarrollo cerrado | Admisión explícita de fuentes, unidades y catálogo dentro del alcance registrado. |
| 2. Cinética de una y dos celdas | Desarrollo cerrado con límites | Acoplamientos conservativos y trayectorias registradas; certificado dinámico de malla incompleto. |
| 3. Infraestructura espacial, bordes y circuito | Desarrollo cerrado con pendientes explícitos | Energía y corriente comunes, empalme de campo 2D–1D, potencial, cargas de reservorio y circuito de la memoria; evidencia estática e instantánea. |
| **3.5. Dominio de confianza físico y numérico** | **Investigación iniciada: fuentes e inventario; rangos finales pendientes** | Establecer qué parámetros, estados, geometrías y ventanas de observación tienen fundamento físico, y qué márgenes numéricos necesitamos para estudiarlos. |
| 4. Núcleo y disipación | No iniciada | Contrastar los cierres físicos dentro del dominio que proponga 3.5 y con las capacidades dinámicas requeridas. |
| 5. Transientes completos y comparación experimental | No iniciada | Comparar el dispositivo y sus observables con el experimento, incluyendo incertidumbre y cadena de lectura. |

La [decisión de cierre de etapa 3](stage3/closure_20260923/closure_decision.json)
no declara completado todo el contrato histórico D.4.3. La [etapa 3.5](stage3_5/README.md)
puede avanzar documental y analíticamente sin una nueva campaña de transientes.

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
