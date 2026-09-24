# Respuesta de carga en frecuencia: revisión independiente

Las 24 evaluaciones terminaron en **21,19 s**, reutilizando espectros. El límite de frecuencia cero recupera las cuatro referencias estáticas con una diferencia máxima de **8,88×10⁻¹⁶**. Las 24 matrices completas tienen certificado remoto; además se recalcularon localmente las normas de tres mapas representativos.

## Resultado útil para el cierre lento

La tabla usa la sonda angular del núcleo radial 65² con eta/Delta=0,01. Cada cambio compara la respuesta espectral completa a esa frecuencia con su propio límite estático. El parámetro nu es omega·tD, con tD=0,441516 ps.

| nu | 1/omega (ps) | Cambio de corriente | Cambio de torque | Alcance |
|---:|---:|---:|---:|---|
| 0.0001 | 4415.2 | 0.002576 % | 0.004220 % | Ancla lenta |
| 0.001 | 441.52 | 0.025764 % | 0.042199 % | Ancla lenta |
| 0.01 | 44.152 | 0.257615 % | 0.421962 % | Exploratorio |
| 0.1 | 4.4152 | 2.551590 % | 4.194585 % | Exploratorio |
| 1 | 0.44152 | 18.905181 % | 33.694234 % | Exploratorio |

En los cuatro grupos, nu=0,001 cambia la corriente un 0,02564–0,02612 % y el torque un 0,04201–0,04368 %. Esto favorece conservar la forma espectral de la carga y eliminarla algebraicamente en un dominio lento justificado. No favorece reducir esa forma a una sola susceptibilidad térmica: su discrepancia angular sigue aproximadamente en 20,55 % para corriente y 77,53 % para torque en el caso de referencia.

**1/omega no es el período de la oscilación.** El período vale 2pi/omega. Tampoco esas escalas representan el tiempo de detección ni la duración de un transiente simulado.

## Signos y balances

Con la convención exp(−i omega t), la matriz es −L_TT−i nu M y el balance comprobado es rT+i nu M hT=0. El almacenamiento M=m Re(g) es positivo en la rama admitida; eta no se añadió como baño o tasa de colisión. En el límite normal, −L_TT es la rigidez difusiva positiva y el signo corresponde a difusión decreciente. Esto no sustituye un análisis de estabilidad de todos los modos del núcleo inhomogéneo.

Los máximos absolutos son: identidad de fuerza/corriente **1,6635×10⁻¹⁶**, balance dinámico integrado **2,2178×10⁻¹⁶** y ecuación de carga por energía **1,8323×10⁻¹⁶**. Son comprobaciones algebraicas y del solver, no precisiones físicas del detector.

La fuerza conserva por separado sus dos coordenadas cartesianas complejas, evitando confundir la fase temporal con la fase del condensado. La fase de la proyección espacial sobre la respuesta estática puede tener signos distintos para corriente y torque: no es un único retardo físico del dispositivo.

## Dónde termina la interpretación física

nu=0,0001 y 0,001 son anclas lentas, no un certificado uniforme de validez AC. A nu=0,01 ya se obtiene hbar·omega/eta=1,134 en el contorno más cercano. A nu=0,1 se obtiene hbar·omega/kBTb=1,922; a nu=1 la frecuencia también supera la escala del gap. Esas filas son pruebas exploratorias del cierre con espectro fijo, no mediciones de la banda del detector ni una justificación para una respuesta de pocos picosegundos.

La dinámica del espectro y el trabajo cuando cambia el condensado siguen pendientes. No se ha conectado un fotón ni un circuito nuevo y no se ha cerrado la etapa 4. El siguiente control debe conservar explícito ese alcance.
