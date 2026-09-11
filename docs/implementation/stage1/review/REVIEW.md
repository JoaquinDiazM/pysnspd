# Revisión independiente de la etapa 1

La implementación conserva la energía y sus derivadas dentro del mismo objeto
numérico. La revisión no encontró un error de factores, unidades o derivación
en los casos comprobados. La procedencia del archivo NbN pudo establecerse;
su contenido todavía no supera la admisión material. La prueba independiente
del regulador muestra que la etapa 1 aún no está cerrada para validación física:
el regulador del piloto altera un 8.61% la respuesta a corrientes pequeñas en
una población con excitaciones de baja energía.

## Controles independientes

- El vacío se comparó con una integral independiente sobre frecuencia imaginaria,
  en 24 estados que incluyen las ramas con y sin gap: error escalado máximo
  `4.15e-9`. La referencia no importa el módulo experimental.
- Se comprobaron 126 combinaciones de amplitud, depareamiento y regulador
  retardado. No aparecieron soluciones rechazadas ni DOS negativa; el residuo
  normalizado máximo de la identidad espectral fue `8.89e-16`.
- La energía no térmica, las dos fuerzas y sus derivadas cruzadas coinciden
  con diferencias independientes del mismo catálogo: error absoluto máximo
  `1.07e-9`, incluidas las derivadas cruzadas. Una derivada unilateral en
  `Gamma=0`, incluidos los extremos de amplitud, dio error menor que `1.07e-9`.
- La conversión de la pendiente a corriente SI se verificó variando directamente
  el gradiente de fase: error relativo no térmico menor que `3.29e-11`.
- Guardar y cargar conserva exactamente los datos, el regulador, las escalas
  físicas y las referencias. Se comprobaron nombres sin extensión, población
  fuera del intervalo de Pauli, NaN, infinito y consultas fuera del soporte.
- El límite normal reproduce exactamente el factor de conteo `4 N0` y anula
  las fuerzas del condensado. Doscientas consultas interiores mantuvieron
  energías positivas y ordenadas.

Los scripts `independent_vacuum_reference.py` y `verify_query_contract.py`
reproducen estos controles en segundos. Sus JSON identifican el código revisado.
Las poblaciones empleadas son sintéticas y no representan un evento del detector.

## Correcciones y límites que importan para continuar

La revisión detectó un defecto de persistencia con nombres sin extensión y pidió
verificar el conteo espectral, las conversiones SI y las poblaciones próximas al
borde inferior del espectro. El nombre de archivo ya se conserva. El catálogo
añadió esas comprobaciones y una prueba con ocupaciones de baja energía.

La interpolación directa de energías perdió monotonía en un ensayo grueso y se
rechazó. La implementación final interpola los logaritmos de incrementos positivos
y reconstruye la energía mediante su suma. El primer interpolante de esos
logaritmos todavía erraba un 14.96% en la pendiente de corriente del caso
`|Delta|/Delta0=0.72`, `Gamma=0`, a regulador y cuadratura idénticos. Las
pendientes espectrales explícitas se incorporaron mediante interpolación Hermite
en Gamma. La revisión verificó la regla de la cadena, el borde y la persistencia
de esta reconstrucción; el esquema del archivo pasó a v2. La aceptación de la
malla final corresponde a los resultados del generador, no a una identidad
interna de derivadas.

Hay otro error independiente que persiste aunque se corrija el interpolante.
Para `p(x)=0.2 exp[-x/(0.25 Delta0)]`, con la misma amplitud y Gamma cero, la
integración adaptativa da una pendiente normalizada `0.737011909` cuando
`eta/Delta0=0.001`. El límite causal exacto es `0.678584013`: la diferencia es
8.61%. Reducir eta a `0.0001` y `0.00001` deja diferencias de 2.99% y 0.986%.
La corriente en Gamma cero es cero; estos porcentajes afectan su pendiente al
aumentar una corriente pequeña. El cálculo separa el regulador de la malla de
campos y de la cuadratura de estados. `low_energy_regulator_reference.py` y su
JSON lo reproducen sin importar el módulo experimental.

El archivo NbN contiene 361 muestras negativas de DOS y cuatro pares de abscisas
duplicadas con ordenadas distintas. Las unidades y la normalización átomo/celda
siguen sin acreditarse. La cabecera local de 48 bytes no pertenece al archivo
numérico publicado. La trazabilidad por hash no autoriza a renormalizar,
recortar o fusionar esos datos. Las pruebas de unidades y energía fonónica usan
espectros sintéticos admitidos expresamente como tales.

La recomendación es continuar con la convergencia del catálogo antes de promoverlo
a la siguiente etapa de validación física. Pueden prepararse ensayos de balance
en una celda sintética con el regulador declarado, pero no acreditan haber cerrado
la etapa 1. No corresponde iniciar predicciones materiales NbN ni transitorios
completos con estos datos. La sensibilidad de núcleo y la admisibilidad espacial
detectadas en 0.4 tampoco quedan resueltas por construir un catálogo uniforme.

La etiqueta anotada remota `v1.0.0` se verificó de manera independiente: apunta
al commit `5ea0cd6a08d2cae414c82944f8230469ca5aa7d6`, anterior a esta implementación.
