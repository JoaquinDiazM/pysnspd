# Etapa 3.5: dominio de confianza físico y numérico

**Estado: investigación iniciada en su base documental y preparada para el estudio
de rangos. No se ha ejecutado una campaña nueva de simulaciones ni se han adoptado
rangos finales del dispositivo.** La secuencia vigente está en
[SECUENCIA_VIGENTE.md](../SECUENCIA_VIGENTE.md).

El objetivo es establecer dónde tiene sentido físico aplicar el modelo y qué
márgenes necesitamos para resolver sus observables. «Dominio de confianza» se
usa aquí como dominio de validez y uso sustentado; no significa un intervalo
estadístico de confianza con un porcentaje asignado.

El experimento principal es el de **Boris Korzh y colaboradores (2020)**, sobre
resolución temporal sub-3 ps. La tesis de Allmaras (2020) complementa la descripción
del dispositivo, los procesos físicos y la lectura. La identificación de una
muestra concreta y sus datos es previa a fijar un conjunto de parámetros.

## Qué se investiga

- Geometría física y dominio simulado: ancho, espesor, longitud 2D, continuaciones,
  distancia del impacto a los planos de borde y región donde observar la hotbelt.
- Material, espectros, colisiones, transporte, condensado y cierres efectivos:
  fuentes, unidades, escalas, incertidumbres, correlaciones y estados admisibles.
- Preparación por fotón: energía retenida, perfil, posición, ancho y origen de
  tiempo, separando la etapa óptica omitida de la evolución modelada.
- Bordes, sustrato, reservorios, circuito de la memoria y cadena de medida.
- Parámetros numéricos: soportes y cuadraturas espectrales, regularizaciones,
  interpolación o aceleración, mallas espaciales, tiempo y presupuestos de error.
- Definición de observables: formación de hotbelt, eventos de fase, latencia,
  señal eléctrica y distribución de tiempos; sus incertidumbres no son intercambiables.

El [inventario](parameter_inventory.json) recoge las familias y sus variables.
El [plan de investigación](research_plan.md) organiza cómo justificar límites;
la [ficha del experimento](reference_experiment.md) distingue los datos publicados
de lo que todavía debe verificarse.

## Geometría y protección de los extremos

Se estudian dominios espaciales 2D y, si la evolución lo justifica, regiones 1D.
La reducción 1D requiere evidencia de variación transversal despreciable en
los campos y distribuciones relevantes durante la ventana de interés. El
empalme de traza uniforme implementado en etapa 3 no aporta por sí solo esa prueba.

**1,5–6 anchos es un ejemplo de intervalo para investigar**, no un requisito
cerrado ni un rango físico verificado. La longitud debe permitir identificar
la hotbelt, mantener una distancia suficiente a los bordes y controlar el coste.
La posición del fotón, la extensión inicial y el tiempo observado afectan ese
compromiso. Una ecuación difusiva tiene colas inmediatas: la condición útil es
limitar su efecto sobre los observables, no exigir que ninguna excitación alcance
matemáticamente el borde.

## Entrega que debe producir esta etapa

Un conjunto trazable de rangos, márgenes y restricciones conjuntas. Cada entrada
debe indicar: fuente; valor o intervalo y unidades; muestra o escenario al que
pertenece; hipótesis; incertidumbre; dependencias con otras variables; observable
afectado; criterio para detectar extrapolación; y estado de la evidencia. También
debe declarar lo que permanece desconocido y qué medición o estudio permitiría
resolverlo. No se asignan tolerancias extremas por defecto.

Los contrastes dinámicos pendientes de etapa 3 seguirán siendo necesarios antes
de las etapas 4–5. Esta investigación prepara su selección física; no los sustituye.
