# Preparación de la siguiente implementación

Esta entrega cierra el contrato que se indique en `closure_admission.json`.
Los ensayos de uso de una y dos celdas no sustituyen el sistema dinámico completo
de D.1. La siguiente implementación debe completar las reacciones electrón-fonón
en celdas con espectro móvil antes de incorporar la malla espacial y el circuito.

## Entradas que ya quedan disponibles

- Catálogo electrónico R2 persistido, con energía y sus derivadas en una misma
  representación; su archivo y módulo se verifican por hash. La prueba nueva
  `catalog_recheck.json` vuelve a contrastarlo con las referencias independientes.
- Operadores y pruebas de celdas en `pysnspd/experimental/cell_validation.py`.
  Los resultados documentan el contrato de cada operador y sus límites; un tiempo
  de ensayo sintético no debe interpretarse como tiempo material de NbN.
- Receta explícita y evidencia de preprocesado fonónico en esta entrega. El
  dictamen material delimita qué convención está sustentada por la fuente y qué
  usos siguen condicionados. El archivo original permanece identificado por hash.
- Una entrada fonónica Debye analítica para cerrar primero las identidades del
  algoritmo sin depender de la normalización absoluta de una tabla material.

## Trabajo concreto, en orden

1. Implementar dispersión, recombinación y ruptura de pares de D.14-D.19 como una
   lista compartida de eventos: la misma tasa actualiza electrones y fonones.
   Integrar el evento en los soportes reales del catálogo y comprobar el factor
   de dos de recombinación. Una reacción que envíe energía fuera del soporte debe
   detectarse, no perderse ni convertirse en calor añadido.
2. Acoplar esos eventos con el transporte a energía fija y la evolución de la
   amplitud. La reconstrucción de `p(x)` incorpora el movimiento del espectro;
   no se debe añadir otra vez la deriva de D.20. Comparar las fuerzas de energía
   con las consultadas en el catálogo durante la trayectoria.
3. Sustituir, cuando corresponda, la movilidad de ensayo por D.11-D.13 y separar
   sus parámetros de `tau_kin`. Incorporar la fuente D.18, escape D.19 y sus
   registros de energía compartidos. Volver a verificar por separado cada
   término antes de ejecutar todas las contribuciones juntas.
4. Sólo después de aceptar esas celdas, preparar el ensayo espacial D.4, punto 3:
   condición D.36, interfaces 2D/1D, reservorios, Poisson y circuito D.28-D.33.
   Un fallo de la condición espacial detiene ese ensayo; refinar el tiempo no
   convierte una rigidez negativa en una física admitida.

## Decisiones numéricas para empezar

La nueva implementación debe reutilizar la actualización exponencial del
intercambio lineal entre dos celdas. Su forma como combinación convexa evita
restar dos números casi iguales cuando una ocupación es muy pequeña. A campos
fijos, la relajación BGK también dispone de una actualización exponencial exacta.
El paso global se elige por los términos acoplados y su error de separación; no
por resolver innecesariamente estos dos intercambios elementales.

Para las colisiones, preparar primero la lista de transiciones y sus pesos
conservativos. El producto entre la energía de cada sector y los cambios de una
misma transición debe sumar cero antes de evolucionar una trayectoria. Separar
el almacenamiento de esa geometría de las tasas dependientes de las ocupaciones
permite reutilizarla. La cuadratura fonónica debe resolver bordes y picos de la
entrada, sin confundir las miles de muestras del archivo original con una
obligación de evaluar cada muestra en cada paso temporal.

El prototipo de transporte conserva su energía reconstruida y mide por separado
el sesgo respecto de R2 y el defecto del número de cuasipartículas. La siguiente
etapa debe integrar esa distinción en el balance cinético completo, o sustituir
la representación por otra que cumpla simultáneamente los momentos requeridos.
No debe ocultar el defecto con una corrección posterior de las poblaciones.

## Criterios que debe conservar la siguiente etapa

- Equilibrio de Fermi-Dirac/Bose-Einstein: cada par directo/inverso debe cancelar
  con la misma discretización, también entre celdas de gaps distintos.
- Conservación local del evento y global de la trayectoria, con fuentes, trabajo
  espectral y escape contabilizados una sola vez.
- `0 <= p <= 1`, `n >= 0`, rechazo de salidas de soporte y ausencia de recortes de
  poblaciones. Registrar los casos rechazados y los mínimos, no sólo promedios.
- Refinamientos independientes en tiempo, representación de energía y campos.
  Fijar las tolerancias antes de obtener los resultados de esa implementación.
- Comparación con referencias independientes y permanencia de todos los casos
  que hayan fallado en las iteraciones anteriores.
- Identificar por separado la aceptación numérica, los parámetros efectivos y
  la validez material. Una tabla saneada con normalización todavía condicional
  no acredita tasas absolutas ni latencias de un detector.

## Ejecución y entrega

Trabajar directamente en `main`, manteniendo inmutable `v1.0.0` y conservando la
evidencia histórica. Los nuevos operadores siguen en `pysnspd.experimental` hasta
la promoción explícita del sistema completo. El cuaderno
`/home/jdiaz/GEMINGA_COMMANDS.md` contiene la reproducción ligera de esta entrega.

Antes de un cálculo nuevo, estimar tiempo y memoria a partir de un ensayo pequeño
que tenga utilidad diagnóstica propia. Todo trabajo previsto de más de cinco
minutos se entrega al usuario mediante ese cuaderno. Un trabajo que agote su
límite de tiempo permanece incompleto; no se reinicia por fragmentos para eludir
la política. Cada siguiente informe debe presentar resultados y decisión de
continuación, con la teoría limitada a las referencias D/B necesarias.
