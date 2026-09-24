# Ruta física práctica después del ensayo térmico no lineal

Revisión del 24 de septiembre de 2026. Esta nota decide qué trabajo hace falta
para conectar los bloques existentes. No sustituye los datos ni los veredictos
de la corrida anterior y no modifica producción.

**Recomendación: aceptar el alcance térmico medido y trasladar ahora el núcleo
común a la malla dual que utilizará el dispositivo.** No repetir la trayectoria
de 1 ps para resolver una pequeña señal cruzada al final de su relajación. La
precisión de esa señal se declara como límite. Las identidades de carga y
energía siguen siendo necesarias; una precisión relativa universal de cada
componente que tiende a cero no debe bloquear la conexión del modelo.

## Qué está publicado y se puede reutilizar

| Bloque | Evidencia consultada | Consecuencia para esta implementación |
|---|---|---|
| Cinética electrónica y fonónica, trabajo de un gap variable | Vodolazov (2017), ecuaciones 1–5 y 9–11 | Usar el término de trabajo y los integrales existentes con sus convenciones, sin crear una termodinámica alternativa. La ecuación 1 contiene explícitamente el término proporcional a la velocidad de la amplitud y a la derivada energética de la distribución. |
| Separación entre espectro y ocupación; modos de carga y energía | Belzig y colaboradores, secciones 2.3–2.7 | Reutilizar el oráculo espectral y el operador cinético ya implementados. El artículo anuncia su restricción a equilibrio y estados estacionarios fuera del equilibrio: no acredita por sí solo nuestro cierre temporal completo. |
| Código de cinética y datos materiales asociado al trabajo de 2025 | Repositorio de los autores `qnngroup/proj-KE-solver` | Su README identifica `solver.m`, las tablas DFPT y el uso de pyTDGL. Es una referencia directa para conservar la separación entre cinética y evolución espacial, no un motivo para reconstruir la malla. |
| Geometría dual y operadores espaciales | Código local derivado de pyTDGL y operadores FV de la memoria | Mantener nodos, caras, masas, contactos y orientación de aristas. El nuevo cierre espectral debe recibir esos objetos. |

Fuentes primarias: [Vodolazov, texto completo](https://arxiv.org/pdf/1611.06060),
[Belzig y colaboradores, texto completo](https://arxiv.org/html/cond-mat/9812297v2)
y [README del código de los autores de 2025](https://github.com/qnngroup/proj-KE-solver/blob/main/README.md).
La interpretación de qué reutilizar es la decisión de ingeniería de esta
revisión; no se atribuye a esos autores nuestro ensamblaje concreto.

## La transferencia a la malla dual requiere un adaptador

`ThermalGraph` ya acepta un grafo general. `rectangular_graph` construye uno
particular; las ecuaciones espectrales no dependen de que sea cartesiano. Las
entradas geométricas necesarias existen en `pysnspd.mesh`:

| Entrada del núcleo experimental | Geometría existente |
|---|---|
| Área adimensional del nodo | Área de su celda dual dividida por ℓ0² |
| Conductancia geométrica de una arista | Longitud de la cara dual dividida por longitud de la arista primal |
| Coordenadas adimensionales | Coordenadas físicas divididas por ℓ0 |
| Orientación y contactos | Los índices de aristas y etiquetas físicas de la malla existente |

Debe utilizarse una sola pareja compatible de áreas y caras: `Mesh.areas` con
su `EdgeMesh`, o los campos compatibles del `FVOperators` de producción. No
mezclar áreas de una construcción y caras de otra. Las longitudes y áreas se
transforman de unidades una sola vez. Las corrientes y fuerzas continúan siendo
derivadas de la misma acción de grafo.

La malla dual evita prolongar una implementación geométrica auxiliar. No se
afirma que cure un error temporal ni que aporte un término físico ausente.
La comparación útil es entre soluciones y observables sobre la geometría que
se conservará, no entre nombres de métodos de malla.

## Un siguiente ensayo concreto, sin otra cadena de diagnósticos

El próximo ensayo ejecutable debe tener **una geometría dual, una preparación
suave y dos trayectorias: referencia y referencia perturbada**. Así se mide el
cambio causado por la perturbación sin confundirlo con la deriva de la referencia.

1. Usar el rectángulo provisional de 80 por 160 nm preparado en
   [dual_mesh](dual_mesh/README.md), con referencia térmica uniforme y contactos
   de extremo explícitos. El material y cierre térmico se conservan. Es un
   **nuevo control suave**, no una continuación del vórtice con su borde radial
   anterior ni una geometría de Korzh ya admitida. Si se compara directamente
   una malla cartesiana con ésta, ambas deben resolver este mismo control.
   No declarar estacionaria una referencia cuya fuerza no sea cero.
2. Usar una perturbación suave de amplitud y fase, compacta en el interior y
   resuelta por varias aristas. No normalizarla dividiendo por un gap que se
   anula en el centro de un núcleo. Se pueden reutilizar las direcciones suaves
   ya declaradas; no hace falta introducir un depósito fotónico ni fijar ahora
   su ancho desconocido.
3. Conservar la movilidad KWT y el cierre normal de este control térmico.
   Reutilizar el integrador ya implementado si su coste es adecuado; comparar
   métodos alternativos sólo si un fallo o el coste lo exige. No rehacer la
   prueba temporal cartesiana por el cambio de constructor geométrico.
4. Registrar amplitud espacial, diferencia compleja respecto a la trayectoria
   base, corriente, potencia KWT y balance integrado. Si se necesita una
   comparación de resolución, variar una sola resolución sobre esos observables;
   no abrir un barrido de radios, energías, tolerancias y perturbaciones a la vez.

Este ensayo comprueba la transferencia al soporte espacial final. No necesita
una precisión relativa prefijada para todo torque residual ni una nueva prueba
radial del continuo. Antes de ejecutarlo se declara el cambio físico de interés
y se compara el error con esa escala y con la señal absoluta, incluyendo el
error de interpolación inicial. Si una señal queda bajo el error medido, se
informa como no resuelta; no se transforma en un fallo de todo el desarrollo.

Un defecto de conservación, de signo, de unidades, una solución no finita o
un error de los observables principales que impida interpretar el ensayo sí
obligan a corregir. La ejecución larga se entrega al usuario, con avance/ETA y
el presupuesto compartido de recursos ya establecido.

## Cómo unir después los bloques ya construidos

El módulo `longitudinal_reciprocal.py` es una pieza útil, pero su nombre no debe
ocultar su alcance: es lineal, tiene gap real, fase constante y corriente nula.
Admite distribuciones con forma energética no térmica y conserva una identidad
de disponibilidad cuadrática. **Disponibilidad significa energía libre de la
perturbación; no es la energía interna total del dispositivo.** No sirve como
cierre implícito de carga, calor y circuito en un detector polarizado.

Su primer uso debe hacerse en la misma malla dual, con una perturbación suave
y pequeña, conservando las variables de población y sin reajustar su forma a
una temperatura. El coeficiente recíproco se contrasta con el término de
amplitud de Vodolazov ya citado. No se requiere resolver una familia extensa de
casos singulares antes de esta unión. El estado uniforme proporciona además
un espectro analítico y permite probar el ensamblaje sin repetir la resolución
completa de un núcleo para cada frecuencia.

La conexión con etapa 2 debe conservar sus eventos electrón–fonón y sus momentos
conservativos, pero evaluarlos con el espectro y las unidades comunes. No se
puede insertar sin transformación la antigua coordenada de conteo local cuando
la DOS depende de los vecinos: el transporte a energía física común y el
trabajo al cambiar el gap deben permanecer explícitos. Ésta es una adaptación
concreta del ensamblaje, no una invitación a rederivar toda la teoría.

Para la siguiente trayectoria polarizada se conectan los dos modos cinéticos,
los contactos y las tres variables del circuito de la memoria. El intercambio
con reservorios, la disipación KWT y el calor de Joule se contabilizan una sola
vez. La prueba mínima de esa conexión es un transiente espacial suave acoplado
con su balance integrado y salida eléctrica identificada. No hace falta exigir
primero la recuperación completa del detector.

## Qué se pospone y qué no puede omitirse antes del fotón

Se posponen el refinamiento de colas casi nulas, otra certificación puntual de
la DOS, nuevos barridos de núcleos extremos y cualquier prueba de recuperación
de nanosegundos. Tampoco se hace de este control térmico una validación de
latencia, hotbelt o jitter.

Antes de interpretar un fotón siguen haciendo falta la transferencia fotónica
admitida, la normalización material de las tasas utilizadas, el trabajo cinético
compatible con el gap espacial, los bordes del dispositivo y el circuito
completo. Son entradas físicas de la predicción, no tolerancias superpuestas.
La incertidumbre de cascada y el ancho inicial continúan explícitos conforme al
cierre 3.5. No se resuelven imponiendo un pulso suave de prueba.

En los informes, cada figura debe indicar la variable, unidad, dominio,
preparación, instante y normalización. Si se resta la referencia, debe decirlo
en la figura o su pie. «Perturbación angular del gap» no se presenta como
voltaje; una norma global no se presenta como valor local; energía libre no se
presenta como energía interna. Cuando dos curvas se superponen, añadir una
curva de diferencia con unidad y escala explícitas es más útil que cambiar
sólo sus colores.
