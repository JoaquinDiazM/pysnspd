# Referencia experimental para preparar la etapa 3.5

Estado al 23 de septiembre de 2026: **investigación inicial; ningún nuevo transiente ejecutado ni dominio físico admitido**. El cierre del desarrollo espacial de la etapa 3 permite preparar esta comparación. No demuestra todavía transporte cinético en una interfaz 2D–1D, formación de un hotbelt ni reproducción de un experimento. El registro numérico y de procedencia está en [source_register.json](source_register.json).

## Fuentes y versiones

- **K20:** Korzh et al., *Demonstration of sub-3 ps temporal resolution with a superconducting nanowire single-photon detector*, Nature Photonics **14**, 250–255 (2020), [DOI y registro editorial](https://doi.org/10.1038/s41566-020-0589-x). Se leyó la [versión aceptada y suplemento, 37 páginas](https://eprints.lancs.ac.uk/id/eprint/140252/3/Binder1.pdf), identificada como tal por el [repositorio institucional de Lancaster](https://eprints.lancs.ac.uk/id/eprint/140252/). Las páginas siguientes son páginas del PDF, empezando en 1: artículo 1–26; suplemento 27–37. No son las páginas de la revista.
- **A20:** J. P. Allmaras, *Modeling and Development of Superconducting Nanowire Single-Photon Detectors*, tesis doctoral, Caltech (2020), [registro](https://thesis.caltech.edu/13748/), [PDF original](https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf), [DOI](https://doi.org/10.7907/wgak-vs11). Se leyó la copia local ya utilizada por modelo 0.4, con SHA-256 registrado. En el capítulo 3, página impresa = página PDF − 12. Se verificaron visualmente las páginas impresas 81 y 95, donde la extracción textual degrada símbolos.
- **K18:** el [preprint de 2018](https://arxiv.org/abs/1804.06839) es antecedente, no sustituto de K20: informa 2.7 ± 0.2 ps y 4.6 ± 0.2 ps, frente a 2.6 ± 0.2 ps y 4.3 ± 0.2 ps del resumen publicado. No combinar sus números, circuitos detallados o instrumental sin identificación de versión.

## Qué se comparará

**Jitter FWHM** es el ancho a media altura de la distribución de tiempos registrados —la IRF—, no la duración de un transiente individual. **Latencia** es el retardo de respuesta; **latencia relativa** compara dos energías ópticas con la misma cadena y calibración. En K20 se obtiene la posición del máximo ajustado de la IRF. La latencia absoluta no se identifica simplemente con ese máximo. **Resolución temporal de sistema** incluye la cadena de adquisición; los récords del título no son una medida aislada del tiempo intrínseco del condensado. [K20, PDF 3–5, 12, 25–26; A20, pp. 64–65 (PDF 76–77).]

Para un transiente determinista se compararán primero formas de señal y retardos definidos con un umbral y una cadena de lectura explícitos. Una comparación de FWHM exige además una población de eventos y fuentes de dispersión justificadas. Normalizar una curva de conteos no la convierte en eficiencia absoluta.

## Casos que deben mantenerse separados

| Caso | Geometría y condiciones | Uso y localización primaria |
|---|---|---|
| K20, familia NbN | Espesor nominal 7 nm; longitud activa 5 μm; anchos 60/80/100/120 nm; baño 0.9 K | Fig. 2–3 y métodos, PDF 9–11, 23–24. No identifica un único dispositivo con todos los récords. |
| K20, comparación de latencia | Ancho 80 nm; 775/1550 nm; 21.5 y 15.5 μA | Fig. 1, PDF 22. A 21.5 μA: diferencia 1550−775 de 4.2 ± 0.4 ps, PDF 5. |
| K20, dependencia espectral | Ancho 120 nm; 273–1550 nm | Fig. 4, PDF 25–26. Adquisiciones a 40 y 80 GS/s no son intercambiables. |
| K20, simulación 1D | Ancho 80 nm; hotbelt inicial 40 nm | Tabla suplementaria 1, PDF 33. El hotbelt es una condición modelada; no mide el diámetro del foco ni demuestra homogeneización transversal. |
| A20, experimento adicional | Anchos 60/80/100/120 nm; baños 1 y 4 K; 1550/1064/775/532 nm; polarizaciones TE/TM | Cap. 3, pp. 65–71 (PDF 77–83). La figura 3.3 corresponde a 100 nm y 1 K. |
| A20, modelos distintos | 1D a 2 K, hotbelt 40 nm; ajuste 2D posterior a 4 K | pp. 95 y 98–101 (PDF 107, 110–113). No asignar automáticamente 0.9 K a estas simulaciones. |

K20 también estudia WSi y una demostración de escaneo de distancia. Quedan fuera de la referencia NbN inicial: no transferirles parámetros ni tratar su IRF de 6.2 ps como el récord del detector. La tesis contiene otros dispositivos en capítulos posteriores; tampoco se incorporan a este caso.

## Material, circuito y deposición: tres procedencias distintas

| Procedencia | Valores localizados | Restricción de uso |
|---|---|---|
| Fabricación K20 | NbN sobre Si/300 nm SiO₂; resistencia de hoja inicial 340 Ω/□ a temperatura ambiente; RRR 0.8; Tc 8.65 K | Métodos, PDF 9–10. No sustituir por la resistencia ajustada de la simulación. |
| Lectura K20 | Inductor añadido de 1.5 mm, 96 nH estimados con 64 pH/□; línea de 1 μm; amplificador CITLF1 a 4 K, 50 dB y 1.5 GHz nominales; osciloscopio 80 GS/s y filtro de 6 GHz | PDF 9–12, 23. Los 6 GHz del procesamiento no son el ancho de banda nominal del amplificador. |
| Tabla de simulación K20 | D=0.5 cm²/s; R□=608 Ω/□; parámetro fonónico γ=60; τee(Tc)=6 ps; τep(Tc)=24.7 ps; τesc=20 ps; fracción retenida 0.667; σF(1550)=92 meV y σno uniforme(1550)=40 meV | Suplemento, PDF 31–33. Son entradas/ajustes del modelo histórico, no un catálogo microscópico experimental admitido. |
| Modelo 2D anterior de A20 | D=0.5 cm²/s; R□=600 Ω/□; τesc=9.4 ps; γ calculado 23.8 y aproximado por 20 | p. 81 (PDF 93); esta parametrización tampoco es la tabla de K20. |

La adenda circuital vigente de `modelo_v0_4` conserva el circuito de la memoria. Para reproducir K20 habrá que justificar su correspondencia con el montaje experimental; **96 nH no reemplaza automáticamente un parámetro del circuito de la memoria**. La parte inductiva ya resuelta espacialmente se separará de la externa, y la señal se comparará en el puerto y después del filtrado declarados.

La óptica de K20 usa fuentes y conversiones distintas según la longitud de onda; la latencia 775/1550 nm exige sincronización y calibración del recorrido óptico. A20, figura 3.2, describe fuentes 1550/1064 nm y generación de segundo armónico 775/532 nm. No se fija aquí un pulso térmico universal. A20 pp. 80 y 101 (PDF 92, 113) distingue la inicialización térmica del proceso de conversión inicial y advierte del posible doble conteo de escape fonónico. Esto exige registrar qué energía entra en la cinética y qué pérdidas ya están incluidas.

## Dominio de confianza que falta construir

Los siguientes son **objetivos de investigación y criterios por prerregistrar**, no tolerancias aprobadas ni resultados:

1. Elegir un caso experimental identificable: material, muestra, ancho, baño, corriente, longitud de onda, polarización y adquisición. Resolver antes la discrepancia entre los 33 μA usados en la estimación de corriente de ruptura y los 35 μA del extremo de la familia en K20 (PDF 13 y 24); no inferir que son el mismo punto.
2. Seguir un transiente 2D con distribuciones electrónicas y fonónicas, amplitud, fase, corriente y energía. Mostrar la formación —o ausencia— de una región excitada que cubra el ancho. Un gráfico de temperatura promedio no demuestra un hotbelt.
3. Admitir la continuación 1D sólo donde las variaciones transversales relevantes para **los observables elegidos** sean pequeñas durante la ventana temporal usada, y comprobarlo contra una extensión 2D. Si quedan vórtices, asimetrías de deposición o corrientes transversales relevantes, conservar 2D. A20 p. 98 (PDF 110) no respalda una equivalencia física general entre sus ajustes 1D y 2D.
4. Variar la distancia hasta interfaces y terminales, manteniendo el experimento físico, y comprobar que energía, retardo y señal no cambien por la posición artificial del corte. Distinguir convergencia espacial, temporal, espectral y tamaño del dominio. La cifra ilustrativa **L₂D/W = 1.5–6 no se adopta como intervalo de confianza**.
5. Cerrar conservación y remapeo a energía física común en la interfaz cinética antes de etapas 4–5. Comparar tanto campos como flujos, incluidas poblaciones no térmicas, y contabilizar el circuito/lectura. La validación estática anterior no acredita estos pasos.

## Pendientes de procedencia y acceso

No se ha obtenido un conjunto público de eventos crudos que permita repetir el ajuste de la IRF; el registro editorial indica disponibilidad de datos a petición del autor correspondiente. Tampoco se han digitalizado curvas ni atribuido incertidumbre cero a valores sin barras. Quedan por precisar: identificador de cada dispositivo, corriente de cada récord, geometría completa de los tapers, foco/absorción espacial, transferencia compleja y ruido de lectura, umbral exacto por adquisición, covarianzas del ajuste y parámetros del sustrato de la muestra elegida. La tabla suplementaria fue leída textualmente; la captura visual remota falló, por lo que su tipografía de símbolos debe corroborarse antes de transcribir ecuaciones desde ella. D y las unidades de tiempo se contrastaron con A20 p. 95.

Esta preparación permite formular el primer experimento numérico trazable. No afirma que el modelo actual reproduzca los 2.6 ps, ni asigna al modelo un rango de validez a partir de una sola geometría o de los parámetros ajustados de otro modelo.
