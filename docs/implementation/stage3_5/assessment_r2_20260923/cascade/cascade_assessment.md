# Cascada: reparto y perfil fonónico recuperados de Allmaras

23 de septiembre de 2026. **El ancho permanece abierto.** Este diagnóstico completa una parte documental de la caracterización: separa electrones y fonones en figuras históricas y comprueba una restricción de forma. No simula una cascada ni traslada el resultado a Korzh de 775/1550 nm a 0.9 K.

## Resultado útil

En los dos instantes de A20 que pueden contrastarse entre las figuras 2.6(a) y 2.8(a), la energía fonónica domina, pero queda una fracción electrónica apreciable. Al normalizar el perfil fonónico por **su propia energía**, sus radios siguen sin ser compatibles simultáneamente con una sola gaussiana 2D dentro del presupuesto de lectura:

| Tiempo después del inicio A20 | $E_e/E_\gamma$ | $E_{\rm ph}/E_\gamma$ | $R_{50,\rm ph}$ | $R_{90,\rm ph}$ | $R_{90,\rm ph}/R_{50,\rm ph}$ |
|---:|---:|---:|---:|---:|---:|
| 0.0935 ps | 6.56–7.82% | 92.18–93.44% | 1.51–1.62 nm | 3.30–3.63 nm | 2.039–2.406 |
| 0.1870 ps | 7.07–8.97% | 91.03–92.93% | 1.55–1.67 nm | 3.43–3.88 nm | 2.056–2.510 |

Los intervalos combinan lectura del gráfico y energía no mostrada fuera del radio usado; **no son intervalos de confianza física**. Sus extremos están correlacionados y no deben muestrearse como parámetros independientes.

Toda gaussiana 2D isotrópica tiene $R_{90}/R_{50}=\sqrt{\ln10/\ln2}=1.823$, cualquiera que sea su anchura. Este valor queda por debajo de ambos intervalos fonónicos. La diferencia que antes se observó para energía total, por tanto, no se explica únicamente por mezclar la cola electrónica con fonones. La nube fonónica de este caso histórico también tiene un núcleo y una cola que una gaussiana única no reproduce simultáneamente.

Esta restricción no descarta cualquier reducción gaussiana útil. Obliga a declarar qué observable se preserva y qué error introduce en el otro, y a estudiar si ese error modifica la latencia relativa o la formación del hotbelt. No se obtiene un ancho admisible eligiendo el percentil que mejor convenga.

## Fuente y reloj

Se usó el PDF local de [Allmaras 2020](https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf), figuras 2.6(a) y 2.8(a), páginas PDF 34–35, páginas impresas 22–23. El escenario del capítulo 2 usa 1 eV, dos pares electrón–hueco iniciales dentro de 1 nm, $D=0.5$ cm²/s y baño $T_c/2=4.325$ K. La comparación radial omite ee, escape y difusión fonónica; usa espectro electrónico normal y coordenadas cilíndricas. El radio inicial de 1 nm pertenece a ese cálculo, no se adopta para el proyecto.

El eje se convierte con el tiempo **cinético** $\tau_0=1.87$ ns del capítulo: $5\times10^{-5}\tau_0=0.0935$ ps y $10^{-4}\tau_0=0.187$ ps. No se usa el tiempo GL del modelo posterior. La figura 2.8 también permite muestras visibles a $3\times10^{-5}\tau_0=0.0561$ ps, pero la figura 2.6 no ofrece ese mismo instante: no se interpola entre tiempos para inventar una comparación.

La figura 2.6 muestra energías acumuladas radialmente divididas por la energía del fotón; distingue electrones mediante curvas sólidas y fonones mediante curvas discontinuas. La figura 2.8 muestra energías dentro de radios fijos y añade la curva punteada de energía total. Estos denominadores se conservaron antes de calcular percentiles de cada subsistema.

## Extracción y controles

Las figuras son imágenes RGB originales incrustadas en el PDF, no curvas vectoriales. Se extrajeron sus bytes sin reescalar ni suavizar: figura 2.6(a), 1750×1380 píxeles; figura 2.8(a), 1750×1313. Las regiones ocupadas por las leyendas se excluyeron expresamente.

La selección usa los colores originales, posición de ejes y vecindades visibles de cada trazo. La interpolación local atraviesa sólo los pequeños huecos del patrón discontinuo; la inversión de las curvas acumuladas limita esos huecos a 28 píxeles. No se reconstruyen curvas ocultas ni se extrapolan extremos. La incertidumbre de lectura incluye grosor, dispersión local y ±2 píxeles de calibración de ejes. Para invertir percentiles se reserva además una envolvente vertical de seis píxeles y un margen horizontal de calibración. Estos son presupuestos declarados de digitalización, no estimaciones estadísticas del cálculo original.

Se recuperaron 12 muestras de la figura 2.8: radios 2, 3, 4 y 5 nm, en tiempos de 0.0561, 0.0935 y 0.187 ps, con las tres curvas separadas. Sus sumas electrónica + fonónica difieren de la curva total en como máximo $0.00120E_\gamma$, dentro de las envolventes de lectura. Las cuatro comparaciones independientes entre figuras 2.6 y 2.8, a radio 5 nm y los dos tiempos comunes, difieren en menos de $0.00040E_\gamma$ y sus envolventes se solapan.

## Cómo se obtuvo el denominador fonónico sin borrar la cola

Se leen $e_R=E_e(r<R)/E_\gamma$ y $p_R=E_{\rm ph}(r<R)/E_\gamma$ a $R=9.5$ nm. Bajo los supuestos del cálculo original —excesos no negativos, energía total inicial igual a $E_\gamma$, sin escape ni trabajo externo—,

$$
e_R\le\frac{E_e}{E_\gamma}\le1-p_R,\qquad
p_R\le\frac{E_{\rm ph}}{E_\gamma}\le1-e_R.
$$

Se propagan las envolventes de lectura al aplicar estas cotas. No se supone que 9.5 nm contenga toda la energía ni que el último píxel sea una meseta exacta.

Para un percentil fonónico $q$ se busca el radio donde la energía fonónica acumulada por fotón alcanza $q\,E_{\rm ph}/E_\gamma$. Como ese denominador está acotado, se obtienen intervalos de radios, no un radio artificialmente exacto. Este procedimiento difiere de buscar directamente 0.5 o 0.9 de la energía del fotón: esa búsqueda produciría otros radios.

Las cotas anteriores **no** pueden aplicarse sin cambios a un caso con escape, pérdidas laterales o trabajo de corriente. Allí se necesita el balance completo antes de inferir el denominador.

## Qué permite decidir y qué sigue faltando

Se han cerrado dos ambigüedades de la lectura de A20:

1. El radio total no era un radio fonónico normalizado. Ahora existen cotas fonónicas propias para dos instantes y verificaciones cruzadas entre figuras.
2. La hipótesis de una gaussiana única no conserva ambos percentiles del perfil fonónico histórico. Una eventual reducción debe evaluar su error de forma.

No basta para cerrar la caracterización física de Korzh. Faltan perfiles a 0.9 K para cada color, reparto y pérdidas al instante de transferencia, profundidad y superficies, espectros espaciales y datos materiales con tasas admitidas. Tampoco se recuperan de estos gráficos un segundo momento completo ni la densidad pico con precisión controlada. La omisión de ee, propagación fonónica y escape limita la transferencia del resultado histórico.

En particular, una entrada puramente fonónica con electrones en el baño no reproduce exactamente este escenario A20 a los tiempos extraídos: queda al menos aproximadamente un 6.6–7.1% de energía electrónica según las cotas inferiores. No se declara que ese residuo sea irrelevante para latencia; se debe evaluar su efecto o escoger una transferencia efectiva con justificación independiente.

La siguiente decisión útil para la preparación fotónica es identificar/admitir los datos que alimentan la cascada y diseñar el cálculo reducido o la recuperación de datos originales que faltan. Este resultado evita pedir una simulación para volver a demostrar la incompatibilidad de forma histórica, pero no sustituye la caracterización por color. **No se adopta ningún ancho.** Después de este diagnóstico, el usuario eligió «Cerrar investigación 3.5 y preparar etapa 4 sin fotón». Se cierra así la investigación de rangos y límites; la caracterización fotónica sigue pendiente, como registra [la decisión de cierre](../closure_decision.json).

## Reproducción

El script [cascade_extract_partition.py](../../../../../sandbox/stage3_5/assessment_r2_20260923/cascade_extract_partition.py) sólo lee el PDF y extrae píxeles. Desde la raíz del repositorio:

    python sandbox/stage3_5/assessment_r2_20260923/cascade_extract_partition.py

Requiere Python, pdfplumber y Pillow, junto al PDF con el hash esperado. La salida [cascade_partition.json](cascade_partition.json) guarda el hash del PDF y de sus imágenes, el hash del script, las calibraciones, puntos locales, cotas, percentiles y controles. No descarga fuentes ni ejecuta física.
