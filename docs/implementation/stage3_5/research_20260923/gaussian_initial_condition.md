# Del disco térmico a una condición inicial gaussiana

Investigación documental del 23 de septiembre de 2026. No se ejecutaron cascadas ni transientes. **La gaussiana fue una elección del proyecto; las fuentes revisadas permiten precisar su significado, pero no convierten 5–20 nm en un rango aceptado.** Tampoco se adopta aquí un intervalo de retención de energía. El registro [gaussian_sources.json](gaussian_sources.json) separa datos publicados, conversiones e hipótesis.

## 1. Qué hizo realmente Zotova–Vodolazov en 2012

El artículo omite la absorción y termalización inicial. Su cálculo comienza con temperatura electrónica uniforme elevada dentro de un círculo —semicírculo en el borde—, conservando energía con capacidad térmica constante. Propone $R_{\rm init}\sim\sqrt{D\tau_{ee}}$ y después evoluciona temperatura, condensado y potencial. **No impone inicialmente un agujero con $\Delta=0$.** El supuesto es $\tau_{ee}\ll\tau_{ep}$; la circularidad es una simplificación isotrópica, no una imagen medida. [ZV12, sección II, pp. 2–3](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevB.85.024509/fulltext).

| Parámetro publicado | Valor |
|---|---:|
| $D$, $\tau_{ee}$, $\sqrt{D\tau_{ee}}$ | 0.45 cm²/s; 7 ps; aproximadamente 18 nm |
| Radios probados / radio presentado | 9 y 18 nm / 9 nm |
| Espesor; $T_c$; baño | 4 nm; 10 K; 5 K |
| Coherencia GL citada | 7.5 nm |

El mismo trabajo **sí usa una gaussiana** en su estimación analítica:

$$
T(r,t)-T_b=\frac{\beta}{4\pi Dt}\exp\!\left(-\frac{r^2}{4Dt}\right),
\qquad \beta=\frac{E_\gamma}{C_vd}.
$$

Es el núcleo de difusión de una fuente puntual, omitiendo pérdidas y calentamiento adicional [ec. (13), p. 5]. Otra aproximación analítica, de London, reemplaza la supresión espacial por un disco efectivo. Son tres construcciones distintas. Los autores advierten el alcance cualitativo a baja temperatura [sección IV, p. 7]. Ninguna determina directamente la anchura fonónica de Korzh.

## 2. Antecedentes: el contorno normal no es toda la nube

El refinamiento de Semenov et al. escribe la concentración electrónica como

$$
C(r,t)=\frac{M(t)}{4\pi Ddt}\exp[-r^2/(4Dt)].
$$

La población total $M(t)$ crece durante la cascada. El radio donde la concentración supera un umbral para destruir superconductividad es distinto del ancho de esta nube. Quasipartículas fuera de ese contorno también reducen la capacidad de transportar corriente; por tanto, ni el diámetro normal ni la longitud de coherencia identifican automáticamente una desviación estándar. [Preprint de los autores, sección III, ec. (2), pp. 4–5](https://arxiv.org/pdf/cond-mat/0410633), antecedente de [EPJB 47, 495 (2005)](https://doi.org/10.1140/epjb/e2005-00351-8).

El registro institucional identifica el antecedente de [Semenov–Gol’tsman–Korneev de 2001](https://elib.dlr.de/18281/), pero no ofrece texto completo. No se atribuyen a ese artículo números ni ecuaciones que aquí sólo se verificaron en 2005/2012.

La circularidad de un **contorno** es compatible con una densidad gaussiana: todos sus contornos de igual densidad son círculos. La diferencia relevante es entre un perfil con borde abrupto y uno que decrece suavemente, y entre energía electrónica, energía fonónica y supresión del condensado.

## 3. Una anchura con significado inequívoco

Se conserva la convención espacial de D.30. Para una deposición isotrópica 2D, lejos de bordes y uniforme a través del espesor $d$:

$$
u(r,t_0)=\frac{E_{\rm in}}{2\pi d\,s^2}
              \exp\!\left(-\frac{r^2}{2s^2}\right),\qquad
\int u\,d^2r\,d=E_{\rm in}.
$$

$u$ es densidad de **energía del subsistema al que se deposita**, no temperatura ni $|\Delta|$. $s$ es la desviación estándar de cada coordenada cartesiana. $t_0$ es el instante de transferencia desde la etapa omitida; no tiene por qué coincidir con la absorción óptica. Las siguientes identidades son derivaciones de esta definición:

$$
\langle r^2\rangle=2s^2,\quad
u(0)=\frac{E_{\rm in}}{2\pi d s^2},\quad
\frac{E(r<R)}{E_{\rm in}}=1-e^{-R^2/(2s^2)}.
$$

De ahí $R_{50}=1.17741s$, $R_{90}=2.14597s$ y $R_{95}=2.44775s$. La anchura a mitad de altura de un corte que pasa por el centro es $2.35482s$; no es un diámetro que encierre el 50% de energía. Reducir $s$ a la mitad cuadruplica la densidad central a energía fija.

**Conversión de un disco uniforme.** Un disco de radio $R$ contiene $\langle r^2\rangle=R^2/2$. Igualar energía total y segundo momento a la gaussiana da

$$
s=\frac R2.
$$

Esta conversión preserva la extensión cuadrática, pero la gaussiana tiene un pico dos veces mayor que el disco. Igualar energía y pico daría, en cambio, $s=R/\sqrt2$ y duplicaría el segundo momento. No existe una conversión que conserve simultáneamente energía, pico y segundo momento. Así, los discos históricos de 9/18 nm equivalen a $s=4.5/9$ nm **sólo bajo la primera regla**, sin transformarse en mediciones gaussianas.

La [comparación pedagógica de perfiles](figures/03_disk_gaussian.png) ilustra estas diferencias; es una comparación analítica de formas, no una simulación de cascada.

Un radio que contiene el 90% permite otra conversión, $s=R_{90}/2.14597$, pero no determina el segundo momento de una distribución no gaussiana: el 10% restante puede transportar una cola amplia. Una esfera uniforme 3D tampoco es un disco: su proyección tiene $\langle x^2\rangle=R^2/5$, no $R^2/4$. Antes de convertir un radio hay que declarar qué distribución y qué momento se están igualando.

## 4. Qué aporta la cascada radial de Allmaras

Allmaras resuelve cinética radial de un metal normal en coordenadas cilíndricas y esféricas. En el caso relevante de su capítulo 2 deposita **1 eV dentro de 1 nm** en dos pares electrón–hueco; usa $D=0.5$ cm²/s, $T_c=8.65$ K, baño $T_c/2=4.325$ K y parámetros materiales que incluyen espesor de 7 nm. Para esa comparación omite colisiones ee, escape y difusión fonónica. La geometría esférica y la cilíndrica son aproximaciones, no una simulación 3D de las dos superficies de la película. [A20, pp. 15–16 y 20–24; PDF 27–28 y 32–36](https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf).

Las figuras 2.6–2.8 separan la nube electrónica extendida de la energía fonónica localizada. La figura 2.7 y su explicación sitúan aproximadamente el **90% de la energía total** dentro de 3–4 nm en el caso cilíndrico. La transferencia rápida ocurre en una escala aproximada $3\times10^{-5}\tau_0\simeq0.056$ ps; el eje de la figura 2.7 llega a $10^{-4}\tau_0\simeq0.187$ ps, usando el $\tau_0=1.87$ ns cinético de ese capítulo. La figura 2.8 continúa hasta aproximadamente 0.374 ps. Estas conversiones no usan el tiempo GL que comparte el símbolo $\tau_0$.

Si se representa esa **energía localizada** mediante una gaussiana 2D y se iguala únicamente $R_{90}$, resulta

$$
s_{90,\mathrm{A20}}\simeq 1.40\text{–}1.86\ \mathrm{nm}.
$$

Es una traducción útil de un cálculo histórico, no una incertidumbre experimental ni un ajuste gaussiano publicado. Tampoco es directamente el radio fonónico normalizado por su propia energía: el resultado citado usa energía total. Identificarlo con una fuente puramente fonónica requiere comprobar la pequeña fracción electrónica residual al instante escogido.

Hay razones para no dar a ese intervalo más precisión física: se parte de un radio electrónico de 1 nm impuesto, el perfil no es necesariamente gaussiano, se omite ee, no se resuelven superficies ni propagación fonónica, y el baño/energía difieren de 0.9 K y 775/1550 nm. A20 pp. 106–107 muestra además que la aproximación de dos temperaturas con difusión puede conservar demasiada energía en electrones y extenderla demasiado frente a la cinética. La aparente insensibilidad de un antiguo TDGL al tamaño inicial no certifica el nuevo modelo.

En particular, $R_{90}\simeq3$–4 nm es menor que el espesor material de 7 nm. Su traslado a D.30 presupone una representación 2D promediada en profundidad: **no demuestra** que la cascada sea ya uniforme entre ambas superficies. El perfil cilíndrico y la esfera radial de A20 no resuelven esa transición con interfaces. La deposición a distinta profundidad y la pérdida al sustrato permanecen parte de la incertidumbre de transferencia.

Simon et al. aporta una condición **espectral** de burbuja fonónica, pero supone energía espacialmente uniforme en $\pi\xi_c^2d$ y omite difusión en la cinética inicial. Para su TDGL usa después un radio isotrópico $\sqrt{Dt}$, indicando la necesidad de tratar la difusión con más rigor. Por tanto, ese artículo no proporciona una medición adicional de $s$ ni resuelve la elección gaussiana. [SI25, sección III y apéndices sobre condición inicial/TDGL](https://arxiv.org/html/2501.13791v3).

## 5. El reloj cambia la anchura, pero no siempre con D electrónico

Para **difusión libre de una densidad conservada** con coeficiente constante, una gaussiana conserva su forma y

$$
s^2(t)=s^2(t_0)+2D(t-t_0).
$$

Con $D=50$ nm²/ps, una fuente puntual puramente electrónica da $s=3.16$ nm a 0.1 ps, 10 nm a 1 ps y 20 nm a 4 ps. Son conversiones del núcleo de difusión, no predicciones de la cascada ni argumentos para escoger 5–20 nm. Aplicar la misma ley a energía inmovilizada en fonones contradice precisamente el mecanismo que A20 identifica.

Una forma útil de ver la limitación es integrar por partes las ecuaciones de energía de dos subsistemas, en un plano infinito, sin pérdidas ni calentamiento externo. Si sólo los electrones difunden con coeficiente $D$ y las reacciones locales intercambian energía sin transportarla,

$$
\frac{d}{dt}\langle r^2\rangle_{\rm total}
 =4D\,\frac{E_e(t)}{E_{\rm total}},\qquad
s_{\rm mom}^2(t)-s_{\rm mom}^2(t_0)
 =2D\int_{t_0}^{t}\frac{E_e(t')}{E_{\rm total}}\,dt'.
$$

Aquí $s_{\rm mom}^2=\langle r^2\rangle_{\rm total}/2$ es una anchura por segundo momento, incluso si el perfil no es gaussiano. Cuando la energía está mayormente en fonones, el ensanchamiento total es menor que $2D(t-t_0)$. Esta identidad es un diagnóstico bajo esos supuestos: DOS espacial, transporte espectral, propagación fonónica, escape, Joule y bordes requieren sus términos adicionales; no se incorpora como sustituto de las ecuaciones del proyecto.

## 6. Recomendación para el próximo ensayo

**No hay todavía un intervalo físico validado de $s$ para la muestra de Korzh.** El primer intervalo estrecho con una procedencia verificable es **$s\simeq1.4$–$1.9$ nm**, exclusivamente como *escenario puente A20 equivalente en $R_{90}$*, con depósito después de la formación temprana de la burbuja, aproximadamente 0.06–0.19 ps en aquel cálculo. No se adopta automáticamente: primero se debe comprobar que trasladar ese perfil total a la fuente fonónica conserva el reparto requerido y que el modelo continuo resuelve una escala subcoherencia sin confundirla con una supresión impuesta de $\Delta$.

Para el objetivo 80 nm/775–1550 nm/0.9 K, falta establecer el perfil energético **al instante de transferencia elegido**: como mínimo $E_e/E_{\rm in}$, $E_{\rm ph}/E_{\rm in}$, $R_{90}$ o segundo momento, y energía que ya escapó. La elección inicial informada es probar la validez de ese escenario localizado; si falla o si $t_0$ es posterior, hay que obtener una anchura compatible con la cascada y el nuevo instante. No basta inflar $s$ por conveniencia de malla.

Un depósito gaussiano bien resuelto puede después producir un hotbelt por transporte y realimentación; poner desde el inicio una nube que ocupa el ancho presupone parte del fenómeno que se quiere observar. En impactos cercanos a un borde, truncar/renormalizar preserva la energía **declarada** pero no calcula la energía perdida antes de $t_0$. La comparación de latencias debe conservar separado ese tiempo omitido del tiempo que empieza a integrar el modelo.

La decisión que puede plantearse ahora es **usar esta referencia compacta condicionada para evaluar sensibilidad y formación del hotbelt, o caracterizar primero la cascada que debe alimentar D.30**. La primera opción permite avanzar con una hipótesis trazable; la segunda reduce la incertidumbre inicial y exige trabajo adicional. Los controles térmicos equivalentes de 4.5/9 nm responden a otra condición inicial histórica: no se unen con 1.4–1.9 nm para fabricar un rango físico continuo.

## Decisión recibida tras revisar estas fuentes

Se eligió **caracterizar primero la cascada y dejar el ancho abierto**. Por tanto,
1,4–1,9 nm sólo permanece como traducción histórica de R90; no es un escenario
seleccionado ni un rango de barrido. La [caracterización requerida](cascade_characterization.md)
especifica el reparto, la extensión y las pérdidas que deben preceder a la elección.
