# Preparación DC de una sección interior de la nanocinta

Actualización de alcance, 24 de septiembre de 2026. Incorpora las observaciones
posteriores al cierre de etapa 4. La siguiente campaña prepara un equilibrio
sin fotón. No adopta una fuente AC ni una inductancia circuital dependiente del
hotspot. El registro trazable está en [source_decisions.json](source_decisions.json).

## Dispositivo y región que representamos

La geometría descrita por Korzh contiene un tramo activo de 5 µm y un inductor
serie de 1.5 mm de longitud y 1 µm de ancho. Este último es suficientemente
ancho para no actuar como detector de un fotón. Sus 96 nH son una **estimación**
con 64 pH/cuadro de películas similares. Se añadió para evitar el enclavamiento
resistivo. La longitud activa corta reduce el jitter geométrico longitudinal.
El suplemento interpreta las diferencias de subida con inductancia fija.
[Korzh et al., métodos, figura 2 y nota suplementaria 1; PDF 10, 23–24 y 28](https://eprints.lancs.ac.uk/id/eprint/140252/3/Binder1.pdf).

Allmaras explica el compromiso: retardar suficientemente el retorno de corriente
para evitar el enclavamiento, conservando una subida rápida del pulso. También
distingue la posición longitudinal, que altera el tiempo de propagación hasta la
lectura, de la posición transversal, que modifica la detección local.
[Allmaras, tesis de 2020, pp. 60–63, PDF 72–75](https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf).

La decisión de modelado es representar una **ventana interior recta**, de ancho
80 nm, alrededor del futuro punto de absorción. Sus dos cortes longitudinales
son cortes de una misma nanocinta superconductora. No representan contactos
metálicos, tapers ni extremos físicos. La coordenada longitudinal de absorción se
mantiene fija. La futura variación de la coordenada transversal permitirá
estudiar el jitter geométrico del ancho; no se agrega una distribución artificial
de tiempos de propagación a lo largo de los 5 µm.

Esto no presupone que la iluminación experimental seleccione exclusivamente esa
ventana. Define la familia de eventos interiores que se busca representar.
El alcance tampoco incluye auto-campo, rugosidad, defectos ni fluctuaciones
materiales adicionales. En esa geometría ideal, el equilibrio esperado tiene
amplitud uniforme y corriente longitudinal uniforme. Un perfil de recuperación
del condensado cerca de los cortes sería un efecto de borde del modelo a
diagnosticar, no una predicción sobre contactos presentes en el dispositivo.

Hay además una comprobación geométrica de la aproximación sin auto-campo: la
inductancia de hoja citada implica una longitud de Pearl
$\Lambda=2L_\square/\mu_0\simeq102$ µm. El ancho de 80 nm es unas 1270 veces
menor. Esta conversión usa el modelo de película delgada; no reemplaza una
medición independiente de penetración magnética.

## El piso inicial incluye el efecto de la corriente

Un equilibrio con corriente puede ser homogéneo en sus magnitudes físicas aunque
la fase varíe. En gauge con potencial vector cero se escribe

$$
\Delta(x)=\Delta_q e^{iqx},\qquad
F_n(x)=F_{nq}e^{iqx},\qquad G_n(x)=G_{nq},
$$
$$
\Gamma_q F_{nq}G_{nq}=\Delta_qG_{nq}-\epsilon_nF_{nq},\qquad
F_{nq}^2+G_{nq}^2=1,\qquad \Gamma_q=\frac{\hbar Dq^2}{2}.
$$

Aquí $q$ es el gradiente de fase, $\epsilon_n=(2n+1)\pi k_BT$ la energía
de Matsubara y $\Gamma_q$ la energía de desemparejamiento por corriente. La
autoconsistencia determina $\Delta_q$ y la corriente selecciona $q$ en la rama
estable conectada con corriente cero. Es la solución uniforme de Usadel de
[Clem y Kogan, sección II.A, ecuaciones (3)–(8) y (11), PDF 2–3](https://arxiv.org/pdf/1207.6421).

Para construir la referencia numérica se resuelven juntos el espectro, el gap y
la corriente con la misma normalización, corte espectral y regularización usados
en la malla. En magnitud, usando $R_\square$ como resistencia de hoja y $w$ como
ancho, la relación de corriente es

$$
|I_q|=\frac{2\pi k_BT\,w}{eR_\square}|q|
\sum_{n\ge0}F_{nq}^{2}.
$$

La orientación del signo se adapta al convenio de puerto ya existente. No se
impone una caída de fase independiente de la longitud al comparar dominios:
se mantiene la misma corriente y se usa $qL$ para cada longitud $L$.

La preparación de borde e interior debe coincidir en:

- $|\Delta|=\Delta_q$, espectro con $\Gamma_q$ y fase coherente con $q$;
- distribución electrónica térmica, $h_L(E)=\tanh[E/(2k_BT)]$ y $h_T(E)=0$;
- población fonónica de Bose a la temperatura de baño;
- potencial eléctrico constante y ausencia de caída de tensión DC;
- corriente normal nula y corriente total igual a la superconductora inicial.

En los bordes laterales se anulan los flujos normales. En los cortes
longitudinales se usa la continuación del equilibrio portador de corriente.
Poner allí un espectro BCS de corriente cero, aun fijando una fase distinta en
cada extremo, no representa ese exterior. Las funciones de Green en el borde
deben incluir también la fase y el desemparejamiento correspondientes.

La solución discreta se relaja y se compara con la referencia uniforme; copiar
una fórmula continua a los nodos no prueba por sí solo estacionariedad discreta.
Se registran amplitud, corriente de sección, potencial, distribuciones y su
deriva temporal. Las diferencias de malla se miden respecto del mismo piso
inicial. Se comparan longitudes y mallas antes de decidir cuánto dominio basta.
Una malla más fina debe reducir o estabilizar el desvío, no sostener un gradiente
longitudinal macroscópico impuesto por los terminales.

Estas condiciones sirven para el ensayo **anterior al fotón**. Durante un futuro
pulso, no se congelará artificialmente la fase de ambos extremos mientras la
corriente circuital cambia. La adaptación de los puertos al circuito, y la
comprobación de que la perturbación no alcanza los cortes durante la ventana
observada, siguen siendo parte del transiente fotónico.

## Inductancia exterior fija y circuito de la memoria

Se conserva la topología de tres variables de la
[adenda CM](../../../modelo_v0_4/actualizaciones/circuito_memoria_20260922.md):

$$
\begin{aligned}
L_b\dot I_b&=V_b-R_bI_b-v_c-R_L(I_b-I_s),\\
L_{k,\rm ext}\dot I_s&=v_c+R_L(I_b-I_s)-V_{\rm dev},\\
C\dot v_c&=I_b-I_s,\qquad V_{\rm out}=R_L(I_b-I_s).
\end{aligned}
$$

La fuente es constante: $V_b=R_bI_{\rm DC}$. El equilibrio circuital es
$(I_b,I_s,v_c)=(I_{\rm DC},I_{\rm DC},0)$, con
$V_{\rm dev}=V_{\rm out}=0$. Se mantienen $R_b=10$ kΩ, $L_b=1$ µH,
$R_L=50$ Ω y $C=100$ pF como parámetros del circuito CM. Adoptar geometría e
inductancia exterior motivadas por Korzh no convierte esos otros componentes
en una caracterización medida de su electrónica.

**Se adopta una inductancia exterior fija durante cada transiente.** La evidencia
geométrica respalda esa reducción; no demuestra que toda respuesta inductiva
microscópica del segmento excitado sea exactamente constante. El voltaje y la
energía del segmento resuelto siguen procediendo de sus ecuaciones. No se añade
un elemento $L_k(t)$ ni una fuente $I_s\dot L_k$ para representar nuevamente esa
misma respuesta.

En el mismo film y a la misma corriente serie, el tramo de 1 µm lleva sólo
$80/1000=0.08$ veces la densidad de corriente del tramo activo. Por eso el
exterior ancho está mucho más lejos del desemparejamiento por corriente, además
de no recibir el depósito fotónico localizado. Es un argumento para conservar
su inductancia de referencia; no exige eliminar la depresión del condensado ni
la dinámica local del segmento estrecho.

La partición recomendada para una longitud resuelta $L$ es

$$
L_{k,\rm ext}^{\rm ref}
=96\,\mathrm{nH}
+(5\,\mathrm{\mu m}-L)\,\mathcal L_{\rm bulk}^{\rm ref},\qquad
\mathcal L_{\rm bulk}^{\rm ref}
=\frac{\hbar}{2e}\left.\frac{dq}{dI}\right|_{I_{\rm DC}}.
$$

$\mathcal L_{\rm bulk}^{\rm ref}$ es la inductancia diferencial por longitud
del equilibrio uniforme elegido. Se evalúa **una sola vez**, usando la misma
rama y el mismo material del dominio espacial; después se guarda como constante.
No es una derivada recalculada durante el pulso. La fórmula excluye los tapers y
conexiones, cuya contribución no está cuantificada, y requiere $0<L<5$ µm.
Si se incluyen continuaciones 1D explícitas, sus longitudes forman parte de $L$.

Como comparación de escala, la aproximación de hoja publicada daría
$64$ pH/cuadro $\times(5\,\mathrm{\mu m}/80\,\mathrm{nm})=4$ nH para todo
el tramo activo: unos 100 nH sumando el inductor añadido, antes de tapers y
conexiones. Para ventanas de 120–480 nm, sólo 0.096–0.384 nH de esa estimación
corresponderían al material resuelto. Son cálculos geométricos, **no nuevas
mediciones ni cotas al cambio de inductancia cerca de un núcleo normal**.

El valor de hoja experimental procede de películas similares, mientras que
$D=0.5$ cm²/s y $R_\square=608$ Ω son la referencia de ajuste seleccionada para
el modelo. Por eso no se impone que ambos produzcan exactamente el mismo
$\mathcal L_{\rm bulk}$. El informe de la campaña debe mostrar las dos
estimaciones y su procedencia, sin ajustar la corriente para hacerlas coincidir.

Los 96 nH son un componente **exterior añadido**: no se les resta la inductancia
de la ventana. Tampoco se suman 100 nH y la ventana otra vez. Si en el futuro se
dispone de un total medido consistente, podrá emplearse la contabilidad
alternativa $L_{\rm ext}=L_{\rm total}-L_{\rm res}^{\rm ref}$, una sola vez.
Este escenario reemplaza los 10 nH históricos para la nueva campaña motivada
por el dispositivo; no reinterpreta los resultados previos.

## Qué podrá quedar admitido antes del fotón

El pase de esta campaña acredita un piso DC homogéneo, sin señal eléctrica
espuria apreciable, con bordes de continuación y partición inductiva explícitos.
No se programa una nueva familia de excitaciones sinusoidales. Los resultados
AC de etapa 4 quedan como diagnóstico histórico de acoplamiento.

La ausencia de potencia en equilibrio no verifica cómo se reparte el calor de
una perturbación finita. Por tanto, el desbalance de potencia de etapa 4 no se
declara resuelto porque este ensayo no lo excite. La futura inyección necesita
el balance de energía no lineal, las tasas materiales admitidas y la transferencia
fotónica ya identificados como pendientes. Ninguno se sustituye por ruido
numérico, relajación de tolerancias o por fijar la inductancia exterior.

La siguiente entrega debe mostrar perfiles $|\Delta(x,y)|/\Delta_q$, corriente
longitudinal, espectro bulk y de borde, deriva sin fotón y $V_{\rm out}$ referido
a su línea basal. Cada gráfico distinguirá magnitud física, referencia restada,
unidades y tamaño del dominio. El objetivo inmediato es eliminar artefactos de
contacto de la preparación, manteniendo el mismo sistema DC que luego recibirá
un único fotón.
