# Circuito de la memoria para el modelo previsto

**Actualización normativa documental, 22 de septiembre de 2026. No implementada.**

El modelo previsto adopta la red de polarización y lectura utilizada en los
resultados de la memoria: fuente de tensión, resistencia e inductancia de
polarización, rama del detector con inductancia serie y lectura acoplada mediante
un condensador. Esta adenda sustituye, para el plan vigente, C.32–C.34 y D.28–D.29,
el balance global D.33 y la parte circuital de la inicialización D.1.9/D.30.
Amplía el estado `I(t)` de D.1 a tres variables circuitales. El resto de las
ecuaciones de A–D conserva su alcance; el transporte y el calor del detector no
se reemplazan por una resistencia prescrita.
En las condiciones de rama y de borde D.26–D.27, la corriente del detector
antes llamada $I$ pasa a ser $I_s$, nunca $I_b$.

Los documentos C/D 0.4 y sus pruebas se conservan como historia. El circuito de
corriente ideal y carga directamente paralela que verificaron corresponde al
esquema conceptual de la memoria, no a su red de resultados. Su pase numérico no
valida automáticamente las ecuaciones de esta adenda. Esta actualización no
modifica producción, no cierra la etapa 2 ni inicia la etapa 3.

## 1. Fuente cotejada y parámetros de referencia

Se utiliza *Multiscale Modeling of the Transient Response of Superconducting
Nanowire Single-Photon Detectors*, J. A. Díaz Monge, copia `memoria_02.pdf`,
173 páginas, procedente de `/home/jdiaz/memoria/main/memoria_02.pdf`.
Su SHA-256 es
`50a75f1bd84f06f32820dbf9477c0fdcaf809fb502485e84249ebb46b286c1e8`,
coincidente con el [manifiesto 0.4](../verificaciones/manifiesto_v0_4.json).

| Referencia de la memoria | Contenido usado |
|:--|:--|
| §2.3.8, figura 2.7 y ec. (2.36), p. 33 / PDF 64 | Circuito conceptual de corriente ideal; antecedente, no circuito seleccionado |
| §4.3.3, ec. (4.15), p. 77 / PDF 108 | Voltaje central de 100 nm y voltaje entre terminales del segmento completo de 360 nm |
| §4.4.1, ec. (4.16), p. 83 / PDF 114 | Tres ecuaciones diferenciales de la red utilizada en los resultados |
| Figura 4.15 y ec. (4.17), p. 84 / PDF 115 | Conexiones, intercambio con el solver e inicialización estacionaria ideal |
| Ec. (4.18), p. 85 / PDF 116; tabla 4.3, p. 54 / PDF 85 | Parámetros circuitales e identificación del estado histórico |
| §3.2.2, p. 40 / PDF 71; ec. (4.1), p. 51 / PDF 82 | Región resuelta, longitud exterior e interpretación de la inductancia total |

Los parámetros históricos son:

| Parámetro | Valor | Interpretación |
|:--|--:|:--|
| $R_{\rm bias}$ | $10\,\mathrm{k\Omega}$ | Resistencia de polarización |
| $L_{\rm bias}$ | $1\,\mathrm{\mu H}$ | Inductancia de polarización |
| $R_{\rm load}$ | $50\,\mathrm{\Omega}$ | Carga de lectura |
| $C_{\rm couple}$ | $100\,\mathrm{pF}$ | Condensador serie de lectura |
| $V_{\rm bias}$ | $0.300\,\mathrm V$ | Tensión de la fuente en los resultados históricos |
| $L_k^{\rm total,mem}$ | $10\,\mathrm{nH}$ | Inductancia **total** de referencia; véase su partición en §3 |

Son entradas del escenario de la memoria, no una nueva caracterización
experimental ni parámetros que se ajusten después de mirar el pulso. La memoria
asocia los 10 nH a una longitud uniforme equivalente de aproximadamente
22.9 µm mediante su estimación de London. Esa longitud no es la longitud del
dominio numérico.

## 2. Variables, conexiones y signo pasivo

Se toma la tierra de la red como referencia eléctrica. $V_d$ es el voltaje del
nodo donde se separan las ramas, y $V_{\rm out}$ el voltaje del nodo superior de
la carga. Se define $v_c=V_d-V_{\rm out}$, positivo en la placa del condensador
del lado del detector. El estado circuital es

$$
\mathbf y_c=(I_b,I_s,v_c).
\tag{CM.1}
$$

$I_b$ fluye desde la fuente por $R_{\rm bias}$ y $L_{\rm bias}$ hasta el nodo
$d$. Desde allí, $I_s$ entra en la rama del detector e $I_{\rm RF}$ en el
condensador y la carga. **$I_s$ es la corriente total del ramal del nanohilo**:
en la película contiene corriente superconductora y normal. No es la integral
de $\mathbf j_s$ por sí sola, y $I_b(t)$ ya no es una corriente impuesta constante.

$$
I_{\rm RF}=I_b-I_s,\qquad
V_{\rm out}=R_{\rm load}I_{\rm RF},\qquad
V_d=v_c+V_{\rm out}.
\tag{CM.2}
$$

$V_{\rm dev}$ es la caída pasiva de la región mesoscópica que recibe $I_s$:
$I_sV_{\rm dev}$ es el trabajo eléctrico que entra en esa región. Para conservar
la orientación izquierda–derecha declarada en C/D, los planos de puerto son
$L,R$ y

$$
V_{\rm dev}=\phi_L-\phi_R\quad(\mathbf A=0),\qquad
\int_{\Gamma_L}\mathbf j_{\rm tot}\cdot\widehat{\mathbf n}\,dS=-I_s,
\qquad
\int_{\Gamma_R}\mathbf j_{\rm tot}\cdot\widehat{\mathbf n}\,dS=I_s.
\tag{CM.3}
$$

Aquí los potenciales de puerto son equipotenciales, como en la reducción 1D
de las continuaciones. Si una sección no es equipotencial, su trabajo se calcula
con $-\int\phi\,\mathbf j\cdot\widehat{\mathbf n}\,dS$; promediar $\phi$ sin
ponderar por la corriente no garantiza la identidad de potencia. Se mantienen
la aproximación sin apantallamiento de C/D y el espesor en las medidas de cara.

Las tres ecuaciones de la memoria quedan, con la misma topología y con la
partición energética de §3,

$$
\begin{aligned}
L_{\rm bias}\dot I_b
 &=V_{\rm bias}-R_{\rm bias}I_b-v_c-R_{\rm load}(I_b-I_s),\\
L_{k,\rm ext}\dot I_s
 &=v_c+R_{\rm load}(I_b-I_s)-V_{\rm dev},\\
C_{\rm couple}\dot v_c&=I_b-I_s.
\end{aligned}
\tag{CM.4}
$$

Los elementos pasivos de esta red son constantes y positivos en el escenario
seleccionado. No se añade una segunda inductancia en paralelo, ni se supone
$V_{\rm out}=V_d$: el condensador sostiene la diferencia $v_c$. Tampoco se
añade una función de transferencia de amplificador sin datos que la definan.

**Cotejo del dibujo.** La figura 4.15 de la memoria dibuja el signo positivo de
la fuente $V_{\rm bias}$ junto a tierra, mientras que (4.16), (4.17) y el código
usan $V_{\rm bias}>0$ como tensión del terminal que alimenta $R_{\rm bias}$.
La convención normativa es la de las ecuaciones: terminal superior de la fuente
positivo respecto de tierra. Se conserva la conexión de la figura, sin copiar
esa discrepancia de polaridad.

### Adaptación explícita al backend histórico

El [solver circuital actual](../../../pysnspd/circuit/readout.py#L318) contiene
las tres ecuaciones (4.16). Su
[observable central](../../../pysnspd/circuit/readout.py#L373) es
$V_{\rm TDGL}^{\rm hist}=\phi_R-\phi_L$. No basta comparar esta fórmula con
CM.3 para decidir que el voltaje histórico tiene un signo incorrecto: hay que
seguir también la corriente.

El [ensamblador de terminales](../../../pysnspd/solver/stationary.py#L270)
entrega `{left: -I_backend, right: +I_backend}`. El
[adaptador Neumann](../../../pysnspd/mesh/device.py#L98) devuelve
$g_\Gamma=-\sum_{\Gamma'\ne\Gamma}I_{\Gamma'}/(\ell'_\Gamma\sigma_n dV_0)$;
al sumar cero las corrientes, resulta $g_\Gamma=I_\Gamma/(\ell'_\Gamma\sigma_n dV_0)$.
El [operador](../../../pysnspd/gtdgl/tdgl_operators.py#L111) introduce $g$ como
derivada normal en
$L_\mu\mu=D\mathbf j_s-Bg$, y el
[núcleo](../../../pysnspd/solver/core.py#L415) calcula
$\mathbf j_n=-G\mu$.
Por tanto, $I_{\rm backend}>0$ corresponde en esa implementación a corriente
física de derecha a izquierda. El comentario del adaptador que identifica su
argumento con una corriente exterior no coincide con este signo algebraico;
la selección automática de signo de un diagnóstico de divergencia no redefine
el campo físico.

La conversión al convenio izquierda–derecha de CM.3 es entonces

$$
I_{\rm backend}=-I_s,\qquad
V_{\rm dev}=-V_{\rm TDGL}^{\rm hist},
\tag{CM.5}
$$

si ambos voltajes se miden entre **los mismos planos**. De forma equivalente,
con la orientación histórica derecha–izquierda se conservan ambos valores
brutos, $I_{\rm backend}=I_s$ y $V_{\rm dev}=V_{\rm TDGL}^{\rm hist}$, y se
invierten los signos de los flujos de CM.3. No se mezclan ambas orientaciones,
no se toma valor absoluto y no se invierte sólo una variable. CM.5 documenta
el adaptador que requeriría el acoplamiento futuro; no modifica el backend.

## 3. Región resuelta e inductancia sin doble conteo

La memoria exporta el voltaje central de 100 nm y conserva por separado el
voltaje de los terminales del segmento completo de 360 nm. Su texto llama
$L_k$ a una inductancia total que incluye longitud no simulada. En el nuevo
modelo, la energía común ya almacena superflujo en el dominio 2D y sus
continuaciones 1D. No se puede sumar de nuevo su inductancia mediante un
elemento denominado «total».

Para el acoplamiento previsto se fijan los puertos en los extremos del dominio
energético resuelto de D.8. $V_{\rm dev}$ abarca ese mismo dominio; sustituye
a la entrada mesoscópica de la memoria, conservando la red exterior CM.4. El
voltaje central de 100 nm permanece como observable separado para comparación
histórica. Exportar sólo ese voltaje y afirmar a la vez que el circuito entrega
todo el trabajo de la región 2D más continuaciones omitiría el trabajo de las
zonas excluidas del puerto.

El elemento de CM.4 es exclusivamente $L_{k,\rm ext}$, correspondiente a la
parte no incluida en $U_{\rm res}$. Antes de un ensayo acoplado se registran
geometría, planos y partición de flujo del estado de referencia. Si se utiliza
la referencia total de 10 nH, su partición diferencial en ese estado es

$$
L_k^{\rm total,mem}
 =L_{k,\rm res}^{\rm diff}(I_{\rm ref},T_b)+L_{k,\rm ext},
\qquad L_{k,\rm ext}>0.
\tag{CM.6}
$$

La respuesta resuelta se obtiene de la rama y del funcional espacial admitidos;
la estimación geométrica de London de la memoria puede documentar una
aproximación, pero no reemplaza esa identificación sin declarar su error.
La topología no cambia por esta partición. **10 nH es la referencia total, no
un segundo almacenamiento de 10 nH que se añade íntegro a los campos.**

$L_{k,\rm ext}$ se fija antes del transitorio. No se resta en cada instante una
«inductancia del hotspot» para forzar que la suma permanezca en 10 nH: ello haría
responder artificialmente la parte exterior al estado local. Si la partición
no es conocida o no da un valor positivo admisible, falta una entrada del
ensayo y no se atribuye todavía una predicción circuital absoluta al modelo.
La comparación con la red histórica de 10 nH puede hacerse como prueba de
circuito con una entrada mesoscópica prescrita, identificada como tal.

## 4. Energía y balance que sustituyen C.34/D.29/D.33

Para los elementos constantes seleccionados,

$$
U_c=\frac12 L_{\rm bias}I_b^2+
\frac12 L_{k,\rm ext}I_s^2+\frac12 C_{\rm couple}v_c^2.
\tag{CM.7}
$$

Multiplicar las tres ecuaciones CM.4 por $I_b$, $I_s$ y $v_c$ y sumarlas
cancela los términos de intercambio $v_c(I_b-I_s)$. El término de la carga
queda $-V_{\rm out}(I_b-I_s)=-R_{\rm load}I_{\rm RF}^2$. Por ello,

$$
\dot U_c+I_sV_{\rm dev}
 +R_{\rm bias}I_b^2+R_{\rm load}I_{\rm RF}^2
 =V_{\rm bias}I_b.
\tag{CM.8}
$$

La fuente aporta $V_{\rm bias}I_b$, no $I_bV_{\rm out}$. Los dos resistores
disipan fuera del dominio de la película; su calor no entra otra vez en la
fuente Joule de B. El almacenamiento capacitivo y el de $L_{\rm bias}$ no son
pérdidas ni pueden omitirse del balance transitorio.

Con $U_{\rm res}=\int_{\mathcal V}(e_\delta+u_{\rm ph})\,dV$, D.32 y los
puertos pasivos proporcionan el balance global, para $t>t_0$,

$$
\begin{aligned}
\frac{d}{dt}(U_{\rm res}+U_c)
={}&V_{\rm bias}I_b-R_{\rm bias}I_b^2-R_{\rm load}I_{\rm RF}^2\\
&-\int_{\mathcal V}P_{\rm esc}\,dV
-\oint_{\partial\mathcal V}
(\mathbf Q_e-\mathbf J_{\rm tr})\cdot\widehat{\mathbf n}\,dS.
\end{aligned}
\tag{CM.9}
$$

Aquí $\int\mathbf j_{\rm tot}\cdot\mathbf E\,dV=I_sV_{\rm dev}$ cancela
el trabajo que sale del circuito. Se conserva el trabajo de los reservorios
de amplitud y fase en $\mathbf J_{\rm tr}$. La disipación $Q_\Delta$ sigue
cancelándose internamente entre condensado y ocupaciones; no se añade como
segunda potencia externa. En $t_0$ se contabiliza sólo el salto fotónico
$E_{\rm in}$ de D.30.

## 5. Estado inicial, interfaz y salidas

Un estado estacionario de la red satisface

$$
I_b^*=I_s^*=I_{\rm SS},\qquad
I_{\rm RF}^*=V_{\rm out}^*=0,\qquad
v_c^*=V_{\rm dev}^*,\qquad
V_{\rm bias}=R_{\rm bias}I_{\rm SS}+V_{\rm dev}^*.
\tag{CM.10}
$$

Se resuelve junto con el equilibrio de los campos, las poblaciones y los
reservorios de la misma energía. Para el fondo superconductivo sin caída
resuelta de D.1.9, $V_{\rm dev}^*=v_c^*=0$ e
$I_{\rm SS}=V_{\rm bias}/R_{\rm bias}$, sujeto a que exista la rama estable
admitida. No se pueden imponer de forma independiente una corriente inicial,
un $V_{\rm bias}$ incompatible y afirmar que la inicialización es estacionaria.

Las cifras históricas de la memoria —corriente solicitada 30 µA, estado guardado
29.3427 µA y voltaje central 50.9939 µV— describen una inicialización de
estacionariedad parcial. No se insertan como equilibrio exacto del nuevo
funcional. Si se estudia una trayectoria histórica se conserva su estado
circuital guardado y se identifica su línea de base.

En el depósito fotónico se mantiene la actualización de fonones de D.30 y
se sustituye su condición circuital por la continuidad de los tres estados:

$$
I_b(t_0^+)=I_b(t_0^-),\qquad
I_s(t_0^+)=I_s(t_0^-),\qquad
v_c(t_0^+)=v_c(t_0^-).
\tag{CM.11}
$$

No se reinicia $I_s$ a una corriente de fuente ideal. Tampoco se deposita el
fotón en el condensador o en las inductancias. El circuito recibe
$V_{\rm dev}(t)$ y devuelve $I_s(t)$ a los bordes; ese intercambio y la
evolución simultánea del espectro deben cumplir CM.9.

Las salidas anteriores a cualquier amplificación son
$I_b,I_s,I_{\rm RF},v_c,V_d,V_{\rm dev},V_{\rm out}$, el voltaje central
de comparación cuando se defina, y el balance CM.9. El pulso de lectura es
$V_{\rm out}$ **después del condensador**, no $V_d$ ni el voltaje de la
película. Cada latencia debe declarar el observable, la línea de base y la
regla de cruce; la actualización del circuito no valida por sí sola una latencia.

La implementación histórica intercambia voltaje y corriente por intervalos y
avanza el circuito mediante RK2. Esa elección numérica no se impone al nuevo
integrador ni acredita su orden acoplado. Antes de admitirlo deberán comprobarse
el punto fijo CM.10, el signo pasivo, la contabilidad CM.8–CM.9, la resolución
temporal del intercambio y la partición de CM.6. Los cálculos costosos siguen
la política de ejecución manual registrada en `GEMINGA_COMMANDS.md`.

**Límite de esta entrega:** queda definido el circuito físico que se pretende
acoplar, con los puertos y la contabilidad requeridos. No se han cambiado
fuentes físicas, corrido un transitorio del detector ni concedido admisión numérica,
material o de producción.
