# Preparación fotónica: qué entra en el modelo y qué falta conocer

Investigación del 23 de septiembre de 2026. Referencia prioritaria: NbN de 80 nm de Korzh, comparación 775/1550 nm y formación de un hotbelt. Se aceptó estudiar primero **desde la transferencia**, conservando explícito el retardo óptico desconocido; también usar el escenario material ajustado de Korzh como referencia. Estas decisiones no fijan la retención ni la anchura fotónica. No se ejecutaron transientes, no se cambiaron parámetros y no se admitió un dominio físico. Procedencia y clasificación de cada cifra: [photon_sources.json](photon_sources.json).

La revisión solicitada de los artículos de Allmaras y Vodolazov de 2019 se incorpora en §3.1. Los intervalos anteriormente sugeridos $\chi\in[1/3,1]$ y $s_\gamma\in[5,20]$ nm **no fueron aceptados y no se adoptan**. La forma gaussiana es una elección del modelo solicitada por el usuario; su anchura necesita una justificación propia.

## 1. Una entrega de energía entre dos descripciones

La **cascada de conversión energética** comienza cuando el fotón absorbido crea excitaciones electrónicas muy energéticas. Estas generan otras excitaciones y fonones. El modelo 0.4 empieza después de una parte de esa cascada: recibe un estado preparado y continúa su evolución. Es análogo a recibir un recipiente con una cantidad de energía conocida: conocer esa cantidad no determina cómo se repartió dentro ni cuánto tardó en llegar.

En B.38–B.39 y D.30 se distingue:

| Dato de transferencia | Significado y unidad | Qué no determina |
|---|---|---|
| $E_\gamma=hc/\lambda$ | Energía del fotón absorbido, J o eV | Probabilidad de absorción o detección |
| $E_{\rm in}$ | Energía excedente retenida en los sectores resueltos en $t_0$ | Reparto espectral y espacial |
| $\chi=E_{\rm in}/E_\gamma$ | Fracción retenida de un escenario determinista, adimensional | Varianza entre eventos |
| $s_\gamma$ | Desviación estándar por coordenada del perfil gaussiano no truncado, m | Diámetro del foco óptico o longitud de hotbelt |
| $\mathbf r_0$ | Centro del perfil efectivo en la película, m | Distribución de posiciones de absorción |
| $t_0-t_{\rm abs}$ | Duración física de la etapa omitida, s | Paso temporal numérico o tiempo de termalización posterior |

Para una preparación sin otro aporte impulsivo,

$$E_{\rm in}=E_\gamma-E_{\rm perdido}^{<t_0},\qquad 0\leq E_{\rm in}\leq E_\gamma.$$

El trabajo eléctrico posterior se contabiliza por separado. El escape fonónico que la cinética resuelve **después** de $t_0$ tampoco se descuenta en $E_{\rm perdido}^{<t_0}$. Esto evita contar una misma pérdida dos veces. Si toda la energía del fotón sigue dentro de electrones y fonones al transferir el estado, $E_{\rm in}=E_\gamma$ aunque su reparto sea complejo. [Modelo 0.4, B.8](../../../modelo_v0_4/B_cinetica_energia_y_temperaturas_v0_4.md), [D.30](../../../modelo_v0_4/D_sintesis_plan_y_verificaciones_v0_4.md).

Para las longitudes nominales de interés, $E_\gamma(775\,\mathrm{nm})=1.599796\,\mathrm{eV}$ y $E_\gamma(1550\,\mathrm{nm})=0.799898\,\mathrm{eV}$. Son conversiones de $hc/\lambda$, no medidas del ancho espectral del láser. Una **hipótesis común** $\chi_{775}=\chi_{1550}$ produce una razón de energías retenidas igual a dos; no garantiza una razón simple entre retardos.

## 2. Estado inicial que ya especifica 0.4

La *burbuja fonónica* significa un exceso localizado de ocupación de modos vibracionales; no es una cavidad física. El candidato conserva electrones inicialmente en Fermi–Dirac al baño, el condensado continuo y los estados circuitales continuos. Deposita el exceso sólo en fonones:

$$
u_\gamma(\mathbf r)=\frac{E_{\rm in}}{d_f}w(\mathbf r),\qquad
w(\mathbf r)=\frac{\exp[-|\mathbf r-\mathbf r_0|^2/(2s_\gamma^2)]}
{\int_{\mathcal W}\exp[-|\mathbf r'-\mathbf r_0|^2/(2s_\gamma^2)]\,d^2r'},
$$

$$
\delta n(\Omega,\mathbf r)=
\frac{u_\gamma(\mathbf r)\,\alpha^2F(\Omega)}{g_{\rm ph}(\Omega)M_1},
\qquad M_1=\int_0^\infty \Omega\alpha^2F(\Omega)\,d\Omega.
$$

Aquí $d_f$ es el espesor, $\mathcal W$ la región de deposición de la película, $\Omega$ la **energía** de un fonón, $g_{\rm ph}$ su densidad de modos por volumen y energía, y $\alpha^2F$ el espectro de acoplamiento electrón–fonón. $n$ es una ocupación por modo. Al multiplicar por $g_{\rm ph}\Omega$ e integrar se obtiene $u_\gamma$; al integrar además el volumen se obtiene $E_{\rm in}$. Esta identidad comprueba la energía depositada, no demuestra que el espectro inicial sea exacto.

Simon et al. motivan un exceso proporcional a $\alpha^2=\alpha^2F/F$ a partir de una excitación electrónica concentrada a alta energía. Lo normalizan a la energía del fotón y dejan inicialmente los electrones en equilibrio. Esta sustitución evita resolver hasta energías ópticas. Su cálculo usa un volumen cilíndrico uniforme: no determina un $s_\gamma$ gaussiano universal. En 0.4 la normalización usa la energía **retenida**, que puede ser menor. [S25, §III y §VII.1, después de (8)](https://arxiv.org/html/2501.13791v3).

Al implementar la preparación se deberá usar el catálogo admitido y sus unidades, con $g_{\rm ph}\geq0$, $\alpha^2F\geq0$, soporte compatible y $0<M_1<\infty$. Donde ambos espectros son cero no se evaluará un cociente $0/0$. El mismo soporte y la misma cuadratura deben producir el normalizador y el balance; la energía omitida por un corte numérico no se reinterpretará como escape físico. El recorte de soporte NbN ya documentado exige recalcular $M_1$, no conservar una amplitud previa. [Cierre del catálogo](../../stage1_closure/material_findings.md).

El archivo vigente `pysnspd/excitation/photon.py` pertenece a la implementación térmica anterior: actualiza una temperatura fonónica mediante una tabla de energía. Su valor nominal de 10 nm es procedencia histórica, no evidencia de que exista ya un inyector cinético D.30 conectado al transiente mixto.

## 3. Qué aportan las fuentes y qué cifras no son transferibles

**Vodolazov 2017.** Compara preparaciones electrónicas, fonónicas y térmicas a energía igual bajo hipótesis específicas; al estudiar la relajación inicial fija el condensado y simplifica la expansión espacial. Usa $V_{\rm init}=\pi\xi^2d$ y también examina una meseta de ocupación fonónica. Es evidencia a favor de investigar una preparación efectiva, no una medida de su radio o duración en el dispositivo de Korzh. Una gaussiana en **energía** de sus ecuaciones iniciales tampoco es una anchura espacial. [V17, §III, (16)–(22)](https://arxiv.org/pdf/1611.06060).

**Allmaras 2020.** En el modelo Debye, la meseta fonónica es aproximadamente constante en ocupación hasta el corte: no es una temperatura. La tesis contrasta esta sustitución con una cascada electrónica y encuentra diferencias tempranas y dependientes de la difusión. Sus distribuciones radiales localizadas se obtienen con hipótesis materiales y geométricas concretas; un radio que contiene cierta fracción de energía no equivale a $s_\gamma$. [A20, §§2.3.1–2.3.4, pp. 19–26; PDF 31–38](https://thesis.caltech.edu/13748/).

La posterior *burbuja modificada* de A20 introduce, en un modelo térmico reducido, una fuente electrónica $E\exp(-t/\tau_{DC})/(V_{HS}\tau_{DC})$. Allí $\tau_{DC}\simeq1.4$ ps; luego se ensaya ralentizarla aproximadamente cuatro veces para aproximar mejor la transferencia cinética. **No es una medición de $t_0$.** Añadir esa fuente sobre la burbuja B.39 depositaría o transferiría de nuevo energía que la nueva cinética ya sigue. La tesis muestra que cambiar la preparación cambia los parámetros ajustados: por ejemplo, en un caso 2D, $\chi$ pasa de 0.92 con inicio térmico a 0.61 con burbuja, sin ajuste satisfactorio completo. Estos números no forman un intervalo admisible para 0.4. [A20, §2.4.6 y pp. 101–105; PDF 113–117](https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf).

**Korzh 2020.** La tabla del modelo 1D de 80 nm contiene longitud de hotbelt 40 nm y fracción retenida 0.667; esta última se ajusta conjuntamente con dispersión energética y relajación electrón–electrón. El suplemento separa $\sigma_F(1550)=92$ meV y una contribución no uniforme de 40 meV. Son parámetros de aquel modelo, no mediciones independientes de transferencia. El propio ajuste gaussiano simple presenta colas energéticas imposibles. [K20, nota suplementaria 3 y tabla 1, PDF 31–33](https://eprints.lancs.ac.uk/id/eprint/140252/3/Binder1.pdf).

### 3.1. Artículos de 2019: energía depositada y valor de referencia

Se cotejaron **A19v2**, preprint de *Intrinsic Timing Jitter and Latency…* asociado a PRApplied **11**, 034062, y **V19v1**, preprint de *Minimal Timing Jitter…* asociado a PRApplied **11**, 014016. Sus versiones son de agosto y julio de 2018, respectivamente. La edición publicada de A19 requiere acceso editorial; los valores siguientes se atribuyen al preprint, sin presumir que su numeración o parámetros coincidan con la revisión publicada.

| Fuente y cierre | Retención utilizada | Convención y localizador |
|---|---|---|
| A19v2, hotbelt simplificado | $\chi=0.43$ | Ajuste de conteo; su justificación incluye pérdidas durante la latencia y reparto electrónico. §II.D, figura 3, PDF 7. |
| A19v2, TDGL generalizada 1D | $\chi=0.16,0.39,0.52$ para $\tau_{ee}(T_c)=0,5,10$ ps | Ajuste PCR1550; hotbelt de 80 nm, baño de 2 K. §III.B–C, PDF 10–11. |
| A20, revisión 1D de la tesis | $\chi=0.37,0.65,0.79$ para $\tau_{ee}(T_c)=0,5,10$ ps | Hotbelt de 40 nm, baño de 2 K y nuevo ajuste de las fluctuaciones. pp. 95–97 / PDF 107–109. |
| K20, caso prioritario | $\chi=0.667$, con $\tau_{ee}(T_c)=6$ ps | Ajuste conjunto; tabla suplementaria 1 y figuras suplementarias 3–4, PDF 32–35. |
| V19v1, modelo térmico 2D | No fija una $\chi$ óptica estándar | Ejemplos de energía retenida media 1.5 y 2.5 eV, dispersión $0.1\bar E$; §IV, PDF 4–5. |

A19 define $E$ mediante $E_e(T_i)+E_{ph}(T_i)-E_e(T_b)-E_{ph}(T_b)=E$: es exceso **electrónico más fonónico**, no sólo energía de QP. La preparación TDGL eleva ambos sectores a una temperatura común. El $0.43$ del cierre simplificado no debe reinterpretarse como pérdida exclusivamente anterior a $t_0$. No se confirmó $\chi=0.49$ en las versiones consultadas. [A19v2, (5), §II.D y §III.A–C](https://arxiv.org/pdf/1805.00130v2), [A20, pp. 95–97](https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf).

V19 distingue explícitamente energía del fotón $E_\nu$ y energía entregada a electrones **y** fonones $E<E_\nu$. Su $t=0$ inicializa calentamiento térmico simultáneo; no calcula la duración óptica de la cascada omitida. Los ejemplos de $\bar E$ no permiten inferir $\bar E/E_\nu$ sin especificar $E_\nu$. El límite de jitter estudiado tampoco calibra retención. [V19v1, §§II y IV](https://arxiv.org/pdf/1807.07709v1).

**Conclusión para el escenario elegido:** proponemos $\chi_{\rm ref}=0.667$ como estándar operativo del **escenario Korzh**, porque corresponde al caso experimental prioritario. No se obtiene promediando ajustes de otros cierres. En K20 se define como retención tras conversión energética; su ilustración electrónica no autoriza depositar esa energía dos veces. [K20, métodos, PDF 14; suplemento, PDF 32–33](https://eprints.lancs.ac.uk/id/eprint/140252/3/Binder1.pdf).

En D.30 esta propuesta significaría $E_{\rm in}=0.667hc/\lambda$: aproximadamente **0.533532 eV a 1550 nm** y **1.067064 eV a 775 nm**. Ambas cifras son un supuesto trazable para la energía total entregada al nuevo modelo, no una calibración de su partición fonónica. El cambio de cierre puede alterar la respuesta aun con igual $E_{\rm in}$; eso se diagnostica, no se oculta reajustando la fracción después de ver una sola curva. La propuesta queda pendiente de adopción; no hay rango de retención aprobado.

| Magnitud | Qué puede fijarse por esta investigación | Qué permanece sin inferir |
|---|---|---|
| $E_\gamma$ | Conversión para cada longitud nominal | Energía retenida en el traspaso |
| $\chi$ | Definición y límites energéticos condicionales | Valor común o dependencia con el color |
| $s_\gamma$ | Convención matemática; 10 nm histórico | Anchura calibrada de la cascada |
| $t_0$ | Origen de la descripción resuelta | Retardo físico respecto de la absorción |
| Espectro inicial | B.39 como cierre vigente | Exactitud del reparto postcascada en esta muestra |
| Posición | Coordenada de un evento controlado | Ley óptica de posiciones y eventos de tapers |

Ninguna fuente consultada identifica simultáneamente $(E_{\rm in},s_\gamma,t_0)$ para el experimento prioritario. Los límites energéticos no son intervalos estadísticos de confianza.

## 4. Hotbelt y dos relojes distintos

Un **hotbelt** es una región excitada que abarca el ancho del hilo. Imponerla al inicio y observar su respuesta es una pregunta válida, distinta de investigar si un depósito localizado llega a formarla. Los 40 nm del caso 1D son longitud longitudinal de la región inicial; no hay conversión única a una gaussiana 2D. Igualar volumen, valor máximo o segundo momento conduce a convenios distintos. Tampoco una temperatura media transversal demuestra homogeneidad de poblaciones, condensado y corriente.

Para el objetivo de **formación** proponemos conservar una preparación localizada 2D y comprobar después la evolución transversal; un hotbelt impuesto serviría como control declarado, no como sustituto de esa prueba. El centro inicial puede fijarse en un evento controlado interior. No se presume que represente el promedio experimental: tal promedio exige posición de absorción, polarización y lectura. Cerca de un borde, renormalizar la gaussiana conserva $E_{\rm in}$, pero no calcula pérdidas previas ni la deformación real de la cascada.

Es lícito llamar $t=0$ al instante de transferencia. Para relacionarlo con el reloj óptico conviene escribir, con $\Delta$ indicando **1550 menos 775 nm**,

$$\Delta t_{\rm registrado}
=\Delta(t_0-t_{\rm abs})+\Delta\tau_{\rm resuelta}+\Delta t_{\rm lectura/recorrido}.$$

Esta descomposición es contabilidad temporal; los términos pueden estar correlacionados entre eventos. Un desplazamiento común de $t_0$ cancela en la diferencia; un desplazamiento dependiente del color no. El máximo de una distribución de tiempos tampoco se obtiene, en general, sumando los máximos de distribuciones de estos términos. K20 compara máximos ajustados de la respuesta temporal con una cadena calibrada. Un retardo determinista desde $t_0$ permite contrastar mecanismos, pero no reproduce por sí solo esa estimación experimental. [K20, métodos y figura 1](https://doi.org/10.1038/s41566-020-0589-x).

## 5. Fano permanece fuera del candidato determinista

El factor de Fano suele describir la varianza del número de excitaciones, $F_N=\mathrm{Var}(N)/\langle N\rangle$. En modelos efectivos también se parametriza la dispersión de energía electrónica disponible. Debe identificarse el sector antes de trasladar una varianza: con energía total fija pueden fluctuar electrones y fonones en sentidos opuestos. La dispersión de energía perdida al sustrato es otra contribución. [Kozorezov et al. 2017, §II, (11)–(12)](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=920691).

B.38a sólo podría representar la varianza de la **energía total retenida** si se justifican pérdidas fluctuantes o sectores excluidos. No se activa aquí. Una extensión futura deberá especificar soporte físico $[0,E_\gamma]$, momentos y correlaciones con posición, anchura y tiempo; una gaussiana sin límites no basta. Un único transiente determinista no predice jitter FWHM. Tampoco corresponde identificar variación de un parámetro de escenario con fluctuación aleatoria de eventos.

## 6. Decisión de alcance y protocolo pendiente

**Decisión aceptada:** estudiar primero desde la transferencia y mantener visible el retardo de cascada desconocido. La elección adicional del escenario material ajustado de Korzh ($D=0.5\,\mathrm{cm^2/s}$ y $R_\square=608\,\Omega/\square$) no lo convierte en una caracterización independiente ni adopta automáticamente todos sus parámetros fotónicos.

**Protocolo recomendado dentro de ese alcance:** un par determinista con forma espectral, posición y convención espacial comunes entre colores; retención común como hipótesis explícita; reloj de transferencia común. Esto aísla la respuesta a la energía retenida y permite estudiar la formación del hotbelt. Se registrarán $E_{\rm in}$ y los datos de transferencia completos. No se adopta aquí un valor nuevo de $\chi$, $s_\gamma$ o duración de cascada, ni se ajustarán para acercar una sola latencia a K20.

La revisión de 2019–2020 permite concretar la recomendación a **$\chi_{\rm ref}=0.667$ para el escenario Korzh**, con la interpretación y energías de §3.1; no constituye un estándar universal ni una decisión numérica ya aceptada. Para $s_\gamma$, 10 nm sigue siendo únicamente procedencia de implementación. Se retira la sugerencia de un intervalo 5–20 nm: el estudio específico de la preparación espacial deberá justificar una anchura para la gaussiana elegida, sin convertir un círculo o un hotbelt histórico en una medida de $s_\gamma$.

Antes de ejecutar ese par habrá que completar el transiente acoplado y la interfaz cinética, verificar energía de inyección en cuadratura, continuidad de los estados que no reciben impulso, soporte y convergencia. La elección numérica de un escenario inicial debe venir acompañada de su justificación y de una sensibilidad posterior prerregistrada. Esta investigación no inicia un barrido ni sustituye esos controles.

No hace falta volver a preguntar por el alcance ya aceptado. Una futura latencia óptica cuantitativa requerirá caracterizar o modelar pérdidas, extensión y duración de cascada por color, junto con la correspondencia de lectura. Las fuentes actuales no completan esos datos.

La única aclaración conceptual útil si “hotbelt” queda ambiguo es: **¿se desea comprobar su formación desde un depósito localizado, o imponerlo al inicio para estudiar la respuesta posterior?** La recomendación corresponde a comprobar su formación. Para preparar la primera corrida bastará aprobar un escenario explícito de transferencia con sus limitaciones; no hace falta pedir números supuestamente medidos que las fuentes no proporcionan.

## 7. Alcance de la comprobación de fuentes

Se consultaron fuentes primarias y el contrato local de 0.4. A20 se leyó también en la copia PDF local de 259 páginas: se comprobaron visualmente la meseta (p. 26/PDF 38) y la fuente modificada (p. 102/PDF 114), y textualmente pp. 95–97 para la comparación nueva. La tabla K20 y los preprints A19v2/V19v1 se cotejaron textualmente; sus capturas remotas no estuvieron disponibles. El PDF publicado de A19 devolvió acceso no autorizado; no se intentó eludirlo. Los hashes locales, localizadores y valores no adoptados están en el registro JSON. No se digitalizaron curvas ni se atribuyeron incertidumbres nuevas.
