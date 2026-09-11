---
title: "D. Síntesis, plan y verificaciones"
subtitle: "Desarrollo del modelo pySNSPD · Documento 4 de 5"
date: "8 de septiembre de 2026 · Revisión 0.3"
lang: es
---

# D.0. Qué cambia y para qué

La revisión 0.3 organiza la propuesta en tres tareas físicas: **A determina el espectro y la energía estática; B describe las ocupaciones y los intercambios de energía; C acopla el condensado, la corriente y el circuito.** El objetivo sigue siendo mejorar el modelamiento de SNSPD. La construcción microscópica es una herramienta para evaluar y ampliar los cierres del detector; su utilidad final debe medirse en los observables.

Esta iteración añade un quinto documento, E, que permite trabajar los conceptos necesarios mediante figuras, ejemplos resueltos y ejercicios con respuesta editable. Sus cinco temas comienzan abiertos. El avance del cuaderno y la validación física del modelo son registros distintos: resolver un ejercicio no valida un transitorio, y obtener un residuo pequeño no demuestra que un concepto se haya comprendido.

| Documento | Cambio principal | Pregunta de lectura |
|:--|:--|:--|
| A | Derivación de la referencia normal; alcance estático; significado de $\lambda$ y del ajuste de Simon | ¿Qué energía se varía y qué información microscópica contiene? |
| B | Notación FD; burbuja, retención y Fano; traslado de la energía adiabática | ¿Qué se distribuye, qué se conserva y qué puede escapar? |
| C | Ecuación gTDGL de partida y forma estable; conexión explícita con la fuerza propuesta | ¿Qué parte se conserva y qué cierre se propone ampliar? |
| D | Mapa de cambios, evidencia y próximos pasos | ¿Qué afirmaciones están verificadas y con qué alcance? |
| E | Cinco unidades pedagógicas con evaluación pendiente | ¿Se puede explicar y usar cada concepto sin confundir sus objetos? |

# D.1. Acuerdos físicos de esta versión

## D.1.1. Una corriente justificada y una ampliación por comprobar

La corriente Usadel incorporada al cierre modificado tiene una justificación de continuidad. Vodolazov desarrolla el término correctivo en la página 11 de [V]; la memoria explica la cancelación correspondiente en su anexo C [M]. La revisión 0.3 conserva esa justificación. Comparar corrientes críticas de aproximaciones distintas no basta para declarar inválida una de ellas.

El funcional común añade una pregunta: si se calculan la fuerza de amplitud y la corriente desde una misma energía, ¿mejora la consistencia del almacenamiento y su acoplamiento con las ocupaciones? A construye la parte estática a partir de [F]. C plantea cómo usarla sin confundir esa construcción con una derivación microscópica completa de la movilidad temporal.

En una rama uniforme, suave y térmica, las identidades que se quieren respetar son

$$
\begin{gathered}
X_{|\Delta|}^{\rm FD}
=\left.\frac{\partial f_e^{\rm FD}}{\partial|\Delta|}\right|_{T,\mathbf q},
\qquad
\left.\frac{\partial f_e^{\rm FD}}{\partial\mathbf q}\right|_{T,|\Delta|}
=\frac{\hbar}{2e}\mathbf j_s,\\
\frac{\partial X_{|\Delta|}^{\rm FD}}{\partial q}
=\frac{\hbar}{2e}\frac{\partial j_s}{\partial|\Delta|}.
\end{gathered}\tag{D.1}
$$

La segunda línea es igualdad de derivadas mixtas, a dirección de corriente fija. No sustituye la ecuación espacial de conservación de carga. A su vez, ni esta identidad ni la continuidad fijan por sí solas la disipación del condensado.

## D.1.2. Tres objetos que la notación debe separar

| Símbolo | Objeto | Qué permite calcular |
|:--|:--|:--|
| $f_{\rm FD}(E,T)$ | Probabilidad de ocupación, adimensional | Poblaciones térmicas y factores de bloqueo |
| $f_e^{\rm FD}$ | Densidad de energía libre electrónica | Fuerzas termodinámicas y entropía |
| $u_e$ | Densidad de energía interna electrónica | Almacenamiento y balance de energía |

A deriva la contribución normal antes de añadir la diferencia superconductora. B incorpora la extensión adiabática fuera de equilibrio:

$$
\begin{gathered}
u_e=U_{\rm vac}(|\Delta|,q)+u_{\rm qp},\\
u_{\rm qp}=4N_0\int_0^\infty E\rho(E;|\Delta|,q)f(E)\,dE.
\end{gathered}\tag{D.2}
$$

El factor cuatro corresponde a la DOS normal por espín y al conteo de excitaciones positivas de A–B. El término de vacío y el de excitaciones son partes de una misma contabilidad electrónica. Los fonones y el circuito exterior tienen sus propios almacenamientos. La dependencia de $U_{\rm vac}$ con el superflujo no debe volver a añadirse como otra energía local idéntica.

Una temperatura equivalente es una etiqueta de energía a espectro fijado. Si $u_{\rm qp}=\mathcal U(T_E;|\Delta|,q)$ y $\partial_T\mathcal U>0$, se puede invertir ese mapa. **La inversión no recupera una distribución no térmica completa.** B.5 muestra estados de igual energía con fuerzas QP distintas; C.4 explica cómo avanzar variables energéticas y consultar una temperatura cuando hace falta un coeficiente térmico.

## D.1.3. Material, cascada y condensado: dónde entra cada dato

La DOS fonónica cuenta modos. El espectro $\alpha^2F$ incorpora además el peso de interacción electrón–fonón; su integral ponderada define el acoplamiento adimensional $\lambda$. A.5 y E.5 explican esta diferencia. El apéndice de Simon consultado en A reescala resultados de acoplamiento débil por $2.1/1.76$ para ajustar el gap; no resuelve el sistema completo de Eliashberg con funciones dependientes de frecuencia [S]. Este ajuste resulta útil como referencia material, pero no permite atribuirle por sí solo todas las propiedades de un funcional de acoplamiento fuerte.

La condición inicial de burbuja resume una cascada previa. B.3 y B.8 separan absorción, reparto entre sectores y escape al sustrato, usando los alcances concretos de Vodolazov y Allmaras [V, AT, AP]. En un sistema electrónico y fonónico aislado puede fluctuar el número de QP mientras la energía total retenida permanece fija. Por ello una varianza de conteo Fano no se transforma automáticamente en una varianza de la energía total inyectada.

Un depósito compacto tampoco demuestra que cada evento tenga un único máximo espacial. La figura B.5 compara perfiles construidos de igual integral para enseñar esa diferencia. No calcula una tasa de eventos con dos lóbulos ni la distancia entre ellos en NbN. Una predicción de ese tipo requiere un modelo de cascada y transporte con geometría e interacciones especificadas.

# D.2. Mapa de lectura y traslado de ecuaciones

La reorganización evita que A deba enseñar a la vez espectro, cinética y dinámica. Se conservan las etiquetas de A que no se trasladan; los saltos de numeración son deliberados para facilitar comparar revisiones.

| Ubicación en 0.2 | Ubicación en 0.3 | Contenido |
|:--|:--|:--|
| A.1–A.14 | A.1–A.14 | Funcional estático y derivadas; se añaden A.7a–A.7c para la referencia normal |
| Sección A.4; A.15–A.24 | Sección B.11; B.43–B.52 | Energía adiabática, coordenadas espectrales y fuerzas a ocupaciones fijadas |
| A.16a | B.44a | Identidad auxiliar del traslado anterior |
| Sección A.5; A.25–A.31 | Sección C.12; C.35–C.41 | Corriente, continuación uniforme y comparación de cierres |
| A.25a–A.25c | C.35a–C.35c | Desarrollo de fase y continuidad |
| Sección A.6; A.32–A.34 | Sección A.5; A.32–A.34 | Material, $\alpha^2F$, $\lambda$ y acoplamiento fuerte |
| Figura A.2 | Figura C.4 | Comparación uniforme de corriente |
| Apertura de C | C.0a–C.0f | gTDGL de partida, forma estable y puente hacia C.9 |

La ecuación final consultada en la memoria está en el **anexo C, ecuación (C.3), página impresa 121, página 152 del PDF**. Los coeficientes usados en C proceden de sus ecuaciones (2.29)–(2.31), página impresa 31. Se indica la localización efectiva de la copia consultada para evitar ambigüedad al buscar anexos. C coteja también la ecuación (3.24) de Allmaras, página impresa 94 [AT].

La forma estable de C conserva el término temporal $\partial_t|\Delta|^2\,\Delta$ y la multiplicación del lado derecho por su factor de relajación. La propuesta se conecta después con la fuerza de amplitud. La equivalencia algebraica de esa escritura y la continuidad no certifican cualquier discretización cerca de $|\Delta|=0$; ese tratamiento sigue siendo un punto específico de implementación.

Para estudiar el conjunto: E.2–E.3 prepara A; E.4–E.5 prepara B; después C permite seguir el acoplamiento al detector. E.1 aporta el contexto de campos, partículas y excitaciones. También puede leerse A–C de forma continua y consultar E sólo cuando aparezca una dificultad concreta.

# D.3. Qué se ejecutó en esta revisión

Se reejecutaron en Geminga los tres cálculos ligeros heredados de 0.2 y se añadieron las figuras y comprobaciones del cuaderno. Se utilizó el entorno científico existente, con Python 3.10.20, NumPy 2.2.6, SciPy 1.15.3 y Matplotlib 3.10.8, limitando los hilos de cálculo a uno. La carpeta de trabajo es `sandbox/model_v0_3`, dentro de `/home/jdiaz/pysnspd`.

| Comprobación | Resultado registrado en 0.3 | Lo que verifica |
|:--|:--|:--|
| A: derivadas de la energía uniforme | Error relativo máximo $1.34\times10^{-10}$ en tres puntos | Fuerza y corriente frente a derivadas numéricas, en esos puntos |
| A/C: corriente máxima Usadel | $38.8503\,\mu$A; estable entre 400 y 3200 términos Matsubara | Convergencia del problema uniforme de acoplamiento débil con los parámetros declarados |
| B: 3000 eventos de cada reacción | Residuo relativo de intercambio $3.96\times10^{-20}$ | Signos y conservación en eventos sintéticos |
| B: invertir energía térmica | Error absoluto máximo $2.22\times10^{-16}$ en unidades de la prueba | Inversión a espectro fijo y capacidad positiva en el intervalo ensayado |
| B: estados de igual energía | Razón de fuerzas QP $25.1873$ entre dos paquetes elegidos | Energía sola no fija todos los momentos de la distribución |
| C: celda aislada ilustrativa | Deriva de energía $6.75\times10^{-13}$; diferencia entre coordenadas $1.14\times10^{-12}$ o menor, en unidades adimensionales | Balance y equivalencia numérica del modelo de celda declarado |
| C: circuito con resistencia prescrita | Residuo integrado de balance $1.67\times10^{-5}$ o menor, en unidades adimensionales | Contabilidad de potencia para esa entrada temporal y cuadratura |
| E: ejemplos resueltos | BdG: $E=\pm5$, pesos $0.8/0.2$; curvatura reducida $1$; autovalores de rigidez adimensionales $1,3$ | Álgebra de los ejemplos, sin resolver las actividades evaluadas |

Los números coincidentes con 0.2 son **reproducción de pruebas anteriores**, no nuevas predicciones materiales. En la prueba instantánea de B, el residuo de la regla de la cadena es $6.94\times10^{-18}$; su pequeñez se interpreta a precisión numérica, no como cero matemático inferido de una ejecución.

Los parámetros de la comparación uniforme están en `verificaciones/A_resumen.json`; B y C declaran sus unidades y modelos en `B_verificaciones.json` y `C_checks.json`. Las tablas CSV permiten rehacer los gráficos. La celda polinómica, el circuito prescrito y los paquetes sintéticos se identifican como tales en los pies de figura. No se ejecutó un transitorio espacial del detector, una cascada estocástica ni un cálculo DFT/DFPT nuevo.

## D.3.1. Figuras y controles de entrega

El conjunto incorpora 17 figuras distintas: dos en A, cinco en B, cuatro en C y seis en E. Siete son recursos nuevos para B/E; la antigua comparación A.2 acompaña ahora C.12. Cada figura se conserva en PNG y PDF vectorial. Los cinco documentos tienen fuentes Markdown editables y LaTeX compilable.

La comprobación de entrega verifica que todas las etiquetas de ecuación aparezcan en los PDF, que las figuras estén incrustadas y que no haya avisos de símbolos ausentes o desbordes. La revisión visual y el conteo final de páginas se registran en `verificaciones/QA_estructura.json` y `verificaciones/QA_visual.json`. El manifiesto de la entrega identifica los archivos y sus hashes; el paquete conserva la estructura del repositorio para que las rutas relativas de las imágenes sigan funcionando.

# D.4. Próximos pasos del modelo

1. **Cerrar las convenciones estáticas.** Contrastar la derivación de A con la elección material: normalización de DOS, gap, rama espectral y régimen de corriente. Decidir después si el primer catálogo utiliza acoplamiento débil, un ajuste fenomenológico identificado o una construcción de acoplamiento fuerte más completa.
2. **Elegir el estado retenido de B.** Comparar una distribución energética resuelta con cierres por momentos bajo las mismas entradas materiales. Medir errores en potencia electrón–fonón y fuerza de amplitud, además del error en energía. Una temperatura equivalente puede seguir usándose para lectura o evaluación de coeficientes.
3. **Probar el acoplamiento local de C.** Sustituir la celda pedagógica por un cierre microscópico explícito, verificar balance aislado y disipación, y precisar la movilidad. Mantener una prueba de equivalencia de la escritura temporal antes de introducir cambios espaciales.
4. **Integrar espacio y circuito.** Revisar continuidad, fronteras, paso temporal, extensión cerca de ceros del condensado y contabilidad de energía con el circuito. Sólo entonces comparar transitorios, latencias y señales con la referencia anterior y con datos.
5. **Evaluar la condición inicial cuando haya un modelo de cascada.** Separar retención media, fluctuaciones de reparto, escape y forma espacial. La hipótesis de varios lóbulos debe contrastarse como distribución de eventos; no basta cambiar el ancho de una única burbuja.

Estos pasos son propuestas de desarrollo, no cambios realizados en el solver. La revisión local de código utilizada como contexto es `391bb301a6c36f4438754a51a10a7fe2aade0555`. Los módulos de producción no forman parte de esta actualización documental.

# D.5. Seguimiento del aprendizaje

E dispone de cinco temas con identificadores estables E01 a E05. Cada uno incluye propósito, figura pedagógica, desarrollo, ejemplo resuelto, tres actividades y criterio de comprensión. Las respuestas se escriben en su Markdown, debajo de cada enunciado.

| Tema | Estado de 0.3 | Evidencia pendiente |
|:--|:--|:--|
| E01. Modelo Estándar | ABIERTO | Respuestas y explicación de campos frente a ocupaciones |
| E02. Nambu y espín | ABIERTO | Cálculos y distinción entre índices y estados físicos |
| E03. Cálculo funcional | ABIERTO | Variación, bordes y uso correcto de la envolvente |
| E04. Electrón y hueco | ABIERTO | Contabilidad de carga y energía; vínculo con superconductividad |
| E05. DFPT y DOS | ABIERTO | Modos, normalización espectral y pesos de interacción |

El estado cambia a **CERRADO / APROBADO** con al menos 8 de 10 puntos y el criterio esencial del tema satisfecho. En la siguiente revisión se añadirá la pauta correspondiente a las respuestas recibidas, con corrección razonada y nueva oportunidad de resolver lo pendiente. Un tema sin respuestas permanece abierto y sin nota. Si se detecta después una confusión esencial, el registro permite reabrirlo indicando el motivo; el historial no se borra.

# D.6. Fuentes y trazabilidad

[M] *Multiscale Modeling of the Transient Response of Superconducting Nanowire Single-Photon Detectors*, copia `memoria_02.pdf` de `/home/jdiaz/memoria/main`, 173 páginas. Se consultaron las ecuaciones (2.29)–(2.31), (C.3) y el desarrollo de continuidad del anexo C. Las páginas impresas y las páginas del PDF se distinguen en C.

[V] D. Yu. Vodolazov, *Single-Photon Detection by a Dirty Current-Carrying Superconducting Strip Based on the Kinetic-Equation Approach*, Phys. Rev. Applied **7**, 034014 (2017). [Preprint consultado](https://arxiv.org/pdf/1611.06060): ecuaciones (16)–(22) para condiciones iniciales y página 11, ecuación (36), para el término correctivo.

[F] P. Virtanen, A. Vargunin y M. Silaev, *Quasiclassical free energy of superconductors: Disorder-driven first-order phase transition in superconductor/ferromagnetic-insulator bilayers*, Phys. Rev. B **101**, 094507 (2020), ecuación (20). [Artículo](https://doi.org/10.1103/PhysRevB.101.094507) y [preprint v1](https://arxiv.org/pdf/1909.00992v1), donde el título es *Quasiclassical expressions for the free energy of superconducting systems*. A explicita la reducción utilizada.

[AT] J. P. Allmaras, *Modeling and Development of Superconducting Nanowire Single-Photon Detectors*, tesis doctoral, Caltech (2020). [PDF oficial](https://thesis.caltech.edu/13748/08/Allmaras_Thesis_Final.pdf). B localiza las secciones de cascada, escape y condición inicial; C utiliza la ecuación (3.24) y su discusión numérica.

[AP] J. P. Allmaras et al., *Intrinsic Timing Jitter and Latency in Superconducting Nanowire Single-Photon Detectors*, Phys. Rev. Applied **11**, 034062 (2019). [Artículo](https://doi.org/10.1103/PhysRevApplied.11.034062) y [preprint v2](https://arxiv.org/pdf/1805.00130v2). Las referencias específicas de fluctuaciones Fano y escape están desarrolladas en B.

[S] Simon et al., [preprint arXiv:2501.13791v3](https://arxiv.org/pdf/2501.13791v3), sección VII.3, ecuaciones (11)–(14). La numeración citada corresponde a esta versión del preprint. A distingue el ajuste del gap de una solución completa de Eliashberg.

Los documentos 0.1 y 0.2, las observaciones recibidas y los PDF de las fuentes se conservan como procedencia. El paquete 0.3 distribuye los documentos y figuras originales de esta revisión; no redistribuye los PDF de terceros. Todos los archivos de la entrega y las comprobaciones se copian a Geminga bajo las mismas carpetas de versión.
