# Etapa 3: resolver el ensayo espacial admitido con bordes y circuito

**Estado: PREPARADA_NO_INICIADA. Implementación posterior a esta entrega.**

La etapa 2 se cierra **como etapa de desarrollo**, según la
[decisión de cierre](../stage2/closure_20260922/closure_decision.json), para continuar
con la evidencia disponible. Su certificado estricto de convergencia dinámica
entre mallas permanece incompleto. Son dos decisiones distintas: el cierre
habilita la siguiente secuencia de trabajo y no convierte los resultados
faltantes o fallidos en aprobados. La preparación actual no contiene un solver
espacial nuevo ni resultados de etapa 3.

El [dictamen de etapa 2](../stage2/stage2_admission.json) y la
[revisión de evidencia](../stage2/practical_review_20260922/saved_results_audit.json)
conservan los pases temporales propios del método protegido en una y dos celdas,
los casos especiales y el fallo auxiliar del par fonónico de 2049 nodos:
`4.565195781e-5` frente a `2.5e-5`. No se relajan esas cifras ni se atribuye
convergencia a las mallas no comparadas. La
[revisión del alcance](../stage2/practical_review_20260922/validation_scope_review.md)
explica por qué ese certificado incompleto no debe ser un bloqueo general
para trabajo que tiene otras dependencias.

El valor histórico `stage2.numerical_admission = false` describe la falta de
ese certificado; **no bloquea por sí solo el desarrollo de etapa 3**. Cada ensayo
espacial debe declarar lo que quiere medir y sus controles antes de ejecutarse.
La precisión de una futura trayectoria se juzgará en ese alcance, sin heredar
por nombre el pase de una celda ni ocultar la incertidumbre espectral pendiente.

Las trayectorias excitadas de etapa 2 mantuvieron `Gamma = 0`; el equilibrio de
dos celdas usó valores fijos `Gamma = (0, 0.1)`. No evolucionaron fase, potencial,
funcional espacial ni circuito. Su balance con amplitud móvil no verifica el
trabajo debido a una variación de `Gamma(q)` ni el trabajo de un puerto eléctrico.
Estos controles pertenecen al ensayo nuevo y no exigen repetir la etapa cerrada.

## Modelo que se lleva a esta etapa

El contrato continuo es [D.4.3 del modelo 0.4](../../modelo_v0_4/D_sintesis_plan_y_verificaciones_v0_4.md),
con la [adenda del circuito de la memoria](../../modelo_v0_4/actualizaciones/circuito_memoria_20260922.md).
Esta última sustituye C.32–C.34, D.28–D.29, D.33 y la inicialización circuital
de D.30. El [contrato preparado](entry_contract.json) registra las dependencias
por ensayo; sus propuestas de precisión todavía no son tolerancias vinculantes.

- Energía, fuerzas cartesianas y corriente proceden del mismo funcional D.8–D.10.
  Se conserva el catálogo R2 y se declara la representación ocupacional usada.
- Se comprueba D.36: la respuesta a perturbaciones espaciales cortas debe tener
  signo positivo certificado en el dominio admitido. Refinar tiempo o recortar
  un valor propio negativo no corrige una inestabilidad de las ecuaciones.
- Las poblaciones electrónicas y fonónicas son variables distintas. En las
  interfaces se compara `f(E)` a energía física común; no se iguala `p(x)` por
  índice cuando cambia el espectro. El trabajo del espectro móvil se conserva.
- La movilidad KWT, la relajación cinética y el escape mantienen sus papeles.
  `Q_Delta` se entrega una sola vez a las ocupaciones.
- La red exterior es la de la memoria: estado `(Ib, Is, vc)`, fuente de tensión,
  resistencia e inductancia de polarización y capacitor de lectura. `Is` es la
  corriente total del detector; el voltaje medido es `Vout` en la carga.
- La convención es pasiva: para corriente izquierda–derecha,
  `Vdev = phi_L - phi_R` e `Is*Vdev` es trabajo que entra en el dominio resuelto.
  El adaptador histórico debe transformar juntos los signos de corriente y
  voltaje, según CM.5. La inductancia serie representa sólo la región exterior;
  los 10 nH totales de referencia de la memoria no se suman otra vez a la
  energía de superflujo resuelta.

El parámetro efectivo de núcleo sigue siendo `delta = 0.10 Delta0`. Resolver
mejor la malla no lo valida físicamente. Una entrada Debye sintética permite
verificar el algoritmo; la procedencia incompleta de NbN sigue excluyendo
atribuirle tasas o latencias absolutas de un dispositivo.

## Secuencia de implementación posterior

| Paso | Trabajo y condiciones | Resultado que permitirá continuar |
|:--|:--|:--|
| **3A. Funcional espacial estático** | Tira pequeña, poblaciones congeladas, estado uniforme y gradientes débiles. Derivar fuerzas y corriente de una energía discreta, evaluar soporte y D.36. | Comparación independiente de energía/derivadas/corriente y dominio estable identificado. Sin tiempo ni señal. |
| **3B. Bordes y reservorios** | Sustituir el cierre periódico de diagnóstico por un dominio abierto: paredes aislantes, continuaciones 1D e interfaces con el sector 2D, reservorios D.24–D.27. Empezar sin corriente y luego con corriente débil admitida. | Flujos de cara conservativos, trabajo de reservorios contabilizado y región central suficientemente independiente de la longitud de las continuaciones. |
| **3C. Conservación de carga y potencial** | Resolver D.22 con la misma corriente superconductora y la corriente normal. Fijar una referencia de potencial, corriente de terminales y puertos de trabajo. Comprobar el trabajo espectral en recorridos controlados con flujo variable. | Continuidad de corriente, límites uniforme y normal verificables, concordancia entre trabajo volumétrico y de puerto y variación energética consistente con `Gamma(q)`. No se añade una dinámica de acumulación de carga ausente del modelo. |
| **3D. Circuito de la memoria** | Acoplar las tres ODE CM.4 con signos pasivos, identificar `Lk,ext` antes del ensayo y preparar el punto fijo CM.10. Reutilizar primero la referencia con resistencia prescrita. | Evolución y almacenamiento de `Ib`, `Is`, `vc`; balance CM.8–CM.9 y resolución del intercambio circuito–película. |
| **3E. Dinámica débil y depósito local** | Perturbación pequeña sin fotón; después, un depósito localizado sintético de D.30/CM.11. Mantener todos los estados dentro del soporte y del dominio D.36. | Trayectorias exploratorias con error medido en los observables registrados, balances y límites físicos; comparación espacial/temporal/espectral pertinente. No implica una predicción material ni validación de umbral o latencia. |

El ensayo final de esta etapa tendrá **terminales, reservorios y circuito**.
La periodicidad del primer diagnóstico no representa los bordes físicos del
SNSPD. Antes de pasar del piloto longitudinal a 2D se registrarán la geometría
y las direcciones adicionales en las que debe comprobarse D.36: una prueba
1D no certifica automáticamente una película 2D.

## Primer ensayo utilizable: 3A

Se propone una sección longitudinal uniforme de longitud **360 nm**, ancho
**120 nm** y espesor **7 nm**, dimensiones del segmento de referencia de la
memoria. Su uso como tira 1D con **8/16/32 celdas** es una propuesta de diagnóstico,
no una reproducción del dispositivo ni una declaración de resolución suficiente.
La sección fija las medidas de volumen y cara; las tres discretizaciones
representarán exactamente la misma longitud física.

Para aislar el funcional se propone periodicidad de amplitud y poblaciones con
una torsión de fase declarada: `Delta(L) = exp(i*q0*L)*Delta(0)`. Primero se usa
`q0 = 0`; después se registran valores débiles dentro del dominio certificado.
No se presupone una corriente permitida sólo porque esté dentro de la tabla.
Las familias de entrada serán:

1. Estado uniforme, comparado con la energía y corriente uniformes del catálogo.
2. Modulación suave de amplitud y, separadamente, de fase, manteniendo fijo el
   perfil físico al refinar la malla.
3. Los mismos campos con una población no térmica física y de soporte explícito.

Se fijan amplitud de fondo, tamaño de las modulaciones, `q0`, temperatura de
preparación y perfiles ocupacionales **antes de evaluar errores**. Una población
preparada como Fermi–Dirac permanece fija en las variaciones de la energía
interna: no se vuelve a termalizar después de perturbar el campo. Si se prueba
además un funcional libre a temperatura fija, se registra como comprobación
diferente con su derivada correspondiente.

El control negativo reproduce D.36 con `|Delta|/Delta0 = 0.60`, `q*ell0 = 1`
y `delta/Delta0 = 0.10`: debe detectarse y rechazarse para evolución. Se conserva
su resultado; no se integra ni se corrige su rigidez.

Los observables iniciales son energía total, fuerzas cartesianas, corriente de
cara, errores de variación, cambio bajo fase global y margen del menor valor
propio del símbolo principal. La referencia combina el límite uniforme R2 y
una variación/linealización independiente de **las mismas ecuaciones**. Se
separa el error de tabla y cuadratura del espacial. Las variaciones de fase y
potencial vector comprueban el convenio de calibre aunque el escenario físico
se especialice a `A = 0`.

Las variaciones de fase deben ejercitar explícitamente
`Gamma = hbar*D*|q_delta|^2/2`, con el argumento regularizado de C/D. Se compara
la corriente obtenida de la **energía espacial completa**, incluyendo sus
términos de gradiente, con su variación conjugada. En 3C se comprobará además
la variación de energía sobre un recorrido controlado de flujo a ocupaciones
fijas, incluyendo el desplazamiento de niveles espectrales. No se sustituye
ese control por el balance previo a `Gamma` constante ni por la corriente
espectral aislada, que omitiría parte del funcional.

Salidas previstas: registro previo, tablas de errores por caso y malla, mapa del
dominio probado y sus rechazos, figuras de las respuestas, hashes y coste.
La primera decisión será si ese funcional admite el ensayo de bordes 3B, no
si ya reproduce un pulso de detección.

## Precisión y comparaciones proporcionales

No se fija aquí un porcentaje universal nuevo. Antes de cada ensayo se declara
la diferencia física que interesa resolver y el error numérico máximo que
permitiría distinguirla. El registro distribuirá ese presupuesto entre espacio,
tiempo, espectro, tabla y referencia cuando esas fuentes intervengan.

| Tipo de control | Presupuesto y comparación que se deben registrar |
|:--|:--|
| Identidades algebraicas y signos | Escala física, incertidumbre de la referencia y sensibilidad al paso de diferenciación; suficientemente precisos para detectar factores y energía duplicada. No equivalen a precisión experimental. |
| Respuesta estática espacial | Error en energía, fuerza y corriente frente a referencias independientes; tres mallas iniciales y comportamiento de refinamiento explicable dentro de la incertidumbre de referencia. |
| Bordes, potencial y circuito | Residuo de flujo y trabajo normalizado por el intercambio ejercitado, punto fijo y comparación analítica o lineal independiente cuando exista. No normalizar sólo por una gran energía de fondo. |
| Trayectorias débiles o depósito local | Amplitud, energía transferida, corriente y voltajes seleccionados a tiempos comunes. Registrar norma y precisión útil para cada observable; cerca de cero usar una escala absoluta declarada. |
| Futuros umbrales o latencias | Sólo si se incorporan expresamente al alcance: observable, línea de base, regla de cruce y estabilidad de la clasificación frente a la resolución. No se infieren del piloto 3A. |

Las cifras que se elijan serán propuestas hasta congelar el registro previo.
No se escogerán para convertir el par fonónico ya fallido en PASS. Su evidencia
puede reutilizarse para estimar incertidumbre, pero la convergencia no medida
seguirá declarada. Los fallos de unidades, signos, balance físico, finitud,
Pauli/positividad, soporte o D.36 detienen el cálculo afectado. Un presupuesto
de precisión incumplido conserva su dictamen y limita la conclusión; no
bloquea automáticamente un trabajo independiente de ese observable.

## Reutilización y ejecución

La memoria aporta geometría y topología, convenciones históricas y límites
comparables. El [diagnóstico de circuito con entrada prescrita](../stage2/practical_review_20260922/circuit_diagnostic.json)
aporta una referencia para CM.4 y CM.8; usa parámetros sintéticos comparables,
no una película resuelta ni la partición identificada de los 10 nH totales.
Las trazas históricas de la memoria no son una referencia exacta de la nueva
cinética o energía común.

Los 320 pasos usados como referencia práctica de etapa 2 corresponden a sus
ensayos en el intervalo normalizado `0–2`. **No definen el paso temporal de la
PDE espacial**. Éste y el intervalo de intercambio circuital se seleccionan con
las escalas y comparaciones del ensayo nuevo. Cambiar un informe no exige
recalcular física; cambiar el operador o la representación exige las
comprobaciones que afecte.

Antes de una trayectoria se mide el coste de un paso. Un diagnóstico incierto
usa límite de 240 s; un cálculo previsto por encima de cinco minutos se deja
con propósito, salidas y recursos en `/home/jdiaz/GEMINGA_COMMANDS.md` para
que lo ejecute el usuario. Esta preparación no lanza cálculos ni deja un lote
espacial activo. La validación de núcleo corresponde a D.4.4 y los transitorios
completos y cualquier promoción a producción, a D.4.5.
