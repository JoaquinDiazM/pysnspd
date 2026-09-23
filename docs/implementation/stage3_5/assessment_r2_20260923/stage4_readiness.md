# Entrada propuesta a etapa 4: núcleo y disipación

Revisión independiente del 23 de septiembre de 2026. Es una recomendación de alcance, no una nueva admisión física, una implementación ni una orden de cálculo. Se revisaron las 127 entradas de las 15 familias del [ledger](../research_20260923/range_ledger.json), la [secuencia vigente](../../SECUENCIA_VIGENTE.md), D.4 y los módulos experimentales. La trazabilidad y clasificación están en [stage4_readiness.json](stage4_readiness.json).

**Se puede preparar y comenzar la parte diagnóstica de etapa 4 sin elegir un ancho fotónico.** Lo primero es contrastar el funcional, su estabilidad y la respuesta disipativa en estados controlados sin fotón. Para interpretar un transiente espacial posterior deben completarse las capacidades dinámicas que utilice ese ensayo. Los datos actuales no admiten una simulación física del hotbelt de Korzh ni su latencia óptica.

## 1. Qué cierre de 3.5 sería honesto

El [contrato de 3.5](../entry_contract.json) permite terminar la revisión documental con un registro de restricciones, incertidumbres explícitas y una decisión sobre lo que podrán afirmar las etapas 4–5. No exige obtener una medición de cada parámetro ni construir una caja universal de confianza. Esa entrega puede cerrarse con alcance si la síntesis final registra los resultados y la continuación necesaria.

Debe conservarse la decisión **«Caracterizar primero la cascada; dejar el ancho abierto»**. La caracterización documental disponible identifica qué parte de la transferencia no se puede inferir: un percentil de energía total de otra muestra no determina el perfil fonónico de Korzh. Esta conclusión es útil, pero no equivale a haber caracterizado su cascada completa. El [plan de cascada](../research_20260923/cascade_plan.json) sigue siendo una dependencia explícita de cualquier preparación fotónica física. No se adoptan 5–20 nm ni la equivalencia histórica 1.4–1.9 nm.

Dictamen recomendado: **revisión de entradas y límites de 3.5 susceptible de cierre; dominio físico completo y transferencia fotónica no admitidos; etapa 4 preparada para diagnósticos que no usan esa transferencia**. El cierre efectivo corresponde a la síntesis y decisión del trabajo principal. No debe describirse como «todas las incertidumbres resueltas» ni «cascada validada».

## 2. Qué hay realmente en el código

| Componente | Capacidad existente | Límite relevante para etapa 4 |
|---|---|---|
| `energy_catalog.py` | Energía y derivadas coherentes; consulta causal; normal exacto separado | Rechaza $0<|\Delta|<|\Delta|_{\min}$ del catálogo. El punto $\Delta=0$ no completa ese hueco. |
| `spatial_nodal.py`, `spatial_open.py`, `mixed_spatial.py` | Fuerza y corriente derivadas de la misma energía; operadores espaciales y diagnósticos D.36 | La evidencia es estática/instantánea. El empalme mixto sólo identifica el modo de campo transversal uniforme; no transporta automáticamente poblaciones a igual energía. |
| `cell_closures.py` | Matriz KWT positiva y disipación coherente | `KWTMobility.coefficients` fija **literalmente 0.50 y 2.47 ps**; todavía no recibe estos tiempos como parámetros. |
| `coupled_cells.py` | Cinética y BGK de una/dos celdas, con ledgers energéticos | No integra el funcional espacial, la fase, el potencial ni el circuito. `tau_kin=0.7` es el escenario sintético, expresado en tiempo normalizado; sólo vale 0.7 ps cuando $t_{ref}=1$ ps. |
| `spatial_dynamics.py` | RHS instantáneo KWT, reparto de calor y balance común | No proporciona integrador espacial ni kernel de colisiones. Su identificación de cuadraturas es de un único nodo por registro, no interpolación general. |
| `reservoir_spatial_dynamics.py` | Carga de corriente prescrita y radio terminal fijo, con trabajo de reacción | No implementa $a_b(I_s(t))$, su derivada, traza normal exterior ni transiente. |
| `reservoir_transport.py` | RHS de intercambio electrónico con reservorio FD a igual energía | No es el empalme completo entre dos regiones dinámicas 2D–1D ni una garantía de soporte para cualquier paso temporal. |
| `electrical_ports.py` | Tres ODE CM y su potencia | Debe acoplarse a la misma corriente/puerto y a un integrador común antes de atribuir una señal al detector dinámico. |

Estos límites ya están declarados en los propios módulos y en el [cierre de etapa 3](../../stage3/closure_20260923/closure_decision.json). Por tanto, no es necesario repetir toda su campaña estática para comenzar a trabajar; tampoco basta reutilizar sus instantáneas como prueba temporal.

## 3. Secuencia mínima útil

### 4A. Diagnóstico del cierre, sin fotón

Preparar estados uniformes y perturbaciones suaves de amplitud/fase dentro del soporte admitido. Se pueden comparar energía, fuerza, corriente, símbolo principal y respuesta KWT instantánea. Las poblaciones deben declararse: vacío, FD a un baño identificado o una distribución de control. Una perturbación espacial de prueba no representa una absorción óptica ni el ancho que el usuario dejó abierto.

Esta parte requiere escalas materiales coherentes, consultas soportadas y derivadas del mismo funcional. **No requiere** fracción retenida, ancho de depósito, tasas fonónicas absolutas, calibración de amplificador ni una distribución de eventos. Si una población se mantiene fija, el resultado es una variación a población fija; no se llama equilibrio térmico a esa derivada.

La primera comparación temporal puede ser puramente algebraica: coeficientes y autovalores de movilidad radial/fase para las parejas heredada, Allmaras y K20. No es una latencia ni un transiente. Antes de ejecutar una evolución con una pareja distinta, parametrizar los tiempos KWT de manera explícita, conservar por defecto el escenario heredado y registrar los parámetros realmente usados. Cambiar un YAML actualmente no modifica esa fórmula experimental.

### 4B. Completar el transiente débil que el ensayo necesite

Un ensayo espacial dinámico acoplado debe evolucionar campo, poblaciones, potencial y, si se incluye, circuito con el mismo ledger energético. Debe tener comparación temporal de observables comunes e integrar trabajo y pérdidas. Se parte de una perturbación controlada suave, sin fuente fotónica ni pretensión de hotbelt.

Los requisitos se activan **por dependencia**:

- Si hay reservorios que siguen la corriente del circuito, hacen falta su rama estable $a_b(I_s)$, evolución y trabajo correspondientes. Mantener el radio fijo puede seguir siendo un control restringido, pero no completa D.27.
- Si se usa el empalme cinético 2D–1D, hacen falta remapeo a igual energía, flujos conservativos y propiedad única de la población fonónica compartida. Un ensayo enteramente 2D no necesita fingir que ese empalme ya existe.
- Si se usan terminales exteriores, la traza normal y el balance se calculan allí. Una cara interior adyacente no los reemplaza.
- Un ensayo periódico o aislado puede verificar un subconjunto del algoritmo; se identifica como tal y no se presenta como dispositivo con lectura CM.

No se exige solucionar todas estas variantes antes del primer control. Se exige resolver las que ese control utiliza. Los componentes aún no implementados permanecen pendientes, sin trasladarlos a etapa 5 si ya afectan al resultado de etapa 4.

### 4C. Núcleo, disipación y primeros eventos

Después del ensayo débil correspondiente, avanzar hacia supresión de amplitud, perfiles de núcleo y eventos de fase. D.4.4 pide confrontar fuerzas, perfiles y barreras con una referencia espacial adecuada; la conservación de energía no sustituye esa comparación. Una referencia analítica de un régimen límite o una comparación publicada con condiciones identificadas puede ser útil. No se exige equivalencia microscópica completa de toda la materia condensada.

Antes de afirmar que se ha resuelto un núcleo hay que cubrir la amplitud continua hasta el régimen utilizado, incluidos estados cercanos a cero. El catálogo positivo actual empieza en aproximadamente $0.08\Delta_0$; no se puede saltar ese intervalo por interpolación al punto normal. La semilla aceleradora $|\Delta|/\Delta_0\in[0.985,1.005]$, $\Gamma/\Delta_0\in[0.003,0.009]$ tampoco cubre un núcleo. Usar consulta directa o una extensión verificada es una necesidad de soporte, no una mejora opcional de precisión.

Los eventos pueden inducirse mediante un estado o forzamiento de control registrado para estudiar el mecanismo. Eso no convierte su energía o longitud en las de un fotón. Una interpretación de formación de hotbelt **por 775/1550 nm** seguirá esperando una transferencia coherente con la caracterización de cascada.

## 4. Parámetros que no deben resolverse por conveniencia

**Movilidad KWT.** Los pares $(\tau_{ee}(T_c),\tau_{ep}(T_c))$ de 0.50/2.47, 5/24.7 y 6/24.7 ps son escenarios distintos. La elección material $D=0.5$ cm²/s, $R_\square=608$ Ω no selecciona automáticamente el tercero. La matriz KWT multiplica la movilidad radial por $1/R$ y la tangencial por $R$, con $R=\sqrt{1+4|\Delta|^2\tau_\psi^2/\hbar^2}$: cambiar esos tiempos no equivale a cambiar globalmente el reloj. Se recomienda empezar por el contraste algebraico y después por el menor ensayo dinámico que distinga las respuestas. No ajustar tiempos para forzar 4.2 ps.

**BGK.** $\tau_{kin}$ controla redistribución hacia una FD de la misma energía y es independiente de $\tau_\psi$. Para probar mecanismos se puede declarar BGK desactivado o el valor sintético histórico; `infinity` ya desactiva ese término en el ensayo de celdas. Ninguna opción es una identificación cinética de NbN. Si la conclusión sobre disipación depende de BGK, su identificación o sensibilidad debe formar parte de etapa 4; no se posterga hasta después de anunciar una predicción física. No es imprescindible implementar un operador ee microscópico completo para un resultado condicional, pero sí dejar clara la dependencia no identificada.

**$\delta_{reg}$.** $0.10\Delta_0$ es la terminación efectiva del núcleo, no un epsilon de máquina. Las alternativas históricas 0.05/0.10/0.20 pueden organizar un contraste del cierre con procedencia, manteniendo la misma energía, fuerzas y corriente. No son un intervalo admitido ni una extrapolación obligatoria a cero. `delta_regularizer_bar` está fijado como atributo en el funcional actual; cualquier parametrización debe ser explícita y coherente en todas sus derivadas. El refinamiento de malla prueba la discretización para un $\delta$ dado, no valida ese $\delta$.

**D.36.** El signo positivo del símbolo principal, resuelto frente a su incertidumbre numérica, es indispensable para continuar una trayectoria espacial del candidato. No se exige un margen positivo arbitrario adicional. Si el signo se vuelve negativo o no resoluble, se detiene esa trayectoria y se informa la frontera de uso. Evaluar el contraejemplo como diagnóstico sigue siendo válido; integrarlo con una rigidez recortada no lo es. D.36 es necesario, no una garantía completa de estabilidad cinética.

**Circuito y latencia.** Se mantiene CM: $R_b=10$ kΩ, $L_b=1$ μH, $R_L=50$ Ω, $C=100$ pF y referencia total 10 nH. Material y longitud nuevos obligan a recalcular la partición resuelta/exterior en la referencia, con $L_{ext}>0$; no a sumar los 96 nH del experimento. Los puntos 15.5/21.5 μA son referencias de interés, sujetos a una rama inicial admisible. La precisión propuesta de 0.1 ps concierne a la diferencia de cruces cuando exista tal observable. No obliga a todas las fuerzas o poblaciones a ese presupuesto, ni identifica máximos IRF. No hace falta conocer el amplificador para estudiar el núcleo; sí para afirmar una reproducción cuantitativa de su señal medida.

## 5. Qué es indispensable y qué puede esperar

| Requisito | Momento en que pasa a ser indispensable | Qué puede esperar |
|---|---|---|
| Unidades y escalas conjuntas; identidad de fuerza/corriente/energía | Todo diagnóstico del nuevo escenario | Una calibración independiente de cada escala no es condición para un ensayo identificado como modelo. |
| Soporte y causalidad del catálogo | Toda consulta; cobertura cerca de cero antes de núcleos | Extender zonas que el ensayo no visita. |
| Signo D.36 | Toda trayectoria espacial, en los estados que recorre | Certificar una caja global jamás utilizada. |
| Soporte Pauli/fonónico y balance integrado | Toda dinámica que use poblaciones | Repetir pruebas ajenas al camino de código o imponer precisión relativa sobre cantidades despreciables. |
| Normalización material de DOS/acoplamiento y corte energético | Antes de atribuir tasas e intercambio a NbN | Los controles estáticos o cinéticos sintéticos con escala declarada. El corte debe cubrir los procesos retenidos: más nodos no arreglan un corte insuficiente. |
| Movilidad y BGK | Elección explícita en toda evolución; contraste antes de atribuir una mejora a la disipación | Una identificación universal de esos cierres; nunca ocultar su dependencia en un resultado de muestra. |
| Reservorios, potencial, circuito e interfaces | Antes del ensayo que use cada uno | El empalme 2D–1D si se mantiene íntegramente 2D; la lectura si el control no la utiliza. |
| Contraste de $\delta$ y referencia de núcleo | Antes de admitir física de núcleo o barreras/eventos | No puede aplazarse a etapa 5 si etapa 4 afirma haberlo resuelto. |
| Cascada, reparto, perfil y reloj de transferencia | Antes del fotón físico y de inferir su hotbelt/latencia | Los ensayos no fotónicos 4A–4B pueden avanzar sin seleccionar ancho. |
| Transferencia instrumental, umbral y estadística de eventos | Comparación cuantitativa de señal experimental, IRF o jitter en etapa 5 | Diagnósticos del dispositivo y cruces del modelo claramente identificados. |

Las comparaciones deben fijarse antes de ver sus resultados y responder a la pregunta del ensayo. Separar error espacial, temporal y espectral es útil; pedir una matriz exhaustiva de todas las 127 entradas no lo es. Si la aparición de un evento cambia al refinar, una pequeña diferencia entre dos tiempos aislados no certifica la misma dinámica. No se proponen nuevos umbrales retrospectivos ni una repetición automática de lotes antiguos.

## 6. Producto inmediato recomendado

Preparar un registro 4A con el escenario material elegido, $\delta=0.10\Delta_0$ como referencia heredada, poblaciones y geometrías de control, consultas dentro del soporte, observables y criterios acordes al diagnóstico. Incluir la comparación algebraica KWT prevista y un plan explícito para parametrizar sus tiempos. La revisión del ancho continúa por la vía de cascada; no se solicita un valor arbitrario para desbloquear ese registro.

La salida debe separar tres dictámenes: **identidades numéricas verificadas**, **sensibilidad a cierres efectivos** y **admisión física del escenario**. Puede haber progreso útil en los dos primeros sin aprobar el tercero. Cuando corresponda preparar dinámica, primero un piloto acotado y estimación de coste; cualquier cálculo previsto por encima de cinco minutos queda para ejecución manual con comando documentado. Esta revisión no crea ni ejecuta ese cálculo.
