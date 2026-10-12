# Modelo previsto y ruta de implementación vigente

## Decisión actual de desarrollo — 11 de octubre de 2026

La corrida acoplada microscópica de primer gradiente terminó y se rechaza como arranque temporal: al reducir el paso no convergen los incrementos del gap en frío y se pierde soporte Pauli. La auditoría se conserva en /home/jdiaz/scratch/stage51_failure_audit_20261011/audit.json. No se repiten por defecto matrices densas con pasos menores.

Con la autorización actual de adoptar hipótesis útiles y ejecutar corridas largas, se desarrolla un cierre **neutral efectivo**: estados BCS locales en coordenadas nativas fijas, energía variacional común con rigidez dirty radial, corriente derivada del gauge, movilidad KWT calibrable, eventos electrón–fonón y transporte conservativos, potencial cuasiestático y circuito CM completo. La rigidez radial a temperatura de referencia fija es una hipótesis mesoscópica explícita; no equivale al espectro Usadel espacial ni a la expansión Moyal.

Se preservan la malla Delaunay–Voronoi y los parámetros físicos del escenario. El avance temporal implícito se introduce por rigidez y debe acreditar su propia convergencia frente a refinamiento; no se atribuye esa evidencia al Euler histórico. Cortes superconductores rígidos fijan el módulo del nuevo DC y evolucionan su fase mediante Josephson; sus distribuciones térmicas son reservorios cuyo flujo se cuenta explícitamente. El alcance corto se comprobará con longitud del dominio.

Estado: **implementación y validación en curso; etapa 5 sin cerrar todavía**. Las nuevas referencias DC satisfacen corriente, fuerza y tensión basal en cuatro geometrías. No se reutiliza la antigua admisión Usadel como prueba del nuevo cierre. Los nominales Allmaras 2020 tau_ee(Tc)=5 ps y tau_ep(Tc)=24,7 ps siguen siendo calibradores libres. Material acústico, retención 0,667 y anchura fotónica exploratoria conservan su condición de hipótesis sin calibración experimental.


Actualizado el 8 de octubre de 2026. El tag `v1.0.0` y el solver de producción
conservan la implementación de la memoria. Lo siguiente identifica el modelo
previsto y sus pruebas experimentales posteriores.

La base es el [sistema continuo 0.4](../modelo_v0_4/D_sintesis_plan_y_verificaciones_v0_4.md),
con esta sustitución normativa solicitada por el usuario:

- **Circuito:** [red de polarización y lectura de la memoria](../modelo_v0_4/actualizaciones/circuito_memoria_20260922.md).
  Sus tres variables son las corrientes de polarización y del detector y el
  voltaje del capacitor. Sustituye el circuito ideal simplificado de C/D,
  sus balances y su inicialización. La adenda contiene las ecuaciones completas,
  las convenciones de signo y la partición de inductancia que evita doble conteo.
- **Preparación y escenario del detector:** el
  [contrato DC interior del 24 de septiembre](stage5/prephoton_dc_20260924/physical_scope.md)
  sustituye los contactos de corriente cero para la siguiente ventana espacial
  y adopta una inductancia exterior fija motivada por el dispositivo de Korzh.
  Conserva la topología CM, con fuente DC, sin $L_k(t)$ ni jitter longitudinal.
- **Implementación:** el [cierre de desarrollo de etapa 2](stage2/closure_20260922/closure_decision.json),
  autorizado por el usuario, acepta la evidencia disponible para continuar con
  la secuencia espacial, sus bordes y el circuito. No cambia retrospectivamente
  las tolerancias ni convierte fallos en pases.

La etapa 2 queda cerrada como desarrollo. El
[informe final](stage2/closure_20260922/Informe_cierre_etapa_2.md)
([PDF](../../output/pdf/implementation/Informe_cierre_etapa_2.pdf)) reúne las
13 trayectorias del último lote, los pases temporales de una/dos celdas y las
limitaciones medidas. Su [certificado numérico estricto](stage2/stage2_admission.json)
sigue incompleto: se conserva el fallo auxiliar y no se atribuye convergencia
a las mallas dinámicas pendientes.

La etapa 3 se cerró como desarrollo con pendientes dinámicos explícitos, y la
investigación 3.5 se cerró dentro de su alcance autorizado. La
[etapa 4 sin fotón](stage4/closure_20260924/README.md) se cerró como desarrollo: los controles de campo
prescrito permiten separar sensibilidad de malla y de cierre constitutivo.
El cierre local regularizado falló el contraste radial y se conserva como
antecedente rechazado para esa dinámica. La nueva representación espacial
recupera dicho contraste y produce sus propios núcleos térmicos. La
[revisión autoconsistente](stage4/self_consistent_review_20260924/README.md)
acepta los cuatro casos y obtiene su espectro a energías reales y su respuesta
cinética lineal de carga y energía.
Siguen pendientes las capacidades dinámicas antes de interpretar el detector.
La transferencia fotónica, su ancho y la normalización de tasas NbN no se
admiten mediante este avance estático.

## Etapa 5.1 activa - 9 de octubre de 2026

La corrida conjunta de 4940 columnas termino. La SVD truncada estaba dominada por escala; el QR equilibrado recupera rango 4939/4939, pero tasas y voltajes son inadmisibles y el residuo original no pasa. No se cierra el modelo. Se implemento un primer paso Rosenbrock-Euler conjunto con rigidez fisica, M+dt Lh, puertos/cortes moviles y CM; incluye escape fononico BE con parametro efectivo libre de 15 ps heredado D.19. Un piloto real de 112 energias y 256 modos termino en 99 s. La unica corrida pendiente compara 1 y 0.5 fs y reconstruye R0/Mats en los estados finitos, con cero pasos aceptados de antemano. Comando: /home/jdiaz/GEMINGA_COMMANDS_PRESENT.md; plan versionable: sandbox/stage5_photon/coupled_first_step_plan.json. El proveedor Born general, la anchura fisica en R0 y el limite temporal siguen pendientes; 5.2 requiere energia completa (no solo observador neutral) y 5.3 admision. No hay nuevos informes ni transiente/reset completo.

## Decisión tras el RHS fotónico inicial completo del 8 de octubre

El lote de 224 consultas termina correctamente en 45,8 s, pero no admite cerrar la etapa 5 ni su publicación como modelo concluido. El trabajo K1 es el 92,59/89,40 % del K0 para 775/1550 nm; el defecto independiente es 14,271/8,725 nW. Las derivadas Moyal dominan la corrección, y faltan neutralidad/gauge, los puertos de corriente y el balance temporal conjunto. La [decisión vigente](stage5/completed_photon_rhs_20261008/closure_decision.json) conserva las hipótesis de material/fuente y márgenes 1 %/2 %. No repetir DC, UV ni los controles del actuador; no corregir con clipping, corriente óhmica añadida o calor igual al residuo.

La siguiente iteración es de implementación: cierre microscópico inicial G0+G1=0, Gamma_res=0 como candidato explícito, con KWT sólo como precondicionador numérico y no una segunda disipación física. No se presume bien condicionado ni admitido; requiere cortes superconductores dinámicos, neutralidad y puerto CM del mismo sistema, energía al mismo orden y fuentes consistentes. No hay un photon-run del nuevo backend que pueda ejecutarse rigurosamente todavía. Las pipelines nuevas revisan y plotean los datos iniciales guardados, con cero pasos físicos. Se conserva provisionalmente la ventana de 200 ps sin gatillo y 5 ps extra, umbral positivo 0,1 mV, sin reset.

El [informe de esta revisión](../../output/pdf/implementation/Informe_revision_RHS_fotonico_completo_20261008.pdf) reúne cuatro páginas de plots físicos y la lectura SNSPD. Las pipelines `pipelines/04_initial_photon_rhs_review.py` y `plot_pipelines/04_plot_initial_photon_rhs_review.py` verificaron los datos guardados y producen diez figuras, sin evolución temporal; el [recibo](stage5/completed_photon_rhs_20261008/pipeline_review_receipt.json) registra unidades, comprobaciones y recursos. No hay cálculo pesado pendiente del usuario.

## Límite anterior tras las cinco campañas (histórico)

Los controles DC, transporte y cuadratura del saldo examinado están cerrados en su alcance. No se solicitan nuevos barridos estáticos para perseguir más cifras. Se mantiene el material acústico efectivo, el depósito gaussiano y la fuerza térmica finita F256 de la referencia DC; la cola analítica queda como observador separado. La implementación fotónica general aún no está completa: el siguiente ensayo usa el RHS del estado inicial real de 775/1550 nm y mide los defectos de neutralidad, fuerza y trabajo, sin simular una trayectoria o un reset total. La [decisión actual](stage5/final_review_20261008/closure_decision.json) fija márgenes y limita las siguientes iteraciones a defectos determinantes del acoplamiento.

El piloto inicial de 775 nm midió ||R1||F/||R0||F=1,218 cerca del gap (máximo nodal 4,469). El lote posterior de 224 consultas ya terminó y su revisión vigente está al comienzo de este documento; no repetir ese comando. La razón local no se interpreta como error de latencia y no hubo avance temporal.

## Continuación anterior del 8 de octubre (histórica; sustituida por la revisión superior)

El usuario ya ejecutó las cinco auditorías RHS solicitadas. Su revisión detecta falta de convergencia de cuadratura en el borde de~0,804 meV, aunque la extensión UV 6 a 12 meV es estable para ese ensayo. No se declara completado el balance por residuos proyectados pequeños. Hay dos repeticiones de transporte necesarias tras reparar el soporte del paso; no se repiten los casos h5 que pasaron ni la preparación DC admitida.

La [revisión actual de etapa 5](stage5/README.md) documenta los nuevos momentos de Nambu, contraste estacionario, regulador causal, material efectivo, fuente y snapshots. Material acústico efectivo y ancho gaussiano exploratorio de 7,549834435 nm fueron autorizados el 8 de octubre; [el contrato físico](stage5/photon_readiness_20261008/source_and_observation_contract.md) separa la aproximación de una calibración NbN. Esa autorización sustituye la decisión de mantener el ancho abierto sólo para esta familia exploratoria; no admite una cascada óptica calculada.

El modelo de producción conserva su versión. El nuevo backend fotónico acoplado todavía requiere cierre de trabajo, neutralidad y fuentes sobre un único RHS. Las siguientes campañas del usuario son controles para resolver ese cierre; no son comandos de reset total ni corridas fotónicas de producción reetiquetadas.

## Preparación vigente antes de la inyección

La campaña de preparación DC ya fue ejecutada y admitida. La
[revisión de resultados](stage5/dc_admission_20261005/README.md) verifica cuatro
referencias y cinco holds de 10 ps, con comparaciones de malla, longitud y paso.
No debe repetirse por instrucciones históricas. Su región representa una sección interior recta
de la cinta de 80 nm. Su estado inicial debe tener $|\Delta|$ uniforme con la
depresión correspondiente a la corriente, fase $qx$, espectro autoconsistente,
distribuciones térmicas y potencial uniforme. Los extremos numéricos continúan
ese mismo estado; no son contactos metálicos ni reservorios de corriente cero.
Se conservan la malla dual Delaunay-Voronoi y el avance KWT Euler heredado.

La [preparación DC](stage5/README.md) compara longitudes a igual corriente y
verifica que el circuito no genere señal basal espuria. Su fuente es constante,
$V_b=R_bI_{\rm DC}$, y mantiene los tres estados CM. Para una ventana de longitud
$L$ se fija una sola vez

$$L_{k,\rm ext}^{\rm ref}=96\,\mathrm{nH}
+(5\,\mathrm{\mu m}-L)\frac{\hbar}{2e}
\left.\frac{dq}{dI}\right|_{I_{\rm DC}}.$$

Los 96 nH son una estimación publicada del inductor exterior añadido, no el
total de la rama. La contribución del resto activo usa el mismo equilibrio
material que la región resuelta; no se recalcula durante el pulso. El escenario
histórico de 10 nH y las respuestas AC de etapa 4 no se reinterpretan. El
[registro de fuentes y decisiones](stage5/prephoton_dc_20260924/source_decisions.json)
separa valores publicados, conversiones y contribuciones no cuantificadas.

La coordenada longitudinal de absorción permanecerá fija. La futura familia
de posiciones estudia el ancho; se excluye el jitter geométrico longitudinal.
No se requieren nuevas fuentes AC ni perturbaciones artificiales intermedias.
El ensayo sin fotón verifica la preparación y conservación DC, no la cinética
general de una población excitada. En particular, ausencia de calentamiento DC
no resuelve el balance no lineal de energía pendiente antes del fotón.

## Ventana inicial de observación del mismo dispositivo

La prioridad acordada es el retardo relativo de 775 frente a 1550 nm en la cinta
de 80 nm de Korzh (2020). Los futuros transientes se seguirán hasta un cruce
confirmado del gatillo en Vout más un margen declarado; no se exigirá primero
la recuperación completa durante nanosegundos. La
[política de horizonte](stage4/followup_20260923/horizon_policy.md) define el
observable, el tratamiento de corridas sin cruce y las decisiones aún abiertas.

**Sólo cambia el horizonte de integración/observación.** Ecuaciones, material,
depósito, geometría, malla, bordes, circuito completo y política de precisión
permanecen iguales entre una corrida corta y su extensión. No se eliminan
variables lentas ni se acortan constantes de tiempo. Se conserva el estado
final para continuar. Umbral y margen se especificarán para el ensayo; no se
adoptan automáticamente los valores predeterminados de producción.

El modo heredado `latency` espera un máximo confirmado y todavía no implementa
la regla gatillo más margen. Esta adenda registra la política futura; no cambia
el solver de producción. El reloj sigue comenzando en la transferencia al modelo,
con el retardo óptico previo desconocido por separado.

Antes de una predicción dinámica se fija qué diferencia del dispositivo se
quiere resolver y se registra su presupuesto numérico. Las identidades de
energía, unidades, dominios y signos siguen siendo controles necesarios; una
precisión numérica alta no acredita por sí sola el modelo físico del detector.

Los C/D y PDF históricos se mantienen íntegros. Cuando se solicite su siguiente
edición integrada, deberá incorporar esta adenda; no debe recuperarse el circuito
simplificado como contrato vigente por leer un PDF anterior.
La [revisión previa del alcance](stage2/practical_review_20260922/validation_scope_review.md)
se conserva como antecedente de esta decisión; su preparación limitada a 3A
no sustituye la secuencia completa vigente.

## Referencia térmica espacial regular

La [entrega espacial](stage4/spatial_energy_20260924/README.md) implementa
un oráculo térmico con campos espectrales complejos, sin completado δ ni
gradiente K0 añadido. Recupera el contraste radial y deriva fuerza y corriente
de la misma energía. Los cuatro casos de la campaña autoconsistente cumplen el criterio declarado.
No sustituye el cierre de poblaciones no térmicas: la rama retardada y el
operador cinético congelado ya están implementados, pero falta el balance
dinámico de trabajo espectral, transporte y disipación. Producción sigue intacta.


La [decisión actual](stage4/self_consistent_review_20260924/acceptance_decision.json)
cierra únicamente el problema estático térmico al nivel medido, sin otra
refinación. Su continuación retardada usa coordenadas complejas independientes
y el mismo entorno radial de borde. Es una transformación de la acción
espectral, no una sustitución de distribuciones por temperatura. El
[puente físico pendiente](stage4/self_consistent_review_20260924/physics_review.md)
ya obtiene incrementos de fuerza y corriente con una identidad común de
conversión de carga. La [revisión de momentos y frecuencia](stage4/moment_review_20260924/README.md)
resuelve la discrepancia de la proyección a un potencial por nodo: 20,55 % en
corriente y 77,53 % en torque de fase, frente a variaciones numéricas de 2,87 %
y 3,64 %. No se admite esa compresión energética para este ensayo. La eliminación
algebraica espectral completa conserva el detalle energético y aproxima bien la
respuesta lenta medida; no acredita todavía el régimen ultrarrápido del detector.
Tampoco se rechaza por este control congelado toda la fase y el potencial de la
memoria. El operador temporal térmico ya se contrastó con 2304 raíces espectrales;
la [trayectoria afín completada](stage4/time_review_20260924/README.md) mantiene
la deriva de referencia y pasa el contraste temporal. El posterior contraste
no lineal también pasa tras corregir cancelación de redondeo en Newton. Sigue
la conexión con la malla dual y el paso KWT heredados. La
[trayectoria no lineal completada](stage4/practical_time_review_20260924/README.md)
permite cerrar ese ensayo como desarrollo con el límite explícito del torque
secundario tardío; se conserva su certificado temporal incompleto. No se pide
otra repetición cartesiana ni otro integrador para resolver esa cola. Falta unir trabajo y transporte no térmicos antes de evaluar
un fotón; el circuito completo y el horizonte acordado siguen vigentes.

## Decisión posterior al ensayo Euler dual

El [resultado actual](stage4/final_kwt_20260924/README.md) admite el paso KWT
heredado sobre la malla dual: 1 ps en 153,81 s, sondas de amplitud/fase y medio
paso dentro del margen práctico del 2 %. La pequeña discrepancia del ensayo
anterior también se acepta al margen revisado y deja de ser trabajo pendiente.

Queda explícita una distinción que no resuelve la elección del integrador:
el sistema 0.4 usaba un espectro local adiabático, mientras la referencia
actual incluye respuesta espectral espacial. Sus sectores térmico y
longitudinal débil están comprobados, pero la unión no térmica general no se
obtiene sumando sin más sus fuerzas y balances. La
[investigación solicitada](stage4/final_kwt_20260924/quasiclassical_assessment.md)
separa aproximación cuasiclásica, límite difusivo, espectro adiabático y
relajación de las poblaciones. Recomienda mantener Usadel difusivo y la
evolución de las distribuciones, sin atribuir a la neutralidad eléctrica
un ajuste instantáneo del desequilibrio electrón–hueco. El espectro adiabático
sigue siendo una aproximación independiente: depende de la amplitud y rapidez
local del cambio, no sólo de la duración total del ensayo. Esta investigación
no modifica todavía las ecuaciones ejecutadas ni sustituye el cierre pendiente
de trabajo espectral y calor. Posteriormente se preparó el
[lote polarizado de unión dinámica](stage4/coupled_closure_20260924/README.md),
que terminó con dos referencias y cinco respuestas. Conserva la evolución de ambas
distribuciones y usa los desplazamientos energéticos exactos de la respuesta
armónica. Mide el posible solapamiento con la relajación KWT antes de adoptar
un cierre radial. El calor no lineal y la preparación fotónica siguen fuera
de la aceptación de este control débil.

## Cierre de desarrollo de etapa 4 y puerta de entrada a etapa 5

La [decisión vigente](stage4/closure_20260924/closure_decision.json) admite
el control reactivo condicionado, no la disipación absoluta ni el contrato
dinámico general. El refinamiento cambia 1,589 % la admitancia compleja, pero
45,66 % su parte real respecto al refinado. El calor radial candidato supera
la potencia de puerto en 31,37 %; falta el balance independiente con trabajo
DC y reservorios. La continuidad de corriente también exige revisar cortes
interiores, no sólo terminales y residuos proyectados.

No se modifican ecuaciones físicas para forzar esta aceptación. La
[entrada a etapa 5](stage5/README.md) concreta las interfaces aún necesarias
antes de simular el fotón. La transferencia, ancho gaussiano y tasas NbN siguen
sin admitirse. La evidencia débil actual usa movilidad KWT heredada, no una
calibración completa del experimento Korzh.


## Continuación posterior a la preparación DC admitida

La etapa 5 completa sigue abierta. El pase DC conserva poblaciones térmicas
prescritas y cierre óhmico heredado; no completa calor no lineal ni tasas
absolutas NbN. La [revisión vigente](stage5/dc_admission_20261005/README.md)
identifica implementación disponible, pendientes físicos y evidencia del servidor.
La [continuación del 5 de octubre](stage5/continuation_20261005/README.md)
prepara energía interna independiente y controlador gatillo más margen, sin
atribuir a esas utilidades la integración fotónica general todavía pendiente.
No se propone una recuperación eléctrica completa como ensayo inicial.

El techo solicitado sin gatillo es 200 ps desde la transferencia; un cruce
confirmado permite completar su margen posterior. Se preserva la dirección
global DC mediante parada con diagnóstico ante inversión, sin modificar corrientes
ni circuito. El rango provisional de estudio de Vout es 0,1–0,5 mV; no se
atribuye a un umbral experimental calibrado.


## RHS espacial y balance independiente — actualización del 5 de octubre

Se implementaron el almacenamiento disperso y RHS temporal de ambos modos
Usadel, el sector cinético nodal conservativo y la derivada independiente de
energía interna. La evaluación de 256 modos sobre la referencia DC215 real
verifica la derivada con diferencias finitas: error relativo 8,6425e-7.
La validación integrada terminó con 145 pruebas y 13 subpruebas. Se corrigió
además la pérdida de parte imaginaria del Jacobiano retardado con gap de tipo real;
la preparación DC usa el sector térmico y no debe repetirse por esa corrección.

El informe principal es el [PDF con plots físicos](../../output/pdf/implementation/Informe_etapa 5_RHS_balance.pdf).
La [derivación](stage5/rhs_balance_20261005/rhs_energy_derivation.md) explica
energía, trabajo, fuente y alcance. `/home/jdiaz/GEMINGA_COMMANDS_PRESENT.md`
contiene cinco auditorías manuales: espectro móvil, convergencia de energía,
transporte temporal, cuadratura a igual corte y extensión del corte UV.

No se admite todavía el detector no térmico completo: el pase de proyecciones
cinéticas no acredita la self-energy anómala, neutralidad, fuerza y trabajo
consistentes ni tasas NbN. La consulta real conserva un residuo Nambu sin esas
fuentes. Se mide ese defecto y se separan los errores numéricos antes de
acoplar un fotón; no se compensa mediante calor residual ni se repite DC.
Se mantienen el circuito, sentido DC y horizonte de 200 ps sin gatillo más
margen tras el cruce confirmado. No se preparó un reset eléctrico completo.
