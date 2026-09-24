# Modelo previsto y ruta de implementación vigente

Actualizado el 24 de septiembre de 2026. El tag `v1.0.0` y el solver de producción
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

## Preparación vigente antes de la inyección

La siguiente campaña, todavía pendiente, representa una sección interior recta
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
