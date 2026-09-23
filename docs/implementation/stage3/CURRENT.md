# Estado vigente de la etapa 3

Actualizado el 23 de septiembre de 2026. **Desarrollo cerrado con los límites
autorizados por el usuario; investigación de etapa 3.5 abierta.**

La [decisión y guía de cierre](closure_20260923/README.md) reúne la evidencia
espacial estática y los balances instantáneos de material, potencial, KWT,
calentamiento, circuito de tres estados y reservorio con carga prescrita.
El lote manual de seis instantáneas terminó en **355,651 s**. En el perfil
perturbado, el cambio entre las mallas intermedia y fina es **0,043961%** en
calentamiento total del condensado y **0,96757%** en velocidad material máxima.
Son sensibilidades observadas, no un certificado nuevo de convergencia.
La regresión focalizada aprobó **243 pruebas y 46 subpruebas en 10,29 s**.

El [informe final](closure_20260923/Informe_cierre_desarrollo_etapa_3_y_apertura_3_5.md)
([PDF](../../../output/pdf/implementation/Informe_cierre_desarrollo_etapa_3_y_apertura_3_5.pdf))
distingue este cierre de desarrollo de completar el contrato original:
la admisión temporal, D.27 completa, la interfaz cinética y la promoción a
producción siguen pendientes. Se conservan los indicadores originales falsos
y los resultados históricos, incluido el alcance instantáneo de las pruebas.

Antes de las etapas 4–5 siguen siendo necesarios el radio de reservorio
dependiente de Is y su trabajo, el transporte conservativo a energía común
en la unión 2D–1D, el control de la frontera exterior y una trayectoria débil
acoplada con balances integrados. El empalme actual admite únicamente la traza
transversal uniforme del campo; no se ha demostrado formación de un hotbelt,
latencia ni señal de detección.

La [secuencia vigente](../SECUENCIA_VIGENTE.md) abre la
[investigación 3.5](../stage3_5/README.md). Su
[inventario](../stage3_5/parameter_inventory.json) contiene 127 entradas en
15 familias; el [plan](../stage3_5/research_plan.md) y las
[fuentes Korzh/Allmaras](../stage3_5/reference_experiment.md) separan material,
escenarios, decisiones numéricas, cierres y lectura. No se ejecutan barridos ni se adopta
un rango de confianza. L2D/W = 1,5–6 es un ejemplo para investigar; la suficiencia
de una continuación 1D depende de la evolución transversal relevante y queda
para el ensayo temporal. Un dominio físico de confianza no es un intervalo
estadístico.

**No hay cálculo largo pendiente. No repetir el lote completado.** El chequeo
actual es `sandbox/stage3_spatial/closure_20260923/verify_delivery.py`.
Los comandos útiles están en [la libreta](../../GEMINGA_COMMANDS.md), cuya versión
anterior se conserva en [el archivo de cierre](closure_20260923/GEMINGA_COMMANDS_before_closure.md).
Todo futuro cálculo de más de cinco minutos se entregará al usuario en esa
libreta y mediante un comando exacto copiable en el chat, con propósito,
salidas, recursos y duración estimada. No se inicia automáticamente.

Las iteraciones [de empalme](coupled_20260923/README.md),
[de bordes](ports_20260923/README.md), [nodal](nodal_20260923/README.md)
e [inicial](iteration_20260923/README.md) permanecen como historia.
`README.md` y `entry_contract.json` de esta carpeta conservan la propuesta
del 22 de septiembre; sus estados de preparación no describen el cierre actual.
El tag `v1.0.0` y los módulos de producción no cambian.
