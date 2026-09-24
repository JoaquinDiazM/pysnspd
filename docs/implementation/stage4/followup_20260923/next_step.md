# Siguiente paso: revisar el cierre con una referencia reutilizable

La referencia Usadel ya detectó un desacuerdo de fuerza mucho mayor que la
variación numérica observada. El paso útil es examinar **qué respuesta espacial
está perdiendo el cierre**, antes de ejecutar un transiente costoso de núcleo.
El [contrato estructurado](next_step_contract.json) concreta este desarrollo
sin seleccionar todavía un modelo físico nuevo.

Primero se hará reutilizable el cálculo espacial térmico actual, con perfiles,
temperatura, geometría, bordes y corte espectral explícitos. Se comenzará por
el vórtice radial prescrito y los controles uniforme, normal y GL compatibles.
La geometría radial no reemplaza silenciosamente los perfiles 2D anteriores.

Después se comparará el trabajo de las fuerzas frente a variaciones de
amplitud localizadas en el núcleo y en su exterior. Cada candidato conservará
una energía común para sus derivadas. Los cambios se comprobarán con
variaciones independientes y los límites aplicables; pasar esas identidades
no bastará para declarar fidelidad física. Las perturbaciones de fase no
uniformes esperarán una referencia que pueda resolverlas.

No se ajustará δ para hacer coincidir una curva aislada. Se prepararán opciones
de revisión con sus consecuencias físicas y coste. Si elegir entre ellas
requiere ampliar la representación aprobada del dispositivo, se consultará
cuando existan alternativas concretas; esa decisión no impide preparar ahora
el diagnóstico. No se añaden tolerancias universales ni nuevas corridas largas
obligatorias.

Se conservan los contratos válidos de cinética de etapa 2 y el circuito de la
memoria. Un cambio futuro de energía exigirá revisar sus acoplamientos, no
descartar automáticamente todo lo anterior. Siguen pendientes el transiente
espacial débil, los bordes y demás dependencias que utilice, y la admisión
física del núcleo. Producción permanece sin cambios.

La [política de horizonte](horizon_policy.md) autorizada afecta sólo la ventana
de futuros transientes fotónicos: mantiene las ecuaciones y el circuito hasta
la parada. Todavía no está implementada y no restringe este contraste térmico
sin fotón.
