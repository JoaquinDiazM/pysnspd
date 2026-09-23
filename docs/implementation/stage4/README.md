# Etapa 4 iniciada: núcleo y disipación sin fotón

La [primera implementación 4A](start_20260923/README.md) está lista para su campaña
manual. [Estado actual](start_20260923/implementation_status.json) y
[plan registrado](start_20260923/campaign_plan.json). El contrato de entrada
original se conserva como fotografía de la preparación; su “no iniciada” es histórico.

Ahora los tiempos KWT y el parámetro de núcleo son explícitos, con los valores
heredados por defecto. Se corrigieron tres denominadores fijos que habrían roto
la correspondencia entre energía, fuerza y corriente al variar el núcleo.
El nuevo rectángulo usa todos sus nodos 2D, incluidos los bordes, sin continuación
1D ni identificación transversal.

El piloto y las pruebas de rutas modificadas pasaron. Faltan los 40 controles
locales y seis estados 2D que ejecutará el usuario con los comandos de
[/home/jdiaz/GEMINGA_COMMANDS.md](../../GEMINGA_COMMANDS.md).
Son estados estáticos y respuestas instantáneas. No hay todavía un transiente
espacial, núcleo físico admitido, fotón ni circuito acoplado en esta campaña.

Después de recoger esos resultados se decidirá si el núcleo y la disipación
permiten continuar, requieren otra resolución o justifican reformular el cierre.
Los pendientes dinámicos de etapa 3 se resolverán antes del ensayo que los use.
La memoria sigue definiendo el circuito futuro; producción y v1.0.0 permanecen iguales.
