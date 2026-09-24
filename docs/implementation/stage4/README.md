# Etapa 4: control térmico completado y conexión con la malla final

La [revisión actual](practical_time_review_20260924/README.md) recupera las dos
trayectorias no lineales completas hasta 1 ps. La interrupción ocurrió después,
al fallar una comparación de una señal cruzada casi extinguida; el original
se conserva como certificado incompleto. Los campos principales concuerdan
al nivel medido y el balance instantáneo térmico se mantiene.

Se cierra este ensayo térmico como desarrollo con un límite explícito en la
cola de torque de fase de la sonda de amplitud. **No se repite la campaña.**
La etapa 4 permanece abierta para el acoplamiento no térmico y los puertos.
El [informe](../../../output/pdf/implementation/Informe_etapa_4_revision_practica.pdf)
identifica exactamente campos, normas, unidades, sustracciones y denominadores.

Las [fuentes y ruta numérica](practical_time_review_20260924/published_methods.md)
indican reutilizar el backend Delaunay–Voronoi y el paso KWT de la memoria.
La nueva acción ya es un grafo general: el adaptador convierte áreas y caras
sin reconstruir Usadel. La adaptación temporal requiere convertir fuerza,
movilidad, tiempo y potencial; la corriente también se toma de esa acción.
Una publicación respalda cada método, pero no sustituye comprobar las unidades
y transferencias entre bloques al conectarlos.

La [ruta física](practical_time_review_20260924/physical_route.md) conserva el
trabajo espectral no térmico como tarea constitutiva. Los
[criterios de implementación y figuras](../REPORTES_Y_CRITERIOS.md) priorizan
el sistema final y las decisiones útiles, sin nuevas baterías extremas por defecto.
El [cuaderno](../../GEMINGA_COMMANDS.md) ya no pide otra corrida larga.

## Evidencia anterior conservada

- [Trayectoria afín, corrección de Newton y preparación ETD2](time_review_20260924/README.md).
- [Momentos de carga, respuesta armónica y operador térmico](moment_review_20260924/README.md).
- [Núcleos autoconsistentes y espectros](self_consistent_review_20260924/README.md).
- [Energía espacial](spatial_energy_20260924/README.md).
- [Diagnóstico del cierre local](followup_20260923/README.md).

Producción y v1.0.0 siguen intactos. El circuito es el de tres estados de la
memoria. Para el futuro fotón se conserva el mismo sistema hasta V_out más
margen; transferencia fotónica y tasas materiales siguen abiertas. No hay
predicción de hotbelt o retardo experimental ni inicio de etapa 5.
