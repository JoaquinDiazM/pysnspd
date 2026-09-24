# Etapa 3 iniciada: funcional espacial y campaña estática

23 de septiembre de 2026. Estado: **piloto aprobado; campaña espacial pendiente
de ejecución del usuario**. La etapa 3 completa permanece abierta.

El bloque 3A calcula energía, fuerza cartesiana y corriente conjugada desde el
mismo funcional discreto en una tira periódica. Incluye la dependencia del
espectro respecto del flujo regularizado y el diagnóstico local de estabilidad
D.36. El nuevo módulo es experimental; aún no integra bordes físicos, potencial
eléctrico, circuito ni dinámica espacial.

## Resultados obtenidos

- 74 pruebas del funcional, las barras de progreso y la reutilización del piloto
  aprobadas en 4,37 s en Geminga.
- Piloto de ocho celdas aprobado en 33,9 s. Diferencias absolutas respecto de
  variaciones independientes de energía: fuerza 7,08e-12 y corriente 3,42e-13.
- Invariancia de fase global y referencia helicoidal verificadas.
- El control negativo conserva el valor propio negativo esperado (-0,30662)
  y se rechaza para evolución.

El [informe con figuras](Informe_inicio_etapa_3_20260923.md)
([PDF](https://github.com/JoaquinDiazM/pysnspd/blob/c80c0f8612d8fd5a2b390938ca5ed4b9d58d3caa/output/pdf/implementation/Informe_inicio_etapa_3_20260923.pdf))
distingue estos resultados de las verificaciones pendientes. La corriente del
perfil prescrito no tiene por qué ser constante: no es todavía una solución
estacionaria con conservación de corriente total.

## Próximo cálculo

La [prerregistración](registration.json) fija seis casos y mallas de 8, 16 y
32 celdas. El contraste espectral 630/1260 se limita a puntos representativos
de los casos ocupados de 32 celdas. El objetivo espacial del 1 % se aplica a
las respuestas registradas; las referencias incluyen una estimación de error
y orden observado, sin imponer una tolerancia temporal de transitorio completo.

Duración estimada: unos 16 minutos (orientativamente 12-25), un proceso CPU y
menos de 1 GB de RAM estimado. El agente **no lanzó** este lote. El bloque activo
de `/home/jdiaz/GEMINGA_COMMANDS.md` incluye la orden, salidas y recursos.
La terminal muestra barras por caso y lote, tiempo transcurrido y ETA aproximada.
`progress.jsonl` conserva el avance. El piloto se reutiliza solo si coinciden
fuentes, condiciones, malla y recibos de integridad.

Al terminar, la revisión decidirá qué ensayos pueden pasar a **3B: bordes y
reservorios**. Después siguen 3C (carga y potencial), 3D (las tres ecuaciones
circuitales de la memoria) y 3E (dinámica débil admitida). No se sustituye el
circuito de la memoria por una fuente ideal.

## Evidencia y reproducción

- [Pruebas](checks/result.json), [resumen del piloto](pilot/summary.json) y
  [entrega pendiente al usuario](handoff.json).
- Fuentes: `pysnspd/experimental/spatial_functional.py` y
  `sandbox/stage3_spatial/`; las versiones quedan identificadas por SHA-256.
- `verify_delivery.py` comprueba integridad y correspondencia del piloto;
  no ejecuta física ni convierte una entrega íntegra en validación completa.
- El cierre de etapa 2 y la propuesta original de etapa 3 se conservan intactos.
  La nueva implementación no activa el modelo en producción ni cambia `v1.0.0`.
