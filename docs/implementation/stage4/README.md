# Etapa 4: núcleo aceptado y reducción cinética en evaluación

Los cuatro núcleos autoconsistentes cumplen el criterio anunciado. La
[revisión actual](self_consistent_review_20260924/README.md) acepta el problema
térmico estático y obtiene su espectro a energías reales y su respuesta
cinética de carga y energía. La etapa 4 permanece abierta hasta
resolver esa unión física; la etapa 5 y producción aún no se activan.

La corrida del usuario terminó en 8,08 min. Los residuos RMS finales son
0,084–0,098 %. El refinamiento de malla cambia la densidad de corriente
0,0798 %; duplicar las frecuencias la cambia 0,6151 %. No se solicita otra
relajación térmica. El informe muestra perfiles, diferencias e historia de
convergencia, y distingue el residuo del núcleo del máximo global.

El nuevo módulo experimental `retarded_spatial_usadel.py` continúa la misma
acción mediante dos campos complejos independientes. Se conservan el entorno
radial de borde y los condensados calculados. Pasaron 72 consultas en 88,48 s.
El operador cinético añade 144 respuestas en 5,13 s, con una identidad conjunta
de fuerza y corriente. La comparación con un potencial por nodo tiene una
cuadratura todavía gruesa: la siguiente prueba concentra energías donde falta
resolución antes de juzgar la reducción física. Sigue pendiente el balance
dinámico de energía al mover el gap.

El [comando vigente](../../GEMINGA_COMMANDS.md) identifica la siguiente prueba,
sus salidas y recursos. Los cálculos de más de cinco minutos siguen siendo
ejecuciones manuales del usuario; cada comando largo se copia también en el chat.
El presupuesto compartido llega como máximo a 28 de 32 hilos en Geminga,
dejando dos núcleos físicos completos libres.

## Evidencia anterior conservada

- [Energía espacial y piloto autoconsistente](spatial_energy_20260924/README.md):
  recupera la fuerza radial con 0,180 % de diferencia al mismo corte.
- [Diagnóstico del cierre local](followup_20260923/README.md): contrasta su
  fuerza con Usadel y documenta por qué no bastaba refinar esa malla.
- [Controles iniciales](review_20260923/Informe_avance_etapa_4A.md): conserva sus valores,
  conclusiones y límites originales.

La infraestructura de la memoria y su circuito completo se mantienen. La
ventana futura hasta el gatillo en Vout más margen cambia sólo el tiempo
observado del mismo dispositivo. Los pendientes cinéticos, de transferencia
fotónica y tasas materiales continúan explícitos en la secuencia vigente.
