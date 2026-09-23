# Cierre de desarrollo de etapa 3 y apertura de 3.5

La decisión explícita «Cerrar desarrollo e iniciar 3.5» cierra el desarrollo
espacial disponible. No certifica el contrato físico y temporal completo de la
etapa original. La [decisión](closure_decision.json) mantiene pendientes la
dinámica débil acoplada, el reservorio dinámico, toda D.27 y el transporte cinético
en la unión. Producción y `v1.0.0` permanecen sin cambios.

El [informe final](Informe_cierre_desarrollo_etapa_3_y_apertura_3_5.md)
([PDF](../../../../output/pdf/implementation/Informe_cierre_desarrollo_etapa_3_y_apertura_3_5.pdf))
resume los resultados con gráficos. La [auditoría independiente](full_reservoir_review.md)
y su [registro](full_reservoir_review.json) revisan las seis instantáneas guardadas,
las fuentes y los balances; no ejecutan nueva física. El lote del usuario duró
355,65 s. Las [pruebas focalizadas](checks.json) pasan: 243 pruebas y 46 subpruebas.

La [secuencia vigente](../../SECUENCIA_VIGENTE.md) abre la
[investigación 3.5](../../stage3_5/README.md): rangos, márgenes y límites físicos y
numéricos antes de las etapas 4–5. Ningún rango final queda adoptado en esta entrega.
El ejemplo 1,5–6 anchos se estudiará; la reducción 1D exige evidencia transversal
durante un futuro transiente completo.

Verificación ligera, sin consultar de nuevo el modelo físico:

```bash
python sandbox/stage3_spatial/closure_20260923/verify_delivery.py
```

`delivery_manifest.json` fija esta entrega y el verificador recorre la cadena
anterior utilizando las copias exactas de las entradas reemplazadas. El directorio
`raw/` conserva los 61 archivos recuperados. Los flags originales y los fallos
anteriores se preservan; el cierre de desarrollo no reescribe sus dictámenes.

La libreta previa está archivada íntegramente en
[GEMINGA_COMMANDS_before_closure.md](GEMINGA_COMMANDS_before_closure.md).
No hay cálculos largos pendientes para iniciar 3.5.
