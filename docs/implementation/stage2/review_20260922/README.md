# Etapa 2: revisión del 22 de septiembre

**No cerrada.** El nuevo lote RK4 repitió la convergencia temporal aprobada en la
malla candidata y el fallo de positividad en la malla electrónica fina. La
recuperación SSP40 completó su trayectoria, pero su error energético fue
3,65224e-6 frente al límite 1e-7. El limitador no actuó en ese ensayo.

El [informe de revisión](../../../../output/pdf/implementation/Informe_revision_etapa_2_20260922.pdf)
presenta cuatro páginas de resultados y plots, sin emitir un certificado de
cierre. La [versión editable](Informe_revision_etapa_2_20260922.md) conserva
el mismo contenido.

- `raw/`: 73 archivos recuperados, con hashes verificados contra Geminga.
- `import_inventory.json`: rutas originales, tamaños y hashes de la copia.
- `independent_audit.json`: 118 verificaciones de procedencia y resultados;
  dictamen de cierre NO.
- `diagnosed_ssp.json`: diagnóstico por posprocesamiento, sin dinámica nueva.
- `figures/` y `report_content.json`: plots y contenido reproducible del informe.
- `GEMINGA_COMMANDS_before_cleanup.md`: cuaderno anterior íntegro. Sus antiguos
  lanzamientos de etapa 2 no representan una cola de trabajo pendiente.

El único siguiente cálculo preparado es `manual_time_plan.json`: una celda
SSP con 160/320/640 pasos frente a 1280, estimado en 26 minutos. Se ejecuta
manualmente desde el cuaderno activo de Geminga; no se lanzó en esta revisión.
Su aprobación no sustituye los controles posteriores de dos celdas, mallas,
fronteras de población, campos y soporte. Las tolerancias originales, el
catálogo R2 y el tag v1.0.0 permanecen intactos.

El archivo `selected_continuum.json` vuelto a ejecutar sólo difiere en su tiempo
de ejecución. La nueva copia se conserva en `raw/diagnostics/`; la entrega
histórica mantiene sus bytes originales.
