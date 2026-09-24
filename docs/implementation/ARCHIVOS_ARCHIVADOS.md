# Archivos archivados al cierre de la etapa 4

La entrega vigente conserva los informes de cierre de las etapas 1, 2, 3 y 3.5,
el informe final de la etapa 4 y los datos físicos que sustentan sus conclusiones.
Se retiran del árbol actual **19 PDF intermedios y 15 trayectorias NPZ repetidas**:
128 494 498 bytes (122,54 MiB). Las trayectorias retiradas son copias exactas,
verificadas por SHA256, de archivos que permanecen en el repositorio.

El [manifiesto de archivo](archive_20260924/manifest.json) identifica cada ruta,
tamaño, hash y copia canónica, cuando corresponde. Todos los archivos también
pueden recuperarse del commit
[`c80c0f8612d8fd5a2b390938ca5ed4b9d58d3caa`](https://github.com/JoaquinDiazM/pysnspd/tree/c80c0f8612d8fd5a2b390938ca5ed4b9d58d3caa).
Los enlaces de los informes históricos apuntan a esa versión. Se conservan sus
fuentes Markdown y los generadores; los snapshots llamados
`previous_delivery_exact` permanecen sin modificar.

Esta limpieza reduce los archivos del directorio de trabajo. No reescribe el
historial de Git, no reduce sus objetos históricos y no modifica `v1.0.0`.
Tampoco retira las figuras únicas, el catálogo electrónico, la malla dual usada
por los ensayos actuales ni los datos de regresión de Newton. Los archivos del
cuaderno E y de la reunión de Virtanen quedan fuera de esta limpieza.

## Recuperación para auditorías anteriores

Algunos verificadores históricos esperan un PDF o una trayectoria en su ruta
original. Antes de ejecutarlos, restaurar los archivos del manifiesto desde la
raíz del repositorio:

```bash
python sandbox/maintenance/restore_archived_artifacts.py --all
```

Para listar las rutas o recuperar sólo una de ellas:

```bash
python sandbox/maintenance/restore_archived_artifacts.py --list
python sandbox/maintenance/restore_archived_artifacts.py --path output/pdf/implementation/Informe_etapa_4_Euler_y_acoplamiento.pdf
```

El restaurador utiliza primero la copia canónica y después el objeto local de
Git. Verifica SHA256 y tamaño, sólo admite rutas del manifiesto dentro del
repositorio y rechaza sobrescribir contenido diferente. No descarga archivos ni
necesita volver a ejecutar una simulación. Una clonación superficial debe incluir
el commit de archivo antes de recuperar los PDF.

Una auditoría integral de una entrega antigua que además compruebe hashes de la
documentación debe ejecutarse en una copia de su commit correspondiente. Restaurar
estos artefactos no revierte los enlaces documentales actualizados por la limpieza.

En particular, `audit_one_cell.py` de `time_pass_20260922` y `audit_two_cell.py`
de `closure_prep_20260922` requieren restauración porque leen directorios de
trabajo anteriores. Los datos canónicos de cierre permanecen en
`stage2/practical_review_20260922/raw` y `stage2/review_20260922/raw`.

La copia de resguardo adicional en Geminga está en
`/home/jdiaz/scratch/pysnspd_archive_20260924/superseded_artifacts.tar.gz`,
verificada en [su recibo](archive_20260924/geminga_archive_receipt.json).
En Windows un PDF antiguo abierto por otra aplicación puede permanecer como
copia local ignorada aunque ya no esté en el índice de Git; el manifiesto
identifica los 34 artefactos retirados de la versión publicada.
